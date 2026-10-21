"""Terminal adapter and primer removers extracted from the thesis prototype notebooks.

Matching logic and default thresholds are unchanged from bit-vector.ipynb, sassy.ipynb,
semi-global.ipynb and adapter-removal.ipynb; only module-level globals became explicit settings.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import multiprocessing as mp
import re
import time
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from contextlib import ExitStack
from dataclasses import dataclass, replace
from itertools import islice
from pathlib import Path
from typing import Any, Protocol

from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

REPOSITORY_ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class FastqRecord:
    title: str
    sequence: str
    qualities: str


@dataclass(frozen=True)
class Primer:
    name: str
    sequence: str


@dataclass(frozen=True)
class EndCall:
    """The primer or adapter construct accepted at one read end."""

    name: str
    aligned_bases: int
    cost: int | None = None
    detail: str = ""


@dataclass(frozen=True)
class TrimResult:
    record: FastqRecord
    original_length: int
    left_removed: int
    right_removed: int
    left: EndCall | None
    right: EndCall | None


class Trimmer(Protocol):
    name: str

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]: ...


def read_fastq(path: Path | str) -> Iterator[FastqRecord]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for title, sequence, qualities in FastqGeneralIterator(handle):
            yield FastqRecord(title, sequence.upper(), qualities)


def read_primers(path: Path | str) -> list[Primer]:
    primers = [Primer(record.id, str(record.seq).upper()) for record in SeqIO.parse(path, "fasta")]
    if not primers:
        raise ValueError(f"No primer records found in {path}")
    return primers


def split_primers(primers: Sequence[Primer]) -> tuple[tuple[Primer, ...], tuple[Primer, ...]]:
    left_primers, right_primers = [], []
    for primer in primers:
        if "_LEFT" in primer.name:
            left_primers.append(primer)
        elif "_RIGHT" in primer.name:
            right_primers.append(primer)
        else:
            raise ValueError(f"Primer name has no LEFT/RIGHT orientation: {primer.name}")
    if not left_primers or not right_primers:
        raise ValueError("Expected at least one left and one right primer")
    return tuple(left_primers), tuple(right_primers)


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def _result(
    record: FastqRecord, left_cut: int, right_cut: int, left: EndCall | None, right: EndCall | None
) -> TrimResult:
    """Build a result from absolute cut coordinates; the retained read is [left_cut, right_cut)."""
    length = len(record.sequence)
    trimmed = FastqRecord(record.title, record.sequence[left_cut:right_cut], record.qualities[left_cut:right_cut])
    return TrimResult(trimmed, length, left_cut, length - right_cut, left, right)


# --------------------------------------------------------------------------------------
# Myers bit-vector (bitvector_longest)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class BitvectorSettings:
    max_error_rate: float = 0.22
    min_aligned_primer_fraction: float = 0.45
    min_aligned_primer_bases: int = 11
    overhang_penalty: float = 0.35
    offset_penalty: float = 0.0
    max_terminal_slack: int = 180
    max_terminal_slack_fraction: float = 0.35
    seed_length: int = 7
    window_fraction: float = 0.35
    boundary_slack_threshold: int = 30
    embedded_max_error_rate: float = 0.15
    embedded_min_aligned_primer_fraction: float = 0.70
    embedded_min_aligned_primer_bases: int = 16
    truncated_max_error_rate: float = 0.10
    pre_seed_indel_allowance: int = 2
    exact_truncated_min_primer_bases: int = 10
    enable_tandem_trim: bool = True
    tandem_max_gap: int = 40
    tandem_max_overlap: int = 12
    tandem_max_error_rate: float = 0.15
    tandem_min_aligned_primer_fraction: float = 0.60
    tandem_min_aligned_primer_bases: int = 18
    tandem_min_retained_bases: int = 50
    enable_cognate_rescue: bool = True
    cognate_truncated_min_primer_bases: int = 8
    cognate_truncated_max_error_rate: float = 0.22
    cognate_terminal_noise_bases: int = 5
    cognate_terminal_noise_min_primer_bases: int = 10


@dataclass(frozen=True)
class PreparedPrimer:
    primer: Primer
    orientation: str
    sequence: str
    eq_table: dict[str, int]


@dataclass(frozen=True)
class BitvectorMatch:
    primer: Primer
    orientation: str
    start: int
    end: int
    cost: int
    aligned_bases: int


def build_eq_table(pattern: str) -> dict[str, int]:
    eq = {base: 0 for base in "ACGTN"}
    for index, symbol in enumerate(pattern.upper()):
        for base in ambiguous_dna_values.get(symbol, ""):
            eq[base] |= 1 << index
    eq["N"] = (1 << len(pattern)) - 1
    return eq


def prepare_primers(
    primers: tuple[Primer, ...], orientation: str, reverse_for_prefix_search: bool
) -> tuple[PreparedPrimer, ...]:
    prepared = []
    for primer in primers:
        sequence = reverse_complement(primer.sequence) if orientation == "reverse_complement" else primer.sequence
        if reverse_for_prefix_search:
            sequence = sequence[::-1]
        prepared.append(PreparedPrimer(primer, orientation, sequence, build_eq_table(sequence)))
    return tuple(prepared)


def build_seed_index(
    prepared_primers: tuple[PreparedPrimer, ...], seed_length: int
) -> dict[str, tuple[tuple[int, int], ...]]:
    index = defaultdict(list)
    for primer_index, prepared in enumerate(prepared_primers):
        for offset in range(len(prepared.sequence) - seed_length + 1):
            index[prepared.sequence[offset : offset + seed_length]].append((primer_index, offset))
    return {seed: tuple(locations) for seed, locations in index.items()}


def myers_global_scores(text: str, pattern: str, eq_table: dict[str, int] | None = None) -> list[int]:
    """Return Levenshtein distances from pattern to every text prefix."""
    if not pattern:
        raise ValueError("A primer pattern cannot be empty")

    eq_table = eq_table if eq_table is not None else build_eq_table(pattern)
    high_bit = 1 << (len(pattern) - 1)
    positive = ~0
    negative = 0
    score = len(pattern)
    scores = []

    for base in text:
        matched = eq_table.get(base, 0)
        vertical = matched | negative
        horizontal = (((matched & positive) + positive) ^ positive) | matched
        positive_horizontal = negative | ~(horizontal | positive)
        negative_horizontal = positive & horizontal

        if positive_horizontal & high_bit:
            score += 1
        elif negative_horizontal & high_bit:
            score -= 1

        positive_horizontal = (positive_horizontal << 1) | 1
        negative_horizontal <<= 1
        positive = negative_horizontal | ~(vertical | positive_horizontal)
        negative = positive_horizontal & vertical
        scores.append(score)
    return scores


def bitvector_match_rank(match: BitvectorMatch) -> tuple[float, int, int]:
    return (-(match.aligned_bases / len(match.primer.sequence)), match.cost, match.start)


class BitvectorLongestRemover:
    """Seeded Myers matching; equal-cost, equal-start candidates resolve to the longest primer coverage."""

    name = "bitvector_longest"

    def __init__(self, primer_fasta: Path | str, settings: BitvectorSettings = BitvectorSettings()) -> None:
        self.settings = settings
        left_primers, right_primers = split_primers(read_primers(primer_fasta))
        self.left_prepared = prepare_primers(left_primers + right_primers, "forward", False)
        self.right_prepared = prepare_primers(right_primers + left_primers, "reverse_complement", True)
        self.left_seed_index = build_seed_index(self.left_prepared, settings.seed_length)
        self.right_seed_index = build_seed_index(self.right_prepared, settings.seed_length)
        self.maximum_primer_length = max(len(prepared.sequence) for prepared in self.left_prepared)
        family_names = defaultdict(set)
        for prepared in self.left_prepared:
            family_names[self.primer_group_key(prepared.primer.name)].add(prepared.primer.name)
        self.primer_family_names = {group: frozenset(names) for group, names in family_names.items()}

    def minimum_aligned_bases(self, primer: Primer, settings: BitvectorSettings | None = None) -> int:
        settings = settings or self.settings
        return max(
            settings.min_aligned_primer_bases,
            math.ceil(len(primer.sequence) * settings.min_aligned_primer_fraction),
        )

    def effective_terminal_slack(self, read_length: int, settings: BitvectorSettings | None = None) -> int:
        settings = settings or self.settings
        if settings.max_terminal_slack <= settings.tandem_max_gap:
            return settings.max_terminal_slack
        return min(settings.max_terminal_slack, max(50, int(read_length * settings.max_terminal_slack_fraction)))

    def terminal_window_length(self, read_length: int) -> int:
        slack = self.effective_terminal_slack(read_length)
        error_allowance = math.ceil(self.maximum_primer_length * self.settings.max_error_rate)
        minimum_window = self.maximum_primer_length + error_allowance + slack
        return min(read_length, max(int(read_length * self.settings.window_fraction), minimum_window))

    @staticmethod
    def primer_group_key(primer_name: str) -> str:
        match = re.match(r"^(.*?_\d+)_(?:LEFT|RIGHT)(?:_alt\d+)?$", primer_name)
        return match.group(1) if match else primer_name

    def seed_supported_candidates(
        self,
        text: str,
        prepared_primers: tuple[PreparedPrimer, ...],
        seed_index: dict[str, tuple[tuple[int, int], ...]],
        settings: BitvectorSettings | None = None,
        allowed_primers: frozenset[str] | None = None,
    ) -> set[tuple[int, int, int]]:
        settings = settings or self.settings
        seed_length = settings.seed_length
        candidates = set()
        slack = settings.max_terminal_slack
        # The tandem search uses its own settings and keeps the stricter behaviour.
        primary_pass = settings == self.settings
        cognate_mode = allowed_primers is not None
        for text_offset in range(len(text) - seed_length + 1):
            seed = text[text_offset : text_offset + seed_length]
            for primer_index, primer_offset in seed_index.get(seed, ()):
                prepared = prepared_primers[primer_index]
                if cognate_mode and prepared.primer.name not in allowed_primers:
                    continue
                shift = text_offset - primer_offset
                minimum_overlap = self.minimum_aligned_bases(prepared.primer, settings)
                if cognate_mode:
                    min_truncated = min(minimum_overlap, settings.cognate_truncated_min_primer_bases)
                elif primary_pass and settings.exact_truncated_min_primer_bases > 0:
                    min_truncated = min(minimum_overlap, settings.exact_truncated_min_primer_bases)
                else:
                    min_truncated = minimum_overlap

                if -(len(prepared.sequence) - min_truncated) <= shift <= slack:
                    if shift < 0:
                        candidates.add((primer_index, -shift, 0))
                    else:
                        candidates.add((primer_index, 0, shift))
                # A slightly negative shift is usually an indel among the first primer bases, which leaves only
                # seeds on a shifted diagonal. Also try the whole primer from the read start so the alignment can
                # absorb the indel instead of treating the read as starting inside the primer.
                if (primary_pass or cognate_mode) and shift < 0 and -shift <= math.ceil(
                    len(prepared.sequence) * settings.max_error_rate
                ):
                    candidates.add((primer_index, 0, 0))
                # An indel in the 5' primer prefix before the seed shifts text_offset relative to primer_offset.
                # Propose nearby start positions so Myers does not include leading non-primer bases or clip the 5' end.
                if (
                    (primary_pass or cognate_mode)
                    and settings.pre_seed_indel_allowance > 0
                    and shift >= 0
                    and primer_offset > 0
                ):
                    for delta in range(-settings.pre_seed_indel_allowance, settings.pre_seed_indel_allowance + 1):
                        if delta != 0:
                            shifted_start = shift + delta
                            if 0 <= shifted_start <= min(slack, text_offset):
                                candidates.add((primer_index, 0, shifted_start))
                # For cognate rescue only: allow a few noisy terminal bases before an exact 3' primer suffix.
                if cognate_mode and shift < 0 and settings.cognate_terminal_noise_bases > 0:
                    for noise_bases in range(1, settings.cognate_terminal_noise_bases + 1):
                        suffix_start = noise_bases - shift
                        if (
                            len(prepared.sequence) - suffix_start >= settings.cognate_terminal_noise_min_primer_bases
                            and noise_bases <= text_offset
                        ):
                            candidates.add((primer_index, suffix_start, noise_bases))
        return candidates

    def candidate_matches(
        self,
        text: str,
        prepared: PreparedPrimer,
        suffix_start: int,
        text_start: int,
        settings: BitvectorSettings | None = None,
        cognate_rescue: bool = False,
    ) -> list[BitvectorMatch]:
        settings = settings or self.settings
        suffix = prepared.sequence[suffix_start:]
        overhang_cost = math.ceil(suffix_start * settings.overhang_penalty)
        offset_cost = math.ceil(text_start * settings.offset_penalty)

        if suffix_start > 0 and not cognate_rescue:
            std_min_aligned = self.minimum_aligned_bases(prepared.primer, settings)
            if len(suffix) < std_min_aligned:
                # Below the standard truncated minimum, accept only an error-free 3' primer suffix at the read start.
                if (
                    settings == self.settings
                    and text_start == 0
                    and len(suffix) >= settings.exact_truncated_min_primer_bases
                    and text.startswith(suffix)
                ):
                    return [
                        BitvectorMatch(
                            prepared.primer,
                            prepared.orientation,
                            0,
                            len(suffix) - 1,
                            overhang_cost,
                            len(suffix),
                        )
                    ]
                return []

        if cognate_rescue and suffix_start > 0 and text_start > 0:
            if (
                len(suffix) >= settings.cognate_terminal_noise_min_primer_bases
                and text[text_start:].startswith(suffix)
            ):
                return [
                    BitvectorMatch(
                        prepared.primer,
                        prepared.orientation,
                        text_start,
                        text_start + len(suffix) - 1,
                        overhang_cost + text_start,
                        len(suffix),
                    )
                ]
            return []

        eq_table = prepared.eq_table if suffix_start == 0 else build_eq_table(suffix)

        if cognate_rescue:
            eff_max_error_rate = settings.cognate_truncated_max_error_rate if suffix_start > 0 else settings.max_error_rate
            min_aligned = (
                settings.cognate_truncated_min_primer_bases
                if suffix_start > 0
                else self.minimum_aligned_bases(prepared.primer, settings)
            )
        elif suffix_start > 0:
            eff_max_error_rate = settings.truncated_max_error_rate
            min_aligned = self.minimum_aligned_bases(prepared.primer, settings)
        elif text_start > settings.boundary_slack_threshold:
            eff_max_error_rate = settings.embedded_max_error_rate
            min_aligned = max(
                settings.embedded_min_aligned_primer_bases,
                math.ceil(len(prepared.primer.sequence) * settings.embedded_min_aligned_primer_fraction),
            )
        else:
            eff_max_error_rate = settings.max_error_rate
            min_aligned = self.minimum_aligned_bases(prepared.primer, settings)

        maximum_errors = math.ceil(len(suffix) * eff_max_error_rate)
        maximum_length = min(len(text) - text_start, len(suffix) + maximum_errors)
        candidates = []

        for relative_end, edit_distance in enumerate(
            myers_global_scores(text[text_start : text_start + maximum_length], suffix, eq_table)
        ):
            aligned_bases = relative_end + 1
            observed_error_budget = math.ceil(aligned_bases * eff_max_error_rate)
            if aligned_bases >= min_aligned and edit_distance <= observed_error_budget:
                candidates.append(
                    BitvectorMatch(
                        prepared.primer,
                        prepared.orientation,
                        text_start,
                        text_start + relative_end,
                        edit_distance + overhang_cost + offset_cost,
                        aligned_bases,
                    )
                )
        return candidates

    def best_terminal_match(
        self,
        text: str,
        prepared_primers: tuple[PreparedPrimer, ...],
        seed_index: dict[str, tuple[tuple[int, int], ...]],
        settings: BitvectorSettings | None = None,
        excluded_primers: frozenset[str] = frozenset(),
    ) -> BitvectorMatch | None:
        settings = settings or self.settings
        candidates = [
            candidate
            for primer_index, suffix_start, text_start in self.seed_supported_candidates(
                text, prepared_primers, seed_index, settings
            )
            for candidate in self.candidate_matches(
                text, prepared_primers[primer_index], suffix_start, text_start, settings
            )
            if prepared_primers[primer_index].primer.name not in excluded_primers
        ]
        if not candidates:
            return None
        # Remaining ties share cost, start and end, so they imply the same cut whatever the primer identity.
        return min(candidates, key=lambda candidate: (candidate.cost, candidate.start, -candidate.end))

    def cognate_terminal_match(
        self,
        text: str,
        prepared_primers: tuple[PreparedPrimer, ...],
        seed_index: dict[str, tuple[tuple[int, int], ...]],
        anchor: BitvectorMatch,
    ) -> BitvectorMatch | None:
        group = self.primer_group_key(anchor.primer.name)
        family = self.primer_family_names.get(group, frozenset()) - {anchor.primer.name}
        if not family:
            return None
        candidates = [
            candidate
            for primer_index, suffix_start, text_start in self.seed_supported_candidates(
                text, prepared_primers, seed_index, self.settings, allowed_primers=family
            )
            for candidate in self.candidate_matches(
                text, prepared_primers[primer_index], suffix_start, text_start, self.settings, cognate_rescue=True
            )
        ]
        if not candidates:
            return None
        return min(candidates, key=lambda candidate: (candidate.cost, candidate.start, -candidate.end))

    def tandem_terminal_match(
        self,
        text: str,
        prepared_primers: tuple[PreparedPrimer, ...],
        seed_index: dict[str, tuple[tuple[int, int], ...]],
        primary: BitvectorMatch,
        overlap_bases: int = 0,
    ) -> BitvectorMatch | None:
        settings = self.settings
        error_allowance = math.ceil(self.maximum_primer_length * settings.tandem_max_error_rate)
        window_length = self.maximum_primer_length + settings.tandem_max_gap + overlap_bases + error_allowance
        excluded_primers = frozenset({primary.primer.name})
        tandem_settings = replace(
            settings,
            max_error_rate=settings.tandem_max_error_rate,
            min_aligned_primer_fraction=settings.tandem_min_aligned_primer_fraction,
            min_aligned_primer_bases=settings.tandem_min_aligned_primer_bases,
            max_terminal_slack=settings.tandem_max_gap + overlap_bases,
            exact_truncated_min_primer_bases=0,
            pre_seed_indel_allowance=0,
        )
        match = self.best_terminal_match(
            text[:window_length],
            prepared_primers,
            seed_index,
            settings=tandem_settings,
            excluded_primers=excluded_primers,
        )
        if match is not None and match.end < overlap_bases:
            return None
        return match

    def trim_record(self, record: FastqRecord) -> TrimResult:
        sequence = record.sequence
        if not sequence:
            return _result(record, 0, 0, None, None)

        window = self.terminal_window_length(len(sequence))
        left_window = sequence[:window]
        right_window = sequence[-window:][::-1]
        left_match = self.best_terminal_match(left_window, self.left_prepared, self.left_seed_index)
        right_match = self.best_terminal_match(right_window, self.right_prepared, self.right_seed_index)

        if self.settings.enable_cognate_rescue:
            if left_match is None and right_match is not None:
                left_match = self.cognate_terminal_match(
                    left_window, self.left_prepared, self.left_seed_index, right_match
                )
            elif right_match is None and left_match is not None:
                right_match = self.cognate_terminal_match(
                    right_window, self.right_prepared, self.right_seed_index, left_match
                )

        left_cut = left_match.end + 1 if left_match else 0
        right_cut = len(sequence) - (right_match.end + 1) if right_match else len(sequence)

        if left_cut > right_cut:
            if bitvector_match_rank(left_match) <= bitvector_match_rank(right_match):
                right_match = None
                right_cut = len(sequence)
            else:
                left_match = None
                left_cut = 0

        left_tandem = right_tandem = None
        if self.settings.enable_tandem_trim:
            if left_match:
                max_overlap = min(self.settings.tandem_max_overlap, left_match.aligned_bases - 1)
                slice_start = max(left_match.start + 1, left_cut - max_overlap)
                overlap_bases = left_cut - slice_start
                left_tandem = self.tandem_terminal_match(
                    sequence[slice_start:right_cut],
                    self.left_prepared,
                    self.left_seed_index,
                    left_match,
                    overlap_bases=overlap_bases,
                )
                if left_tandem:
                    new_left_cut = slice_start + left_tandem.end + 1
                    if new_left_cut <= right_cut - self.settings.tandem_min_retained_bases:
                        left_cut = new_left_cut
                    else:
                        left_tandem = None

            if right_match:
                right_end_offset = right_match.end + 1
                max_overlap = min(self.settings.tandem_max_overlap, right_match.aligned_bases - 1)
                slice_start = max(right_match.start + 1, right_end_offset - max_overlap)
                overlap_bases = right_end_offset - slice_start
                right_remaining = sequence[left_cut : len(sequence) - slice_start][::-1]
                right_tandem = self.tandem_terminal_match(
                    right_remaining,
                    self.right_prepared,
                    self.right_seed_index,
                    right_match,
                    overlap_bases=overlap_bases,
                )
                if right_tandem:
                    new_right_cut = len(sequence) - (slice_start + right_tandem.end + 1)
                    if new_right_cut >= left_cut + self.settings.tandem_min_retained_bases:
                        right_cut = new_right_cut
                    else:
                        right_tandem = None

        def call(match: BitvectorMatch | None, tandem: BitvectorMatch | None) -> EndCall | None:
            if match is None:
                return None
            detail = match.orientation
            if tandem:
                detail += (
                    f";tandem={tandem.primer.name}:{tandem.cost}:"
                    f"{tandem.aligned_bases}@{tandem.start}-{tandem.end}"
                )
            return EndCall(match.primer.name, match.aligned_bases, match.cost, detail)

        return _result(record, left_cut, right_cut, call(left_match, left_tandem), call(right_match, right_tandem))

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]:
        return [self.trim_record(record) for record in records]


# --------------------------------------------------------------------------------------
# Sassy (sassy)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SassySettings:
    max_error_rate: float = 0.30
    min_aligned_primer_fraction: float = 0.27
    min_aligned_primer_bases: int = 8
    overhang_penalty: float = 0.5
    window_fraction: float = 0.5


@dataclass(frozen=True)
class SassyMatch:
    primer: Primer
    start: int
    end: int
    cost: int
    cigar: str


def cigar_error_count(cigar: str) -> int:
    return sum(int(length) for length, operation in re.findall(r"(\d+)([XID])", cigar) if operation in "XID")


def is_boundary_overhang(match: Any, primer_length: int, edge_length: int, terminal: str) -> bool:
    if terminal == "left":
        return match.text_start == 0 and match.pattern_start > 0 and match.pattern_end == primer_length
    return match.text_end == edge_length and match.pattern_start == 0 and match.pattern_end < primer_length


def sassy_match_rank(match: SassyMatch) -> tuple[float, int, int]:
    observed_bases = match.end - match.start
    return (-(observed_bases / len(match.primer.sequence)), cigar_error_count(match.cigar), match.cost)


class SassyRemover:
    name = "sassy"

    def __init__(self, primer_fasta: Path | str, settings: SassySettings = SassySettings()) -> None:
        import sassy

        self.settings = settings
        self.searcher = sassy.Searcher("iupac", rc=False, alpha=settings.overhang_penalty)

        left_primers, original_right_primers = split_primers(read_primers(primer_fasta))
        right_primers = tuple(Primer(p.name, reverse_complement(p.sequence)) for p in original_right_primers)
        reverse_left_primers = tuple(Primer(p.name, reverse_complement(p.sequence)) for p in left_primers)
        reverse_right_primers = tuple(Primer(p.name, reverse_complement(p.sequence)) for p in right_primers)
        self.left_groups = self.prepare_primer_groups(left_primers + reverse_right_primers)
        self.right_groups = self.prepare_primer_groups(right_primers + reverse_left_primers)

    @staticmethod
    def prepare_primer_groups(
        primers: tuple[Primer, ...],
    ) -> tuple[tuple[tuple[Primer, ...], tuple[bytes, ...]], ...]:
        groups: dict[int, list[Primer]] = {}
        for primer in primers:
            groups.setdefault(len(primer.sequence), []).append(primer)
        return tuple(
            (tuple(group), tuple(primer.sequence.encode("ascii") for primer in group)) for group in groups.values()
        )

    def minimum_aligned_bases(self, primer: Primer) -> int:
        return max(
            self.settings.min_aligned_primer_bases,
            math.ceil(len(primer.sequence) * self.settings.min_aligned_primer_fraction),
        )

    def maximum_match_cost(self, primer: Primer, aligned_bases: int) -> int:
        observed_error_budget = math.ceil(aligned_bases * self.settings.max_error_rate)
        overhang_bases = max(0, len(primer.sequence) - aligned_bases)
        return observed_error_budget + math.ceil(overhang_bases * self.settings.overhang_penalty)

    def full_match_budget(self, primers: tuple[Primer, ...]) -> int:
        return math.ceil(len(primers[0].sequence) * self.settings.max_error_rate)

    def overhang_match_budget(self, primers: tuple[Primer, ...]) -> int:
        primer = primers[0]
        return self.maximum_match_cost(primer, self.minimum_aligned_bases(primer))

    def is_acceptable(self, primer: Primer, match: Any) -> bool:
        aligned_bases = match.text_end - match.text_start
        # The CIGAR regex is the expensive check, so it runs last.
        return (
            aligned_bases >= self.minimum_aligned_bases(primer)
            and match.cost <= self.maximum_match_cost(primer, aligned_bases)
            and cigar_error_count(match.cigar) <= math.ceil(aligned_bases * self.settings.max_error_rate)
        )

    def terminal_candidates(
        self,
        windows: Sequence[str],
        primer_groups: tuple[tuple[tuple[Primer, ...], tuple[bytes, ...]], ...],
        terminal: str,
    ) -> list[list[SassyMatch]]:
        encoded_windows = [window.encode("ascii") for window in windows]
        candidates: list[list[SassyMatch]] = [[] for _ in windows]

        for primers, patterns in primer_groups:
            primer_length = len(primers[0].sequence)

            for match in self.searcher.search_many(
                patterns, encoded_windows, k=self.full_match_budget(primers), threads=1, mode="batch_patterns"
            ):
                primer = primers[match.pattern_idx]
                if self.is_acceptable(primer, match):
                    candidates[match.text_idx].append(
                        SassyMatch(primer, match.text_start, match.text_end, match.cost, match.cigar)
                    )

            # The larger overhang budget is only affordable on the read boundary, not the whole window.
            edge_length = primer_length + math.ceil(primer_length * self.settings.max_error_rate)
            if terminal == "left":
                edges = [window[:edge_length] for window in encoded_windows]
            else:
                edges = [window[-edge_length:] for window in encoded_windows]

            for match in self.searcher.search_many(
                patterns, edges, k=self.overhang_match_budget(primers), threads=1, mode="batch_patterns"
            ):
                edge = edges[match.text_idx]
                if not is_boundary_overhang(match, primer_length, len(edge), terminal):
                    continue
                primer = primers[match.pattern_idx]
                if self.is_acceptable(primer, match):
                    offset = 0 if terminal == "left" else len(encoded_windows[match.text_idx]) - len(edge)
                    candidates[match.text_idx].append(
                        SassyMatch(primer, match.text_start + offset, match.text_end + offset, match.cost, match.cigar)
                    )

        return candidates

    @staticmethod
    def select_terminal_match(candidates: list[SassyMatch], window_length: int, terminal: str) -> SassyMatch | None:
        if terminal == "left":
            placement = lambda match: match.start
        elif terminal == "right":
            placement = lambda match: window_length - match.end
        else:
            raise ValueError(f"Unknown terminal: {terminal}")

        candidates = [match for match in candidates if placement(match) <= window_length]
        if not candidates:
            return None
        return min(candidates, key=lambda match: (match.cost, placement(match)))

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]:
        window = lambda record: int(len(record.sequence) * self.settings.window_fraction)
        left_windows = [record.sequence[: window(record)] for record in records]
        right_window_starts = [max(0, len(record.sequence) - window(record)) for record in records]
        right_windows = [record.sequence[start:] for record, start in zip(records, right_window_starts)]

        left_candidates = self.terminal_candidates(left_windows, self.left_groups, "left")
        right_candidates = self.terminal_candidates(right_windows, self.right_groups, "right")

        def call(match: SassyMatch | None) -> EndCall | None:
            if match is None:
                return None
            return EndCall(match.primer.name, match.end - match.start, match.cost, match.cigar)

        results = []
        for record, left_window, right_window_start, right_window, left_found, right_found in zip(
            records, left_windows, right_window_starts, right_windows, left_candidates, right_candidates
        ):
            left_match = self.select_terminal_match(left_found, len(left_window), "left")
            right_match = self.select_terminal_match(right_found, len(right_window), "right")

            left_cut = left_match.end if left_match else 0
            right_cut = right_window_start + right_match.start if right_match else len(record.sequence)

            if left_cut > right_cut:
                if sassy_match_rank(left_match) <= sassy_match_rank(right_match):
                    right_match = None
                    right_cut = len(record.sequence)
                else:
                    left_match = None
                    left_cut = 0

            results.append(_result(record, left_cut, right_cut, call(left_match), call(right_match)))
        return results


# --------------------------------------------------------------------------------------
# Parasail semi-global (semiglobal)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SemiglobalSettings:
    match_score: int = 3
    mismatch_score: int = -2
    gap_open: int = 3
    gap_extend: int = 1
    min_overlap: int = 8
    min_identity: float = 0.70
    min_score_fraction: float = 0.65
    window_fraction: float = 0.5


@dataclass(frozen=True)
class TerminalAlignment:
    primer_name: str
    primer_sequence: str
    score: int
    overlap_bases: int
    matches: int
    identity: float
    cut: int
    cigar: str


class SemiglobalRemover:
    name = "semiglobal"

    def __init__(self, primer_fasta: Path | str, settings: SemiglobalSettings = SemiglobalSettings()) -> None:
        import parasail

        self.settings = settings
        self.parasail = parasail
        self.primers = [(primer.name, primer.sequence) for primer in read_primers(primer_fasta)]
        self.score_matrix = parasail.matrix_create("ACGTN", settings.match_score, settings.mismatch_score)

    def score_fraction(self, hit: TerminalAlignment) -> float:
        maximum = self.settings.match_score * hit.overlap_bases
        return hit.score / maximum if maximum else 0.0

    def align_terminal_prefix(self, window: str, primer_name: str, primer: str) -> TerminalAlignment:
        """Align the primer to the left edge of a window with free primer-prefix overhang."""
        result = self.parasail.sg_trace_scan_sat(
            primer, window, self.settings.gap_open, self.settings.gap_extend, self.score_matrix
        )

        cigar = result.cigar.decode
        if isinstance(cigar, bytes):
            cigar = cigar.decode("ascii")

        operations = [(int(length), operation) for length, operation in re.findall(r"(\d+)([=XIDM])", cigar)]

        # A leading insertion is the free primer prefix; a trailing deletion is the unused window suffix.
        if operations and operations[0][1] == "I":
            operations = operations[1:]
        if operations and operations[-1][1] == "D":
            operations = operations[:-1]

        matches = sum(length for length, operation in operations if operation in {"=", "M"})
        query_bases = sum(length for length, operation in operations if operation in {"=", "X", "M", "I"})
        reference_bases = sum(length for length, operation in operations if operation in {"=", "X", "M", "D"})

        alignment_columns = sum(length for length, _ in operations)
        overlap_bases = min(query_bases, reference_bases)

        return TerminalAlignment(
            primer_name=primer_name,
            primer_sequence=primer,
            score=result.score,
            overlap_bases=overlap_bases,
            matches=matches,
            identity=matches / alignment_columns if alignment_columns else 0.0,
            cut=result.end_ref + 1,
            cigar=cigar,
        )

    def best_terminal_alignment(self, window: str) -> TerminalAlignment | None:
        if not window:
            return None

        accepted = []
        for primer_name, primer_sequence in self.primers:
            hit = self.align_terminal_prefix(window, primer_name, primer_sequence)
            if (
                hit.overlap_bases >= self.settings.min_overlap
                and hit.identity >= self.settings.min_identity
                and self.score_fraction(hit) >= self.settings.min_score_fraction
                and hit.cut <= len(window)
            ):
                accepted.append(hit)

        return max(accepted, key=lambda hit: (hit.score, hit.overlap_bases, hit.identity), default=None)

    def trim_record(self, record: FastqRecord) -> TrimResult:
        window = int(len(record.sequence) * self.settings.window_fraction)
        left_window = record.sequence[:window]
        # Window 0 makes [-0:] the whole read; kept for output parity.
        right_window_as_prefix = reverse_complement(record.sequence[-window:])

        left_hit = self.best_terminal_alignment(left_window)
        right_hit = self.best_terminal_alignment(right_window_as_prefix)

        left_removed = left_hit.cut if left_hit else 0
        right_removed = right_hit.cut if right_hit else 0
        right_boundary = len(record.sequence) - right_removed

        if left_removed > right_boundary:
            left_hit = None
            right_hit = None
            left_removed = 0
            right_boundary = len(record.sequence)

        def call(hit: TerminalAlignment | None) -> EndCall | None:
            if hit is None:
                return None
            detail = f"score={hit.score};identity={hit.identity:.3f};cigar={hit.cigar}"
            return EndCall(hit.primer_name, hit.overlap_bases, None, detail)

        return _result(record, left_removed, right_boundary, call(left_hit), call(right_hit))

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]:
        return [self.trim_record(record) for record in records]


# --------------------------------------------------------------------------------------
# Native ligation adapter and barcode removal (adapters)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class AdapterSettings:
    max_error_rate: float = 0.10
    terminal_window: int = 160
    terminal_slack: int = 50
    left_partial_adapter_length: int = 12
    left_partial_terminal_slack: int = 5
    left_partial_max_edits: int = 5
    right_partial_terminal_slack: int = 5
    right_partial_max_edits: int = 5
    overhang_penalty: float = 0.5


@dataclass(frozen=True)
class AdapterMatch:
    barcode: str
    start: int
    end: int
    cost: int
    cigar: str


class AdapterRemover:
    name = "adapters"

    def __init__(
        self,
        barcodes_tsv: Path | str = REPOSITORY_ROOT / "adapters.tsv",
        ligation_adapters_json: Path | str = REPOSITORY_ROOT / "ligation_adapters.json",
        barcode_flanks_json: Path | str = REPOSITORY_ROOT / "barcode_flanks.json",
        settings: AdapterSettings = AdapterSettings(),
    ) -> None:
        import sassy

        self.settings = settings
        self.searcher = sassy.Searcher("iupac", rc=False, alpha=settings.overhang_penalty)

        ligation_adapters = json.loads(Path(ligation_adapters_json).read_text())
        ligation_adapter_top = ligation_adapters["ligation_adapter"]["top_strand"]
        ligation_adapter_bottom = ligation_adapters["ligation_adapter"]["bottom_strand"]

        barcode_flanks = json.loads(Path(barcode_flanks_json).read_text())
        forward_flank_left = barcode_flanks["forward_flanks"]["left"]
        forward_flank_right = barcode_flanks["forward_flanks"]["right"]
        reverse_flank_left = barcode_flanks["reverse_flanks"]["left"]
        reverse_flank_right = barcode_flanks["reverse_flanks"]["right"]

        with Path(barcodes_tsv).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise ValueError(f"No barcode records found in {barcodes_tsv}")

        barcode_names = [row["barcode"] for row in rows]
        forward_barcodes = [row["Forward sequence"].upper() for row in rows]
        reverse_barcodes = [row["Reverse sequence"].upper() for row in rows]
        if len(set(barcode_names)) != len(barcode_names):
            raise ValueError("Barcode names must be unique")
        if any(len(barcode) != 24 for barcode in forward_barcodes + reverse_barcodes):
            raise ValueError("Every native barcode sequence must be 24 bases")

        left_flank_length = len(forward_flank_left) + len(forward_flank_right)
        left_partial_barcode_length = 24 + len(forward_flank_right)
        right_partial_barcode_length = 24 + len(reverse_flank_left)

        self.left_constructs = {
            name: ligation_adapter_top + forward_flank_left + forward + forward_flank_right
            for name, forward in zip(barcode_names, forward_barcodes)
        }
        self.left_partial_adapter_constructs = {
            name: construct[-(settings.left_partial_adapter_length + left_flank_length + 24) :]
            for name, construct in self.left_constructs.items()
        }
        self.left_partial_barcode_constructs = {
            name: (forward + forward_flank_right)[:left_partial_barcode_length]
            for name, forward in zip(barcode_names, forward_barcodes)
        }
        self.right_barcode_constructs = {
            name: reverse_flank_left + reverse + reverse_flank_right
            for name, reverse in zip(barcode_names, reverse_barcodes)
        }
        self.right_partial_barcode_constructs = {
            name: construct[:right_partial_barcode_length] for name, construct in self.right_barcode_constructs.items()
        }
        self.right_adapter_constructs = {
            name: self.right_barcode_constructs[name] + ligation_adapter_bottom for name in barcode_names
        }

    def max_edits(self, constructs: dict[str, str]) -> int:
        lengths = {len(construct) for construct in constructs.values()}
        if len(lengths) != 1:
            raise ValueError("Sassy batch patterns must have a common length")
        return math.ceil(lengths.pop() * self.settings.max_error_rate)

    def best_terminal_match(
        self,
        window: str,
        constructs: dict[str, str],
        terminal: str,
        terminal_slack: int | None = None,
        max_allowed_edits: int | None = None,
        require_unique_barcode: bool = False,
    ) -> AdapterMatch | None:
        terminal_slack = self.settings.terminal_slack if terminal_slack is None else terminal_slack
        barcode_names = tuple(constructs)
        patterns = [constructs[barcode].encode("ascii") for barcode in barcode_names]
        cost = max_allowed_edits if max_allowed_edits is not None else self.max_edits(constructs)
        matches = self.searcher.search_many(
            patterns, [window.encode("ascii")], k=cost, threads=1, mode="batch_patterns"
        )

        if terminal == "left":
            matches = [match for match in matches if match.text_start <= terminal_slack]
            placement = lambda match: match.text_start
        elif terminal == "right":
            matches = [match for match in matches if len(window) - match.text_end <= terminal_slack]
            placement = lambda match: len(window) - match.text_end
        else:
            raise ValueError(f"Unknown terminal: {terminal}")

        if not matches:
            return None

        if require_unique_barcode and len({barcode_names[match.pattern_idx] for match in matches}) > 1:
            return None

        match = min(matches, key=lambda candidate: (candidate.cost, placement(candidate)))
        return AdapterMatch(
            barcode=barcode_names[match.pattern_idx],
            start=match.text_start,
            end=match.text_end,
            cost=match.cost,
            cigar=match.cigar,
        )

    def best_left_terminal_match(self, window: str) -> AdapterMatch | None:
        settings = self.settings
        return (
            self.best_terminal_match(window, self.left_constructs, "left")
            or self.best_terminal_match(window, self.left_partial_adapter_constructs, "left")
            or self.best_terminal_match(
                window,
                self.left_partial_barcode_constructs,
                "left",
                settings.left_partial_terminal_slack,
                settings.left_partial_max_edits,
                True,
            )
        )

    def best_right_terminal_match(self, window: str) -> AdapterMatch | None:
        settings = self.settings
        return (
            self.best_terminal_match(window, self.right_adapter_constructs, "right")
            or self.best_terminal_match(window, self.right_barcode_constructs, "right")
            or self.best_terminal_match(
                window,
                self.right_partial_barcode_constructs,
                "right",
                settings.right_partial_terminal_slack,
                settings.right_partial_max_edits,
                True,
            )
        )

    def trim_record(self, record: FastqRecord) -> TrimResult:
        window_size = self.settings.terminal_window
        left_window = record.sequence[:window_size]
        right_window_start = max(0, len(record.sequence) - window_size)
        right_window = record.sequence[right_window_start:]

        left_match = self.best_left_terminal_match(left_window)
        right_match = self.best_right_terminal_match(right_window)

        left_cut = left_match.end if left_match else 0
        right_cut = right_window_start + right_match.start if right_match else len(record.sequence)
        if left_cut > right_cut:
            left_match = None
            right_match = None
            left_cut = 0
            right_cut = len(record.sequence)

        def call(match: AdapterMatch | None) -> EndCall | None:
            if match is None:
                return None
            return EndCall(match.barcode, match.end - match.start, match.cost, match.cigar)

        return _result(record, left_cut, right_cut, call(left_match), call(right_match))

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]:
        return [self.trim_record(record) for record in records]


# --------------------------------------------------------------------------------------
# FASTQ runner and command line
# --------------------------------------------------------------------------------------

METHODS = ("adapters", "bitvector_longest", "sassy", "semiglobal")
PRIMER_METHODS = METHODS[1:]
CALL_COLUMNS = ("read_id", "original_length", "retained_start", "retained_end")
END_COLUMNS = ("name", "cost", "aligned_bases", "detail")

_ACTIVE_TRIMMER: Trimmer | None = None

# Debugging aid: methods listed here are run twice in a row (the second pass trims the output of the first).
# Add or remove method names, e.g. TWO_PASS_METHODS = {"bitvector_longest"}; an empty set keeps single-pass behaviour.
TWO_PASS_METHODS: set[str] = set()


class TwoPassTrimmer:
    """Runs a trimmer a second time on its own output; cuts and calls from both passes are combined."""

    def __init__(self, trimmer: Trimmer) -> None:
        self.trimmer = trimmer
        self.name = trimmer.name

    @staticmethod
    def _merge_call(first: EndCall | None, second: EndCall | None) -> EndCall | None:
        if second is None:
            return first
        if first is None:
            return replace(second, detail=f"{second.detail};pass2")
        return replace(first, detail=f"{first.detail};pass2={second.name}:{second.aligned_bases}")

    def trim_records(self, records: Sequence[FastqRecord]) -> list[TrimResult]:
        first_pass = self.trimmer.trim_records(records)
        second_pass = self.trimmer.trim_records([result.record for result in first_pass])
        return [
            TrimResult(
                second.record,
                first.original_length,
                first.left_removed + second.left_removed,
                first.right_removed + second.right_removed,
                self._merge_call(first.left, second.left),
                self._merge_call(first.right, second.right),
            )
            for first, second in zip(first_pass, second_pass)
        ]


def create_trimmer(method: str, primer_fasta: Path | str | None = None) -> Trimmer:
    trimmer = _create_single_trimmer(method, primer_fasta)
    return TwoPassTrimmer(trimmer) if method in TWO_PASS_METHODS else trimmer


def _create_single_trimmer(method: str, primer_fasta: Path | str | None = None) -> Trimmer:
    if method == "adapters":
        return AdapterRemover()
    if primer_fasta is None:
        raise ValueError(f"Method {method!r} needs a primer FASTA")
    if method == "bitvector_longest":
        return BitvectorLongestRemover(primer_fasta)
    if method == "sassy":
        return SassyRemover(primer_fasta)
    if method == "semiglobal":
        return SemiglobalRemover(primer_fasta)
    raise ValueError(f"Unknown method {method!r}; choose from {METHODS}")


def _batched(items: Iterable[FastqRecord], size: int) -> Iterator[list[FastqRecord]]:
    iterator = iter(items)
    while batch := list(islice(iterator, size)):
        yield batch


def _trim_batch(batch: list[FastqRecord]) -> list[TrimResult]:
    return _ACTIVE_TRIMMER.trim_records(batch)


def _end_fields(call: EndCall | None) -> list[object]:
    if call is None:
        return ["", "", "", ""]
    return [call.name, "" if call.cost is None else call.cost, call.aligned_bases, call.detail]


def trim_fastq(
    trimmer: Trimmer,
    input_path: Path | str,
    output_path: Path | str,
    *,
    calls_path: Path | str | None = None,
    workers: int = 1,
    batch_size: int = 2000,
) -> dict[str, float | int]:
    """Trim a FASTQ in input order and return counts and the wall-clock time of reading, trimming and writing."""
    global _ACTIVE_TRIMMER
    _ACTIVE_TRIMMER = trimmer
    totals = {"total_reads": 0, "bases_removed": 0, "left_calls": 0, "right_calls": 0}

    started = time.perf_counter()
    with ExitStack() as stack:
        output = stack.enter_context(open(output_path, "w"))
        writer = None
        if calls_path is not None:
            calls_handle = stack.enter_context(open(calls_path, "w", newline=""))
            writer = csv.writer(calls_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [*CALL_COLUMNS, *(f"{side}_{column}" for side in ("left", "right") for column in END_COLUMNS)]
            )

        batches = _batched(read_fastq(input_path), batch_size)
        if workers > 1:
            pool = stack.enter_context(mp.get_context("fork").Pool(workers))
            trimmed_batches = pool.imap(_trim_batch, batches)
        else:
            trimmed_batches = map(_trim_batch, batches)

        for trimmed_batch in trimmed_batches:
            for result in trimmed_batch:
                record = result.record
                if len(record.sequence) != len(record.qualities) or len(record.sequence) > result.original_length:
                    raise ValueError(f"Inconsistent trim for read {record.title.split()[:1]}")
                output.write(f"@{record.title}\n{record.sequence}\n+\n{record.qualities}\n")

                totals["total_reads"] += 1
                totals["bases_removed"] += result.original_length - len(record.sequence)
                totals["left_calls"] += result.left is not None
                totals["right_calls"] += result.right is not None
                if writer is not None:
                    writer.writerow(
                        [
                            (record.title.split() or [""])[0],
                            result.original_length,
                            result.left_removed,
                            result.original_length - result.right_removed,
                            *_end_fields(result.left),
                            *_end_fields(result.right),
                        ]
                    )

    return {**totals, "elapsed": time.perf_counter() - started}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Trim terminal adapters or primers from a FASTQ file.")
    parser.add_argument("method", choices=METHODS)
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("-p", "--primers", type=Path, help="primer FASTA, required for primer methods")
    parser.add_argument("--calls", type=Path, help="write the per-read calls as TSV")
    parser.add_argument("--stats", type=Path, help="write total_reads,bases_removed,elapsed as CSV")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=2000)
    args = parser.parse_args(argv)

    if args.method in PRIMER_METHODS and args.primers is None:
        parser.error(f"{args.method} needs --primers")

    trimmer = create_trimmer(args.method, args.primers)
    stats = trim_fastq(
        trimmer,
        args.input,
        args.output,
        calls_path=args.calls,
        workers=args.workers,
        batch_size=args.batch_size,
    )
    if args.stats is not None:
        args.stats.write_text(
            f"total_reads,bases_removed,elapsed\n{stats['total_reads']},{stats['bases_removed']},{stats['elapsed']}\n"
        )
    print(
        f"{args.method}: {stats['total_reads']:,} reads, left={stats['left_calls']:,}, "
        f"right={stats['right_calls']:,}, bases removed={stats['bases_removed']:,}, "
        f"{stats['elapsed']:,.1f} s with {args.workers} worker(s)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
