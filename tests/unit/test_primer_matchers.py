"""Regression tests for terminal and tandem bit-vector primer matching."""

import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_matchers import (  # noqa: E402
    BitvectorLongestRemover,
    FastqRecord,
    reverse_complement,
)


LEFT_PRIMER = "AACCGGTTAACCGGTT"
RIGHT_PRIMER = "TTGCAACGTTGCAACG"
INSERT = "GATTACAAGCTGATCGTACCGATGCTAGCTACGATCGTACGATGCTAGCATCGATGCTAGCATCG"


def make_remover(directory: str, primers: list[tuple[str, str]]) -> BitvectorLongestRemover:
    primer_fasta = Path(directory) / "primers.fasta"
    primer_fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in primers))
    return BitvectorLongestRemover(primer_fasta)


def test_noisy_terminal_primer_beats_exact_internal_match():
    decoy = "GGTCAAGTCCGATGCA"
    primers = [("scheme_1_LEFT", LEFT_PRIMER), ("scheme_2_LEFT", decoy), ("scheme_1_RIGHT", RIGHT_PRIMER)]
    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        noisy_terminal = list(LEFT_PRIMER)
        noisy_terminal[2] = "A" if noisy_terminal[2] != "A" else "C"
        noisy_terminal[11] = "A" if noisy_terminal[11] != "A" else "C"
        # Decoy placed 60 nt into the insert: terminal primer at start 0 beats internal decoy
        sequence = "".join(noisy_terminal) + "G" * 60 + decoy + INSERT + reverse_complement(RIGHT_PRIMER)
        window = remover.terminal_window_length(len(sequence))

        match = remover.best_terminal_match(
            sequence[:window], remover.left_prepared, remover.left_seed_index
        )

        assert match is not None
        assert match.primer.name == "scheme_1_LEFT"
        assert match.start == 0


def test_embedded_primer_behind_adapter_noise_is_trimmed():
    primers = [("scheme_1_LEFT", LEFT_PRIMER), ("scheme_1_RIGHT", RIGHT_PRIMER)]
    adapter_noise = "A" * 42  # 42 nt of unremoved adapter/barcode noise
    sequence = adapter_noise + LEFT_PRIMER + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("embedded", sequence, "I" * len(sequence)))

    assert result.record.sequence == INSERT
    assert result.left_removed == len(adapter_noise) + len(LEFT_PRIMER)
    assert result.left is not None
    assert result.left.name == "scheme_1_LEFT"


def test_match_beyond_terminal_slack_is_rejected():
    primers = [("scheme_1_LEFT", LEFT_PRIMER), ("scheme_1_RIGHT", RIGHT_PRIMER)]
    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        sequence = "G" * (remover.settings.max_terminal_slack + 1) + LEFT_PRIMER + INSERT
        window = remover.terminal_window_length(len(sequence))

        match = remover.best_terminal_match(
            sequence[:window], remover.left_prepared, remover.left_seed_index
        )

        assert match is None


def test_primer_with_small_terminal_offset_is_trimmed():
    primers = [("scheme_1_LEFT", LEFT_PRIMER), ("scheme_1_RIGHT", RIGHT_PRIMER)]
    adapter_remnant = "GATTACCA"
    sequence = adapter_remnant + LEFT_PRIMER + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("shifted", sequence, "I" * len(sequence)))

    assert result.record.sequence == INSERT
    assert result.left_removed == len(adapter_remnant) + len(LEFT_PRIMER)
    assert result.left is not None
    assert result.left.cost is not None and result.left.cost >= 0


def test_left_tandem_partial_primer_is_removed_after_primary_match():
    second_primer = LEFT_PRIMER[-6:] + "ATGCAATGCCATGCAATG"
    primers = [
        ("scheme_1_LEFT", LEFT_PRIMER),
        ("scheme_2_LEFT", second_primer),
        ("scheme_1_RIGHT", RIGHT_PRIMER),
    ]
    sequence = LEFT_PRIMER + second_primer[6:] + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("left-tandem", sequence, "I" * len(sequence)))

    assert result.record.sequence == INSERT
    assert result.left_removed == len(LEFT_PRIMER) + len(second_primer) - 6
    assert result.left is not None
    assert result.left.name == "scheme_1_LEFT"
    assert "tandem=scheme_2_LEFT" in result.left.detail


def test_right_tandem_partial_primer_is_removed_after_primary_match():
    second_primer = RIGHT_PRIMER[-6:] + "ATGCAATGCCATGCAATG"
    primers = [
        ("scheme_1_LEFT", LEFT_PRIMER),
        ("scheme_1_RIGHT", RIGHT_PRIMER),
        ("scheme_2_RIGHT", second_primer),
    ]
    sequence = LEFT_PRIMER + INSERT + reverse_complement(second_primer[6:]) + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("right-tandem", sequence, "I" * len(sequence)))

    assert result.record.sequence == INSERT
    assert result.right_removed == len(RIGHT_PRIMER) + len(second_primer) - 6
    assert result.right is not None
    assert result.right.name == "scheme_1_RIGHT"
    assert "tandem=scheme_2_RIGHT" in result.right.detail


def test_same_amplicon_alternative_primer_can_be_treated_as_tandem():
    alternative = LEFT_PRIMER[-6:] + "ATGCAATGCCATGCAATG"
    primers = [
        ("scheme_1_LEFT", LEFT_PRIMER),
        ("scheme_1_LEFT_alt1", alternative),
        ("scheme_1_RIGHT", RIGHT_PRIMER),
    ]
    sequence = LEFT_PRIMER + alternative[6:] + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("same-amplicon-alt", sequence, "I" * len(sequence)))

    assert result.record.sequence == INSERT
    assert result.left is not None
    assert "tandem=scheme_1_LEFT_alt1" in result.left.detail


def test_tandem_scan_rejects_a_15nt_partial_match():
    second_primer = LEFT_PRIMER[-9:] + "ATGCAATGCCATGCA"
    primers = [
        ("scheme_1_LEFT", LEFT_PRIMER),
        ("scheme_2_LEFT", second_primer),
        ("scheme_1_RIGHT", RIGHT_PRIMER),
    ]
    sequence = LEFT_PRIMER + second_primer[9:] + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("short-tandem", sequence, "I" * len(sequence)))

    assert result.record.sequence.startswith(second_primer[9:] + INSERT)
    assert result.left is not None
    assert "tandem=" not in result.left.detail


def test_tandem_scan_does_not_reach_past_adjacent_gap():
    second_primer = LEFT_PRIMER[-6:] + "ATGCAATGCCATGCAATG"
    primers = [
        ("scheme_1_LEFT", LEFT_PRIMER),
        ("scheme_2_LEFT", second_primer),
        ("scheme_1_RIGHT", RIGHT_PRIMER),
    ]
    gap = 50  # Past tandem_max_gap = 40
    sequence = LEFT_PRIMER + "G" * gap + second_primer + INSERT + reverse_complement(RIGHT_PRIMER)

    with TemporaryDirectory() as directory:
        remover = make_remover(directory, primers)
        result = remover.trim_record(FastqRecord("non-adjacent", sequence, "I" * len(sequence)))

    assert result.record.sequence.startswith("G" * gap + second_primer)
    assert result.left is not None
    assert "tandem=" not in result.left.detail


def test_real_read_regressions():
    primer_fasta = Path("manual_test/ESIB_EQA_2026_SARS1.primers.fasta")
    if not primer_fasta.is_file():
        return
    remover = BitvectorLongestRemover(primer_fasta)

    qc1 = Path("manual_test/processed/benchmark2/pipeline/shared/EQA_01/02.qc.fastq")
    qc2 = Path("manual_test/processed/benchmark2/pipeline/shared/EQA_02/02.qc.fastq")
    if not qc1.is_file() or not qc2.is_file():
        return

    import primer_matchers as pm

    target_cases = [
        (qc2, "383ce6a0-e377-4eff-a500-3123214a012c", 352, "2_LEFT 15nt overlap"),
        (qc2, "d7b963e4-8632-490c-a413-e658fb74d3e9", 344, "2_RIGHT + 153nt noise 2_LEFT"),
        (qc2, "eb2d6597-45a0-43c8-9d1f-44006aa01d0e", 352, "7_LEFT 15nt overlap"),
        (qc2, "dd3e7f32-16d2-4a48-97d9-f505be2ee07b", 194, "10_RIGHT + 10_RIGHT_alt1"),
        (qc2, "066362ac-f465-40aa-bbad-62f029a55e06", 346, "16_LEFT 19nt overlap"),
        (qc2, "b876b5f2-0505-45f0-b07f-d754dcaf0840", 465, "78_RIGHT 13nt overlap"),
        (qc2, "c9da0731-d8ac-4503-ac28-2f965ad3c76d", 354, "leading noise + 20_LEFT"),
        (qc2, "8ac18e84-2ec4-4d8e-a698-b856aeb27785", 357, "48_LEFT 21S+9nt overlap"),
        (qc2, "1d2b15da-fdf8-48c3-a55e-b2c0b85c7c3a", 349, "17_RIGHT + 42nt noise 17_LEFT"),
        (qc1, "055eb27c-55da-468c-9d72-823b4c87fdff", 360, "38_RIGHT + 38_LEFT (no overtrim)"),
        (qc1, "322b28d6-9939-4b43-9090-1070c6c4ea9f", 365, "78_RIGHT + tandem 39_RIGHT + 78_LEFT"),
    ]

    for qc_path, qname, exp_len, desc in target_cases:
        r = next((rec for rec in pm.read_fastq(qc_path) if qname in rec.title), None)
        assert r is not None, f"Read {qname} not found in {qc_path}"
        res = remover.trim_record(r)
        assert len(res.record.sequence) == exp_len, (
            f"Read {qname} ({desc}): expected {exp_len} nt, got {len(res.record.sequence)} nt "
            f"(left_removed={res.left_removed}, right_removed={res.right_removed})"
        )


if __name__ == "__main__":
    test_noisy_terminal_primer_beats_exact_internal_match()
    test_embedded_primer_behind_adapter_noise_is_trimmed()
    test_match_beyond_terminal_slack_is_rejected()
    test_primer_with_small_terminal_offset_is_trimmed()
    test_left_tandem_partial_primer_is_removed_after_primary_match()
    test_right_tandem_partial_primer_is_removed_after_primary_match()
    test_same_amplicon_alternative_primer_can_be_treated_as_tandem()
    test_tandem_scan_rejects_a_15nt_partial_match()
    test_tandem_scan_does_not_reach_past_adjacent_gap()
    test_real_read_regressions()
    print("Bit-vector primer matcher tests passed.")