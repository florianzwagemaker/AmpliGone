"""Entry point for future within-read primer scanning and trimming."""

from dataclasses import dataclass

import pandas as pd

from AmpliGone.log import log


@dataclass(frozen=True)
class SequenceTrimResult:
    """The sequence and quality values returned by sequence-aware trimming."""

    sequence: str
    qualities: str


def trim_coordinate_unassociated_read(
    sequence: str, qualities: str, primer_df: pd.DataFrame | None = None
) -> SequenceTrimResult:
    """Run sequence-aware primer trimming for a coordinate-unassociated read.

    The actual primer search method and trimming logic will be initialized from this point.

    ``primer_df`` is passed through unchanged so sequence-aware trimming can use
    the full primer information while it is being prototyped.
    """
    return SequenceTrimResult(sequence=sequence, qualities=qualities)