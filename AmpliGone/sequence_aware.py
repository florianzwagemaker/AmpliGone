"""Entry point for future within-read primer scanning and trimming."""

from dataclasses import dataclass
from AmpliGone.log import log


@dataclass(frozen=True)
class SequenceTrimResult:
    """The sequence and quality values returned by sequence-aware trimming."""

    sequence: str
    qualities: str


def trim_coordinate_unassociated_read(
    sequence: str, qualities: str
) -> SequenceTrimResult:
    """Run sequence-aware primer trimming for a coordinate-unassociated read.

    the actual primer search method and trimming logic will be initialized from this point.
    """
    return SequenceTrimResult(sequence=sequence, qualities=qualities)