"""Counting helpers that attribute consensus differences to primer removal against the no-trim control."""

from __future__ import annotations

from collections import Counter

import numpy as np


def transition_matrix(reference_codes, method_codes, n_codes: int) -> np.ndarray:
    """Count positions by (reference class, method class); positions with a negative code are skipped."""
    reference_codes = np.asarray(reference_codes, dtype=np.int64)
    method_codes = np.asarray(method_codes, dtype=np.int64)
    keep = (reference_codes >= 0) & (method_codes >= 0)
    flat = reference_codes[keep] * n_codes + method_codes[keep]
    return np.bincount(flat, minlength=n_codes * n_codes).reshape(n_codes, n_codes)


def attribute_transitions(matrix, groups, correct: str = "correct", bad=("wrong", "missing")) -> dict:
    """Split a reference-by-method transition matrix into introduced, fixed, shared and mixed counts."""
    matrix = np.asarray(matrix)
    introduced, fixed, shared, mixed = Counter(), Counter(), Counter(), 0
    for i, j in zip(*np.nonzero(matrix)):
        count = int(matrix[i, j])
        group_reference, group_method = groups[i], groups[j]
        if i == j:
            if group_reference in bad:
                shared[group_reference] += count
        elif group_reference == correct and group_method != correct:
            introduced[group_method] += count
        elif group_method == correct and group_reference != correct:
            fixed[group_reference] += count
        elif group_reference != correct and group_method != correct:
            mixed += count
    return {"introduced": dict(introduced), "fixed": dict(fixed), "shared": dict(shared), "mixed": mixed}


def match_events(reference_positions, reference_classes, method_positions, method_classes, tolerance: int) -> dict:
    """Count method events shared with the reference, introduced by the method, and reference events it fixed."""

    def matched(positions, classes, other_positions, other_classes):
        other_positions = np.asarray(other_positions, dtype=np.int64)
        other_classes = np.asarray(other_classes, dtype=object)
        hits = np.zeros(len(positions), dtype=bool)
        for k, (position, name) in enumerate(zip(positions, classes)):
            hits[k] = bool(np.any((other_classes == name) & (np.abs(other_positions - position) <= tolerance)))
        return hits

    method_hit = matched(method_positions, method_classes, reference_positions, reference_classes)
    reference_hit = matched(reference_positions, reference_classes, method_positions, method_classes)
    return {"shared": int(method_hit.sum()), "introduced": int((~method_hit).sum()), "fixed": int((~reference_hit).sum())}


def bootstrap_mean(values, resamples: int = 5000, seed: int = 0, level: float = 0.95) -> dict:
    """Mean of per-sample values with a percentile bootstrap interval over samples."""
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, len(values), size=(resamples, len(values)))].mean(axis=1)
    tail = (1 - level) / 2
    low, high = np.quantile(draws, [tail, 1 - tail])
    return {"mean": float(values.mean()), "low": float(low), "high": float(high)}


def break_even_price(wrong_a: float, missing_a: float, wrong_b: float, missing_b: float) -> float:
    """Missing positions worth one wrong base at which methods a and b score equally; nan if the ranking never flips."""
    wrong_gap = wrong_b - wrong_a
    missing_gap = missing_a - missing_b
    if wrong_gap == 0 or missing_gap == 0 or np.sign(wrong_gap) != np.sign(missing_gap):
        return float("nan")
    return float(missing_gap / wrong_gap)
