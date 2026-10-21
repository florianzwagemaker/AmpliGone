"""Unit tests for the `consensus_attribution` helpers."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from consensus_attribution import (  # noqa: E402
    attribute_transitions,
    bootstrap_mean,
    break_even_price,
    match_events,
    transition_matrix,
)

GROUPS = ["correct", "wrong", "missing"]


def test_transition_matrix_skips_negative_codes():
    matrix = transition_matrix([0, 0, 1, -1], [0, 1, 1, 1], n_codes=2)
    assert matrix.tolist() == [[1, 1], [0, 1]]


def test_attribute_transitions_splits_introduced_fixed_and_shared():
    matrix = [
        [5, 2, 1],  # no-trim correct: method wrong 2 and missing 1 are introduced
        [3, 4, 0],  # no-trim wrong: method correct 3 is fixed, method wrong 4 is shared
        [2, 0, 6],  # no-trim missing: method correct 2 is fixed, method missing 6 is shared
    ]
    parts = attribute_transitions(matrix, GROUPS)
    assert parts["introduced"] == {"wrong": 2, "missing": 1}
    assert parts["fixed"] == {"wrong": 3, "missing": 2}
    assert parts["shared"] == {"wrong": 4, "missing": 6}
    assert parts["mixed"] == 0


def test_attribute_transitions_counts_cross_bad_changes_as_mixed():
    matrix = [
        [5, 0, 0],
        [0, 4, 7],  # wrong in no-trim, missing in method: neither is correct
        [0, 0, 6],
    ]
    parts = attribute_transitions(matrix, GROUPS)
    assert parts["mixed"] == 7
    assert parts["introduced"] == {}
    assert parts["fixed"] == {}


def test_match_events_within_tolerance_is_shared():
    parts = match_events(
        reference_positions=[100, 500],
        reference_classes=["false_deletion", "false_insertion"],
        method_positions=[110, 900],
        method_classes=["false_deletion", "false_deletion"],
        tolerance=25,
    )
    assert parts == {"shared": 1, "introduced": 1, "fixed": 1}


def test_match_events_requires_the_same_class():
    parts = match_events([100], ["false_deletion"], [100], ["false_insertion"], tolerance=25)
    assert parts == {"shared": 0, "introduced": 1, "fixed": 1}


def test_bootstrap_mean_of_constant_values_is_exact():
    interval = bootstrap_mean([2.0] * 8, resamples=200, seed=1)
    assert interval == {"mean": 2.0, "low": 2.0, "high": 2.0}


def test_bootstrap_mean_interval_brackets_the_mean():
    interval = bootstrap_mean([1.0, 3.0], resamples=1000, seed=0)
    assert interval["mean"] == 2.0
    assert 1.0 <= interval["low"] < interval["high"] <= 3.0


def test_bootstrap_mean_is_reproducible_with_a_seed():
    values = np.arange(8, dtype=float)
    assert bootstrap_mean(values, seed=7) == bootstrap_mean(values, seed=7)


def test_break_even_price_makes_both_methods_score_equally():
    price = break_even_price(wrong_a=0, missing_a=10, wrong_b=2, missing_b=0)
    assert price == 5.0
    assert 0 + 10 / price == 2 + 0 / price


def test_break_even_price_is_nan_when_ranking_never_flips():
    assert np.isnan(break_even_price(wrong_a=0, missing_a=0, wrong_b=2, missing_b=5))
    assert np.isnan(break_even_price(wrong_a=1, missing_a=3, wrong_b=1, missing_b=0))
