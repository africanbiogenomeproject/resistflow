"""
Tests for the genetic drift model (drift.py).

Planned tests
-------------
test_drift_returns_valid_frequency
    apply_drift() should always return a float in [0, 1].

test_absorbing_states
    If p = 0.0, result should remain 0.0 (allele already lost).
    If p = 1.0, result should remain 1.0 (allele already fixed).

test_large_N_low_variance
    With a very large effective population size (e.g. N=1_000_000),
    drift per generation should be negligible — the returned frequency
    should stay very close to p across many replicates.

test_small_N_high_variance
    With a very small effective population size (e.g. N=10),
    there should be substantial variance in outcome across replicates.

test_output_is_multiple_of_1_over_2N
    The returned frequency must be a multiple of 1/(2N), since it results
    from integer binomial sampling. This verifies the binomial mechanism.
"""

import pytest
import numpy as np
from resistflow.simulation.drift import apply_drift


@pytest.mark.parametrize("p, expected", [
    (1.0, 1.0),
    (0.0, 0.0),

])
def test_absorbing_states(p, expected):
    assert apply_drift(p, N=  10_000) == expected

def test_large_N_low_variance():
    np.random.seed(42)
    result = apply_drift(p=0.5, N=1_000_000) 
    assert  result == pytest.approx(0.5, abs=0.01)

