"""
Tests for the Wright-Fisher simulation engine.

Tests are added progressively alongside implementation.
Each function we implement in simulation/ gets a corresponding test here.

Planned tests
-------------
test_selection_no_pressure
    With wAA = wAa = waa = 1.0, allele frequency after selection == p0.

test_selection_strong_pressure
    With waa = 0.0, resistance allele should increase each generation.

test_drift_fixation
    With very small N, all replicates should fix at 0 or 1 eventually.

test_drift_large_N
    With very large N, drift should be negligible (frequency stays near p0).

test_simulate_output_shape
    simulate() should return n_replicates trajectories of length n_generations + 1.

test_simulate_baseline_scenario
    Under baseline scenario, frequency should perform a random walk around p0.
"""

import pytest
import numpy as np
from resistflow.simulation.wright_fisher import simulate

def test_simulate_output_shape():
    trajectories = simulate(p0 = 0.5, 
                      N = 10_000, 
                      n_generations= 10,
                      w_AA = 1.0, w_Aa = 0.5, w_aa = 0.0, 
                      n_replicates = 3)
    
    assert len(trajectories) == 3
    for t in trajectories:
        assert len(t) == 11

def test_baseline_scenario():
    np.random.seed(42)
    trajectories = simulate(p0 = 0.5, 
                      N = 1_000_000, 
                      n_generations= 10,
                      w_AA = 1.0, w_Aa = 1.0, w_aa = 1.0, 
                      n_replicates = 3)
    final_freqs = [t[-1] for t in trajectories]
    assert sum(final_freqs) / len(final_freqs) == pytest.approx(0.5, abs=0.05)
