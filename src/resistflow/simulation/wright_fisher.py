"""
Wright-Fisher simulation engine.

Implements the core generation loop:
  1. Selection  - apply_selection() shifts allele frequency based on genotype fitness
  2. Drift      - apply_drift() adds stochasticity via binomial sampling

The engine is intentionally decoupled from the selection and drift models:
both are imported as separate modules so either can be swapped out independently
(e.g. frequency-dependent selection, alternative demographic models) without
touching this file.

Main entry point: simulate()
"""

from .selection import apply_selection
from .drift import apply_drift


def simulate(
    p0,
    N,
    n_generations,
    w_AA,
    w_Aa,
    w_aa,
    n_replicates=1,
):
    """
    Run a Wright-Fisher simulation.

    Parameters
    ----------
    p0 : float
        Initial frequency of the resistance allele A [0, 1].
    N : int
        Effective population size (Ne).
    n_generations : int
        Number of generations to simulate.
    w_AA : float
        Fitness of homozygous resistant genotype (AA).
    w_Aa : float
        Fitness of heterozygous genotype (Aa).
    w_aa : float
        Fitness of homozygous susceptible genotype (aa).
    n_replicates : int
        Number of independent simulation runs (default: 1).
        Multiple replicates capture the stochasticity of drift.
    Returns
    -------
    list of list of float
        Allele frequency trajectories, one list per replicate.
        Each inner list has length n_generations + 1 (includes p0).
    """

    if not 0 <= p0 <= 1:
        raise ValueError(f"p0 must be in [0, 1], got {p0}.")
    if not isinstance(N, int) or N <= 0:
        raise ValueError(f"N must be a positive integer, got {N}.")
    if n_generations <= 0:
        raise ValueError(f"n_generations must be a positive integer, got {n_generations}.")
    if n_replicates <= 0:
        raise ValueError(f"n_replicates must be a positive integer, got {n_replicates}.")
    if w_AA < 0 or w_Aa < 0 or w_aa < 0:
        raise ValueError("Fitness values must be non-negative.")
    
    allele_frequencies = []

    for i in range(n_replicates):
        p = p0
        replicates = [p]
        for j in range(n_generations):
            p_prime = apply_selection(p, w_AA, w_Aa, w_aa)
            p_new = apply_drift(p_prime, N)
            replicates.append(p_new)
            p = p_new
        allele_frequencies.append(replicates)

    return allele_frequencies
    