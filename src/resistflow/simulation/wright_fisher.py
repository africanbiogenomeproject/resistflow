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

from resistflow.simulation.selection import apply_selection
from resistflow.simulation.drift import apply_drift


def simulate(
    p0,
    N,
    n_generations,
    w_AA,
    w_Aa,
    w_aa,
    n_replicates=1,
    generations_per_year=10,
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
    generations_per_year : int
        Biological generations per year for An. gambiae (default: 10).
        Used to label time axis in years rather than generations.

    Returns
    -------
    list of list of float
        Allele frequency trajectories, one list per replicate.
        Each inner list has length n_generations + 1 (includes p0).
    """
    raise NotImplementedError("Not yet implemented.")