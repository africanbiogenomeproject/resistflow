"""
Genetic drift models.

Currently implements Wright-Fisher binomial sampling:
the next generation's allele count is drawn from Binomial(2N, p),
where p is the post-selection frequency. Dividing by 2N gives
the new frequency.

    p_new = Binomial(2N, p) / 2N

Separating this module from wright_fisher.py allows future substitution
with alternative demographic models (e.g. Moran process, age-structured
populations) without modifying the core simulation loop.
"""

import numpy as np


def apply_drift(p, N):
    """
    Apply one generation of genetic drift via binomial sampling.

    Parameters
    ----------
    p : float
        Current frequency of the resistance allele A (post-selection) [0, 1].
    N : int
        Effective population size (Ne). Larger N → less drift.

    Returns
    -------
    float
        Post-drift allele frequency of A.
    """
    if not 0 <= p <= 1:
        raise ValueError(f"p must be in [0, 1], got {p}.")
    
    p_new = np.random.binomial(2 * N, p) / (2 * N)

    return p_new
