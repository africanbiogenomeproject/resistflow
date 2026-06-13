"""
Selection models.

Currently implements simple genotypic selection under Hardy-Weinberg assumptions:
allele frequency shifts deterministically each generation based on relative
genotype fitness values (wAA, wAa, waa).

The formula applied each generation:

    q  = 1 - p
    w̄  = p²·wAA + 2pq·wAa + q²·waa        (mean population fitness)
    p' = (p²·wAA + pq·wAa) / w̄             (post-selection frequency of A)

Separating this module from wright_fisher.py allows future extensions:
  - Stage-specific selection (larval vs. adult)
  - Frequency-dependent selection
  - Epistasis between loci (multi-locus extension)
without modifying the core simulation loop.
"""


def apply_selection(p, w_AA, w_Aa, w_aa):
    """
    Apply one generation of genotypic selection to allele frequency p.

    Parameters
    ----------
    p : float
        Current frequency of the resistance allele A [0, 1].
    w_AA : float
        Fitness of homozygous resistant genotype (AA).
        Under no insecticide pressure, set equal to w_Aa and w_aa.
    w_Aa : float
        Fitness of heterozygous genotype (Aa).
    w_aa : float
        Fitness of homozygous susceptible genotype (aa).

    Returns
    -------
    float
        Post-selection allele frequency of A.

    Raises
    ------
    ValueError
        If mean fitness w̄ is zero (all genotypes have zero fitness).
    """
    if not  0<= p <= 1:
        raise ValueError(f"p must be in [0, 1], got {p}.")
    if w_AA < 0 or w_Aa < 0 or w_aa < 0:
        raise ValueError("Fitness value must be non-negative.")
    q = 1 - p
    w_bar = p**2 * w_AA + 2 * p * q * w_Aa + q**2 * w_aa
    if w_bar == 0:
        raise ValueError(
            "Mean fitness is zero - all genotype fitness values are zero."
        )
    p_prime = (p**2 * w_AA + p * q * w_Aa) / w_bar
    return p_prime