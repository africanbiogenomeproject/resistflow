"""
MalariaGEN data integration.

Fetches allele frequency data from the MalariaGEN API (malariagen_data package)
to initialise simulations from real observations rather than synthetic starting points.

The simulation engine (wright_fisher.simulate) accepts p0 as a plain float.
This module is the only place that knows about MalariaGEN - the engine itself
stays independent of the data source.

Target for MVP
--------------
- Variant:    Vgsc-1014F (kdr resistance mutation, well-documented trajectory)
- Population: Burkina Faso
- Baseline:   2012 observed frequency (~0.24)
- Validation: 2020 observed frequency (~0.68)

References
----------
MalariaGEN Vector Observatory: https://www.malariagen.net/data
malariagen_data package:        https://malariagen.github.io/vector-data/
Ag3 release documentation:      https://malariagen.github.io/vector-data/ag3/
"""


def fetch_allele_frequency(variant, country, year):
    """
    Fetch observed allele frequency for a resistance variant from MalariaGEN.

    Parameters
    ----------
    variant : str
        Variant identifier (e.g. 'Vgsc-1014F').
    country : str
        Country name matching MalariaGEN metadata (e.g. 'Burkina Faso').
    year : int
        Observation year.

    Returns
    -------
    float
        Observed allele frequency of the resistance allele [0, 1].

    Raises
    ------
    NotImplementedError
        Until MalariaGEN integration is implemented.
    ImportError
        If malariagen_data is not installed (install with: poetry install --extras malariagen).
    """
    raise NotImplementedError(
        "MalariaGEN integration not yet implemented. "
        "Install optional dependencies with: poetry install --extras malariagen"
    )