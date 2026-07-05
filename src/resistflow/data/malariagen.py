"""
MalariaGEN data integration.

Fetches allele frequency data from the MalariaGEN API (malariagen_data package)
to initialise simulations from real observations rather than synthetic starting points.

The simulation engine (wright_fisher.simulate) accepts p0 as a plain float.
This module is the only place that knows about MalariaGEN - the engine itself
stays independent of the data source.

Target for MVP
--------------
- Variant:    Vgsc L995F (kdr, MalariaGEN numbering; equivalent to Vgsc-1014F in literature)
- Population: Burkina Faso, An. gambiae
- Baseline:   2004 observed frequency (Ag3 release 3.0)
- Validation: 2012 observed frequency (Ag3 release 3.0)

References
----------
MalariaGEN Vector Observatory: https://www.malariagen.net/data
malariagen_data package:        https://malariagen.github.io/vector-data/
Ag3 release documentation:      https://malariagen.github.io/vector-data/ag3/
"""
import malariagen_data
import pandas as pd

def fetch_allele_frequency(
    aa_change,
    country,
    year,
    transcript="AGAP004707-RD",
    species=None,
    sample_sets="3.0",
    min_cohort_size=10,
):
    """
    Fetch observed allele frequency for a resistance variant from MalariaGEN Ag3.

    Parameters
    ----------
    aa_change : str
        Amino acid change identifier using MalariaGEN numbering (e.g. 'L995F').
        Note: L995F in MalariaGEN corresponds to Vgsc-1014F in the conventional
        kdr literature.
    country : str
        Country name matching MalariaGEN metadata (e.g. 'Burkina Faso').
    year : int
        Observation year.
    transcript : str
        Gene transcript identifier (default: 'AGAP004707-RD' for Vgsc).
    species : str or None, optional
        Filter by species: 'gambiae' or 'coluzzii'. If None, averages across
        all available species cohorts for that year.
    sample_sets : str, optional
        MalariaGEN sample set release (default: '3.0').
    min_cohort_size : int, optional
        Minimum number of samples required for a cohort to be included (default: 10).

    Returns
    -------
    float
        Observed allele frequency of the resistance allele [0, 1].
        Averaged across admin1 regions if multiple cohorts exist for that year.

    Raises
    ------
    ValueError
        If the variant is not found, no data exists for the requested country/year,
        or all frequency values are NaN.
    """
    ag3 = malariagen_data.Ag3()

    df = ag3.aa_allele_frequencies(
        transcript=transcript,
        cohorts="admin1_year",
        sample_query=f"country == '{country}'",
        sample_sets=sample_sets,
        min_cohort_size=min_cohort_size,
    )

    row = df[(df["aa_pos"] == int(''.join(filter(str.isdigit, aa_change)))) &
             (df["alt_aa"] == aa_change[-1])]

    if row.empty:
        raise ValueError(f"Variant {aa_change} not found for {country}.")

    row = row.iloc[0]

    freq_cols = [col for col in df.columns
                 if col.startswith("frq_") and str(year) in col]

    if species is not None:
        abbrev = {"gambiae": "gamb", "coluzzii": "colu"}.get(species, species)
        freq_cols = [col for col in freq_cols if abbrev in col]

    if not freq_cols:
        raise ValueError(f"No data found for {country} in {year}.")

    freqs = [row[col] for col in freq_cols if not pd.isna(row[col])]

    if not freqs:
        raise ValueError(f"All frequency values are NaN for {country} in {year}.")

    return float(sum(freqs) / len(freqs))

