"""
Local data loading.

Allows ResistFlow to be used with a researcher's own allele frequency
observations rather than only MalariaGEN Ag3 data. The simulation engine
and parameter estimator take plain numbers and have no knowledge of where
those numbers came from, so a local file is a first-class input path
rather than a workaround.

Expected file format (.csv, .tsv, or .xlsx)
-------------------------------------------
Required columns:
    year       - sampling year (whole years)
    frequency  - observed allele frequency [0, 1]

Optional columns:
    count      - number of chromosomes carrying the variant
    nobs       - total chromosomes sampled

Supplying count and nobs enables a principled fitting target when an
observation is at fixation (see analysis.parameter_estimation).

Column names are matched case-insensitively against a set of common
aliases. The earliest year is taken as the baseline (p0).

Known limitation
----------------
Time is assumed to be in whole years. Dates and sub-year resolution are
not currently supported. If your data needs this, please open an issue.
"""

from pathlib import Path
import pandas as pd

COLUMN_ALIASES = {
    "year": ["year", "sampling_year", "time"],
    "frequency": ["frequency", "freq", "allele_frequency", "af"],
    "count": ["count", "allele_count", "ac"],
    "nobs": ["nobs", "n", "an"],
}


def _resolve_column(df, canonical):
    """Find which column in df corresponds to a canonical name, or None."""
    lowered = {col.lower().strip(): col for col in df.columns}
    for alias in COLUMN_ALIASES[canonical]:
        if alias in lowered:
            return lowered[alias]
    return None


def load_local_data(path):
    """
    Load observed allele frequency data from a local file.

    Parameters
    ----------
    path : str or Path
        Path to a .csv, .tsv, or .xlsx file. See module docstring for
        the expected format.

    Returns
    -------
    dict
        - p0 : float — baseline frequency (earliest year)
        - observed_years : list of float — years elapsed since baseline
        - observed_frequencies : list of float — frequencies at those years
        - n_chromosomes : int or None — nobs for the final observation
        - p0_count : int or None — count for the baseline observation
        - p0_nobs : int or None — nobs for the baseline observation
        - baseline_year : int — the calendar year used as baseline

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the file type is unsupported, a required column is missing,
        or fewer than two observations are present.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")

    suffix = path.suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t")
    elif suffix in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        raise ValueError(
            f"Unsupported file type '{suffix}'. Expected .csv, .tsv, or .xlsx."
        )

    year_col = _resolve_column(df, "year")
    freq_col = _resolve_column(df, "frequency")

    if year_col is None or freq_col is None:
        missing = []
        if year_col is None:
            missing.append(f"year (accepted: {', '.join(COLUMN_ALIASES['year'])})")
        if freq_col is None:
            missing.append(f"frequency (accepted: {', '.join(COLUMN_ALIASES['frequency'])})")
        raise ValueError(
            "Missing required column(s): " + "; ".join(missing) +
            f". Found columns: {list(df.columns)}"
        )

    count_col = _resolve_column(df, "count")
    nobs_col = _resolve_column(df, "nobs")

    df = df.sort_values(year_col).reset_index(drop=True)

    if len(df) < 2:
        raise ValueError(
            f"At least two observations are required; found {len(df)}."
        )

    baseline_year = df.loc[0, year_col]
    p0 = float(df.loc[0, freq_col])

    rest = df.iloc[1:]
    observed_years = [float(y - baseline_year) for y in rest[year_col]]
    observed_frequencies = [float(f) for f in rest[freq_col]]

    n_chromosomes = int(df.loc[len(df) - 1, nobs_col]) if nobs_col else None
    p0_count = int(df.loc[0, count_col]) if count_col else None
    p0_nobs = int(df.loc[0, nobs_col]) if nobs_col else None

    return {
        "p0": p0,
        "observed_years": observed_years,
        "observed_frequencies": observed_frequencies,
        "n_chromosomes": n_chromosomes,
        "p0_count": p0_count,
        "p0_nobs": p0_nobs,
        "baseline_year": int(baseline_year),
    }