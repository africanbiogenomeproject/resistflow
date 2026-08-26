# ResistFlow

A forward-time simulation framework for exploring insecticide resistance evolution in African malaria vectors.

**Status:** Early proof of concept
**License:** MIT
**Lead developer:** Mohamed Laarej
**Affiliation:** AfricaBP Open Institute 2026

---

## What it does

ResistFlow simulates how insecticide-resistance alleles spread through *Anopheles* mosquito populations under different intervention assumptions. It initialises from real observed allele frequencies — either from the MalariaGEN Ag3 API or from your own data — and estimates population-specific fitness parameters automatically.

> This is a proof of concept, not a forecasting system. Results are exploratory scenario comparisons, not predictions. See [Scope and limitations](#scope-and-limitations).

---

## Installation

Requires Python 3.10–3.12 and [Poetry](https://python-poetry.org/).

```bash
git clone git@github.com:africanbiogenomeproject/resistflow.git
cd resistflow
poetry install
```

That's all you need for simulation, scenario comparison, parameter estimation, and local data. MalariaGEN API access is optional and covered in [section 3](#3-use-malariagen-data) below.

---

## Quickstart

Everything you need is importable from the top-level package:

```python
from resistflow import (
    simulate,
    get_scenario,
    get_or_estimate_fitness_parameters,
    load_local_data,
    plot_trajectories,
    plot_scenario_comparison,
)
```

### 1. Compare intervention scenarios

The simplest thing you can run — no external data required:

```python
from resistflow import simulate, get_scenario, plot_scenario_comparison

results = {}
for name in ["baseline", "continuous_pressure", "reduced_coverage"]:
    params = get_scenario(name)
    results[name] = simulate(
        p0=0.24,              # starting resistance allele frequency
        N=10_000,             # effective population size
        n_generations=80,     # 80 generations = 8 years at 10 gen/year
        w_AA=params["w_AA"],
        w_Aa=params["w_Aa"],
        w_aa=params["w_aa"],
        n_replicates=20,
    )

plot_scenario_comparison(results, generations_per_year=10)
```

### 2. Use your own data

Put your observations in a CSV (see [Input format](#input-format) below):

```python
from resistflow import (
    load_local_data,
    get_or_estimate_fitness_parameters,
    simulate,
    plot_trajectories,
)

data = load_local_data("examples/example_data.csv")

params = get_or_estimate_fitness_parameters(
    p0=data["p0"],
    observed_years=data["observed_years"],
    observed_frequencies=data["observed_frequencies"],
    country="MyPopulation",
    aa_change="L995F",
    species="gambiae",
    n_chromosomes=data["n_chromosomes"],
)
print(params)

trajectories = simulate(
    p0=data["p0"],
    N=10_000,
    n_generations=80,
    w_AA=params["w_AA"],
    w_Aa=params["w_Aa"],
    w_aa=params["w_aa"],
    n_replicates=20,
)

plot_trajectories(
    trajectories,
    generations_per_year=10,
    observed_points={
        "year": [0] + data["observed_years"],
        "frequency": [data["p0"]] + data["observed_frequencies"],
    },
)
```
Paths are relative to the repository root. If you're running from notebooks/, use ../examples/example_data.csv.

### 3. Use MalariaGEN data

**Before you start:** MalariaGEN data are no longer accessible anonymously. You need to request access once — this is free and open to everyone.

1. **Request access** by filling out [this form](https://docs.google.com/forms/d/e/1FAIpQLSdACxDwFhK0sOH9hX8QSZjjmuLXSSaXaSN7fW6CwCbQbUQQ7w/viewform). You'll need an email address associated with a Google account (a personal Gmail account, or a work address if your institution uses Google Workspace). Requests are granted subject to verification and agreement to reasonable use. See [Changes to MalariaGEN cloud data access](https://www.malariagen.net/article/changes-to-malariagen-cloud-data-access/) for background.

2. **Install with the extra:**

   ```bash
   poetry install --extras malariagen
   ```

   `malariagen-data` is already declared in `pyproject.toml` as an optional extra — there is no need to run `poetry add`. If you see a message saying the package already exists, use the command above instead.

3. **Authenticate.** Inside Google Colab this happens automatically. Anywhere else, install the [Google Cloud CLI](https://cloud.google.com/sdk/docs/install) and authenticate with the account you registered.

Then:

```python
from resistflow.data.malariagen import fetch_allele_frequency

p0 = fetch_allele_frequency(
    aa_change="L995F",
    country="Cameroon",
    year=2005,
    species="gambiae",
)
```

The first call is slow — the package fetches and caches sample metadata before returning anything.

Pass the fetched frequencies into `get_or_estimate_fitness_parameters` exactly as in the local-data example above.

> **Publication note:** some Ag3 releases carry publication embargoes, and using project data for publication may require prior permission from the originating partner studies. Check the [terms of use](https://malariagen.github.io/vector-data/ag3/ag3.0.html) for the release you are using before publishing, or email support@malariagen.net.
>
> If you would rather not depend on API access at all, `load_local_data` accepts your own observations — see above.

### Full worked example

`notebooks/01_demo.ipynb` runs everything end to end — scenario comparison, MalariaGEN data fetching, parameter estimation, and both validation cases.

---

## Input format

For `load_local_data`. Accepts `.csv`, `.tsv`, and `.xlsx`.

**Required columns:** `year`, `frequency`
**Optional columns:** `count`, `nobs`

```csv
year,frequency,count,nobs
2004,0.076923,2,26
2012,1.0,198,198
```

- `count` and `nobs` are chromosome counts. Mosquitoes are diploid, so 13 mosquitoes = 26 chromosomes. Supplying them lets the parameter estimator derive a principled fitting target when an observation is at fixation.
- Column names are matched case-insensitively. Accepted alternatives: `sampling_year` or `time` for year; `freq`, `allele_frequency` or `af` for frequency; `allele_count` or `ac` for count; `n` or `an` for nobs.
- Rows are sorted by year internally, so file order does not matter. The earliest year becomes the baseline.
- **Limitation:** time is assumed to be in whole years. Dates and sub-year resolution are not supported. If you need this, please [open an issue](https://github.com/africanbiogenomeproject/resistflow/issues).

---

## Key functions

| Function | Purpose |
|---|---|
| `simulate` | Run a Wright-Fisher simulation, returns trajectories |
| `get_or_estimate_fitness_parameters` | Fit fitness parameters to observed data (cached) |
| `load_local_data` | Load observations from CSV/TSV/Excel |
| `fetch_allele_frequency` | Fetch an observation from MalariaGEN Ag3 |
| `get_scenario` | Retrieve a predefined intervention scenario |
| `plot_trajectories` | Plot one scenario, optionally with observed points |
| `plot_scenario_comparison` | Plot multiple scenarios together |

Estimated parameters are cached at `~/.resistflow/params_cache.json`. Pass `force_reestimate=True` to recompute.

---

## Scope and limitations

**Currently supported**

- Single-locus Wright-Fisher simulation (selection and genetic drift)
- Configurable intervention scenarios
- Automatic fitness parameter estimation with identifiability detection
- MalariaGEN Ag3 and local file input

**Not supported**

- Multi-locus simulation or epistasis
- Spatial structure or gene flow between populations
- Fitness costs in the absence of insecticide
- Formal uncertainty quantification (point estimates only)
- Forward prediction beyond the observed window — the framework currently fits historical data rather than forecasting

Where observations cannot uniquely determine the fitness parameters — for example a baseline and a fixation endpoint with nothing in between — the estimator detects this and returns a bounded estimate with `identifiable: False`, rather than a misleading point value.

---

## Motivation

Genomic surveillance initiatives such as [MalariaGEN](https://www.malariagen.net) provide increasingly rich data for observing resistance-associated variation across Africa. These data describe where resistance is, with growing resolution. They are less readily used to explore where it may be going.

ResistFlow aims to complement that work with a lightweight, transparent tool for scenario exploration — moving from observation toward asking "what if" questions about intervention strategy, using the data surveillance programmes already collect.

---

## Project structure

```
resistflow/
├── examples/
│   └── example_data.csv          # Example input for load_local_data
├── doc/
│   ├── mvp_scope.md              # Original MVP scope
│   └── implementation_notes.md   # What was built, and what changed
├── notebooks/
│   └── 01_demo.ipynb             # Full worked example
├── src/resistflow/
│   ├── simulation/               # Wright-Fisher engine, selection, drift
│   ├── analysis/                 # Parameter estimation
│   ├── scenarios/                # Intervention presets
│   ├── data/                     # MalariaGEN and local data loading
│   └── visualization/            # Trajectory plots
└── tests/
```

---

## Roadmap

1. **Bayesian parameter estimation** — replace point estimates with posterior distributions, enabling uncertainty-aware forward prediction
2. **Multi-locus simulation** — epistatic interactions between co-occurring variants
3. **Spatial structure** — gene flow between geographic populations
4. **Additional validation cases** — further countries and resistance variants

---

## Contributing

Issues and feature requests are welcome, particularly around input formats and data sources — the current design reflects a narrow set of use cases and will benefit from wider input.

---

## Manuscript context

ResistFlow is a contribution to the **AfricaBP Open Institute 2026 regional workshop manuscript**, as a case study in connecting genomic surveillance data with simulation-based hypothesis generation.

---

## License

MIT. See [LICENSE](LICENSE).