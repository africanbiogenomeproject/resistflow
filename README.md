# ResistFlow

A forward-time simulation framework for exploring insecticide resistance evolution in African malaria vectors.

**Status:** Early proof of concept - MVP in active development  
**License:** MIT  
**Lead developer:** Mohamed Laarej  
**Affiliation:** AfricaBP Open Institute 2026

---

## Overview

ResistFlow is an open-source proof-of-concept framework for exploring how insecticide-resistance alleles may spread through African *Anopheles* malaria vector populations under different intervention assumptions.

The project builds on the foundation created by genomic surveillance initiatives such as [MalariaGEN](https://www.malariagen.net), which provide increasingly rich data for observing resistance-associated variation across Africa. ResistFlow aims to complement this work by offering a lightweight, transparent tool for exploring how resistance alleles might spread under different intervention assumptions - moving from observation toward scenario exploration.

ResistFlow addresses this by combining a transparent Wright–Fisher simulation engine with configurable intervention scenarios and, where feasible, real allele-frequency data from MalariaGEN.

> **Note:** This is an early proof-of-concept. It is not a finished forecasting system or a policy recommendation tool. Results should be interpreted as exploratory scenario comparisons, not predictions.

---

## What ResistFlow does

- Simulates resistance allele-frequency trajectories over time using a Wright–Fisher model
- Models natural selection (genotypic fitness) and genetic drift (finite population stochasticity)
- Compares configurable intervention scenarios: baseline, continuous insecticide pressure, reduced coverage, rotation
- Supports initialization from MalariaGEN-derived allele frequencies (planned)

## What ResistFlow does not do (yet)

- Multi-locus simulation or epistasis
- Spatial structure or gene flow between populations
- GWAS integration or candidate variant prioritization
- Selection signature analysis
- Operational policy recommendations or validated forecasting

---

## Installation

Requires Python 3.10+ and [Poetry](https://python-poetry.org/).

```bash
git clone git@github.com:africanbiogenomeproject/resistflow.git
cd resistflow
poetry install
```

To include MalariaGEN data integration (optional, heavy dependencies):

```bash
poetry install --extras malariagen
```

---

## Project structure

```
resistflow/
│
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
│
├── doc/
│   └── mvp_scope.md              # MVP scope and design decisions
│
├── notebooks/
│   └── 01_demo.ipynb             # Scenario comparison walkthrough (planned)
│
├── src/
│   └── resistflow/
│       ├── __init__.py
│       │
│       ├── simulation/
│       │   ├── __init__.py
│       │   ├── wright_fisher.py  # Core generation loop (selection → drift)
│       │   ├── selection.py      # Genotypic selection model (pluggable)
│       │   └── drift.py          # Binomial drift (pluggable)
│       │
│       ├── scenarios/
│       │   ├── __init__.py
│       │   └── presets.py        # Named intervention configurations
│       │
│       ├── data/
│       │   ├── __init__.py
│       │   └── malariagen.py     # MalariaGEN API integration (planned)
│       │
│       └── visualization/
│           ├── __init__.py
│           └── plots.py          # Trajectory and scenario plots
│
└── tests/
    ├── __init__.py
    ├── test_wright_fisher.py
    ├── test_selection.py
    └── test_drift.py
```

---

## Roadmap

**Phase 1 - MVP (current)**
- Wright–Fisher simulation engine (selection + drift)
- Configurable intervention scenarios
- Example notebook with scenario comparison
- MalariaGEN-backed initialization if feasible within timeline

**Phase 2 - Future development**
- GWAS integration and selection signature analysis
- Temporal allele-frequency change analysis
- Candidate variant risk scoring
- Multi-locus simulation
- Spatial structure and gene flow

---

## Manuscript context

ResistFlow is being developed as a contribution to the **AfricaBP Open Institute 2026 regional workshop manuscript**. It represents a case study in connecting genomic surveillance data with simulation-based hypothesis generation, and an example of lightweight computational biology tooling built for African malaria vector research contexts.

---

## License

MIT License. See [LICENSE](LICENSE) for details.