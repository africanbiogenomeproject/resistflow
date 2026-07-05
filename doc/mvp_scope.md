# ResistFlow - MVP Scope

**Document purpose:** Non-README reference for the AfricaBP manuscript core writing team  
**Status:** Agreed scope (confirmed with ThankGod Ebenezer, June 2026)

---

## Project summary

ResistFlow is an open-source proof-of-concept framework for exploring how insecticide-resistance alleles may spread through African *Anopheles* malaria vector populations under different intervention assumptions.

The MVP is not a finished forecasting system or a policy recommendation tool. It is a focused, reproducible proof of concept designed to support the AfricaBP Open Institute 2026 manuscript and provide a foundation for future development.

---

## Core research question

Can a simple forward-time Wright–Fisher simulation framework, initialised from real or MalariaGEN-derived allele-frequency observations, be used to explore how different insecticide intervention scenarios may influence the spread of a known resistance allele in African malaria vector populations?

---

## MVP objective

The MVP will deliver:

1. A Wright–Fisher simulation engine modelling selection and genetic drift
2. Configurable intervention scenario comparison
3. At least one MalariaGEN-backed example if feasible within timeline
4. A demonstration notebook, scenario comparison plots, and documentation

---

## In scope for MVP

- Single-locus Wright–Fisher simulation (selection + drift)
- Genotypic fitness parameters (wAA, wAa, waa)
- Effective population size as a tunable parameter
- Generations per year parameter (biological time conversion)
- Fitness costs in absence of insecticide (essential for rotation scenarios)
- Intervention scenarios: baseline, continuous pressure, reduced coverage, rotation
- One MalariaGEN-backed initialisation (target: Vgsc-1014F, Burkina Faso, 2012 baseline)
- Scenario comparison plots
- Demonstration notebook

---

## Out of scope for MVP

- Multi-locus simulation and epistasis
- Spatial structure or gene flow
- GWAS integration
- Selection signature analysis (H12, Tajima's D, iHS)
- Candidate variant prioritisation
- Detailed ecological or stage-structured modelling
- Operational policy recommendations
- Validated forecasting across multiple countries or species

---

## Fallback plan

If MalariaGEN integration takes longer than expected, the fallback deliverable is:

- Working Wright–Fisher engine
- Scenario comparison with synthetic starting data
- Clear documentation of the planned MalariaGEN integration path

---

## Manuscript contribution

ResistFlow will contribute to the AfricaBP Open Institute 2026 manuscript:

- An open-source proof-of-concept framework
- A reproducible computational workflow
- A case study connecting genomic surveillance data with simulation-based hypothesis generation
- A roadmap for community-driven future development

---

## Target validation case

| Parameter | Value |
|-----------|-------|
| Variant | Vgsc-1014F (kdr resistance, pyrethroid target-site) |
| Species | *Anopheles gambiae* |
| Country | Burkina Faso |
| Baseline year | 2012 |
| Baseline frequency | ~0.24 |
| Validation year | 2020 |
| Observed 2020 frequency | ~0.68 |
| Implied selection coefficient | s ≈ 0.15–0.20 |

---

## Software

- Language: Python 3.10+
- Package management: Poetry
- Core dependencies: numpy, matplotlib, pandas
- Optional: malariagen-data (MalariaGEN API integration)
- License: MIT
- Repository: https://github.com/africanbiogenomeproject/resistflow