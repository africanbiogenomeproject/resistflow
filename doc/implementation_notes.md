# ResistFlow - Implementation Notes

**Document purpose:** Records what was actually built, what changed from the original MVP scope, and why. Intended as a reference for manuscript writing and future development.

**Status:** Post-MVP implementation, July 2026

---

## 1. What changed from the original MVP scope

### 1.1 Variant naming

The original scope referred to the target variant as **Vgsc-1014F**, following conventional kdr literature numbering (based on *Musca domestica* codon coordinates).

The MalariaGEN Ag3 API uses a different coordinate system based on transcript AGAP004707-RD (AgamP4.12). In this system, the same mutation is called **L995F**.

These are the same mutation. Clarkson et al. (2021) explicitly documents this mapping. All ResistFlow code uses L995F internally when interfacing with the API.

### 1.2 Validation window - Burkina Faso

The original plan was to validate using:
- Baseline: 2012 observed frequency (~24%)
- Validation: 2020 observed frequency (~68%)

When the Ag3 API was queried, the actual data showed:

| Year | Cohort | L995F frequency |
|---|---|---|
| 2004 | BF-07 gambiae | 7.7% |
| 2012 | BF-09 gambiae | 100% |
| 2014 | BF-09 gambiae | 100% |

By 2012, the gambiae population in Burkina Faso had already reached fixation. The sweep occurred entirely within the 2004–2012 window - before dense sampling began. The 2020 data point is not available in Ag3 release 3.0.

**Decision:** Pivot to 2004 (baseline) → 2012 (validation). Both points fully from the API. The original 24% figure cited in the proposal was from older literature about an earlier time point, not from Ag3.

**Implication:** Only two data points are available for Burkina Faso. The optimizer matches the endpoint (7.7% → 100%) but the trajectory shape is unconstrained. This is documented honestly as an endpoint validation, not a trajectory validation.

### 1.3 Fitness parameters - presets

The original presets used placeholder fitness values. These were revised twice:

**First revision:** Approximate values estimated from general knowledge of pyrethroid selection pressure on Vgsc L995F: s ≈ 0.175, h ≈ 0.9. These were working estimates, not taken from a specific source.

**Second revision (current):** After reading Lynd et al. (2010), which explicitly reports maximum likelihood estimates from field temporal data, the values were updated to match the published empirical estimates:
- Selection coefficient: **s = 0.16**
- Dominance coefficient: **h = 0.25**

h = 0.25 indicates L995F is largely recessive - heterozygotes gain relatively little fitness benefit under insecticide pressure. This is the opposite of the h = 0.9 assumption in the first revision.

Derived fitness values for `continuous_pressure` preset:
- w_AA = 1.000 (homozygous resistant)
- w_Aa = 1 − (1 − h) × s = 1 − 0.75 × 0.16 = **0.880**
- w_aa = 1 − s = 1 − 0.16 = **0.840**

`reduced_coverage` values are a 50/50 weighted average of `continuous_pressure` and `baseline`, modeling 50% population exposure.

**Reference:** Lynd A, Weetman D, Barbosa S, et al. Field, genetic, and modeling approaches show strong positive selection acting upon an insecticide resistance mutation in *Anopheles gambiae* s.s. *Molecular Biology and Evolution* 27(5):1117–1125, 2010.

---

## 2. What was added beyond the original MVP scope

### 2.1 Automatic parameter estimation

The original MVP did not include parameter estimation - fitness values were intended to be set manually from literature.

During validation, it became clear that:
- Manual values from literature may not apply to specific populations
- Cameroon and Burkina Faso showed different trajectories under the same parameters
- Manually searching for better parameters is scientifically indefensible

**Added:** `analysis/parameter_estimation.py` - Nelder-Mead optimization (scipy) that automatically finds w_Aa and w_aa values minimizing the difference between simulated mean trajectory and observed data points.

The function `get_or_estimate_fitness_parameters()` wraps the optimizer with a JSON cache at `~/.resistflow/params_cache.json`. First call runs the optimizer (slow); subsequent calls load from cache instantly.

**Initial guess** for the optimizer uses Lynd et al. (2010) values (w_Aa=0.880, w_aa=0.840) for consistency with the preset parameterization.

### 2.2 Cameroon validation case

The original MVP specified only Burkina Faso as the validation case.

During Ag3 data exploration, Cameroon was identified as having three gambiae time points with intermediate values:
- 2005: 20.2%
- 2009: 52.3%
- 2013: 93.3%

Three time points with an intermediate observation at year 4 allow genuine trajectory validation - not just endpoint matching. Cameroon was added as a second validation case.

**Auto-estimated fitness parameters:**
- Cameroon: w_Aa = 0.926, w_aa = 0.911 (weaker selection pressure)
- Burkina Faso: w_Aa = 0.805, w_aa = 0.713 (stronger selection pressure)

The regional difference is biologically meaningful and consistent with known differences in pyrethroid use intensity across West and Central Africa.

---

## 3. Validation results summary

### 3.1 Burkina Faso - endpoint validation

| Parameter | Value |
|---|---|
| Species | *An. gambiae* |
| Data source | MalariaGEN Ag3 release 3.0 |
| Baseline | 2004, p₀ = 0.077 |
| Validation point | 2012, observed = 1.000 |
| Auto-estimated w_AA | 1.000 |
| Auto-estimated w_Aa | 0.805 |
| Auto-estimated w_aa | 0.713 |
| Optimization error | 0.000 (endpoint matched) |

**Limitation:** Only two data points. Trajectory shape is unconstrained - the optimizer matches the endpoint but cannot distinguish between trajectories that pass through different intermediate frequencies.

### 3.2 Cameroon - trajectory validation

| Parameter | Value |
|---|---|
| Species | *An. gambiae* |
| Data source | MalariaGEN Ag3 release 3.0 |
| Baseline | 2005, p₀ = 0.202 |
| Intermediate point | 2009, observed = 0.523 |
| Final point | 2013, observed = 0.933 |
| Auto-estimated w_AA | 1.000 |
| Auto-estimated w_Aa | 0.926 |
| Auto-estimated w_aa | 0.911 |
| Optimization error | 0.000022 (near-perfect fit) |

**Result:** Simulated mean passes close to both the intermediate (year 4) and final (year 8) observed points. This constitutes a genuine trajectory validation.

---

## 4. Honest limitations

- **Single-locus model:** Does not capture epistasis with other resistance loci (e.g. N1570Y co-occurring with L995F in Burkina Faso, documented in Clarkson et al. 2021)
- **Constant selection pressure:** Does not model seasonal variation in insecticide use, changes in intervention over time, or fitness costs in absence of insecticide
- **No spatial structure:** Populations treated as panmictic; no gene flow or migration between regions
- **Parameter estimation without uncertainty:** Point estimates only. Confidence intervals on w_Aa and w_aa are not provided. Bayesian inference is the appropriate next step.
- **Burkina Faso data gap:** The sweep in Burkina Faso gambiae occurred before dense Ag3 sampling. The 2004–2012 window has only two observations. Dense temporal data from this period would substantially improve parameter constraints.

---

## 5. Future development priorities

1. **Bayesian parameter estimation** - replace Nelder-Mead point estimates with posterior distributions over fitness parameters using ABC or PyMC. Requires more temporal data points per population.
2. **Fitness costs** - model reduced fitness of resistant genotypes in absence of insecticide pressure. Essential for rotation scenario modeling.
3. **Multi-locus simulation** - extend engine to handle multiple resistance loci simultaneously, including epistatic interactions.
4. **Spatial structure** - metapopulation model with gene flow between geographic regions.
5. **Additional validation cases** - expand to other countries with sufficient Ag3 temporal coverage (Ghana, Mali) and other resistance variants (Ace-1, CYP6P3).
6. **MalariaGEN API integration** - propose ResistFlow as a community tool via GitHub issue after validation results are published.

---

## 6. References

Clarkson CS, Miles A, Harding NJ, et al. The genetic architecture of target-site resistance to pyrethroid insecticides in the African malaria vectors *Anopheles gambiae* and *Anopheles coluzzii*. *Molecular Ecology* 30:5303–5317, 2021. https://doi.org/10.1111/mec.15845

Lynd A, Weetman D, Barbosa S, et al. Field, genetic, and modeling approaches show strong positive selection acting upon an insecticide resistance mutation in *Anopheles gambiae* s.s. *Molecular Biology and Evolution* 27(5):1117–1125, 2010. https://doi.org/10.1093/molbev/msq002

MalariaGEN Vector Observatory. Ag3 release 3.0. https://www.malariagen.net/data