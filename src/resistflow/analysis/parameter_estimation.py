"""
Parameter estimation for the Wright-Fisher simulation.

Estimates fitness parameters (w_Aa, w_aa) from observed temporal
allele frequency data using least squares optimization.

The optimizer finds parameter values that minimize the difference
between the simulated mean trajectory and observed data points.
This removes the need for manual parameter tuning and makes the
framework generalizable across populations.

Method: Nelder-Mead (derivative-free) optimization via scipy.
Suitable for stochastic objective functions.

Limitations
-----------
- Requires at least 3 time points with intermediate values
- Assumes constant selection pressure across the trajectory
- Single-locus model: does not account for epistasis or spatial structure
- Point estimates only: no uncertainty quantification (Bayesian
  inference is the recommended next step for uncertainty-aware estimation)
"""

import numpy as np
from scipy.optimize import minimize
from resistflow.simulation.wright_fisher import simulate
import json
import os
from pathlib import Path

DEFAULT_CACHE_PATH = Path.home() / ".resistflow" / "params_cache.json"

def _load_cache(cache_path):
    if cache_path.exists():
        with open(cache_path, "r") as f:
            return json.load(f)
    return {}


def _save_cache(cache, cache_path):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2)


def estimate_fitness_parameters(
    p0,
    observed_years,
    observed_frequencies,
    N=10_000,
    generations_per_year=10,
    n_replicates=50,
    w_AA=1.0,
):
    """
    Estimate fitness parameters that best fit observed allele frequency data.

    Parameters
    ----------
    p0 : float
        Initial allele frequency at year 0.
    observed_years : list of float
        Years at which frequencies were observed, relative to p0.
        Example: [4, 8] for observations 4 and 8 years after baseline.
    observed_frequencies : list of float
        Observed allele frequencies at each year in observed_years.
    N : int
        Effective population size (default: 10_000).
    generations_per_year : int
        Biological generations per year (default: 10).
    n_replicates : int
        Number of simulation replicates per objective evaluation (default: 50).
        Higher values reduce stochastic noise but increase runtime.
    w_AA : float
        Fitness of resistant homozygote, fixed at 1.0 (reference genotype).

    Returns
    -------
    dict
        Optimized parameters and diagnostics:
        - w_AA : float (fixed at 1.0)
        - w_Aa : float (optimized)
        - w_aa : float (optimized)
        - optimization_success : bool
        - final_error : float (sum of squared residuals)
    """
    observed_generations = [int(y * generations_per_year) for y in observed_years]
    n_generations = max(observed_generations)

    def objective(params):
        w_Aa, w_aa = params

        # Enforce biological ordering: w_aa <= w_Aa <= w_AA
        if w_aa < 0 or w_Aa < w_aa or w_Aa > w_AA:
            return 1e6  # penalty for biologically invalid parameters

        trajectories = simulate(
            p0=p0,
            N=N,
            n_generations=n_generations,
            w_AA=w_AA,
            w_Aa=w_Aa,
            w_aa=w_aa,
            n_replicates=n_replicates,
        )

        mean_traj = np.mean(trajectories, axis=0)

        simulated_at_obs = [mean_traj[g] for g in observed_generations]

        error = sum(
            (sim - obs) ** 2
            for sim, obs in zip(simulated_at_obs, observed_frequencies)
        )
        return error

    # Initial guess from Lynd et al. (2010) ML estimates for L995F
    # s=0.16, h=0.25 → w_Aa=0.880, w_aa=0.840
    x0 = [0.880, 0.840]

    result = minimize(
        objective,
        x0,
        method="Nelder-Mead",
        options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 200},
    )

    w_Aa_opt, w_aa_opt = result.x

    return {
        "w_AA": w_AA,
        "w_Aa": round(float(w_Aa_opt), 4),
        "w_aa": round(float(w_aa_opt), 4),
        "optimization_success": result.success,
        "final_error": round(float(result.fun), 6),
    }


def get_or_estimate_fitness_parameters(
    p0,
    observed_years,
    observed_frequencies,
    country,
    aa_change,
    species,
    N=10_000,
    generations_per_year=10,
    n_replicates=50,
    w_AA=1.0,
    force_reestimate=False,
    cache_path=DEFAULT_CACHE_PATH,
):
    """
    Return cached fitness parameters for a population case, or estimate
    them if not yet cached.

    Parameters
    ----------
    p0 : float
        Initial allele frequency at year 0.
    observed_years : list of float
        Years at which frequencies were observed, relative to p0.
    observed_frequencies : list of float
        Observed allele frequencies at each year in observed_years.
    country : str
        Country name (used as part of cache key).
    aa_change : str
        Variant identifier (used as part of cache key).
    species : str
        Species name (used as part of cache key).
    N : int
        Effective population size (default: 10_000).
    generations_per_year : int
        Biological generations per year (default: 10).
    n_replicates : int
        Number of simulation replicates per objective evaluation (default: 50).
    w_AA : float
        Fitness of resistant homozygote, fixed at 1.0.
    force_reestimate : bool
        If True, re-run the optimizer even if cached params exist (default: False).
    cache_path : Path
        Path to the JSON cache file (default: ~/.resistflow/params_cache.json).

    Returns
    -------
    dict
        Fitness parameters and diagnostics. See estimate_fitness_parameters().
    """
    cache_path = Path(cache_path)
    case_id = f"{country}_{species}_{aa_change}".replace(" ", "")

    cache = _load_cache(cache_path)

    if case_id in cache and not force_reestimate:
        print(f"Loaded cached parameters for {case_id}.")
        return cache[case_id]

    print(f"Estimating parameters for {case_id} — this may take a minute...")
    params = estimate_fitness_parameters(
        p0=p0,
        observed_years=observed_years,
        observed_frequencies=observed_frequencies,
        N=N,
        generations_per_year=generations_per_year,
        n_replicates=n_replicates,
        w_AA=w_AA,
    )

    cache[case_id] = params
    _save_cache(cache, cache_path)
    print(f"Parameters estimated and cached for {case_id}.")

    return params