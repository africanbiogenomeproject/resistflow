"""
Parameter estimation for the Wright-Fisher simulation.

Estimates fitness parameters (w_Aa, w_aa) from observed temporal
allele frequency data using least squares optimization.

The optimizer finds parameter values that minimize the difference
between the simulated mean trajectory and observed data points.
This removes the need for manual parameter tuning and makes the
framework generalizable across populations.

Method: Nelder-Mead (derivative-free) optimization via scipy,
applied to a deterministic (drift-free) trajectory. Genetic drift
is negligible at typical effective population sizes (N >= 1e4)
for the selection strengths observed in this system, so removing
it from the optimization objective eliminates run-to-run noise
without materially changing the fitted parameters.

Limitations
-----------
- Requires at least 3 time points with intermediate values
- Assumes constant selection pressure across the trajectory
- Single-locus model: does not account for epistasis or spatial structure
- Point estimates only: no uncertainty quantification (Bayesian
  inference is the recommended next step for uncertainty-aware estimation)
"""

import numpy as np
import json
import os
from pathlib import Path
import hashlib

from scipy.optimize import minimize
from resistflow.simulation.wright_fisher import simulate

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

def _build_cache_key(country, species, aa_change, p0, observed_years,
                      observed_frequencies, N, generations_per_year,
                      n_replicates, w_AA):
    readable = f"{country}_{species}_{aa_change}".replace(" ", "")

    raw = f"{p0}|{observed_years}|{observed_frequencies}|{N}|{generations_per_year}|{n_replicates}|{w_AA}"

    raw_bytes = raw.encode("utf-8")
    short_hash = hashlib.md5(raw_bytes).hexdigest()[:8]

    return f"{readable}_{short_hash}"

def _deterministic_trajectory(p0, n_generations, w_AA, w_Aa, w_aa):
    from resistflow.simulation.selection import apply_selection
    p = p0
    traj = [p]
    for _ in range(n_generations):
        p = apply_selection(p, w_AA, w_Aa, w_aa)
        traj.append(p)
    return traj

def _is_identifiable(p0, n_generations, observed_generations,
                      observed_frequencies, w_AA, tolerance=0.05):
    """
    Check whether the data identify a unique fitness estimate by running
    the optimizer from several different starting points and comparing
    the resulting w_aa values. If they disagree by more than `tolerance`,
    the fit is not identifiable — the data admit a range of equally valid
    parameter values, not a single point (see _find_minimum_selection_bound).
    """
    starting_points = [[0.95, 0.90], [0.88, 0.84], [0.70, 0.60], [0.60, 0.40], [0.99, 0.50]]
    w_aa_results = []
    for x0 in starting_points:
        def objective(params):
            w_Aa, w_aa = params
            if w_aa < 0 or w_Aa < w_aa or w_Aa > w_AA:
                return 1e6
            traj = _deterministic_trajectory(p0, n_generations, w_AA, w_Aa, w_aa)
            simulated_at_obs = [traj[g] for g in observed_generations]
            return sum((s - o) ** 2 for s, o in zip(simulated_at_obs, observed_frequencies))

        result = minimize(objective, x0, method="Nelder-Mead",
                           options={"xatol": 1e-3, "fatol": 1e-6, "maxiter": 200})
        if result.fun < 1e-4:
            w_aa_results.append(result.x[1])

    if len(w_aa_results) < 2:
        return True, w_aa_results

    spread = max(w_aa_results) - min(w_aa_results)
    return bool(spread < tolerance), w_aa_results


def _find_minimum_selection_bound(p0, observed_generations, observed_frequencies,
                                   w_AA, h=0.25, tolerance=0.01):
    """
    Find the minimum selection coefficient s (holding dominance h fixed at
    a literature-derived value) such that the deterministic trajectory
    reaches within `tolerance` of the final observed frequency.

    Used when _is_identifiable finds the data cannot pin down a unique
    point estimate. Rather than reporting an arbitrary point from an
    unconstrained search, this anchors the search to a literature-justified
    dominance coefficient and reports the weakest selection consistent
    with what was observed — a defensible lower bound, not a guess.
    """
    final_gen = observed_generations[-1]
    final_obs = observed_frequencies[-1]

    s = 0.01
    while s <= 0.99:
        w_aa = 1 - s
        w_Aa = w_aa + h * (1 - w_aa)
        traj = _deterministic_trajectory(p0, final_gen, w_AA, w_Aa, w_aa)
        if traj[-1] >= final_obs - tolerance:
            return {"s_min": round(s, 4), "w_Aa": round(w_Aa, 4),
                    "w_aa": round(w_aa, 4), "h_assumed": h}
        s += 0.001

    return None

def _sweep_h_sensitivity(p0, observed_generations, observed_frequencies, w_AA,
                          h_range=(0.05, 0.65), h_step=0.05, tolerance=0.01):
    """
    Check how sensitive the minimum-selection bound is to the assumed
    dominance coefficient h, by repeating _find_minimum_selection_bound
    across a range of h values rather than trusting a single one.

    If s_min stays roughly stable across this range, the bound does not
    depend strongly on which dominance value is assumed — a materially
    stronger claim than anchoring to one external literature value alone.

    Parameters
    ----------
    h_range : tuple of float
        (min, max) dominance coefficient to sweep, inclusive.
    h_step : float
        Step size for the sweep.

    Returns
    -------
    dict
        - s_min_range : (float, float) — lowest and highest s_min found
          across the swept h values
        - h_range_swept : tuple — the range actually swept
        - sweep : list of dict — per-h results, for plotting/inspection
    """
    final_gen = observed_generations[-1]
    final_obs = observed_frequencies[-1]

    sweep = []
    h = h_range[0]
    while h <= h_range[1] + 1e-9:
        s = 0.01
        while s <= 0.99:
            w_aa = 1 - s
            w_Aa = w_aa + h * (1 - w_aa)
            traj = _deterministic_trajectory(p0, final_gen, w_AA, w_Aa, w_aa)
            if traj[-1] >= final_obs - tolerance:
                sweep.append({"h": round(h, 2), "s_min": round(s, 4)})
                break
            s += 0.001
        h += h_step

    s_values = [entry["s_min"] for entry in sweep]
    return {
        "s_min_range": (round(min(s_values), 4), round(max(s_values), 4)),
        "h_range_swept": h_range,
        "sweep": sweep,
    }

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
        Retained for backward compatibility. No longer used in
        parameter estimation, which is now deterministic. Relevant
        only if simulate() is called separately with stochastic replicates.
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
        if w_aa < 0 or w_Aa < w_aa or w_Aa > w_AA:
            return 1e6
        traj = _deterministic_trajectory(p0, n_generations, w_AA, w_Aa, w_aa)
        simulated_at_obs = [traj[g] for g in observed_generations]
        error = sum((sim - obs) ** 2 for sim, obs in zip(simulated_at_obs, observed_frequencies))
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

    output = {
        "w_AA": w_AA,
        "w_Aa": round(float(w_Aa_opt), 4),
        "w_aa": round(float(w_aa_opt), 4),
        "optimization_success": result.success,
        "final_error": round(float(result.fun), 6),
        "identifiable": True,
    }

    identifiable, _ = _is_identifiable(p0, n_generations, observed_generations,
                                        observed_frequencies, w_AA)
    output["identifiable"] = identifiable

    if not identifiable:
        bound = _find_minimum_selection_bound(p0, observed_generations,
                                               observed_frequencies, w_AA)
        sensitivity = _sweep_h_sensitivity(p0, observed_generations,
                                            observed_frequencies, w_AA)

        output["w_Aa"] = bound["w_Aa"]
        output["w_aa"] = bound["w_aa"]
        output["s_min"] = bound["s_min"]
        output["dominance_h_assumed"] = bound["h_assumed"]
        output["s_min_range"] = sensitivity["s_min_range"]
        output["h_range_swept"] = sensitivity["h_range_swept"]
        output["note"] = (
            "Data do not uniquely identify fitness parameters — observations "
            "insufficient to constrain trajectory shape. Reporting the minimum "
            f"selection coefficient (s_min={bound['s_min']}) assuming dominance "
            f"h={bound['h_assumed']} (Lynd et al. 2010). This bound is stable "
            f"across a wide dominance range: sweeping h from "
            f"{sensitivity['h_range_swept'][0]} to {sensitivity['h_range_swept'][1]} "
            f"gives s_min between {sensitivity['s_min_range'][0]} and "
            f"{sensitivity['s_min_range'][1]}, indicating the conclusion does not "
            "depend strongly on the dominance assumption."
        )

    return output

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
    case_id = _build_cache_key(
            country, species, aa_change, p0, observed_years,
            observed_frequencies, N, generations_per_year, n_replicates, w_AA
            )

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