"""
Pre-defined intervention scenario presets.

Each scenario is a dictionary of parameters compatible with simulate().
Scenarios differ in their fitness values (wAA, wAa, waa), which represent
the net survival advantage of each genotype under a given intervention.

Included scenarios
------------------
baseline
    No insecticide pressure. All genotypes have equal fitness.
    Allele frequency changes only through drift.

continuous_pressure
    Sustained pyrethroid application at high coverage.
    Fitness values derived from Lynd et al. (2010) maximum likelihood
    estimates for Vgsc L1014F (L995F) in An. gambiae s.s.:
    selection coefficient s = 0.16, dominance coefficient h = 0.25.
    h = 0.25 indicates L995F is largely recessive — heterozygotes
    gain relatively little fitness benefit under insecticide pressure.

reduced_coverage
    Partial insecticide coverage (~50% population exposed).
    Effective fitness values scaled from continuous_pressure
    assuming 50% exposure rate.

rotation
    Planned: alternating insecticide classes on a fixed schedule.
    Requires time-varying selection logic - future extension.

Notes
-----
Fitness values for continuous_pressure are grounded in published
maximum likelihood estimates from field data (Lynd et al. 2010,
Mol Biol Evol 27(5):1117-1125). These are illustrative presets
for scenario comparison. For population-specific validation,
use get_or_estimate_fitness_parameters() in
analysis/parameter_estimation.py.
"""

# NOTE: These scenarios are designed for illustrative comparison only.
# For validation against real MalariaGEN data, use
# get_or_estimate_fitness_parameters() in analysis/parameter_estimation.py
# which automatically estimates population-specific fitness values.


# Scenario definitions
SCENARIOS = {
    "baseline": {
        "description": "No insecticide pressure. Pure genetic drift.",
        "w_AA": 1.000,
        "w_Aa": 1.000,
        "w_aa": 1.000,
    },
    "continuous_pressure": {
        "description": (
            "Sustained pyrethroid application at high coverage. "
            "Fitness values from Lynd et al. (2010) ML estimates for "
            "Vgsc L1014F (L995F) in An. gambiae s.s.: s=0.16, h=0.25. "
            "h=0.25 indicates largely recessive resistance."
        ),
        "w_AA": 1.000,
        "w_Aa": 0.880,  # 1 - (1 - h) * s = 1 - 0.75 * 0.16
        "w_aa": 0.840,  # 1 - s = 1 - 0.16
    },
    "reduced_coverage": {
        "description": (
            "Partial insecticide coverage (~50% population exposed). "
            "Effective fitness values are a 50/50 weighted average of "
            "continuous_pressure and baseline fitness values."
        ),
        "w_AA": 1.000,
        "w_Aa": 0.940,  # 0.5 * 0.880 + 0.5 * 1.0
        "w_aa": 0.920,  # 0.5 * 0.840 + 0.5 * 1.0
    },
    # rotation: planned - requires time-varying selection, future extension
}


def get_scenario(name):
    """
    Retrieve a scenario configuration by name.

    Parameters
    ----------
    name : str
        Scenario name. One of: 'baseline', 'continuous_pressure',
        'reduced_coverage'.

    Returns
    -------
    dict
        Scenario parameters compatible with simulate().

    Raises
    ------
    KeyError
        If the scenario name is not recognised.
    """
    if name not in SCENARIOS:
        available = ", ".join(SCENARIOS.keys())
        raise KeyError(f"Unknown scenario '{name}'. Available: {available}")
    return SCENARIOS[name].copy()