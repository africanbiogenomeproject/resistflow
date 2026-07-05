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
    Strong selection advantage for resistant homozygotes.

reduced_coverage
    Partial insecticide coverage (e.g. 50% of population exposed).
    Intermediate selection pressure — weaker than continuous.

rotation
    Planned: alternating insecticide classes on a fixed schedule.
    Requires time-varying selection logic — future extension.

Notes
-----
Fitness values here are illustrative placeholders based on published
bioassay data for Vgsc-1014F in An. gambiae under pyrethroid exposure.
They will be updated when historical validation is performed in Phase 1.
"""

# Scenario definitions

SCENARIOS = {
    "baseline": {
        "description": "No insecticide pressure. Pure genetic drift.",
        "w_AA": 1.0,
        "w_Aa": 1.0,
        "w_aa": 1.0,
    },
    "continuous_pressure": {
        "description": "Sustained high-coverage pyrethroid application.",
        "w_AA": 1.0,
        "w_Aa": 0.85,
        "w_aa": 0.70,
    },
    "reduced_coverage": {
        "description": "Partial insecticide coverage (~50% population exposed).",
        "w_AA": 1.0,
        "w_Aa": 0.92,
        "w_aa": 0.85,
    },
    # rotation: planned - requires time-varying selection, future extension
}


def get_scenario(name):
    """
    Retrieve a scenario configuration by name.

    Parameters
    ----------
    name : str
        Scenario name. One of: 'baseline', 'continuous_pressure', 'reduced_coverage'.

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