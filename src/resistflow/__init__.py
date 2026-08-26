"""
ResistFlow: a forward-time simulation framework for exploring
insecticide resistance evolution in African malaria vectors.

An early proof-of-concept developed for the AfricaBP Open Institute 2026
manuscript.

Main entry points
-----------------
simulate
    Run a Wright-Fisher simulation and return allele frequency trajectories.
get_or_estimate_fitness_parameters
    Estimate fitness parameters from observed data (cached between runs).
fetch_allele_frequency
    Fetch an observed allele frequency from the MalariaGEN Ag3 API.
load_local_data
    Load observed allele frequencies from a local CSV, TSV, or Excel file.
get_scenario
    Retrieve a predefined intervention scenario.
plot_trajectories, plot_scenario_comparison
    Plot simulation output.
"""

__version__ = "0.1.0"

from resistflow.simulation.wright_fisher import simulate
from resistflow.scenarios.presets import get_scenario, SCENARIOS
from resistflow.analysis.parameter_estimation import (
    estimate_fitness_parameters,
    get_or_estimate_fitness_parameters,
)
from resistflow.data.local import load_local_data
from resistflow.visualization.plots import plot_trajectories, plot_scenario_comparison

__all__ = [
    "simulate",
    "get_scenario",
    "SCENARIOS",
    "estimate_fitness_parameters",
    "get_or_estimate_fitness_parameters",
    "load_local_data",
    "plot_trajectories",
    "plot_scenario_comparison",
]