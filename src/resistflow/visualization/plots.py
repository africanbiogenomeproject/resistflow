"""
Trajectory visualization.

Plots allele-frequency trajectories from simulate() output.
Designed to produce figures suitable for the AfricaBP manuscript.

Planned outputs
---------------
- Trajectory plot: allele frequency over time (years), one line per replicate
- Scenario comparison: multiple scenarios overlaid on the same axes
- Optional: uncertainty bands (mean ± std across replicates)
- Optional: observed MalariaGEN data points overlaid for validation
"""


def plot_trajectories(
    trajectories,
    generations_per_year=10,
    title=None,
    scenario_label=None,
    observed_points=None,
    ax=None,
):
    """
    Plot allele-frequency trajectories over time.

    Parameters
    ----------
    trajectories : list of list of float
        Simulation output from simulate(): one trajectory per replicate.
        Each inner list has length n_generations + 1.
    generations_per_year : int
        Used to convert the x-axis from generations to years (default: 10).
    title : str, optional
        Plot title.
    scenario_label : str, optional
        Label for the legend entry.
    observed_points : dict, optional
        Real MalariaGEN observations for validation overlay.
        Format: {'year': [years], 'frequency': [frequencies]}
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates a new figure if None.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """
    raise NotImplementedError("Not yet implemented.")


def plot_scenario_comparison(results_by_scenario, generations_per_year=10, title=None):
    """
    Plot allele-frequency trajectories for multiple scenarios side by side.

    Parameters
    ----------
    results_by_scenario : dict
        Keys: scenario names (str). Values: trajectory lists from simulate().
    generations_per_year : int
        Used to convert the x-axis from generations to years (default: 10).
    title : str, optional
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object.
    """
    raise NotImplementedError("Not yet implemented.")