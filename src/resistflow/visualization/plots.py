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
import numpy as np
import matplotlib.pyplot as plt

def plot_trajectories(
    trajectories,
    generations_per_year=10,
    title=None,
    scenario_label=None,
    observed_points=None,
    observed_label="observed",
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
        Observed allele frequencies to overlay for validation.
        Format: {'year': [years], 'frequency': [frequencies]}
    observed_label : str, optional
        Legend label for the observed points (default: "observed").
        Set to e.g. "observed (MalariaGEN)" when plotting API-sourced data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates a new figure if None.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object with the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    n_generations = len(trajectories[0]) - 1
    years = [g / generations_per_year for g in range(n_generations + 1)]

    for replicate in trajectories:
        ax.plot(years, replicate, color="steelblue", alpha=0.3, linewidth=1)

    mean_trajectory = np.mean(trajectories, axis=0)
    ax.plot(years, mean_trajectory, color="steelblue", linewidth=2.5,
            label=scenario_label or "mean")

    if observed_points is not None:
        ax.scatter(observed_points["year"], observed_points["frequency"],
                   color="black", zorder=5, label=observed_label)

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Resistance allele frequency")
    ax.set_ylim(0, 1)
    if title:
        ax.set_title(title)
    ax.legend()

    return ax


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
    colors = ["steelblue", "tomato", "seagreen"]

    fig, ax = plt.subplots(figsize=(10, 6))

    for (scenario_name, trajectories), color in zip(results_by_scenario.items(), colors):
        n_generations = len(trajectories[0]) - 1
        years = [g / generations_per_year for g in range(n_generations + 1)]

        for replicate in trajectories:
            ax.plot(years, replicate, color=color, alpha=0.2, linewidth=1)

        mean_trajectory = np.mean(trajectories, axis=0)
        ax.plot(years, mean_trajectory, color=color, linewidth=2.5, label=scenario_name)

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Resistance allele frequency")
    ax.set_ylim(0, 1)
    if title:
        ax.set_title(title)
    ax.legend()

    return fig