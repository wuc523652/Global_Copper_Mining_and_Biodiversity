"""Rendering functions adapted from the existing supplementary figure sources.
All quantitative inputs are supplied by render_supplementary.py; no analysis runs.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE = "#1f77b4"
ORANGE = "#ff7f0e"
AREA_BLUE = "#4C78A8"

def panel_title(ax, letter, title):
    """Separate letters so review markup can highlight renumbering alone."""
    label = ax.text(0.0, 1.055, letter, transform=ax.transAxes, fontsize=10,
                    fontweight="bold", va="bottom", ha="left", clip_on=False)
    ax.set_title(title, fontsize=9, pad=13)
    return label


def style_distance_axis(ax):
    ax.set_xlabel("Buffer distance (km)")
    ax.set_xticks([10, 20, 30, 40, 50])
    ax.set_xlim(8, 52)
    ax.axvline(30, color="0.5", linestyle="--", linewidth=0.8, zorder=1)


def build_plot(core, stability, diagnostics):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "lines.linewidth": 1.4,
        "lines.markersize": 4.5, "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })
    # Compact portrait layout leaves room for the Word caption on one page.
    fig, axes = plt.subplots(3, 2, figsize=(178 / 25.4, 188 / 25.4))
    fig.subplots_adjust(left=0.115, right=0.98, bottom=0.078, top=0.94,
                        wspace=0.42, hspace=0.75)
    letters = []
    x = core["buffer_distance_km"].to_numpy()

    ax = axes[0, 0]
    ax.plot(x, core["affected_species_count"], color=BLUE, marker="o", label="Affected species")
    ax.plot(x, core["priority_species_count"], color=ORANGE, marker="s",
            label="Conservation-priority\nspecies")
    style_distance_axis(ax)
    ax.set_ylabel("Species count")
    ax.set_ylim(-20, 1300)
    ax.set_yticks(range(0, 1201, 200))
    letters.append(panel_title(ax, "a", "Species outcomes"))
    ax.legend(frameon=False, loc="upper left", handlelength=1.8, borderaxespad=0.3)

    ax = axes[0, 1]
    share = core["priority_species_count"] / core["affected_species_count"] * 100
    ax.plot(x, share, color=ORANGE, marker="s")
    style_distance_axis(ax)
    ax.set_ylabel("Conservation-priority species (%)")
    ax.set_ylim(0, 44)
    ax.set_yticks(range(0, 41, 10))
    letters.append(panel_title(ax, "b", "Priority species share"))
    for r, value in zip(x, share):
        offset = (-7, 8) if r == 30 else (0, 8)
        alignment = "right" if r == 30 else "center"
        ax.annotate(f"{value:.2f}", (r, value), xytext=offset,
                    textcoords="offset points", ha=alignment, va="bottom", fontsize=8)
    for idx, increment in enumerate(diagnostics["priority_species_share_change_percentage_points"]):
        mid_x = (x[idx] + x[idx + 1]) / 2
        mid_y = (share.iloc[idx] + share.iloc[idx + 1]) / 2
        ax.annotate(f"+{increment:.2f} pp", (mid_x, mid_y), xytext=(3, -12),
                    textcoords="offset points", ha="center", va="top", fontsize=7.5)

    ax = axes[1, 0]
    efficiency = diagnostics["additional_priority_species_per_100000_km2"].to_numpy()
    positions = np.arange(4)
    ax.bar(positions, efficiency, width=0.66, color=AREA_BLUE, alpha=0.88)
    for position, value in zip(positions, efficiency):
        ax.text(position, value + 0.13, f"{value:.3f}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(positions, ["10–20", "20–30", "30–40", "40–50"])
    ax.set_ylim(0, 8.4)
    ax.set_yticks(range(0, 9, 2))
    ax.set_xlabel("Buffer expansion interval (km)")
    ax.set_ylabel("New priority species per\n100,000 km² of added buffer")
    letters.append(panel_title(ax, "c", "Marginal identification efficiency"))

    ax = axes[1, 1]
    ax.plot(x, core["global_strict_physical_exposure_area_km2"] / 1e6,
            marker="o", color=AREA_BLUE)
    style_distance_axis(ax)
    ax.set_ylabel("Affected habitat area (million km²)")
    ax.set_ylim(0, 6.5)
    ax.set_yticks(range(0, 7))
    letters.append(panel_title(ax, "d", "Affected habitats"))

    ax = axes[2, 0]
    ax.plot(x, stability["priority_set_jaccard_vs_30km"], marker="o", color=BLUE,
            label="Priority-set Jaccard")
    ax.plot(x, stability["investor_rank_correlation_vs_30km"], marker="s", color=ORANGE,
            label="Owner-country or region\ncross-border area-rank\ncorrelation")
    style_distance_axis(ax)
    ax.set_ylabel("Similarity to 30-km result")
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    letters.append(panel_title(ax, "e", "Stability"))
    stability_legend = ax.legend(frameon=True, facecolor="white", framealpha=1,
                                 edgecolor="white", loc="lower right", handlelength=1.8,
                                 borderaxespad=0.3, borderpad=0.1,
                                 labelspacing=0.1, fontsize=8)

    ax = axes[2, 1]
    domestic = core["domestic_equity_associated_area_km2"].to_numpy()
    foreign = core["foreign_equity_associated_area_km2"].to_numpy()
    unattributed = core["unattributed_area_km2"].to_numpy()
    total = domestic + foreign + unattributed
    if not (total > 0).all():
        raise ValueError("Ownership composition requires positive total area.")
    ax.stackplot(x, domestic / total * 100, foreign / total * 100,
                 unattributed / total * 100,
                 labels=["Domestic", "Cross-border", "Unattributed"],
                 colors=["#59A14F", "#E15759", "#BAB0AC"], alpha=0.85)
    style_distance_axis(ax)
    ax.set_ylabel("Share of country or region-attributed area (%)")
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    letters.append(panel_title(ax, "f", "Ownership composition"))
    ax.legend(frameon=True, facecolor="white", framealpha=0.92, edgecolor="0.85",
              loc="lower left", handlelength=1.8, borderaxespad=0.3)

    fig.canvas.draw()
    return fig
