"""Supplementary figure drawing functions; inputs are supplied by the runner."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


DOUBLE_COLUMN_WIDTH_IN = 183.0 / 25.4

def panel_label(ax, label, x=-0.08, y=1.04):
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold", clip_on=False, zorder=100)

BIN_FIELDS = ["area_le50_km2", "area_50_75_km2", "area_75_90_km2", "area_90_99_km2", "area_99_100_km2"]
BIN_LABELS = ["0–50th", "50–75th", "75–90th", "90–99th", "99–100th"]
BIN_COLORS = ["#d9d9d9", "#b7d4ea", "#74add1", "#f4a582", "#b2182b"]

def build_plot(data):
    for field in BIN_FIELDS + ["overlap_area_km2", "pct_ge90", "conflict_mine_count"]:
        data[field] = pd.to_numeric(data[field], errors="coerce").fillna(0)
    data = data[
        data["overlap_area_km2"].gt(0)
        & data["Country"].astype(str).str.casefold().ne("global")
    ].copy()
    data = data.sort_values(["overlap_area_km2", "Country"], ascending=[True, True])
    y = np.arange(len(data))

    height = max(9.0, 0.155 * len(data))
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_COLUMN_WIDTH_IN, height),
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.22},
        sharey=True,
    )
    left = np.zeros(len(data))
    for field, label, color in zip(BIN_FIELDS, BIN_LABELS, BIN_COLORS):
        values = data[field].to_numpy(float)
        ax_a.barh(y, values, left=left, height=0.72, color=color, edgecolor="white", linewidth=0.18, label=label)
        left += values
    ax_a.set_xscale("symlog", linthresh=100)
    ax_a.tick_params(axis="x", labelsize=7.2)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(data["Country"], fontsize=5.0)
    ax_a.set_xlabel("Mine potential influence area (km²; symlog scale)")
    ax_a.set_ylabel("Mine-host country")
    ax_a.grid(axis="x", color="#e4e4e4", linewidth=0.3)
    ax_a.legend(
        title="Rarity-weighted richness percentile",
        loc="lower left",
        bbox_to_anchor=(0, 1.005),
        ncol=3,
        fontsize=5.0,
        title_fontsize=5.2,
        columnspacing=0.7,
    )

    overlap = data["overlap_area_km2"].clip(lower=1)
    mines = data["conflict_mine_count"].clip(lower=0)
    scatter = ax_b.scatter(
        data["pct_ge90"],
        y,
        c=overlap,
        norm=LogNorm(vmin=max(1, float(overlap.min())), vmax=float(overlap.max())),
        cmap="viridis",
        s=np.clip(7 + np.sqrt(mines) * 6, 8, 60),
        edgecolor="white",
        linewidth=0.25,
    )
    ax_b.set_xlim(-2, 102)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_xlabel("Mine potential influence area in the top-decile\nrarity-weighted richness class (%)")
    ax_b.grid(axis="x", color="#e4e4e4", linewidth=0.3)
    cbar = fig.colorbar(
        scatter,
        ax=ax_b,
        orientation="vertical",
        pad=0.035,
        fraction=0.055,
        shrink=0.16,
        aspect=9,
        anchor=(0.0, 0.02),
    )
    cbar.set_label("Affected habitat area (km²; log scale)")
    cbar.ax.tick_params(labelsize=7.2)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=np.sqrt(np.clip(7 + np.sqrt(value) * 6, 8, 60)),
            markerfacecolor="#666666",
            markeredgecolor="white",
            label=str(value),
        )
        for value in [1, 10, 50, 150]
    ]
    ax_b.legend(
        handles=handles,
        title="Mines (n)",
        loc="lower left",
        bbox_to_anchor=(0, 1.005),
        ncol=4,
        fontsize=5.0,
        title_fontsize=5.2,
        columnspacing=0.7,
    )
    panel_label(ax_a, "a", x=-0.31, y=1.02)
    panel_label(ax_b, "b", x=-0.08, y=1.02)
    fig.subplots_adjust(left=0.30, right=0.98, top=0.965, bottom=0.055)
    return fig
