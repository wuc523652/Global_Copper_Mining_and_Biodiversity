"""Supplementary figure drawing functions; inputs are supplied by the runner."""
import matplotlib.pyplot as plt

import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

from matplotlib.ticker import FuncFormatter

DOUBLE_COLUMN_WIDTH_IN = 183.0 / 25.4

def longitude(value, pos=None):
    return f"{abs(int(round(value)))}°{'W' if value < 0 else 'E'}" if value else "0°"

def latitude(value, pos=None):
    return f"{abs(int(round(value)))}°{'S' if value < 0 else 'N'}" if value else "0°"

def build_plot(grid, boundary):
    centres = grid.geometry.centroid
    richness = pd.to_numeric(grid["total_richness"], errors="coerce").fillna(0).clip(0, 30)

    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_WIDTH_IN, 3.85))
    boundary.plot(ax=ax, aspect=None, facecolor="#f1f1f1", edgecolor="none", zorder=0)
    cmap = LinearSegmentedColormap.from_list("threatened_richness", ["#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"], N=256)
    points = ax.scatter(
        centres.x,
        centres.y,
        c=richness,
        s=0.08,
        marker="s",
        cmap=cmap,
        norm=Normalize(0, 30),
        linewidths=0,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    boundary.boundary.plot(ax=ax, aspect=None, color="#777777", linewidth=0.17, zorder=2)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 90)
    ax.set_xticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_yticks([-60, -30, 0, 30, 60, 90])
    ax.xaxis.set_major_formatter(FuncFormatter(longitude))
    ax.yaxis.set_major_formatter(FuncFormatter(latitude))
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, pad=1)
    ax.grid(color="#dedede", linewidth=0.3, zorder=0)
    ax.set_aspect("auto")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#888888")
        spine.set_linewidth(0.5)
    cbar = fig.colorbar(points, ax=ax, orientation="horizontal", pad=0.09, shrink=0.58, aspect=32)
    cbar.set_label("Richness of vertebrates threatened by mining (10-km grid)")
    cbar.set_ticks([0, 5, 10, 15, 20, 25, 30])
    cbar.set_ticklabels(["0", "5", "10", "15", "20", "25", "≥30"])
    fig.subplots_adjust(left=0.055, right=0.99, top=0.93, bottom=0.16)
    fig.set_size_inches(DOUBLE_COLUMN_WIDTH_IN, 4.20)
    ax.set_position([0.055, 0.245, 0.935, 0.668])
    ax.set_aspect("equal", adjustable="box")
    cbar.ax.set_position([0.225, 0.130, 0.540, 0.032])
    return fig
