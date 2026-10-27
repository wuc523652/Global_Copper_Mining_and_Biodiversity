"""Rendering functions adapted from the existing supplementary figure sources.
All quantitative inputs are supplied by render_supplementary.py; no analysis runs.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SUPP_PALETTE = np.array([
    ["#e8e8e8", "#ace4e4", "#5ac8c8"],
    ["#dfb0d6", "#a5add3", "#5698b9"],
    ["#be64ac", "#8c62aa", "#3b4994"],
])


def classify_three(value, breaks):
    if pd.isna(value):
        return -1
    return int(np.searchsorted(breaks, value, side="left"))


def _draw_resource_map(ax, world, stats):
    data = stats.copy()
    data["resource_iso3"] = data["resource_iso3"].astype(str).str.upper().str.strip()
    area_field = "physical_exposure_union_area_km2"
    pct_field = "foreign_share_of_traceable_pct"
    data[area_field] = pd.to_numeric(data[area_field], errors="coerce")
    data[pct_field] = pd.to_numeric(data[pct_field], errors="coerce").clip(0, 100)
    area_breaks = np.array([1000.0, 10000.0])
    pct_breaks = np.asarray([20.0, 50.0], dtype=float)
    mapped = world.merge(data, left_on="iso3", right_on="resource_iso3", how="left")
    mapped["area_class"] = mapped[area_field].apply(lambda value: classify_three(value, area_breaks))
    mapped["pct_class"] = mapped[pct_field].apply(lambda value: classify_three(value, pct_breaks))
    mapped["map_color"] = "#e6e6e6"
    valid = mapped["area_class"].ge(0) & mapped["pct_class"].ge(0) & mapped[area_field].gt(0)
    palette = SUPP_PALETTE
    mapped.loc[valid, "map_color"] = [palette[iy, ix] for iy, ix in zip(mapped.loc[valid, "area_class"].astype(int), mapped.loc[valid, "pct_class"].astype(int))]
    mapped.plot(ax=ax, color=mapped["map_color"], edgecolor="#777777", linewidth=0.14)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 90)
    ax.set_xticks([-120, -60, 0, 60, 120])
    ax.set_yticks([-30, 0, 30, 60])
    ax.set_xticklabels(["120°W", "60°W", "0°", "60°E", "120°E"])
    ax.set_yticklabels(["30°S", "0°", "30°N", "60°N"])
    ax.grid(color="#dedede", linewidth=0.3)
    ax.tick_params(length=0)
    inset = ax.inset_axes([0.025, 0.18, 0.22, 0.38])
    for iy in range(3):
        for ix in range(3):
            inset.add_patch(plt.Rectangle((ix, iy), 1, 1, facecolor=palette[iy, ix], edgecolor="white", linewidth=0.4))
    inset.set_xlim(0, 3)
    inset.set_ylim(0, 3)
    inset.set_aspect("equal")
    inset.set_xticks([0.5, 1.5, 2.5])
    inset.set_xticklabels(["<20%", "20–50%", ">50%"], fontsize=5.7)
    inset.set_yticks([0.5, 1.5, 2.5])
    inset.set_yticklabels(["<10³", "10³–10⁴", ">10⁴"], fontsize=5.7)
    inset.tick_params(length=0, pad=1)
    for spine in inset.spines.values():
        spine.set_visible(False)
    inset.set_title("Affected habitat area", fontsize=5.7, pad=1.5)
    inset.set_xlabel("Share attributed to equity holders\nregistered elsewhere", fontsize=5.7, labelpad=1)
    return int(valid.sum()), area_breaks


def _ranking_panel(ax, data, country_field, value_field, xlabel, color):
    values = pd.to_numeric(data[value_field], errors="coerce").fillna(0)
    plotted = data.assign(_value=values).sort_values("_value", ascending=True)
    y = np.arange(len(plotted))
    ax.barh(y, plotted["_value"], color=color, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(plotted[country_field], fontsize=5.7)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.3)
    ax.set_axisbelow(True)
    return len(plotted)


def build_plot(world, resource, investor):
    """Keep the October 6 layout; include every positive owner-registration row."""
    owner = investor.loc[investor["foreign_project_attributed_area_km2"].gt(0)].copy()
    host = resource.loc[resource["foreign_attributed_area_km2"].gt(0)].copy()
    width = 183.0 / 25.4
    fig = plt.figure(figsize=(width, 8.45))
    ax_a = fig.add_axes([0.15, 0.680, 0.83, 0.83 * width / 2.4 / 8.45])
    ax_b = fig.add_axes([0.15, 0.055, 0.307, 0.550])
    ax_c = fig.add_axes([0.673, 0.055, 0.307, 0.550])
    _draw_resource_map(ax_a, world, resource)
    ax_a.set_aspect("equal", adjustable="box")
    _ranking_panel(
        ax_b, owner, "investor_country", "foreign_project_attributed_area_km2",
        "Equity-attributed affected habitat area outside\n"
        "the owner-registration country or region (km²)", "#4C78A8",
    )
    _ranking_panel(
        ax_c, host, "resource_country", "foreign_attributed_area_km2",
        "Affected habitat area attributed to project\n"
        "equity holders registered elsewhere (km²)", "#B44C63",
    )
    for ax in (ax_b, ax_c):
        ax.tick_params(axis="both", labelsize=5.7, pad=2)
    for ax, letter in zip((ax_a, ax_b, ax_c), "abc"):
        ax.text(-0.10 if ax is ax_a else -0.38, 1.02, letter,
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                fontweight="bold", clip_on=False, zorder=100)
    return fig
