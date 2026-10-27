"""Supplementary figure drawing functions; inputs are supplied by the runner."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from matplotlib.ticker import FixedLocator, NullLocator

DOUBLE_COLUMN_WIDTH_IN = 183.0 / 25.4

def _conflict_limit(area):
    area = np.asarray(area, dtype=float)
    limit = np.empty_like(area)
    narrow = area <= 200
    broad = area >= 250000
    middle = ~(narrow | broad)
    limit[narrow] = 0
    limit[broad] = 70
    position = (np.log10(area[middle]) - np.log10(200)) / (np.log10(250000) - np.log10(200))
    limit[middle] = 70 * position
    return limit

def build_plot(data):
    data["total_habitat_area_km2"] = pd.to_numeric(data["total_habitat_area_km2"], errors="coerce")
    data["Impact_Percent"] = pd.to_numeric(data["Impact_Percent"], errors="coerce")
    data["Priority_ForPlot"] = data["Priority_ForPlot"].astype(str).str.lower().isin(["true", "1", "yes"])
    data = data[data["total_habitat_area_km2"].gt(0) & data["Impact_Percent"].notna()].copy()
    x_min = 10 ** np.floor(np.log10(max(1, data["total_habitat_area_km2"].min() * 0.7)))
    x_max = 10 ** np.ceil(np.log10(max(500000, data["total_habitat_area_km2"].max() * 1.3)))
    curve_x = np.logspace(np.log10(x_min), np.log10(x_max), 500)

    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_WIDTH_IN, 4.9))
    ax.axvspan(x_min, 200, color="#f4d5bc", alpha=0.68, lw=0)
    ax.axvspan(200, 250000, color="#d4eee7", alpha=0.68, lw=0)
    ax.axvspan(250000, x_max, color="#d3ecfa", alpha=0.68, lw=0)
    ordinary = data[~data["Priority_ForPlot"]]
    selected = data[data["Priority_ForPlot"]]
    ax.scatter(ordinary["total_habitat_area_km2"], ordinary["Impact_Percent"], s=12, color="#777777", alpha=0.35, linewidth=0, label="Other affected species")
    ax.scatter(selected["total_habitat_area_km2"], selected["Impact_Percent"], s=18, color="#B44C63", alpha=0.82, linewidth=0, label="Conservation-priority species")
    ax.plot(curve_x, _conflict_limit(curve_x), color="#111111", linewidth=1.4, label="Conservation-priority threshold")
    ax.set_xscale("log")
    ax.tick_params(axis="x", labelsize=7.2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-2, 102)
    ax.set_xlabel("Global habitat area (km²)")
    ax.set_ylabel("Proportion of habitat affected (%)")
    ax.text(0.08, 0.96, "Narrow-habitat", transform=ax.transAxes, ha="center", va="top", fontweight="bold")
    ax.text(0.48, 0.96, "Intermediate-habitat", transform=ax.transAxes, ha="center", va="top", fontweight="bold")
    ax.text(0.88, 0.96, "Widespread", transform=ax.transAxes, ha="center", va="top", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), ncol=2, fontsize=5.5)
    ax.grid(color="#e2e2e2", linewidth=0.3)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.96, bottom=0.24)
    ax.xaxis.set_major_locator(FixedLocator([1, 1e2, 1e4, 1e6, 1e8]))
    ax.xaxis.set_minor_locator(NullLocator())
    return fig
