"""Supplementary figure drawing functions; inputs are supplied by the runner."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd




DOUBLE_COLUMN_WIDTH_IN = 183.0 / 25.4

def panel_label(ax, label, x=-0.08, y=1.04):
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold", clip_on=False, zorder=100)

def build_plot(country, iterations, summary):
    country = country[pd.to_numeric(country["mine_count"], errors="coerce").ge(10)].copy()
    country = country.sort_values("top10_pct_difference_observed_minus_background")
    if len(country) > 30:
        country = country.iloc[-30:].copy()
    y = np.arange(len(country))

    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_COLUMN_WIDTH_IN, 6.2),
        gridspec_kw={"width_ratios": [1.3, 1.0], "wspace": 0.38},
    )
    background = country["background_top10_grid_pct"].to_numpy(float)
    observed = country["observed_top10_mine_pct"].to_numpy(float)
    for yy, x0, x1 in zip(y, background, observed):
        ax_a.plot([x0, x1], [yy, yy], color="#bdbdbd", linewidth=0.75, zorder=1)
    ax_a.scatter(background, y, color="#4C78A8", s=18, label="Available grid cells", zorder=2)
    ax_a.scatter(observed, y, color="#B44C63", s=18, label="Copper mines", zorder=3)
    ax_a.axvline(10, color="#777777", linestyle=":", linewidth=0.7)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(country["Country"], fontsize=5.2)
    ax_a.set_xlabel("Locations in the global\ntop-decile rarity-weighted richness class (%)")
    ax_a.set_ylabel("Mine-host country")
    ax_a.legend(loc="lower right")
    ax_a.grid(axis="x", color="#e5e5e5", linewidth=0.3)

    metric = summary.set_index("metric").loc["mean_rwr_percentile"]
    null_values = pd.to_numeric(iterations["mean_rwr_percentile"], errors="raise")
    observed_value = float(metric["observed"])
    ax_b.hist(null_values, bins=35, color="#b7d4ea", edgecolor="white", linewidth=0.45)
    ax_b.axvspan(float(metric["null_q025"]), float(metric["null_q975"]), color="#777777", alpha=0.14, label="95% null interval")
    ax_b.axvline(observed_value, color="#B44C63", linewidth=1.6, label="Observed")
    ax_b.set_xlabel("Mean rarity-weighted richness percentile")
    ax_b.set_ylabel("Randomization iterations (n)")
    ax_b.text(
        0.97,
        0.96,
        f"Observed = {observed_value:.1f}\nSES = {float(metric['ses']):.2f}\nP = {float(metric['p_upper_one_sided']):.4f}",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
    )
    ax_b.legend(loc="upper left")
    panel_label(ax_a, "a", x=-0.23, y=1.02)
    panel_label(ax_b, "b", x=-0.12, y=1.02)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.96, bottom=0.10)
    return fig
