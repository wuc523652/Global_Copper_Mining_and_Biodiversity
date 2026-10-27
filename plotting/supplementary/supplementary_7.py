"""Rendering functions adapted from the existing supplementary figure sources.
All quantitative inputs are supplied by render_supplementary.py; no analysis runs.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build_plot(
    scenario: pd.DataFrame,
    investor_comparison: pd.DataFrame,
    bilateral_comparison: pd.DataFrame,
):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)

    ax = axes[0, 0]
    labels = ["All mines", "Complete cases"]
    mine_values = scenario["num_mines"].to_numpy()
    area_values = scenario["global_strict_physical_exposure_area_km2"].to_numpy() / 1e6
    x = np.arange(2)
    width = 0.36
    ax.bar(x - width / 2, mine_values / mine_values[0] * 100, width, label="Mine count coverage")
    ax.bar(x + width / 2, area_values / area_values[0] * 100, width, label="Affected habitat coverage")
    ax.set(xticks=x, xticklabels=labels, ylabel="Coverage relative to all mines (%)", ylim=(0, 110), title="a  Complete-case coverage")
    ax.legend(frameon=True, facecolor="white", framealpha=0.92, edgecolor="0.85", loc="upper right")

    ax = axes[0, 1]
    domestic = scenario["domestic_equity_associated_area_km2"].to_numpy()
    foreign = scenario["foreign_equity_associated_area_km2"].to_numpy()
    unattributed = scenario["unattributed_area_km2"].to_numpy()
    total = domestic + foreign + unattributed
    bottoms = np.zeros(2)
    for values, label, color in [
        (domestic / total * 100, "Domestic", "#59A14F"),
        (foreign / total * 100, "Cross-border", "#E15759"),
        (unattributed / total * 100, "Unattributed", "#BAB0AC"),
    ]:
        ax.bar(x, values, bottom=bottoms, label=label, color=color)
        bottoms += values
    ax.set(xticks=x, xticklabels=labels, ylabel="Share of affected habitat area (%)", ylim=(0, 100), title="b  Ownership composition")
    ax.legend(frameon=True, facecolor="white", framealpha=0.92, edgecolor="0.85", loc="upper right")

    ax = axes[1, 0]
    shown = investor_comparison.nsmallest(20, "main_rank")
    ax.scatter(shown["main_rank"], shown["complete_case_rank"], color="#4C78A8")
    limit = max(20, float(shown[["main_rank", "complete_case_rank"]].max().max())) + 1
    ax.plot([1, limit], [1, limit], linestyle="--", color="0.6")
    for row in shown.nsmallest(8, "main_rank").itertuples(index=False):
        ax.annotate(row.investor_iso3, (row.main_rank, row.complete_case_rank), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set(xlabel="All-mine rank", ylabel="Complete-case rank", xlim=(0, limit), ylim=(0, limit), title="c  Ranking of owner-registration\ncountries and regions")

    ax = axes[1, 1]
    shown_flows = bilateral_comparison.nsmallest(20, "main_rank")
    ax.scatter(shown_flows["main_rank"], shown_flows["complete_case_rank"], color="#F28E2B")
    limit = max(20, float(shown_flows[["main_rank", "complete_case_rank"]].max().max())) + 1
    ax.plot([1, limit], [1, limit], linestyle="--", color="0.6")
    for row in shown_flows.nsmallest(6, "main_rank").itertuples(index=False):
        ax.annotate(row.flow_key, (row.main_rank, row.complete_case_rank), xytext=(3, 3), textcoords="offset points", fontsize=6)
    ax.set(xlabel="All-mine rank", ylabel="Complete-case rank", xlim=(0, limit), ylim=(0, limit), title="d  Cross-border responsibility ranking")

    return fig
