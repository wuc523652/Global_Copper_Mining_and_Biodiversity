"""Bivariate map and equity-weighted ownership portfolio panels."""
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from common import require_columns, fix_geographic_aspect
from plot_inputs import map_inputs as map_backend

CATEGORY_SPECS = [
    ("neither_portfolio_share_percent", "Neither criterion", "#d9d9d9"),
    (
        "high_rwr_only_portfolio_share_percent",
        "Regional criterion only",
        "#4C78A8",
    ),
    (
        "priority_screened_species_only_portfolio_share_percent",
        "Species-level criterion only",
        "#F2B134",
    ),
    (
        "dual_ecological_attention_portfolio_share_percent",
        "Both criteria",
        "#B44C63",
    ),
]


# The existing ownership map is retained; alluvial panels have their own module.
MAP_LEGEND_SIDE_MM = 15.0


def draw_bivariate_map(ax, world, investor):
    required = [
        "investor_iso3",
        "total_attributed_area_km2",
        "foreign_project_share_pct",
    ]
    require_columns(investor, required, "owner-registration-country summary")
    data = investor.copy()
    data["investor_iso3"] = data["investor_iso3"].astype(str).str.upper().str.strip()
    data["total_attributed_area_km2"] = pd.to_numeric(
        data["total_attributed_area_km2"], errors="coerce"
    )
    data["foreign_project_share_pct"] = pd.to_numeric(
        data["foreign_project_share_pct"], errors="coerce"
    ).clip(0, 100)

    area_breaks = map_backend().make_area_breaks(
        data["total_attributed_area_km2"], map_backend().MAIN_FIXED_AREA_BREAKS_KM2
    )
    percent_breaks = np.asarray(map_backend().PERCENT_BREAKS, dtype=float)
    mapped = world.merge(data, left_on="iso3", right_on="investor_iso3", how="left")
    mapped["area_class"] = mapped["total_attributed_area_km2"].apply(
        lambda value: map_backend().classify_three(value, area_breaks)
    )
    mapped["percent_class"] = mapped["foreign_project_share_pct"].apply(
        lambda value: map_backend().classify_three(value, percent_breaks)
    )
    mapped["map_color"] = "#e6e6e6"
    valid = (
        mapped["area_class"].ge(0)
        & mapped["percent_class"].ge(0)
        & mapped["total_attributed_area_km2"].gt(0)
    )
    mapped.loc[valid, "map_color"] = [
        map_backend().MAIN_PALETTE[iy, ix]
        for iy, ix in zip(
            mapped.loc[valid, "area_class"].astype(int),
            mapped.loc[valid, "percent_class"].astype(int),
        )
    ]
    antarctica = mapped["iso3"].eq("ATA") | mapped["map_country"].astype(str).str.contains(
        "antarct", case=False, na=False
    )
    mapped.loc[antarctica, "map_color"] = "#f4f4f4"

    simplified = mapped.copy()
    simplified.plot(
        ax=ax,
        color=simplified["map_color"],
        edgecolor="#707070",
        linewidth=0.14,
        zorder=2,
    )
    fix_geographic_aspect(ax, (-180, 180, -60, 90), global_map=True)
    ax.set_xticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_yticks([-60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(["150°W", "100°W", "50°W", "0°", "50°E", "100°E", "150°E"])
    ax.set_yticklabels(["60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, pad=1)
    ax.tick_params(axis="y", pad=1.2)
    ax.grid(color="#dedede", linewidth=0.35, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#b7b7b7")
        spine.set_linewidth(0.45)
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax.apply_aspect()
    position = ax.get_position()
    map_width_mm = position.width * ax.figure.get_figwidth() * 25.4
    map_height_mm = position.height * ax.figure.get_figheight() * 25.4
    inset = ax.inset_axes([14 / map_width_mm, 12 / map_height_mm,
                          MAP_LEGEND_SIDE_MM / map_width_mm, MAP_LEGEND_SIDE_MM / map_height_mm])
    for iy in range(3):
        for ix in range(3):
            inset.add_patch(
                Rectangle(
                    (ix, iy),
                    1,
                    1,
                    facecolor=map_backend().MAIN_PALETTE[iy, ix],
                    edgecolor="white",
                    linewidth=0.45,
                )
            )
    inset.set_xlim(0, 3)
    inset.set_ylim(0, 3)
    inset.set_aspect("equal")
    inset.set_xticks([0.5, 1.5, 2.5])
    inset.set_xticklabels(["<20", "20–50", ">50"], fontsize=5.5)
    inset.set_yticks([0.5, 1.5, 2.5])
    inset.set_yticklabels(["<10³", "10³–10⁴", ">10⁴"], fontsize=5.5)
    inset.tick_params(length=0, pad=0.7)
    for spine in inset.spines.values():
        spine.set_visible(False)
    inset.set_ylabel("Traceable affected\nhabitat area (km²)", fontsize=5.5, rotation=90, rotation_mode="anchor", labelpad=1.0)
    inset.set_xlabel(
        "Share linked to mines outside\nthe owner-registration country (%)",
        fontsize=5.5,
        labelpad=0.7,
    )
    return area_breaks, percent_breaks, int(valid.sum())



def _selection_codes(selection, group):
    selected = selection.loc[selection["selection_group"].eq(group)].copy()
    selected["selection_rank"] = pd.to_numeric(selected["selection_rank"], errors="coerce")
    return selected.sort_values("selection_rank")["country_iso3"].astype(str).tolist()


def draw_portfolios(ax, portfolio, selection):
    order = _selection_codes(selection, "Fig5a investor portfolio")
    data = portfolio[portfolio["investor_iso3"].isin(order)].copy()
    rank = {code: idx for idx, code in enumerate(order)}
    data["_rank"] = data["investor_iso3"].map(rank)
    data = data.sort_values("_rank", ascending=False)
    y = np.arange(len(data))
    left = np.zeros(len(data), dtype=float)
    for field, label, color in CATEGORY_SPECS:
        values = pd.to_numeric(data[field], errors="coerce").fillna(0).to_numpy()
        ax.barh(y, values, left=left, color=color, height=0.68, label=label, edgecolor="white", linewidth=0.25)
        left += values
    mine_counts = pd.to_numeric(data["num_unique_mines"], errors="coerce").fillna(0).astype(int)
    ax.scatter(np.full(len(data), 108.5), y, s=np.clip(np.sqrt(mine_counts) * 3.0, 10, 40), color="#313131", zorder=5)
    for yy, count in zip(y, mine_counts):
        ax.text(112.0, yy, str(count), va="center", ha="left", fontsize=5.0)
    ax.text(110.0, len(data) - 0.12, "Mines (n)", ha="center", va="bottom", fontsize=5.2, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(data["investor_country"].replace({"Democratic Republic of the Congo": "DR Congo"}))
    ax.set_xlim(0, 124)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Share of portfolio by owner-registration country\n(equity-weighted mine equivalents; %)")
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.01),
        ncol=2,
        fontsize=5.0,
        handlelength=1.2,
        columnspacing=0.8,
    )
    return len(data)

