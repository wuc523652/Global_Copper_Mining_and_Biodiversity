# -*- coding: utf-8 -*-
"""Spatial-pattern and randomization panels for Figures 1 and 2."""

from __future__ import annotations

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd

from common import INPUTS, table, panel_label, require_columns, fix_geographic_aspect


# Panel-b inset position in parent-axis coordinates:
# [left, bottom, width, height].
TAXON_INSET_BOUNDS = [0.100, 0.245, 0.145, 0.245]


def mine_land_use(mines):
    area = pd.to_numeric(mines.land_use_km2, errors='raise').to_numpy(float)
    classes = np.select([area <= 1, area <= 10], [0, 1], default=2)
    if not np.isfinite(area).all() or (area < 0).any():
        raise ValueError('Mine land-use areas must be finite and nonnegative')
    return area, classes


def _geo_axis(ax, lon_formatter, lat_formatter) -> None:
    fix_geographic_aspect(ax, (-180, 180, -60, 90), global_map=True)
    ax.set_xticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_yticks([-60, -30, 0, 30, 60, 90])
    ax.xaxis.set_major_formatter(FuncFormatter(lon_formatter))
    ax.yaxis.set_major_formatter(FuncFormatter(lat_formatter))
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, pad=1)
    ax.grid(color="#DEDEDE", linewidth=0.3, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_color("#777777")


def draw_richness_map(ax, inputs):
    rwr, boundary, antarctica, mines = inputs.load_plot_data()
    inputs.safe_gdf_plot(boundary, ax=ax, facecolor="#F2F2F2", edgecolor="none", zorder=0)
    inputs.safe_gdf_plot(antarctica, ax=ax, facecolor="#F2F2F2", edgecolor="none", zorder=0)
    coordinates = rwr[['longitude', 'latitude']].to_numpy(float)
    if not np.isfinite(coordinates).all():
        raise ValueError('RWR display coordinates must be finite')
    for order, code in enumerate(inputs.BIN_ORDER, 1):
        mask = rwr[inputs.CLASS_FIELD].eq(code).to_numpy()
        if mask.any():
            ax.scatter(
                coordinates[mask, 0],
                coordinates[mask, 1],
                s=0.075,
                marker="s",
                color=inputs.BIN_COLOR[code],
                edgecolors="none",
                linewidths=0,
                rasterized=True,
                zorder=order,
            )
    inputs.safe_gdf_plot(boundary.boundary, ax=ax, linewidth=0.18, edgecolor="#8C8C8C", zorder=7)
    land_use_km2, mine_classes = mine_land_use(mines)
    marker_areas = np.array([1.2, 3.0, 6.0])
    inputs.safe_gdf_plot(
        mines,
        ax=ax,
        marker="^",
        markersize=marker_areas[mine_classes],
        facecolor="#111111",
        edgecolor="none",
        alpha=0.72,
        zorder=9,
    )
    _geo_axis(ax, inputs.lon_formatter, inputs.lat_formatter)
    percentile_labels = ["0–50th", "50–75th", "75–90th", "90–99th", "99–100th"]
    handles = [
        Patch(facecolor=inputs.BIN_COLOR[code], edgecolor="none", label=label)
        for code, label in zip(inputs.BIN_ORDER, percentile_labels)
    ]
    richness_legend = ax.legend(
        handles=handles,
        title="Rarity-weighted richness percentile",
        loc="lower center",
        bbox_to_anchor=(0.61, 0.015),
        ncol=5,
        columnspacing=0.7,
        handlelength=1.1,
        title_fontsize=5.7,
        fontsize=5.1,
    )
    ax.add_artist(richness_legend)
    mine_handles = [Line2D([0], [0], marker="^", linestyle="none", markersize=np.sqrt(size),
                           markerfacecolor="#111111", markeredgecolor="none", alpha=0.72, label=label)
                    for size, label in zip(marker_areas, ["0–1", "1–10", ">10"])]
    ax.legend(handles=mine_handles, title="Copper mine\nLand use (km²)",
              loc="lower left", bbox_to_anchor=(0.04, 0.17), ncol=1,
              fontsize=5.5, title_fontsize=5.5, handlelength=1.1, labelspacing=0.5)
    panel_label(ax, "A", x=-0.06, y=1.11)
    return {"rwr_grid_count": len(rwr), "mine_count": len(mines),
            "drawn_rwr_grid_count": len(coordinates)}


def draw_taxon_inset(ax, by_taxon: pd.DataFrame):
    """Add the overlapping-species-by-taxon bar chart to panel b."""
    require_columns(by_taxon, ["taxon", "num_species"], "Figure 1B taxon data")
    taxon = by_taxon[["taxon", "num_species"]].copy()
    taxon["num_species"] = pd.to_numeric(taxon["num_species"], errors="raise")

    inset = ax.inset_axes(TAXON_INSET_BOUNDS, zorder=20)
    x = np.arange(len(taxon))
    inset.bar(
        x,
        taxon["num_species"].to_numpy(float),
        width=0.68,
        color="#3B7FB1",
        edgecolor="none",
        zorder=2,
    )
    inset.set_xticks(x)
    inset.set_xticklabels(
        taxon["taxon"].replace({"Fish": "Freshwater\nfishes"}),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
        fontsize=5.0,
    )
    inset.tick_params(axis="x", length=1.5, width=0.45, pad=0.8)
    inset.tick_params(axis="y", labelsize=5.0, length=1.5, width=0.45, pad=0.8)
    # The five rotated group names already identify the categorical axis.
    # Omitting the redundant x-axis title also keeps the inset completely
    # separate from the range-overlap legend below it.
    inset.set_xlabel("")
    inset.set_ylabel(
        "Affected vertebrate\nspecies (n)",
        fontsize=5.0,
        labelpad=1.2,
    )
    inset.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    inset.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ymax = max(float(taxon["num_species"].max()), 1.0)
    inset.set_ylim(0, ymax * 1.10)
    inset.set_facecolor("white")
    inset.grid(False)
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color("#444444")
        spine.set_linewidth(0.50)
    return inset


def draw_overlap_map(ax, inputs):
    map_gdf, point_gdf, by_taxon = inputs.load_plot_data()
    map_gdf = map_gdf.copy()
    no_overlap = pd.to_numeric(map_gdf[inputs.AREA_FIELD], errors="coerce").fillna(0).le(0)
    map_gdf.loc[no_overlap, "area_color"] = "#F1F1F1"
    map_gdf.plot(
        ax=ax,
        color=map_gdf["area_color"],
        edgecolor="#BDBDBD",
        linewidth=0.16,
        zorder=1,
    )
    labels = [item[2] for item in inputs.SPECIES_COUNT_BINS]
    for label in labels:
        subset = point_gdf[point_gdf["species_count_class"] == label]
        if not subset.empty:
            subset.plot(
                ax=ax,
                color=inputs.get_species_count_color(label),
                markersize=7,
                edgecolor="#222222",
                linewidth=0.25,
                marker="o",
                zorder=5,
            )
    _geo_axis(ax, inputs.lon_formatter, inputs.lat_formatter)
    draw_taxon_inset(ax, by_taxon)

    area_labels = ["0–5,000", "5,000–10,000", "10,000–50,000", "50,000–100,000", ">100,000"]
    area_handles = [
        Patch(facecolor=color, edgecolor="#AAAAAA", linewidth=0.25, label=label)
        for (_, _, _, color), label in zip(inputs.AREA_BINS[1:], area_labels)
    ]
    area_legend = ax.legend(
        handles=area_handles,
        title="Affected habitat area (km²)",
        loc="lower center",
        bbox_to_anchor=(0.65, 0.081),
        ncol=5,
        fontsize=5.0,
        title_fontsize=5.0,
        columnspacing=0.6,
        handlelength=1.0,
        handletextpad=0.3,
        borderaxespad=0,
        borderpad=0,
        labelspacing=0.1,
    )
    ax.add_artist(area_legend)
    species_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markersize=3.8,
            markerfacecolor=color, markeredgecolor="#222222", markeredgewidth=0.25,
            label=label,
        )
        for _, _, label, color in inputs.SPECIES_COUNT_BINS
    ]
    ax.legend(
        handles=species_handles,
        title="Affected vertebrate species (n)",
        loc="lower center",
        bbox_to_anchor=(0.65, 0.002),
        ncol=5,
        fontsize=5.0,
        title_fontsize=5.0,
        columnspacing=0.5,
        handletextpad=0.3,
        borderaxespad=0,
        borderpad=0,
        labelspacing=0.1,
    )
    panel_label(ax, "B", x=-0.06, y=1.11)
    return {
        "country_count": int((map_gdf[inputs.AREA_FIELD] > 0).sum()),
        "species_total": int(by_taxon["num_species"].sum()),
    }


def load_country_area_data(inputs) -> pd.DataFrame:
    data = inputs.country_area_data()
    required = [
        inputs.COUNTRY_FIELD,
        inputs.OVERLAP_AREA_FIELD,
        inputs.TOTAL_AREA_FIELD,
        inputs.PCT_TOP10_FIELD,
        *inputs.AREA_COLS,
    ]
    require_columns(data, required, "Figure 1C source data")
    for field in required[1:]:
        data[field] = pd.to_numeric(data[field], errors="raise")
    data = data[
        data[inputs.COUNTRY_FIELD].astype(str).str.casefold().ne("global")
        & data[inputs.OVERLAP_AREA_FIELD].ge(inputs.MIN_OVERLAP_AREA_KM2)
    ].copy()
    sort_fields = [inputs.OVERLAP_AREA_FIELD]
    ascending = [False]
    if "num_species" in data.columns:
        sort_fields.append("num_species")
        ascending.append(False)
    sort_fields.append(inputs.COUNTRY_FIELD)
    ascending.append(True)
    return data.sort_values(sort_fields, ascending=ascending).reset_index(drop=True)


def draw_country_area(ax, inputs):
    data = load_country_area_data(inputs)
    x = np.arange(len(data))
    bottom = np.zeros(len(data))
    labels = ["0–50th", "50–75th", "75–90th", "90–99th", "99–100th"]
    for field, color, label in zip(inputs.AREA_COLS, inputs.BIN_COLORS, labels):
        values = data[field].to_numpy(float)
        ax.bar(x, values, bottom=bottom, width=0.75, color=color, edgecolor="white",
               linewidth=0.25, label=label, zorder=2)
        bottom += values
    ymax = float(data[inputs.TOTAL_AREA_FIELD].max())
    inputs.set_asinh_area_scale(ax, ymax)
    # Reserve an empty internal header for the legend, above all observations.
    transform = ax.yaxis.get_transform()
    upper = float(transform.inverted().transform(transform.transform(ymax) / 0.78))
    ax.set_ylim(0, upper)
    ax.yaxis.set_major_formatter(FuncFormatter(inputs.format_area))
    ax.set_ylabel("Mine potential influence area\nin richness classes (km²; asinh scale)")
    ax.set_xticks(x)
    ax.set_xticklabels(
        data[inputs.COUNTRY_FIELD].replace({"Democratic Republic of the Congo": "DR Congo"}),
        rotation=65,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.0,
    )
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.35, linestyle="--", zorder=0)
    ax.set_axisbelow(True)

    right = ax.twinx()
    values = data[inputs.PCT_TOP10_FIELD].to_numpy(float)
    right.plot(x, values, color="#3B7F5C", linewidth=1.0, linestyle="--",
               marker="o", markersize=2.6, zorder=5)
    right.set_ylim(0, 128)
    right.set_yticks([0, 25, 50, 75, 100])
    right.set_ylabel("Share of classified mine\npotential influence area in\nthe global top decile (%)",
                     color="#3B7F5C")
    right.tick_params(axis="y", colors="#3B7F5C")
    right.spines["right"].set_visible(True)
    right.spines["right"].set_color("#3B7F5C")
    handles = [Patch(facecolor=c, edgecolor="none", label=l) for c, l in zip(inputs.BIN_COLORS, labels)]
    handles.append(Line2D([0], [0], color="#3B7F5C", linestyle="--", marker="o",
                          markersize=3, label="Share in global top-decile class"))
    # The upper internal band is reserved for the legend; no data are masked.
    ax.legend(
        handles=handles,
        title="Rarity-weighted richness percentile",
        loc="upper left",
        bbox_to_anchor=(0.015, 0.985),
        ncol=6,
        fontsize=5.0,
        title_fontsize=5.4,
        columnspacing=0.75,
        handlelength=1.1,
        handletextpad=0.35,
        borderaxespad=0.0,
        frameon=False,
    )
    panel_label(ax, "C", x=-0.06, y=1.12)
    return {"displayed_country_count": int(len(data)), "internal_legend_header": True,
            "display_headroom_only": "all bars and percentage values unchanged; top 22 percent reserved for legend"}


def draw_randomization(ax, randomization):
    iterations = pd.read_csv(randomization.OUT_ITERATIONS)
    summary = pd.read_csv(randomization.OUT_SUMMARY)
    row = summary.set_index("metric").loc["pct_top10_mines"]
    null = iterations["pct_top10_mines"].to_numpy(float)
    observed = float(row["observed"])
    q025, q975 = float(row["null_q025"]), float(row["null_q975"])
    counts, _, _ = ax.hist(null, bins=30, color="#BFD3E6", edgecolor="white", linewidth=0.40)
    ax.set_ylim(0, counts.max() * 1.45)
    ax.axvspan(
        q025,
        q975,
        color="#6D8FA8",
        alpha=0.14,
        label="95% randomization\ninterval",
    )
    ax.axvline(
        observed,
        color="#CB1C2A",
        linewidth=1.6,
        label=f"Observed = {observed:.2f}%",
    )
    ax.set_xlabel(
        "Mine reference points in grids within the\n"
        "top decile of rarity-weighted richness (%)"
    )
    ax.set_ylabel("Randomization frequency")
    ax.grid(axis="y", color="#E0E0E0", linewidth=0.35, linestyle="--")
    ratio = observed / float(row["null_mean"])
    p_value = float(row["p_upper_one_sided"])
    ax.legend(loc="upper left", frameon=False, fontsize=5.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    panel_label(ax, "A", x=-0.12, y=1.04)
    return {
        "iterations": int(len(iterations)),
        "observed_percent": observed,
        "null_mean_percent": float(row["null_mean"]),
        "observed_expected_ratio": ratio,
        "p_upper_one_sided": p_value,
    }


def load_country_comparison_data(randomization) -> pd.DataFrame:
    """Load and select the country-level descriptive randomization results."""
    country = pd.read_csv(randomization.OUT_COUNTRY_DESC, encoding="utf-8-sig")
    required = [
        randomization.COUNTRY_FIELD,
        "mine_count",
        "observed_top10_mine_pct",
        "background_top10_grid_pct",
        "top10_pct_difference_observed_minus_background",
    ]
    require_columns(country, required, "Figure 2B country randomization data")
    for field in required[1:]:
        country[field] = pd.to_numeric(country[field], errors="raise")

    min_mines = int(getattr(randomization, "COUNTRY_MIN_MINES_FOR_FIG", 10))
    selected = country[country["mine_count"].ge(min_mines)].copy()
    return selected.sort_values(
        "top10_pct_difference_observed_minus_background",
        ascending=False,
    ).reset_index(drop=True)


def draw_country_comparison(ax, randomization):
    """Country land-grid background versus observed mine locations."""
    selected = load_country_comparison_data(randomization)
    if selected.empty:
        min_mines = int(getattr(randomization, "COUNTRY_MIN_MINES_FOR_FIG", 10))
        ax.text(
            0.5,
            0.5,
            f"No country has at least {min_mines} mines",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        panel_label(ax, "B", x=-0.07, y=1.04)
        return {"displayed_country_count": 0, "minimum_mine_count": min_mines}

    x = np.arange(len(selected))
    observed = selected["observed_top10_mine_pct"].to_numpy(float)
    background = selected["background_top10_grid_pct"].to_numpy(float)
    for xpos, obs, bg in zip(x, observed, background):
        ax.plot([xpos, xpos], [bg, obs], color="#B4B4B4", linewidth=0.65, zorder=1)
    ax.scatter(
        x,
        background,
        s=10,
        color="#667788",
        edgecolor="white",
        linewidth=0.30,
        label="Country land-grid background",
        zorder=3,
    )
    ax.scatter(
        x,
        observed,
        s=12,
        marker="D",
        color="#CB1C2A",
        edgecolor="white",
        linewidth=0.30,
        label="Observed mines",
        zorder=4,
    )
    tick_labels = [
        f"{'DR Congo' if name == 'Democratic Republic of the Congo' else name}\n(n={int(count)})"
        for name, count in zip(selected[randomization.COUNTRY_FIELD], selected["mine_count"])
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(
        tick_labels,
        rotation=60,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.0,
    )
    ax.set_ylabel(
        "Locations in grids within the top decile\n"
        "of rarity-weighted richness (%)"
    )
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.grid(axis="y", linestyle="--", linewidth=0.35, alpha=0.45)
    ax.legend(frameon=False, fontsize=5.0, ncol=2, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    panel_label(ax, "B", x=-0.10, y=1.08)
    return {
        "displayed_country_count": int(len(selected)),
        "minimum_mine_count": int(getattr(randomization, "COUNTRY_MIN_MINES_FOR_FIG", 10)),
    }
