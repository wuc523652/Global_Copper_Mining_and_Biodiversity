# -*- coding: utf-8 -*-
"""The six range-size composition panels of Figure 3."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from common import panel_label


GROUP_ORDER = ["All", "Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
RANGE_ORDER = ["Narrow", "Intermediate", "Widespread"]
HABITAT_CLASS_LABELS = {
    "Narrow": "Narrow-habitat",
    "Intermediate": "Intermediate-habitat",
    "Widespread": "Widespread",
}
RANGE_COLORS = {
    "Narrow": "#E28740",
    "Intermediate": "#84CFBB",
    "Widespread": "#7FC7EF",
}
PRIORITIZED_COLOR = "#C95A11"

# Shared bar geometry for all six composition panels.
BAR_WIDTH = 0.62
BAR_OUTLINE_WIDTH = 0.65
UNAFFECTED_HATCH = "////"
UNAFFECTED_FACE_COLOR = "white"
STATUS_LEGEND_COLOR = "#777777"

# Initial shared-legend placement; the spacing pass places it below the panels.
PANEL_A_LEGEND_Y = 1.23


def _facet_strip(ax, title: str) -> None:
    ax.add_patch(
        Rectangle(
            (0, 1.0), 1, 0.12, transform=ax.transAxes, clip_on=False,
            facecolor="#DEDEDE", edgecolor="none", zorder=10,
        )
    )
    ax.text(0.5, 1.06, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=6.1, fontweight="bold", zorder=11)


def _nice_y_limit(max_value: float) -> float:
    """Return a compact, publication-friendly upper y limit."""
    if max_value <= 0:
        return 5.0
    raw = max_value * 1.24
    magnitude = 10 ** np.floor(np.log10(max(raw, 1.0)))
    scaled = raw / magnitude
    for step in [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]:
        if scaled <= step:
            return float(step * magnitude)
    return float(np.ceil(raw / magnitude) * magnitude)


def _draw_composition_facet(ax, data: pd.DataFrame, group: str) -> None:
    sub = data[data["Group"].eq(group)].set_index("Range").reindex(RANGE_ORDER)
    x = np.arange(3)
    bar_left = x - BAR_WIDTH / 2.0
    no_overlap = sub["n_unaffected_species"].fillna(0).to_numpy(float)
    overlap_not_prioritized = sub["n_affected_non_priority_species"].fillna(0).to_numpy(float)
    prioritized = sub["n_priority_species"].fillna(0).to_numpy(float)
    totals = no_overlap + overlap_not_prioritized + prioritized

    # Stack measured counts directly on the same numerical axis.
    affected_bottom = no_overlap
    priority_bottom = no_overlap + overlap_not_prioritized
    visual_tops = totals

    # Every segment uses the same explicit left edge and width.  This keeps the
    # hatched, solid and prioritized outlines perfectly aligned.
    for left, range_class, height in zip(bar_left, RANGE_ORDER, no_overlap):
        ax.bar(
            left,
            height,
            width=BAR_WIDTH,
            align="edge",
            color=UNAFFECTED_FACE_COLOR,
            edgecolor=RANGE_COLORS[range_class],
            linewidth=BAR_OUTLINE_WIDTH,
            hatch=UNAFFECTED_HATCH,
            zorder=2,
        )

    for left, range_class, height, bottom in zip(
        bar_left, RANGE_ORDER, overlap_not_prioritized, affected_bottom
    ):
        ax.bar(
            left,
            height,
            bottom=bottom,
            width=BAR_WIDTH,
            align="edge",
            color=RANGE_COLORS[range_class],
            edgecolor=RANGE_COLORS[range_class],
            linewidth=BAR_OUTLINE_WIDTH,
            zorder=3,
        )

    # Completely skip zero-height prioritized rectangles; otherwise Matplotlib
    # can leave an orange baseline that falsely suggests a non-zero value.
    for left, height, bottom in zip(bar_left, prioritized, priority_bottom):
        if height <= 0:
            continue
        ax.bar(
            left,
            height,
            bottom=bottom,
            width=BAR_WIDTH,
            align="edge",
            color=PRIORITIZED_COLOR,
            edgecolor=PRIORITIZED_COLOR,
            linewidth=BAR_OUTLINE_WIDTH,
            zorder=4,
        )

    ymax = _nice_y_limit(float(visual_tops.max()) if len(visual_tops) else 0.0)
    ax.set_ylim(0, ymax)
    for xpos, total, visual_top, n_prioritized in zip(
        x, totals, visual_tops, prioritized
    ):
        label = f"{n_prioritized / total * 100:.0f}%" if total else "0%"
        ax.text(
            xpos,
            visual_top + ymax * 0.025,
            label,
            ha="center",
            va="bottom",
            fontsize=5.2,
            clip_on=False,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [HABITAT_CLASS_LABELS[name].replace("-habitat", "-\nhabitat") for name in RANGE_ORDER],
        fontsize=5.2,
    )
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.grid(axis="y", color="#D7D7D7", linestyle="--", linewidth=0.3)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.55)
    ax.spines["bottom"].set_linewidth(0.55)
    _facet_strip(ax, "Freshwater fishes" if group == "Fish" else group)


def draw_composition(fig, slot, stats: pd.DataFrame):
    inner = slot.subgridspec(2, 3, wspace=0.28, hspace=0.48)
    axes = []
    for index, group in enumerate(GROUP_ORDER):
        ax = fig.add_subplot(inner[index // 3, index % 3])
        _draw_composition_facet(ax, stats, group)
        axes.append(ax)
    axes[0].set_ylabel("Threatened vertebrate species (n)")
    axes[3].set_ylabel("Threatened vertebrate species (n)")
    range_handles = [
        Patch(facecolor=RANGE_COLORS[name], edgecolor="none", label=HABITAT_CLASS_LABELS[name])
        for name in RANGE_ORDER
    ]
    status_handles = [
        Patch(
            facecolor=UNAFFECTED_FACE_COLOR,
            edgecolor=STATUS_LEGEND_COLOR,
            hatch=UNAFFECTED_HATCH,
            label="No habitat overlap",
        ),
        Patch(
            facecolor=STATUS_LEGEND_COLOR,
            edgecolor="none",
            label="Affected species, not prioritized",
        ),
        Patch(facecolor=PRIORITIZED_COLOR, edgecolor="none",
              label="Conservation-priority species"),
    ]
    # Matplotlib fills a multi-row legend down columns.  Interleaving the two
    # categories therefore places range classes on row 1 and statuses on row 2.
    handles = [handle for pair in zip(range_handles, status_handles) for handle in pair]
    axes[1].legend(
        handles=handles,
        loc="lower center",
            bbox_to_anchor=(0.5, PANEL_A_LEGEND_Y),
        ncol=3,
        fontsize=5.1,
        columnspacing=0.8,
        handlelength=1.1,
    )
    for ax, label in zip(axes, "ABCDEF"):
        panel_label(ax, label, x=-0.20, y=1.15)
    return axes



import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from common import panel_label, fix_geographic_aspect
import geopandas as gpd
import pandas as pd

TAXON_ORDER = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
TAXON_COLORS = {"Amphibians": "#238B45", "Birds": "#9C4E97", "Fish": "#2878B5", "Mammals": "#D95F0E", "Reptiles": "#9A7A38"}
GLOBAL_CONTENT_EXTENT = (-180.0, 180.0, -60.0, 90.0)
SCALEBAR_CANDIDATES_KM = np.array([10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000], dtype=float)

def _axes_box_ratio(ax) -> float:
    """Width/height of the allocated axes box in physical figure units."""
    box = ax.get_position(original=True)
    return float(
        (box.width * ax.figure.get_figwidth())
        / (box.height * ax.figure.get_figheight())
    )


def _fit_geographic_extent(ax, extent) -> tuple[float, float, float, float]:
    """Pad an extent to the axes box while retaining equal ground-distance scale.

    Longitude degrees are converted at the extent's central latitude.  Only
    empty geographic margin is added; the requested source extent is never
    cropped and x/y are never independently stretched.
    """
    xmin, xmax, ymin, ymax = map(float, extent)
    x_span = xmax - xmin
    y_span = ymax - ymin
    if x_span <= 0 or y_span <= 0:
        raise ValueError(f"Invalid geographic extent: {extent}")

    x_center = (xmin + xmax) / 2.0
    y_center = (ymin + ymax) / 2.0
    longitude_scale = max(float(np.cos(np.deg2rad(y_center))), 0.25)
    target_ratio = _axes_box_ratio(ax)
    ground_ratio = x_span * longitude_scale / y_span

    if ground_ratio < target_ratio:
        x_span = target_ratio * y_span / longitude_scale
    else:
        y_span = x_span * longitude_scale / target_ratio

    return (
        x_center - x_span / 2.0,
        x_center + x_span / 2.0,
        y_center - y_span / 2.0,
        y_center + y_span / 2.0,
    )


def _geo_ticks(ax, extent, top=True, tick_extent=None) -> None:
    xmin, xmax, ymin, ymax = map(float, extent)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # One degree latitude and the corresponding longitude distance at the map
    # centre receive the same number of pixels.  Because extent was fitted to
    # the allocated box first, setting this aspect does not alter panel size.
    middle_latitude = (ymin + ymax) / 2.0
    longitude_scale = max(float(np.cos(np.deg2rad(middle_latitude))), 0.25)
    fix_geographic_aspect(ax, extent)

    txmin, txmax, tymin, tymax = (
        map(float, tick_extent) if tick_extent is not None else (xmin, xmax, ymin, ymax)
    )
    x_ticks = np.linspace(txmin, txmax, 4 if txmax - txmin < 80 else 7)
    y_ticks = np.linspace(tymin, tymax, 4)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)

    def lon(value, pos=None):
        if abs(value) < 0.01:
            return "0°"
        return f"{abs(value):.0f}°{'W' if value < 0 else 'E'}"

    def lat(value, pos=None):
        if abs(value) < 0.01:
            return "0°"
        return f"{abs(value):.0f}°{'S' if value < 0 else 'N'}"

    ax.xaxis.set_major_formatter(FuncFormatter(lon))
    ax.yaxis.set_major_formatter(FuncFormatter(lat))
    ax.tick_params(axis="x", top=top, labeltop=top, bottom=not top, labelbottom=not top, pad=1)
    ax.grid(color="#DEDEDE", linewidth=0.25, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_color("#777777")


def _add_scale_bar(ax, extent) -> int:
    """Add an automatically sized kilometre scale bar and return its length."""
    xmin, xmax, ymin, ymax = map(float, extent)
    middle_latitude = (ymin + ymax) / 2.0
    longitude_scale = max(float(np.cos(np.deg2rad(middle_latitude))), 0.25)
    displayed_width_km = (xmax - xmin) * 111.32 * longitude_scale
    target_km = displayed_width_km * 0.20
    eligible = SCALEBAR_CANDIDATES_KM[SCALEBAR_CANDIDATES_KM <= target_km]
    scale_km = float(eligible[-1] if len(eligible) else SCALEBAR_CANDIDATES_KM[0])
    bar_fraction = min(scale_km / displayed_width_km, 0.35)

    x0 = 0.06
    x1 = x0 + bar_fraction
    y0 = 0.070
    cap_half_height = 0.012
    ax.add_patch(
        Rectangle(
            (x0 - 0.014, y0 - 0.024),
            bar_fraction + 0.028,
            0.075,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="none",
            alpha=0.82,
            zorder=19,
        )
    )
    ax.plot([x0, x1], [y0, y0], transform=ax.transAxes,
            color="#222222", linewidth=0.8, solid_capstyle="butt", zorder=20)
    ax.plot([x0, x0], [y0 - cap_half_height, y0 + cap_half_height],
            transform=ax.transAxes, color="#222222", linewidth=0.8, zorder=20)
    ax.plot([x1, x1], [y0 - cap_half_height, y0 + cap_half_height],
            transform=ax.transAxes, color="#222222", linewidth=0.8, zorder=20)
    ax.text(
        (x0 + x1) / 2.0,
        y0 + 0.018,
        f"{scale_km:g} km",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
            fontsize=5.0,
        color="#222222",
        zorder=21,
    )
    return int(scale_km)


def _taxon_gradient(taxon: str):
    return LinearSegmentedColormap.from_list(
        f"{taxon}_gradient", ["#F4F4F4", TAXON_COLORS[taxon]]
    )


def draw_distribution_map(ax, grid: gpd.GeoDataFrame, boundary: gpd.GeoDataFrame,
                 extents: pd.DataFrame, zoom_extents: dict):
    boundary.plot(ax=ax, color="#F1F1F1", edgecolor="none", linewidth=0, zorder=0)
    maxima = {}
    for taxon in TAXON_ORDER:
        subset = grid[grid["dominant_taxon"].eq(taxon)].copy()
        maxima[taxon] = max(int(subset["dominant_taxon_richness"].max()), 1) if len(subset) else 1
        if len(subset):
            cmap = _taxon_gradient(taxon)
            norm = Normalize(1, maxima[taxon])
            subset.plot(
                ax=ax,
                color=[cmap(norm(value)) for value in subset["dominant_taxon_richness"]],
                edgecolor="none",
                linewidth=0,
                rasterized=True,
                zorder=2,
            )
    boundary.boundary.plot(ax=ax, color="#AAAAAA", linewidth=0.18, zorder=3)
    formal_labels = {"Central Andes": "b", "Mesoamerica": "c",
                     "East and Southeast Asia": "d", "Southern Africa": "e"}
    for row in extents.itertuples(index=False):
        xmin, xmax, ymin, ymax = zoom_extents[row.region]
        ax.add_patch(
            Rectangle(
                (xmin, ymin), xmax - xmin, ymax - ymin,
                fill=False, edgecolor="#222222", linewidth=0.65, zorder=5,
            )
        )
        ax.text(xmin + (xmax - xmin) * 0.025, ymax - (ymax - ymin) * 0.035,
                formal_labels[row.region], fontsize=5.4,
                fontweight="bold", ha="left", va="top",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.5})
    # Keep the original map slot intact; extend only the empty geographic margin.
    display_extent = _fit_geographic_extent(ax, GLOBAL_CONTENT_EXTENT)
    _geo_ticks(ax, display_extent, top=True, tick_extent=GLOBAL_CONTENT_EXTENT)
    ax._geographic_scale_contract["kind"] = "global_legacy_center_latitude_scale"
    handles = [Patch(facecolor=TAXON_COLORS[t], edgecolor="none",
                     label="Freshwater fishes" if t == "Fish" else t) for t in TAXON_ORDER]
    ax.legend(handles=handles, title="Dominant taxonomic group", loc="lower center",
              bbox_to_anchor=(0.5, 0.01), ncol=5, fontsize=5.0, title_fontsize=5.3,
              columnspacing=0.6, handlelength=1.0)
    panel_label(ax, "A", x=-0.08, y=1.04)
    return maxima


def _subset_by_extent(frame: gpd.GeoDataFrame, extent):
    xmin, xmax, ymin, ymax = extent
    bounds = frame.geometry.bounds
    return frame.loc[
        (bounds.maxx >= xmin) & (bounds.minx <= xmax)
        & (bounds.maxy >= ymin) & (bounds.miny <= ymax)
    ].copy()


def draw_region_panel(ax, label: str, region_row, display_extent,
                      habitat: gpd.GeoDataFrame, boundary: gpd.GeoDataFrame,
                      species: pd.DataFrame) -> tuple[int, int]:
    source_extent = (region_row.xmin, region_row.xmax, region_row.ymin, region_row.ymax)
    region_boundary = _subset_by_extent(boundary, display_extent)
    region_boundary.plot(ax=ax, color="#F1F1F1", edgecolor="none", linewidth=0, zorder=0)
    ids = set(species.loc[species["region"].eq(region_row.region), "species_id"].astype(str))
    selected = habitat[habitat["species_id"].astype(str).isin(ids)].copy()
    # Select habitat using the original scientifically defined window.  The
    # additional display extent is only neutral padding used to avoid stretch.
    selected = _subset_by_extent(selected, source_extent)
    for taxon in TAXON_ORDER:
        layer = selected[selected["taxon"].eq(taxon)]
        if len(layer):
            layer.plot(ax=ax, color=TAXON_COLORS[taxon], edgecolor=TAXON_COLORS[taxon],
                       linewidth=0.08, alpha=0.70, rasterized=True, zorder=2)
    region_boundary.boundary.plot(ax=ax, color="#999999", linewidth=0.25, zorder=3)
    _geo_ticks(ax, extent=display_extent, top=False)
    scale_bar_km = _add_scale_bar(ax, display_extent)
    count = len(ids)
    ax.set_title(f"{region_row.region} (n = {count})", loc="left", pad=2.5,
                 fontsize=6.2, fontweight="bold")
    panel_label(ax, label, x=-0.12, y=1.08)
    return count, scale_bar_km

