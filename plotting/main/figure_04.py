"""Figure 4: retain the original composite's five-map block in physical units."""

import sys
sys.dont_write_bytecode = True

import geopandas as gpd
from common import INPUTS, ZENODO, new_figure, read_csv, require_expected
import priority_panels as panels


REGIONS = [
    ("B", "Central Andes", 78),
    ("C", "Mesoamerica", 44),
    ("D", "East and Southeast Asia", 29),
    ("E", "Southern Africa", 7),
]


def original_map_slots():
    """Reconstruct the physical map frames of the five-panel composite."""
    from matplotlib.figure import Figure

    reference = Figure(figsize=(183.0 / 25.4, 12.25))
    outer = reference.add_gridspec(
        4, 3, height_ratios=[1.52, 0.86, 0.86, 0.78],
        left=0.065, right=0.985, bottom=0.045, top=0.965,
        wspace=0.22, hspace=0.28,
    )
    maps = outer[1:3, :].subgridspec(
        2, 3, height_ratios=[1.0, 1.0], wspace=0.22, hspace=0.14,
    )
    slots = {"a": maps[0, :2], "b": maps[0, 2], "c": maps[1, 0],
             "d": maps[1, 1], "e": maps[1, 2]}
    boxes = {label: slot.get_position(reference) for label, slot in slots.items()}
    # Translate the block vertically to its standalone canvas. Do not resize it.
    crop_bottom_mm = boxes["c"].y0 * 311.15 - 8.0
    return {label: [box.x0 * 183.0, box.y0 * 311.15 - crop_bottom_mm,
                    box.width * 183.0, box.height * 311.15]
            for label, box in boxes.items()}


def add_original_map_axes(fig):
    bounds = original_map_slots()
    width_mm, height_mm = fig.get_size_inches() * 25.4
    axes = {label: fig.add_axes([x / width_mm, y / height_mm, w / width_mm, h / height_mm])
            for label, (x, y, w, h) in bounds.items()}
    return axes, bounds


def build():
    grid = gpd.read_file(ZENODO / '3Grid_and_Biodiversity/Priority_Species_Nonzero_Grid.gpkg')
    boundary = gpd.read_file(INPUTS / 'main_geography.gpkg', layer='figure_4_land')
    extents = read_csv(INPUTS / 'figure_4_region_extents.csv', columns=['region'], numeric=['xmin', 'xmax', 'ymin', 'ymax'])
    species = read_csv(INPUTS / 'figure_4_region_species.csv', columns=['region', 'species_id'])
    expected_regions = {region for _, region, _ in REGIONS}
    if set(extents["region"]) != expected_regions or extents["region"].duplicated().any():
        raise ValueError("The four established regional windows must occur exactly once")
    require_expected(len(grid), 4771, "Figure 4 occupied grids")
    ranges = gpd.read_file(INPUTS / 'figure_4_habitats.gpkg')
    fig = new_figure(4)
    axes, reference_bounds = add_original_map_axes(fig)
    overview = axes["a"]
    regional_axes = {region: axes[label.lower()] for label, region, _ in REGIONS}
    rows = {row.region: row for row in extents.itertuples(index=False)}
    zoom_extents = {
        region: panels._fit_geographic_extent(regional_axes[region],
                                             (row.xmin, row.xmax, row.ymin, row.ymax))
        for region, row in rows.items()
    }
    maxima = panels.draw_distribution_map(overview, grid, boundary, extents, zoom_extents)
    qa = {"A": {"occupied_grids": len(grid), "dominant_group_maxima": maxima}}
    for label, region, expected in REGIONS:
        count, scale = panels.draw_region_panel(regional_axes[region], label, rows[region],
                                                zoom_extents[region], ranges, boundary, species)
        require_expected(count, expected, f"Figure 4 {region} species")
        qa[label] = {"region": region, "species": count, "scale_bar_km": scale}
    fig.canvas.draw()
    actual_bounds = {}
    for label, ax in axes.items():
        box = ax.get_window_extent()
        actual = [v / fig.dpi * 25.4 for v in (box.x0, box.y0, box.width, box.height)]
        if any(abs(a - b) > 1e-6 for a, b in zip(actual, reference_bounds[label])):
            raise ValueError(f"Figure 4 {label}: original map frame changed: {actual}")
        actual_bounds[label] = actual
    return fig, {"panels": qa, "zoom_display_extents": zoom_extents,
                 "original_map_frame_bounds_mm": reference_bounds,
                 "actual_map_frame_bounds_mm": actual_bounds,
                 "original_frame_preservation_check": "passed; translation only, no resizing",
                 "layout": "original two-row, three-column block; global map spans top-left two columns",
                 "geographic_selection": "original regional windows retained; only empty display margins padded"}


if __name__ == '__main__':
    from common import run
    run(4, build)
