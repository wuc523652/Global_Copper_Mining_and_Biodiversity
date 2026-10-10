"""Tiny synthetic checks of selected production spatial calculations.

Creates all geometry in memory; no IUCN, BirdLife, GADM or mine data are used.
Checks the production buffer, polygon union/intersection and area-based RWR
functions against simple geometric expectations. This is a method demonstration,
not a reproduction of the global spatial analysis or its published results.
Provide --output with a directory outside this package.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time

import geopandas as gpd
import numpy as np
from pyproj import CRS, Transformer
from shapely.geometry import box
from shapely.ops import transform

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
sys.dont_write_bytecode = True


def load(filename, alias):
    spec = importlib.util.spec_from_file_location(alias, ROOT / "analysis" / filename)
    if spec is None or spec.loader is None:
        raise ImportError(filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(output):
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    os.environ["NC_DATA_ROOT"] = str(output / "unused_inputs")
    os.environ.pop("COPPER_BIODIV_STAGING_ROOT", None)

    # Production imports and RWR calculations can create intermediate files.
    # Isolate these in a fresh directory, so no old cache can satisfy a check.
    with tempfile.TemporaryDirectory(prefix="spatial_demo_", dir=output) as temporary:
        os.environ["NC_WORK_ROOT"] = temporary
        buffering = load("15_buffer_distance_sensitivity.py", "spatial_demo_buffer")
        overlap = load("04_calculate_mine_species_overlap.py", "spatial_demo_overlap")
        rwr = load("02_calculate_rarity_weighted_richness.py", "spatial_demo_rwr")

        # A 1 km square centred on (0, 0), expressed first in metres in the
        # local azimuthal equidistant projection, then supplied as WGS84.
        local_crs = CRS.from_proj4(
            "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
        )
        to_wgs84 = Transformer.from_crs(local_crs, "EPSG:4326", always_xy=True)
        to_local = Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)
        footprint = transform(to_wgs84.transform, box(-500, -500, 500, 500))
        buffered, buffer_area = buffering.strict_buffer_wgs84(footprint, 1.0)
        if buffered is None:
            raise AssertionError("Synthetic square produced no buffer")
        projected_buffer = transform(to_local.transform, buffered)
        # Square + four strips + 64-sided rounded corners at resolution=16.
        expected_buffer_area = 1.0 + 4.0 + 32.0 * math.sin(math.pi / 32.0)
        np.testing.assert_allclose(buffer_area, expected_buffer_area, rtol=1e-8)
        np.testing.assert_allclose(projected_buffer.bounds, [-1500, -1500, 1500, 1500], atol=1e-5)
        np.testing.assert_allclose(projected_buffer.area / 1e6, expected_buffer_area, rtol=1e-8)

        # These rectangles are illustrative influence zones, independent of the
        # buffer check above. Coordinates are metres in the study's equal-area CRS.
        zones = [box(0, 0, 2000, 2000), box(1000, 0, 3000, 2000)]
        zone_union = overlap.safe_union(zones)
        ranges = [box(0, 0, 2000, 2000), box(1000, 0, 4000, 2000)]
        intersections = [overlap.safe_intersection(geometry, zone_union) for geometry in ranges]
        if zone_union is None or any(geometry is None for geometry in intersections):
            raise AssertionError("Synthetic rectangles produced empty overlap")
        physical_union = overlap.safe_union(intersections)
        if physical_union is None:
            raise AssertionError("Synthetic exposure union is empty")
        np.testing.assert_allclose(zone_union.area / 1e6, 6.0)
        np.testing.assert_allclose([geometry.area / 1e6 for geometry in intersections], [4.0, 4.0])
        np.testing.assert_allclose(physical_union.area / 1e6, 6.0)
        if overlap.safe_intersection(zones[0], box(2000, 0, 3000, 2000)) is not None:
            raise AssertionError("Boundary-only contact must not count as polygon exposure")

        grid = gpd.GeoDataFrame(
            {"grid_id": [1, 2, 3]},
            geometry=[box(x, 0, x + 1000, 2000) for x in (0, 1000, 2000)],
            crs="ESRI:54009",
        )
        species = gpd.GeoDataFrame(
            {"species_id": ["synthetic_1", "synthetic_2"],
             "species_total_area_m2": [4e6, 6e6]},
            geometry=ranges,
            crs=grid.crs,
        )
        rwr.INTERIM_DIR.mkdir(parents=True, exist_ok=True)
        with contextlib.redirect_stdout(io.StringIO()):
            contributions, species_summary = rwr.compute_species_grid_contributions(
                species, grid, "Mammals", "VU"
            )
        contributions = contributions.sort_values(["grid_id", "species_id"]).reset_index(drop=True)
        expected_pairs = [(1, "synthetic_1"), (2, "synthetic_1"),
                          (2, "synthetic_2"), (3, "synthetic_2")]
        if list(zip(contributions["grid_id"], contributions["species_id"])) != expected_pairs:
            raise AssertionError("Unexpected species-cell intersections, including boundary-only contacts")
        np.testing.assert_allclose(contributions["rwr_contribution"], [0.5, 0.5, 1 / 3, 1 / 3])
        grid_results = rwr.aggregate_rwr(contributions).sort_values("grid_id")
        np.testing.assert_allclose(grid_results["rwr_total"], [0.5, 5 / 6, 1 / 3])
        np.testing.assert_array_equal(grid_results["area_based_richness"], [1, 2, 1])
        np.testing.assert_allclose(
            species_summary.sort_values("species_id")["grid_coverage_ratio"], [1.0, 2 / 3]
        )

    contributions.to_csv(output / "spatial_demo_contributions.csv", index=False)
    grid_results.to_csv(output / "spatial_demo_grid_rwr.csv", index=False)
    summary = {
        "status": "PASS",
        "scope": "Synthetic checks of selected production spatial functions; no real study inputs",
        "full_spatial_reanalysis": False,
        "buffer_square_side_km": 1.0,
        "buffer_radius_km": 1.0,
        "buffer_area_km2": buffer_area,
        "expected_polygonal_buffer_area_km2": expected_buffer_area,
        "influence_zone_summed_area_km2": 8.0,
        "influence_zone_union_area_km2": zone_union.area / 1e6,
        "species_overlap_areas_km2": [geometry.area / 1e6 for geometry in intersections],
        "physical_exposure_union_area_km2": physical_union.area / 1e6,
        "grid_rwr": grid_results["rwr_total"].tolist(),
        "expected_grid_rwr": [0.5, 5 / 6, 1 / 3],
        "rwr_note": "Species 2 extends beyond the grid; the full 6 km2 range remains its denominator",
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }
    (output / "spatial_demo_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Output directory outside this package")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output == ROOT or output.is_relative_to(ROOT):
        parser.error("Choose an output directory outside this package")
    run(output)


if __name__ == "__main__":
    main()
