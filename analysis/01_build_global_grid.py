"""
Build a global 10 km x 10 km land grid in Mollweide equal-area projection.

Input:
  GADM global country polygons -> <data-root>/Global_GADM/GC_GADM.shp (EPSG:4326)

Output:
  <work-root>/01_grid/Global_10km_grid.gpkg
    layer: grid
    fields: grid_id (int64), col (int32), row (int32),
            x_cent (float64, m), y_cent (float64, m)
    geometry: 10 km square polygons in ESRI:54009 (Mollweide)

Land rule: a cell is kept if it intersects any GADM polygon (coastal cells in).
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.strtree import STRtree

GADM_PATH = Path(_nc_data("Global_GADM/GC_GADM.shp"))
OUT_DIR = Path(_nc_project("01_grid"))
OUT_PATH = OUT_DIR / "Global_10km_grid.gpkg"
LAYER = "grid"

MOLLWEIDE_CRS = "ESRI:54009"
CELL_SIZE = 10_000          # meters
ROW_CHUNK = 100             # grid rows processed per chunk
SIMPLIFY_TOL = 500          # meters; far below 10km cell -> kept/dropped decisions unchanged


def log(msg: str, t0: float) -> None:
    print(f"[{time.perf_counter() - t0:7.1f}s] {msg}", flush=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if OUT_PATH.exists():
        OUT_PATH.unlink()

    t0 = time.perf_counter()

    log("Loading GADM country polygons", t0)
    land = gpd.read_file(GADM_PATH, engine="pyogrio")
    land = land[land["NAME_0"] != "Antarctica"].copy()
    log(f"  features={len(land)}, src CRS={land.crs}", t0)

    log(f"Reprojecting to Mollweide ({MOLLWEIDE_CRS})", t0)
    land_moll = land.to_crs(MOLLWEIDE_CRS)

    log(f"Simplifying country polygons (tol={SIMPLIFY_TOL} m, topology preserved)", t0)
    land_geoms = np.asarray(
        land_moll.geometry.simplify(SIMPLIFY_TOL, preserve_topology=True).values
    )

    log("Building STRtree on country polygons", t0)
    tree = STRtree(land_geoms)
    prepared = [shapely.prepared.prep(g) for g in land_geoms]

    minx, miny, maxx, maxy = land_moll.total_bounds
    x0 = np.floor(minx / CELL_SIZE) * CELL_SIZE
    y0 = np.floor(miny / CELL_SIZE) * CELL_SIZE
    x1 = np.ceil(maxx / CELL_SIZE) * CELL_SIZE
    y1 = np.ceil(maxy / CELL_SIZE) * CELL_SIZE
    n_cols = int((x1 - x0) / CELL_SIZE)
    n_rows = int((y1 - y0) / CELL_SIZE)
    log(f"  bounds(m): x[{x0:.0f},{x1:.0f}] y[{y0:.0f},{y1:.0f}]", t0)
    log(f"  fishnet: {n_cols} cols x {n_rows} rows = {n_cols * n_rows:,} candidates", t0)

    cols_idx = np.arange(n_cols, dtype=np.int32)
    xc0 = x0 + cols_idx * CELL_SIZE

    kept_chunks = []
    total_kept = 0
    log(f"Building grid in chunks of {ROW_CHUNK} rows", t0)

    for r_start in range(0, n_rows, ROW_CHUNK):
        r_end = min(r_start + ROW_CHUNK, n_rows)
        rows_idx = np.arange(r_start, r_end, dtype=np.int32)
        yc0 = y0 + rows_idx * CELL_SIZE

        # vectorized fishnet for this row chunk
        xx0 = np.repeat(xc0, len(rows_idx))                 # shape (n_cols * n_rows_chunk,)
        yy0 = np.tile(yc0, n_cols)
        xx1 = xx0 + CELL_SIZE
        yy1 = yy0 + CELL_SIZE
        cells_geom = shapely.box(xx0, yy0, xx1, yy1)

        cell_cols = np.repeat(cols_idx, len(rows_idx))
        cell_rows = np.tile(rows_idx, n_cols)

        # bbox-level prefilter via STRtree
        idx_pairs = tree.query(cells_geom, predicate="intersects")
        # idx_pairs is shape (2, M): [cell_idx, country_idx]
        cell_idx = idx_pairs[0]
        if cell_idx.size == 0:
            elapsed = time.perf_counter() - t0
            pct = (r_end / n_rows) * 100
            print(f"  rows {r_start:>5}-{r_end:<5}/{n_rows} ({pct:5.1f}%) | "
                  f"kept=0 | total={total_kept:>9,} | elapsed={elapsed:6.1f}s",
                  flush=True)
            continue

        unique_cells = np.unique(cell_idx)

        kept = gpd.GeoDataFrame(
            {"col": cell_cols[unique_cells],
             "row": cell_rows[unique_cells]},
            geometry=cells_geom[unique_cells],
            crs=MOLLWEIDE_CRS,
        )
        kept_chunks.append(kept)
        total_kept += len(kept)

        elapsed = time.perf_counter() - t0
        pct = (r_end / n_rows) * 100
        print(f"  rows {r_start:>5}-{r_end:<5}/{n_rows} ({pct:5.1f}%) | "
              f"kept={len(kept):>6,} | total={total_kept:>9,} | "
              f"elapsed={elapsed:6.1f}s", flush=True)

    log("Concatenating kept chunks", t0)
    grid = pd.concat(kept_chunks, ignore_index=True)
    grid = gpd.GeoDataFrame(grid, geometry="geometry", crs=MOLLWEIDE_CRS)

    log("Computing centroids and grid_id", t0)
    centroids = grid.geometry.centroid
    grid["grid_id"] = np.arange(len(grid), dtype=np.int64)
    grid["x_cent"] = centroids.x.to_numpy()
    grid["y_cent"] = centroids.y.to_numpy()
    grid = grid[["grid_id", "col", "row", "x_cent", "y_cent", "geometry"]]

    log(f"Writing {len(grid):,} cells to {OUT_PATH}", t0)
    grid.to_file(OUT_PATH, layer=LAYER, driver="GPKG", engine="pyogrio")

    log("Done.", t0)
    print(f"  total land cells: {len(grid):,}")
    print(f"  output: {OUT_PATH}")


if __name__ == "__main__":
    main()