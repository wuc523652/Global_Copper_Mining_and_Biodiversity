# -*- coding: utf-8 -*-
"""Buffer-distance sensitivity analysis for the copper-biodiversity study.

This script regenerates mine buffers at 10, 20, 30, 40 and 50 km from the
original mapped mine polygons, intersects them with the fixed February 2025
IUCN range polygons, re-applies the frozen priority-screening rule, and then
re-runs the country-balanced ownership attribution. The official 30-km results
are treated as the validation reference and are never modified.

All outputs are written below:
    <work-root>/15_buffer_sensitivity
"""

from __future__ import annotations

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

import importlib.util
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
from shapely.geometry import MultiPolygon
from shapely.ops import transform as geometry_transform
from shapely.ops import unary_union

try:
    from shapely.validation import make_valid
except Exception:  # pragma: no cover - Shapely 1 fallback
    def make_valid(geometry):
        return geometry.buffer(0)


warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Frozen paths and parameters
# ---------------------------------------------------------------------------

BASE_DIR = Path(_nc_project(""))
OUTPUT_ROOT = BASE_DIR / "15_buffer_sensitivity"
BUFFER_DIR = OUTPUT_ROOT / "buffers"
CACHE_DIR = OUTPUT_ROOT / "cache"
COMPARISON_DIR = OUTPUT_ROOT / "comparison"
LOG_DIR = OUTPUT_ROOT / "logs"

MINE_FOOTPRINT = Path(
    _nc_data("Global_Copper_Mine/Global_Copper_Mine/Global_Copper_Mine-allinfo.shp")
)
MINE_MASTER_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Master_30km.csv"
OWNER_LEDGER_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Owner_Share_Ledger.csv"
OWNERSHIP_COVERAGE_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Ownership_Coverage.csv"
IUCN_ROOT = Path(_nc_data("IUCNredlist/202502"))
PRIORITY_BASELINE = BASE_DIR / "07_priority_species" / "Species_PrioritySpecies_Full.csv"
OFFICIAL_SPECIES_OVERLAP = (
    BASE_DIR / "04_mine_species_overlap" / "Species_CopperOverlap_GlobalUnique.csv"
)
OFFICIAL_MINE_EXPOSURE = BASE_DIR / "11_mine_exposure" / "Mine_Exposure_30km.csv"
OFFICIAL_GLOBAL_EXPOSURE = BASE_DIR / "11_mine_exposure" / "Global_Exposure_Union_30km.csv"
OFFICIAL_PRIORITY = BASE_DIR / "07_priority_species" / "Species_PrioritySpecies_Priority.csv"
ATTRIBUTION_MODULE_PATH = (
    Path(_nc_script("12_attribute_exposure_by_ownership.py"))
)

DISTANCES_KM = (10, 20, 30, 40, 50)
PRIMARY_DISTANCE_KM = 30
TAXA = ("Amphibians", "Birds", "Fish", "Mammals", "Reptiles")
STATUSES = ("CR", "EN", "VU")
EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_CRS = "EPSG:4326"
AREA_EPSILON_M2 = 1.0
UNION_BATCH_SIZE = 50
RANGE_NARROW_THRESHOLD_KM2 = 200.0
RANGE_WIDESPREAD_THRESHOLD_KM2 = 250000.0
PRIORITY_EPSILON = 1e-12
AREA_VALIDATION_TOLERANCE_KM2 = 1e-3
RELATIVE_VALIDATION_TOLERANCE = 1e-8
DISPLAY_OMITTED_ISO3 = {"HKG"}


def ensure_inputs() -> None:
    required = [
        MINE_FOOTPRINT,
        MINE_MASTER_CSV,
        OWNER_LEDGER_CSV,
        OWNERSHIP_COVERAGE_CSV,
        PRIORITY_BASELINE,
        OFFICIAL_SPECIES_OVERLAP,
        OFFICIAL_MINE_EXPOSURE,
        OFFICIAL_GLOBAL_EXPOSURE,
        OFFICIAL_PRIORITY,
        ATTRIBUTION_MODULE_PATH,
    ]
    required.extend(IUCN_ROOT / f"{taxon}-{status}" / "data_0.shp" for taxon in TAXA for status in STATUSES)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    for directory in [OUTPUT_ROOT, BUFFER_DIR, CACHE_DIR, COMPARISON_DIR, LOG_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def normalise_id(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)


def clean_id(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def clean_polygon(geometry):
    if geometry is None or geometry.is_empty:
        return None
    try:
        geometry = make_valid(geometry)
    except Exception:
        try:
            geometry = geometry.buffer(0)
        except Exception:
            return None
    if geometry is None or geometry.is_empty:
        return None
    if geometry.geom_type in {"Polygon", "MultiPolygon"}:
        return geometry
    parts = [
        part
        for part in getattr(geometry, "geoms", [])
        if part.geom_type in {"Polygon", "MultiPolygon"} and not part.is_empty
    ]
    if not parts:
        return None
    return clean_polygon(unary_union(parts))


def safe_intersection(left, right):
    if left is None or right is None or left.is_empty or right.is_empty:
        return None
    try:
        result = left.intersection(right)
    except Exception:
        try:
            result = make_valid(left).intersection(make_valid(right))
        except Exception:
            return None
    return clean_polygon(result)


def union_in_batches(geometries, label: str = "geometry"):
    level = [clean_polygon(geometry) for geometry in geometries if geometry is not None and not geometry.is_empty]
    level = [geometry for geometry in level if geometry is not None and not geometry.is_empty]
    if not level:
        return None
    round_number = 0
    while len(level) > 1:
        round_number += 1
        next_level = []
        for start in range(0, len(level), UNION_BATCH_SIZE):
            merged = clean_polygon(unary_union(level[start:start + UNION_BATCH_SIZE]))
            if merged is not None and not merged.is_empty:
                next_level.append(merged)
        print(f"    {label}: union round {round_number}, {len(level):,} -> {len(next_level):,}", flush=True)
        level = next_level
    return level[0]


def local_aeqd(lon: float, lat: float) -> CRS:
    return CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs"
    )


def strict_buffer_wgs84(geometry, radius_km: float, resolution: int = 16):
    """Reproduce the project's per-feature local-AEQD strict-distance buffer."""
    geometry = clean_polygon(geometry)
    if geometry is None:
        return None, np.nan
    representative = geometry.representative_point()
    projection = local_aeqd(float(representative.x), float(representative.y))
    forward = Transformer.from_crs(OUTPUT_CRS, projection, always_xy=True).transform
    inverse = Transformer.from_crs(projection, OUTPUT_CRS, always_xy=True).transform
    projected = clean_polygon(geometry_transform(forward, geometry))
    if projected is None:
        return None, np.nan
    buffered = clean_polygon(projected.buffer(float(radius_km) * 1000.0, resolution=resolution))
    if buffered is None:
        return None, np.nan
    area_km2 = float(buffered.area / 1_000_000.0)
    output = clean_polygon(geometry_transform(inverse, buffered))
    return output, area_km2


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig", dtype={"property_id": str, "species_id": str}, low_memory=False)


def build_buffers(master: pd.DataFrame) -> tuple[dict[int, gpd.GeoDataFrame], pd.DataFrame]:
    print("[1/7] Regenerating strict-distance mine buffers", flush=True)
    source = gpd.read_file(MINE_FOOTPRINT, engine="pyogrio")
    if source.crs is None:
        raise ValueError("Mine footprint layer has no CRS.")
    source = source.to_crs(OUTPUT_CRS)
    source["property_id"] = normalise_id(source["Property_I"])
    if source["property_id"].duplicated().any():
        raise ValueError("Original mine footprint Property_I is not unique.")
    master = master.copy()
    master["property_id"] = normalise_id(master["property_id"])
    if set(source["property_id"]) != set(master["property_id"]):
        raise ValueError("Original mine footprints and ownership master have different Property ID sets.")
    attributes = master[[
        "property_id", "property_name", "resource_country", "resource_iso3",
        "resource_continent", "source_fid", "source_order",
    ]].copy()
    source = source[["property_id", "geometry"]].merge(
        attributes, on="property_id", how="left", validate="one_to_one"
    )
    source = gpd.GeoDataFrame(source, geometry="geometry", crs=OUTPUT_CRS).sort_values(
        "source_order", kind="stable"
    ).reset_index(drop=True)

    outputs: dict[int, gpd.GeoDataFrame] = {}
    qa_rows = []
    for distance in DISTANCES_KM:
        started = time.time()
        geometries = []
        areas = []
        for index, geometry in enumerate(source.geometry, 1):
            if index == 1 or index % 200 == 0 or index == len(source):
                print(f"  {distance} km: mine {index:,}/{len(source):,}", flush=True)
            buffered, area_km2 = strict_buffer_wgs84(geometry, distance)
            if buffered is None:
                raise RuntimeError(f"Buffer generation failed at {distance} km, row {index}.")
            if buffered.geom_type == "Polygon":
                buffered = MultiPolygon([buffered])
            geometries.append(buffered)
            areas.append(area_km2)
        table = source.drop(columns="geometry").copy()
        table["buffer_distance_km"] = distance
        table["local_aeqd_buffer_area_km2"] = areas
        buffers = gpd.GeoDataFrame(table, geometry=geometries, crs=OUTPUT_CRS).to_crs(EQUAL_AREA_CRS)
        buffers["geometry"] = buffers.geometry.map(clean_polygon)
        # The official downstream master recalculates buffer area in Mollweide;
        # retain the local-AEQD area separately but use the same official area field.
        buffers["buffer_area_km2"] = buffers.geometry.area / 1_000_000.0
        if buffers.geometry.isna().any() or buffers["property_id"].duplicated().any():
            raise RuntimeError(f"Invalid or duplicate buffers at {distance} km.")
        output_path = BUFFER_DIR / f"Mine_Buffers_{distance}km.gpkg"
        if output_path.exists():
            output_path.unlink()
        buffers.to_crs(OUTPUT_CRS).to_file(
            output_path, layer="mine_buffers", driver="GPKG", engine="pyogrio", index=False
        )
        outputs[distance] = buffers
        qa_rows.append({
            "buffer_distance_km": distance,
            "num_mines": len(buffers),
            "num_unique_property_ids": buffers["property_id"].nunique(),
            "minimum_buffer_area_km2": buffers["buffer_area_km2"].min(),
            "maximum_buffer_area_km2": buffers["buffer_area_km2"].max(),
            "elapsed_seconds": time.time() - started,
            "status": "PASS",
        })
    qa = pd.DataFrame(qa_rows)
    qa.to_csv(COMPARISON_DIR / "Buffer_Generation_QA.csv", index=False, encoding="utf-8-sig")
    return outputs, qa


def load_species_ranges(baseline: pd.DataFrame) -> gpd.GeoDataFrame:
    cache_path = CACHE_DIR / "IUCN_202502_Species_GlobalRanges_Mollweide.parquet"
    if cache_path.exists():
        print("[2/7] Loading cached IUCN species ranges", flush=True)
        species = gpd.read_parquet(cache_path)
    else:
        print("[2/7] Dissolving fixed February 2025 IUCN polygons by species", flush=True)
        parts = []
        for taxon in TAXA:
            for status in STATUSES:
                path = IUCN_ROOT / f"{taxon}-{status}" / "data_0.shp"
                print(f"  {taxon}-{status}", flush=True)
                layer = gpd.read_file(
                    path, engine="pyogrio", columns=["ID_NO", "SCI_NAME"]
                )
                layer["species_id"] = normalise_id(layer["ID_NO"])
                layer["scientific_name_iucn"] = layer["SCI_NAME"].astype("string").str.strip()
                layer = layer[
                    layer["species_id"].ne("") & layer.geometry.notna() & ~layer.geometry.is_empty
                ].copy()
                layer = layer.to_crs(EQUAL_AREA_CRS)
                layer["geometry"] = layer.geometry.map(clean_polygon)
                layer = layer[layer.geometry.notna() & ~layer.geometry.is_empty].copy()
                rows = []
                geoms = []
                for species_id, group in layer.groupby("species_id", sort=False):
                    merged = union_in_batches(group.geometry.tolist(), label=f"{taxon}-{status}-{species_id}")
                    if merged is None:
                        continue
                    names = group["scientific_name_iucn"].dropna().astype(str)
                    rows.append({
                        "species_id": str(species_id),
                        "scientific_name_iucn": "" if names.empty else names.iloc[0],
                        "taxon": taxon,
                        "status": status,
                    })
                    geoms.append(merged)
                parts.append(gpd.GeoDataFrame(rows, geometry=geoms, crs=EQUAL_AREA_CRS))
        species = pd.concat(parts, ignore_index=True)
        species = gpd.GeoDataFrame(species, geometry="geometry", crs=EQUAL_AREA_CRS)
        species.to_parquet(cache_path, index=False)

    species["species_id"] = normalise_id(species["species_id"])
    key = ["species_id", "taxon", "status"]
    if species.duplicated(key).any():
        raise ValueError("Dissolved IUCN species ranges contain duplicate species/taxon/status keys.")
    baseline = baseline.copy()
    baseline["species_id"] = normalise_id(baseline["species_id"])
    if baseline.duplicated(key).any():
        raise ValueError("Priority baseline contains duplicate species/taxon/status keys.")
    keep = [
        "species_id", "taxon", "status", "scientific_name", "total_habitat_area_km2",
        "Range", "Protection_Target", "Conflict_Target", "Target_Rule_Version",
        "Target_Rule_DOI",
    ]
    species = species.merge(baseline[keep], on=key, how="outer", indicator=True, validate="one_to_one")
    unmatched = species[species["_merge"].ne("both")][key + ["_merge"]].copy()
    unmatched.to_csv(COMPARISON_DIR / "IUCN_Baseline_KeyMatch_QA.csv", index=False, encoding="utf-8-sig")
    if not unmatched.empty:
        raise ValueError(f"IUCN geometry/baseline key mismatch: {len(unmatched)} records.")
    species = species.drop(columns="_merge")
    if len(species) != 5890:
        raise ValueError(f"Expected 5,890 species, found {len(species):,}.")
    print(f"  Fixed species universe: {len(species):,}", flush=True)
    return gpd.GeoDataFrame(species, geometry="geometry", crs=EQUAL_AREA_CRS)


def calculate_priority(overlap: pd.DataFrame) -> pd.DataFrame:
    result = overlap.copy()
    habitat = pd.to_numeric(result["total_habitat_area_km2"], errors="coerce")
    overlap_area = pd.to_numeric(result["total_overlap_area_km2"], errors="coerce").fillna(0.0)
    result["Cal_IAP_raw"] = np.where(habitat > 0, overlap_area / habitat, np.nan)
    result["Cal_IAP"] = result["Cal_IAP_raw"].clip(lower=0.0, upper=1.0)
    result["Impact_Percent"] = result["Cal_IAP"] * 100.0
    result["Is_Copper_Affected"] = overlap_area > 0
    result["Priority"] = np.where(
        result["Cal_IAP"].notna()
        & result["Conflict_Target"].notna()
        & (result["Cal_IAP"] > result["Conflict_Target"] + PRIORITY_EPSILON),
        "Yes",
        "No",
    )
    return result


def calculate_species_overlap(
    species: gpd.GeoDataFrame,
    buffer_union,
    distance: int,
) -> tuple[pd.DataFrame, gpd.GeoDataFrame, object]:
    candidate_positions = list(species.sindex.query(buffer_union, predicate="intersects"))
    rows = []
    geometries = []
    for counter, position in enumerate(candidate_positions, 1):
        if counter == 1 or counter % 500 == 0 or counter == len(candidate_positions):
            print(f"  {distance} km species intersections: {counter:,}/{len(candidate_positions):,}", flush=True)
        record = species.iloc[int(position)]
        intersection = safe_intersection(record.geometry, buffer_union)
        if intersection is None or intersection.area <= AREA_EPSILON_M2:
            continue
        rows.append({
            "species_id": str(record["species_id"]),
            "taxon": record["taxon"],
            "status": record["status"],
            "total_overlap_area_km2": float(intersection.area / 1_000_000.0),
        })
        geometries.append(intersection)
    positive = gpd.GeoDataFrame(rows, geometry=geometries, crs=EQUAL_AREA_CRS)
    key = ["species_id", "taxon", "status"]
    metadata = species.drop(columns="geometry")
    overlap = metadata.merge(positive.drop(columns="geometry"), on=key, how="left", validate="one_to_one")
    overlap["total_overlap_area_km2"] = overlap["total_overlap_area_km2"].fillna(0.0)
    overlap = calculate_priority(overlap)
    physical_union = union_in_batches(positive.geometry.tolist(), label=f"{distance} km species exposure")
    if physical_union is None:
        raise RuntimeError(f"No positive exposure union at {distance} km.")
    return overlap, positive, physical_union


def calculate_mine_species_counts(
    buffers: gpd.GeoDataFrame,
    species: gpd.GeoDataFrame,
    distance: int,
) -> pd.DataFrame:
    candidates = gpd.sjoin(
        buffers[["property_id", "geometry"]],
        species[["species_id", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    species_geometry = species.geometry.to_dict()
    mine_geometry = buffers.set_index("property_id").geometry.to_dict()
    valid_pairs = []
    total = len(candidates)
    for counter, row in enumerate(
        candidates[["property_id", "species_id", "index_right"]].itertuples(index=False), 1
    ):
        if counter == 1 or counter % 5000 == 0 or counter == total:
            print(f"  {distance} km mine-species pairs: {counter:,}/{total:,}", flush=True)
        intersection = safe_intersection(
            mine_geometry[str(row.property_id)], species_geometry[int(row.index_right)]
        )
        if intersection is not None and intersection.area > AREA_EPSILON_M2:
            valid_pairs.append((str(row.property_id), str(row.species_id)))
    pairs = pd.DataFrame(valid_pairs, columns=["property_id", "species_id"]).drop_duplicates()
    counts = pairs.groupby("property_id", as_index=False).agg(num_exposed_species=("species_id", "nunique"))
    counts = buffers[["property_id"]].merge(counts, on="property_id", how="left", validate="one_to_one")
    counts["num_exposed_species"] = counts["num_exposed_species"].fillna(0).astype(int)
    return counts


def calculate_exposure_tables(
    buffers: gpd.GeoDataFrame,
    species: gpd.GeoDataFrame,
    physical_union,
    distance: int,
    distance_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, gpd.GeoDataFrame]:
    counts = calculate_mine_species_counts(buffers, species, distance)
    rows = []
    geometries = []
    for counter, mine in enumerate(buffers.itertuples(index=False), 1):
        if counter == 1 or counter % 200 == 0 or counter == len(buffers):
            print(f"  {distance} km mine exposure: {counter:,}/{len(buffers):,}", flush=True)
        intersection = safe_intersection(mine.geometry, physical_union)
        area_km2 = 0.0 if intersection is None else float(intersection.area / 1_000_000.0)
        rows.append({
            "property_id": str(mine.property_id),
            "property_name": mine.property_name,
            "resource_country": mine.resource_country,
            "resource_iso3": mine.resource_iso3,
            "resource_continent": mine.resource_continent,
            "buffer_area_km2": float(mine.buffer_area_km2),
            "mine_exposure_area_km2": area_km2,
            "exposure_share_of_buffer_pct": (
                area_km2 / float(mine.buffer_area_km2) * 100.0 if mine.buffer_area_km2 > 0 else 0.0
            ),
        })
        geometries.append(intersection)
    mine = pd.DataFrame(rows).merge(counts, on="property_id", how="left", validate="one_to_one")
    mine_gdf = gpd.GeoDataFrame(mine, geometry=geometries, crs=EQUAL_AREA_CRS)
    mine_positive = mine_gdf[mine_gdf["mine_exposure_area_km2"] > 0].copy()

    country_rows = []
    country_geometries = []
    mine_meta = mine.groupby(
        ["resource_iso3", "resource_country", "resource_continent"], as_index=False
    ).agg(
        num_mines=("property_id", "nunique"),
        num_exposed_mines=("mine_exposure_area_km2", lambda values: int((values > 0).sum())),
        mine_exposure_area_sum_km2=("mine_exposure_area_km2", "sum"),
    )
    positive_groups = {
        key: group for key, group in mine_positive.groupby("resource_iso3", sort=False)
    }
    for record in mine_meta.itertuples(index=False):
        group = positive_groups.get(record.resource_iso3)
        country_union = None if group is None else union_in_batches(
            group.geometry.tolist(), label=f"{distance} km {record.resource_iso3}"
        )
        union_area = 0.0 if country_union is None else float(country_union.area / 1_000_000.0)
        country_rows.append({
            "resource_iso3": record.resource_iso3,
            "resource_country": record.resource_country,
            "resource_continent": record.resource_continent,
            "num_mines": int(record.num_mines),
            "num_exposed_mines": int(record.num_exposed_mines),
            "mine_exposure_area_sum_km2": float(record.mine_exposure_area_sum_km2),
            "physical_exposure_union_area_km2": union_area,
            "within_country_cross_mine_duplicate_area_km2": max(
                0.0, float(record.mine_exposure_area_sum_km2) - union_area
            ),
        })
        country_geometries.append(country_union)
    country = pd.DataFrame(country_rows)
    country["within_country_cross_mine_duplicate_pct_of_union"] = np.where(
        country["physical_exposure_union_area_km2"] > 0,
        country["within_country_cross_mine_duplicate_area_km2"]
        / country["physical_exposure_union_area_km2"] * 100.0,
        0.0,
    )
    country_gdf = gpd.GeoDataFrame(country, geometry=country_geometries, crs=EQUAL_AREA_CRS)
    country_positive = country_gdf[country_gdf["physical_exposure_union_area_km2"] > 0].copy()

    strict_global_area = float(physical_union.area / 1_000_000.0)
    country_sum = float(country["physical_exposure_union_area_km2"].sum())
    mine_sum = float(mine["mine_exposure_area_km2"].sum())
    global_table = pd.DataFrame([{
        "global_strict_physical_exposure_area_km2": strict_global_area,
        "country_attributed_union_sum_km2": country_sum,
        "cross_resource_country_duplicate_area_km2": max(0.0, country_sum - strict_global_area),
        "mine_exposure_area_sum_km2": mine_sum,
        "within_country_cross_mine_duplicate_area_km2": max(0.0, mine_sum - country_sum),
        "num_mines": mine["property_id"].nunique(),
        "num_exposed_mines": int((mine["mine_exposure_area_km2"] > 0).sum()),
        "definition": "strict union across species; country attribution follows mine host-country attribute",
    }])

    mine.to_csv(distance_dir / f"Mine_Exposure_{distance}km.csv", index=False, encoding="utf-8-sig")
    country.to_csv(
        distance_dir / f"Country_Exposure_Union_{distance}km.csv", index=False, encoding="utf-8-sig"
    )
    global_table.to_csv(
        distance_dir / f"Global_Exposure_Union_{distance}km.csv", index=False, encoding="utf-8-sig"
    )
    gpkg = distance_dir / f"Exposure_Union_{distance}km.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    gpd.GeoDataFrame(
        [{"scope": "Global", "area_km2": strict_global_area}],
        geometry=[physical_union], crs=EQUAL_AREA_CRS,
    ).to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="global_exposure_union", driver="GPKG", engine="pyogrio", index=False
    )
    mine_positive.to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="mine_exposure_union", driver="GPKG", engine="pyogrio", index=False, mode="a"
    )
    country_positive.to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="country_exposure_union", driver="GPKG", engine="pyogrio", index=False, mode="a"
    )
    return mine, country, global_table, mine_positive


def import_attribution_module():
    spec = importlib.util.spec_from_file_location("hkg_attribution", ATTRIBUTION_MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load attribution module: {ATTRIBUTION_MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure_attribution_outputs(module, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    module.OUTPUT_DIR = output_dir
    module.OUT_MINE_ATTRIBUTION = output_dir / "Mine_InvestorCountry_Attribution.csv"
    module.OUT_MINE_BALANCE = output_dir / "Mine_Attribution_Balance.csv"
    module.OUT_GLOBAL = output_dir / "Global_Responsibility_Summary.csv"
    module.OUT_INVESTOR = output_dir / "Investor_Country_Summary_All.csv"
    module.OUT_RESOURCE = output_dir / "Resource_Country_Summary_All.csv"
    module.OUT_BILATERAL = output_dir / "Bilateral_Responsibility_Flow_All.csv"
    module.OUT_INVESTOR_DISPLAY = output_dir / "Investor_Country_Summary_NamedDisplay.csv"
    module.OUT_RESOURCE_DISPLAY = output_dir / "Resource_Country_Summary_NamedDisplay.csv"
    module.OUT_BILATERAL_DISPLAY = output_dir / "Bilateral_Responsibility_Flow_NamedDisplay.csv"
    module.OUT_DISPLAY_OMISSION_AUDIT = output_dir / "NamedDisplay_Omission_Audit.csv"
    module.OUT_QA = output_dir / "Attribution_Balance_QA.csv"
    module.OUT_METHOD = output_dir / "Attribution_Method_Metadata.txt"


def run_attribution(
    module,
    output_dir: Path,
    master: pd.DataFrame,
    ledger: pd.DataFrame,
    coverage: pd.DataFrame,
    mine_exposure: pd.DataFrame,
    country_exposure: pd.DataFrame,
    global_exposure: pd.DataFrame,
):
    configure_attribution_outputs(module, output_dir)
    module.ATTRIBUTION_AREA_MODE = "country_union_balanced"
    base, factors = module.build_mine_base(master, coverage, mine_exposure, country_exposure)
    attribution, balance, omitted = module.build_investor_attribution(base, ledger)
    global_summary = module.build_global_summary(balance, global_exposure)
    investor = module.build_investor_summary(attribution)
    resource = module.build_resource_summary(balance, country_exposure)
    bilateral = module.build_bilateral_flows(attribution)
    investor_display = investor[~investor["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)].copy()
    resource_display = resource[~resource["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)].copy()
    bilateral_display = bilateral[
        ~bilateral["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
        & ~bilateral["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    investor_display.to_csv(module.OUT_INVESTOR_DISPLAY, index=False, encoding="utf-8-sig")
    resource_display.to_csv(module.OUT_RESOURCE_DISPLAY, index=False, encoding="utf-8-sig")
    bilateral_display.to_csv(module.OUT_BILATERAL_DISPLAY, index=False, encoding="utf-8-sig")
    module.save_qa(base, factors, balance, global_summary, omitted)
    return global_summary, investor, resource, bilateral, balance


def rank_correlation(reference: pd.DataFrame, alternative: pd.DataFrame, key: str, value: str) -> float:
    left = reference[[key, value]].rename(columns={value: "reference"})
    right = alternative[[key, value]].rename(columns={value: "alternative"})
    merged = left.merge(right, on=key, how="inner")
    if len(merged) < 2:
        return np.nan
    return float(merged["reference"].rank(ascending=False).corr(
        merged["alternative"].rank(ascending=False), method="pearson"
    ))


def top_set_overlap(reference: pd.DataFrame, alternative: pd.DataFrame, key: str, value: str, n: int = 10):
    ref = set(reference.nlargest(n, value)[key].astype(str))
    alt = set(alternative.nlargest(n, value)[key].astype(str))
    union = ref | alt
    return len(ref & alt), len(ref & alt) / len(union) if union else np.nan


def validate_primary_distance(
    generated_buffers: gpd.GeoDataFrame,
    generated_species: pd.DataFrame,
    generated_mines: pd.DataFrame,
    generated_global: pd.DataFrame,
    baseline_master: pd.DataFrame,
) -> pd.DataFrame:
    official_species = read_table(OFFICIAL_SPECIES_OVERLAP).rename(
        columns={"global_unique_overlap_area_km2": "official_overlap_area_km2"}
    )
    official_species["species_id"] = normalise_id(official_species["species_id"])
    generated = generated_species.loc[
        generated_species["total_overlap_area_km2"] > 0,
        ["species_id", "taxon", "status", "total_overlap_area_km2"],
    ].copy()
    species_compare = official_species[[
        "species_id", "taxon", "status", "official_overlap_area_km2"
    ]].merge(generated, on=["species_id", "taxon", "status"], how="outer", indicator=True)
    species_compare[["official_overlap_area_km2", "total_overlap_area_km2"]] = species_compare[[
        "official_overlap_area_km2", "total_overlap_area_km2"
    ]].fillna(0.0)
    species_compare["absolute_difference_km2"] = (
        species_compare["official_overlap_area_km2"] - species_compare["total_overlap_area_km2"]
    ).abs()
    species_compare.to_csv(
        COMPARISON_DIR / "Primary30km_SpeciesOverlap_Validation.csv", index=False, encoding="utf-8-sig"
    )

    official_mines = read_table(OFFICIAL_MINE_EXPOSURE).rename(
        columns={"mine_exposure_area_km2": "official_mine_exposure_area_km2"}
    )
    mine_compare = official_mines[["property_id", "official_mine_exposure_area_km2"]].merge(
        generated_mines[["property_id", "mine_exposure_area_km2"]],
        on="property_id", how="outer", indicator=True,
    )
    mine_compare["absolute_difference_km2"] = (
        mine_compare["official_mine_exposure_area_km2"] - mine_compare["mine_exposure_area_km2"]
    ).abs()
    mine_compare.to_csv(
        COMPARISON_DIR / "Primary30km_MineExposure_Validation.csv", index=False, encoding="utf-8-sig"
    )

    official_priority = read_table(OFFICIAL_PRIORITY)
    official_priority_ids = set(normalise_id(official_priority["species_id"]))
    generated_priority_ids = set(generated_species.loc[
        generated_species["Priority"].eq("Yes"), "species_id"
    ].astype(str))
    official_global = float(read_table(OFFICIAL_GLOBAL_EXPOSURE).iloc[0][
        "global_strict_physical_exposure_area_km2"
    ])
    generated_global_area = float(generated_global.iloc[0][
        "global_strict_physical_exposure_area_km2"
    ])
    master_area = baseline_master[["property_id", "buffer_area_km2"]].copy()
    buffer_compare = master_area.merge(
        generated_buffers[["property_id", "buffer_area_km2"]].rename(
            columns={"buffer_area_km2": "generated_buffer_area_km2"}
        ), on="property_id", how="outer", indicator=True,
    )
    buffer_compare["absolute_difference_km2"] = (
        buffer_compare["buffer_area_km2"] - buffer_compare["generated_buffer_area_km2"]
    ).abs()
    buffer_compare.to_csv(
        COMPARISON_DIR / "Primary30km_BufferArea_Validation.csv", index=False, encoding="utf-8-sig"
    )
    checks = [
        {
            "check": "buffer_property_id_exact",
            "value": int(buffer_compare["_merge"].eq("both").sum()),
            "expected": len(baseline_master),
            "status": "PASS" if buffer_compare["_merge"].eq("both").all() else "FAIL",
        },
        {
            "check": "maximum_buffer_area_difference_km2",
            "value": float(buffer_compare["absolute_difference_km2"].max()),
            "expected": f"<={AREA_VALIDATION_TOLERANCE_KM2}",
            "status": "PASS" if buffer_compare["absolute_difference_km2"].max() <= AREA_VALIDATION_TOLERANCE_KM2 else "FAIL",
        },
        {
            "check": "species_key_exact",
            "value": int(species_compare["_merge"].eq("both").sum()),
            "expected": len(official_species),
            "status": "PASS" if species_compare["_merge"].eq("both").all() else "FAIL",
        },
        {
            "check": "maximum_species_overlap_difference_km2",
            "value": float(species_compare["absolute_difference_km2"].max()),
            "expected": f"<={AREA_VALIDATION_TOLERANCE_KM2} or relative<={RELATIVE_VALIDATION_TOLERANCE}",
            "status": "PASS" if (
                species_compare["absolute_difference_km2"].max() <= AREA_VALIDATION_TOLERANCE_KM2
                or species_compare["absolute_difference_km2"].sum()
                / max(official_species["official_overlap_area_km2"].sum(), 1.0)
                <= RELATIVE_VALIDATION_TOLERANCE
            ) else "FAIL",
        },
        {
            "check": "priority_species_set_exact",
            "value": len(generated_priority_ids & official_priority_ids),
            "expected": len(official_priority_ids),
            "status": "PASS" if generated_priority_ids == official_priority_ids else "FAIL",
        },
        {
            "check": "mine_exposure_property_id_exact",
            "value": int(mine_compare["_merge"].eq("both").sum()),
            "expected": len(official_mines),
            "status": "PASS" if mine_compare["_merge"].eq("both").all() else "FAIL",
        },
        {
            "check": "maximum_mine_exposure_difference_km2",
            "value": float(mine_compare["absolute_difference_km2"].max()),
            "expected": f"<={AREA_VALIDATION_TOLERANCE_KM2} or relative<={RELATIVE_VALIDATION_TOLERANCE}",
            "status": "PASS" if (
                mine_compare["absolute_difference_km2"].max() <= AREA_VALIDATION_TOLERANCE_KM2
                or mine_compare["absolute_difference_km2"].sum()
                / max(official_mines["official_mine_exposure_area_km2"].sum(), 1.0)
                <= RELATIVE_VALIDATION_TOLERANCE
            ) else "FAIL",
        },
        {
            "check": "global_physical_exposure_area_difference_km2",
            "value": abs(generated_global_area - official_global),
            "expected": f"<={AREA_VALIDATION_TOLERANCE_KM2} or relative<={RELATIVE_VALIDATION_TOLERANCE}",
            "status": "PASS" if (
                abs(generated_global_area - official_global) <= AREA_VALIDATION_TOLERANCE_KM2
                or abs(generated_global_area - official_global) / official_global <= RELATIVE_VALIDATION_TOLERANCE
            ) else "FAIL",
        },
    ]
    qa = pd.DataFrame(checks)
    qa.to_csv(COMPARISON_DIR / "Primary30km_Reproduction_QA.csv", index=False, encoding="utf-8-sig")
    return qa






def main() -> None:
    started = time.time()
    ensure_inputs()
    master = read_table(MINE_MASTER_CSV)
    ledger = read_table(OWNER_LEDGER_CSV)
    coverage = read_table(OWNERSHIP_COVERAGE_CSV)
    baseline = read_table(PRIORITY_BASELINE)
    for table in [master, ledger, coverage]:
        table["property_id"] = normalise_id(table["property_id"])
    baseline["species_id"] = normalise_id(baseline["species_id"])

    buffers_by_distance, _ = build_buffers(master)
    species = load_species_ranges(baseline)
    attribution_module = import_attribution_module()

    core_rows = []
    distance_results = {}
    print("[3/7] Running spatial, priority and ownership analyses", flush=True)
    for distance in DISTANCES_KM:
        distance_started = time.time()
        print(f"\nDistance: {distance} km", flush=True)
        distance_dir = OUTPUT_ROOT / f"distance_{distance}km"
        distance_dir.mkdir(parents=True, exist_ok=True)
        buffers = buffers_by_distance[distance]
        buffer_union = union_in_batches(buffers.geometry.tolist(), label=f"{distance} km mine buffers")
        if buffer_union is None:
            raise RuntimeError(f"Mine-buffer union is empty at {distance} km.")
        species_result, positive_species, physical_union = calculate_species_overlap(
            species, buffer_union, distance
        )
        species_result["buffer_distance_km"] = distance
        species_result.to_csv(
            distance_dir / f"Species_Overlap_and_Priority_{distance}km.csv",
            index=False, encoding="utf-8-sig",
        )
        priority = species_result[species_result["Priority"].eq("Yes")].copy()
        priority.to_csv(
            distance_dir / f"Priority_Species_{distance}km.csv", index=False, encoding="utf-8-sig"
        )
        mine_exp, country_exp, global_exp, _ = calculate_exposure_tables(
            buffers, species, physical_union, distance, distance_dir
        )
        radius_master = master.drop(columns=["buffer_area_km2"], errors="ignore").merge(
            buffers[["property_id", "buffer_area_km2"]].drop(columns="geometry", errors="ignore"),
            on="property_id", how="left", validate="one_to_one",
        )
        attribution_dir = distance_dir / "ownership_attribution"
        global_owner, investor, resource, bilateral, balance = run_attribution(
            attribution_module, attribution_dir, radius_master, ledger, coverage,
            mine_exp, country_exp, global_exp,
        )
        global_record = global_owner.iloc[0]
        core_rows.append({
            "buffer_distance_km": distance,
            "affected_species_count": int(species_result["Is_Copper_Affected"].sum()),
            "priority_species_count": int(priority.shape[0]),
            "priority_narrow_count": int(priority["Range"].eq("Narrow").sum()),
            "priority_intermediate_count": int(priority["Range"].eq("Intermediate").sum()),
            "priority_widespread_count": int(priority["Range"].eq("Widespread").sum()),
            "global_strict_physical_exposure_area_km2": float(global_exp.iloc[0]["global_strict_physical_exposure_area_km2"]),
            "resource_country_attributed_physical_area_km2": float(global_record["resource_country_attributed_physical_area_km2"]),
            "traceable_equity_associated_area_km2": float(global_record["traceable_equity_associated_area_km2"]),
            "unattributed_area_km2": float(global_record["unattributed_area_km2"]),
            "traceable_rate_pct": float(global_record["traceable_rate_pct"]),
            "domestic_equity_associated_area_km2": float(global_record["domestic_equity_associated_area_km2"]),
            "foreign_equity_associated_area_km2": float(global_record["foreign_equity_associated_area_km2"]),
            "domestic_share_of_traceable_pct": float(global_record["domestic_share_of_traceable_pct"]),
            "foreign_share_of_traceable_pct": float(global_record["foreign_share_of_traceable_pct"]),
            "num_mines": int(global_record["num_mines"]),
            "num_exposed_mines": int(global_record["num_exposed_mines"]),
            "elapsed_seconds": time.time() - distance_started,
        })
        distance_results[distance] = {
            "species": species_result,
            "priority_ids": set(priority["species_id"].astype(str)),
            "mine": mine_exp,
            "global_exposure": global_exp,
            "investor": investor,
            "resource": resource,
            "bilateral": bilateral,
        }

    print("[4/7] Validating regenerated 30-km results", flush=True)
    primary_result = distance_results[PRIMARY_DISTANCE_KM]
    primary_qa = validate_primary_distance(
        buffers_by_distance[PRIMARY_DISTANCE_KM],
        primary_result["species"],
        primary_result["mine"],
        primary_result["global_exposure"],
        master,
    )
    if primary_qa["status"].eq("FAIL").any():
        failed = ", ".join(primary_qa.loc[primary_qa["status"].eq("FAIL"), "check"])
        raise RuntimeError(f"Primary 30-km reproduction failed: {failed}")

    print("[5/7] Comparing sets and country rankings", flush=True)
    primary = distance_results[PRIMARY_DISTANCE_KM]
    stability_rows = []
    for distance in DISTANCES_KM:
        current = distance_results[distance]
        priority_intersection = primary["priority_ids"] & current["priority_ids"]
        priority_union = primary["priority_ids"] | current["priority_ids"]
        investor_intersection, investor_top10_jaccard = top_set_overlap(
            primary["investor"], current["investor"], "investor_iso3",
            "foreign_project_attributed_area_km2", 10,
        )
        resource_intersection, resource_top10_jaccard = top_set_overlap(
            primary["resource"], current["resource"], "resource_iso3",
            "foreign_attributed_area_km2", 10,
        )
        ref_bilateral = primary["bilateral"].copy()
        cur_bilateral = current["bilateral"].copy()
        ref_bilateral["flow_key"] = ref_bilateral["investor_iso3"].astype(str) + "->" + ref_bilateral["resource_iso3"].astype(str)
        cur_bilateral["flow_key"] = cur_bilateral["investor_iso3"].astype(str) + "->" + cur_bilateral["resource_iso3"].astype(str)
        ref_bilateral = ref_bilateral[ref_bilateral["relationship"].eq("Foreign")]
        cur_bilateral = cur_bilateral[cur_bilateral["relationship"].eq("Foreign")]
        bilateral_intersection, bilateral_top10_jaccard = top_set_overlap(
            ref_bilateral, cur_bilateral, "flow_key", "attributed_exposure_area_km2", 10
        )
        stability_rows.append({
            "buffer_distance_km": distance,
            "priority_species_shared_with_30km": len(priority_intersection),
            "priority_species_gained_vs_30km": len(current["priority_ids"] - primary["priority_ids"]),
            "priority_species_lost_vs_30km": len(primary["priority_ids"] - current["priority_ids"]),
            "priority_set_jaccard_vs_30km": len(priority_intersection) / len(priority_union) if priority_union else 1.0,
            "investor_rank_correlation_vs_30km": rank_correlation(
                primary["investor"], current["investor"], "investor_iso3",
                "foreign_project_attributed_area_km2",
            ),
            "investor_top10_shared_vs_30km": investor_intersection,
            "investor_top10_jaccard_vs_30km": investor_top10_jaccard,
            "resource_rank_correlation_vs_30km": rank_correlation(
                primary["resource"], current["resource"], "resource_iso3",
                "foreign_attributed_area_km2",
            ),
            "resource_top10_shared_vs_30km": resource_intersection,
            "resource_top10_jaccard_vs_30km": resource_top10_jaccard,
            "bilateral_rank_correlation_vs_30km": rank_correlation(
                ref_bilateral, cur_bilateral, "flow_key", "attributed_exposure_area_km2",
            ),
            "bilateral_top10_shared_vs_30km": bilateral_intersection,
            "bilateral_top10_jaccard_vs_30km": bilateral_top10_jaccard,
        })
    core = pd.DataFrame(core_rows)
    stability = pd.DataFrame(stability_rows)
    core.to_csv(COMPARISON_DIR / "Table_S_BufferDistance_CoreResults.csv", index=False, encoding="utf-8-sig")
    stability.to_csv(COMPARISON_DIR / "Table_S_BufferDistance_SetAndRankStability.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(COMPARISON_DIR / "Table_S_BufferDistance_Sensitivity.xlsx") as writer:
        core.to_excel(writer, sheet_name="Core results", index=False)
        stability.to_excel(writer, sheet_name="Set and rank stability", index=False)
        primary_qa.to_excel(writer, sheet_name="30km reproduction QA", index=False)


    print("[7/7] Writing reproducibility record", flush=True)
    run_record = {
        "analysis_id": "buffer_distance_sensitivity",
        "status": "completed",
        "distances_km": list(DISTANCES_KM),
        "primary_distance_km": PRIMARY_DISTANCE_KM,
        "iucn_snapshot": "2025-02 (author-provided version description)",
        "mine_buffer_method": "per-feature local azimuthal equidistant projection; strict outward distance",
        "equal_area_crs": EQUAL_AREA_CRS,
        "species_attribute_filtering": "none; downloaded vector ranges used directly",
        "ownership_missing_share_rule": "not normalized; retained as unattributed",
        "hkg_rule": "included in calculations; omitted only from named-country display tables",
        "python": sys.version,
        "platform": platform.platform(),
        "geopandas": gpd.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "elapsed_seconds": time.time() - started,
        "primary_30km_qa_failures": int(primary_qa["status"].eq("FAIL").sum()),
        "output_root": str(OUTPUT_ROOT),
    }
    (LOG_DIR / "BufferDistance_Sensitivity_RunRecord.json").write_text(
        json.dumps(run_record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Completed. Output: {OUTPUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
