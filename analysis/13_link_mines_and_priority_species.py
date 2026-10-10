# -*- coding: utf-8 -*-
"""
13_link_mines_and_priority_species.py

Build the mine-level bridge used by the revised ownership analysis.

Each of the 1,278 copper mines is assigned to exactly one of four mutually
exclusive ecological-attention categories:

1. Neither
2. High RWR only
3. Priority-screened species only
4. Dual ecological attention

"High RWR" means that the mine reference point falls in the global top 10%
of rarity-weighted richness (RWR) grid cells, as fixed by the corrected
randomisation input. "Priority-screened species" means that the 30-km mine
zone has a positive-area intersection with the mapped range of at least one
of the 203 retained species. The 30-km zone is therefore treated as a
potential spatial interaction zone, not as confirmed habitat loss.

No legacy output is modified. All results are written to a new directory.
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd

try:
    from shapely.validation import make_valid
except Exception:  # Shapely < 2
    def make_valid(geometry):
        return geometry.buffer(0)


warnings.filterwarnings("ignore")


BASE_DIR = Path(_nc_project(""))
OUTPUT_DIR = BASE_DIR / "13_mine_species_links"

MINE_MASTER_GPKG = (
    BASE_DIR / "10_mine_ownership" / "Mine_Master_30km.gpkg"
)
MINE_MASTER_LAYER = "mine_master"
RWR_INPUT = (
    BASE_DIR
    / "06_randomization"
    / "Mine_Randomization_Input_QA.csv"
)
PRIORITY_FILE = (
    BASE_DIR / "07_priority_species" / "Species_PrioritySpecies_Priority.csv"
)
SPECIES_OVERLAP_GPKG = (
    BASE_DIR
    / "04_mine_species_overlap"
    / "Species_CopperOverlap_GlobalUnique.gpkg"
)
SPECIES_OVERLAP_LAYER = "global_unique_species_overlap"

OUT_MINE_CSV = OUTPUT_DIR / "Mine_EcologicalAttention_30km.csv"
OUT_MINE_GPKG = OUTPUT_DIR / "Mine_EcologicalAttention_30km.gpkg"
OUT_LINKS = OUTPUT_DIR / "Mine_PriorityScreenedSpecies_Link_AllMines.csv"
OUT_GLOBAL = OUTPUT_DIR / "EcologicalAttention_Global_Summary.csv"
OUT_QA = OUTPUT_DIR / "Mine_EcologicalAttention_QA.csv"

EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_CRS = "EPSG:4326"
EXPECTED_MINES = 1278
EXPECTED_PRIORITY_SPECIES = 203
MIN_INTERSECTION_AREA_KM2 = 1e-8

CATEGORY_ORDER = [
    "Neither",
    "High RWR only",
    "Priority-screened species only",
    "Dual ecological attention",
]


def normalise_id(series: pd.Series) -> pd.Series:
    """Return stable string identifiers and remove spreadsheet-style '.0'."""
    values = series.astype("string").str.strip()
    values = values.str.replace(r"\.0$", "", regex=True)
    return values


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input not found: {path}")


def clean_polygon(geometry):
    if geometry is None or geometry.is_empty:
        return None
    try:
        geometry = make_valid(geometry)
    except Exception:
        geometry = geometry.buffer(0)
    if geometry is None or geometry.is_empty:
        return None
    if geometry.geom_type not in {"Polygon", "MultiPolygon"}:
        parts = [
            part
            for part in getattr(geometry, "geoms", [])
            if part.geom_type in {"Polygon", "MultiPolygon"} and not part.is_empty
        ]
        if not parts:
            return None
        from shapely.ops import unary_union

        geometry = unary_union(parts)
    return geometry


def qa_row(check, status, value, expected, detail):
    return {
        "check": check,
        "status": status,
        "value": value,
        "expected": expected,
        "detail": detail,
    }


def pass_fail(condition: bool) -> str:
    return "PASS" if bool(condition) else "FAIL"


def main() -> None:
    for path in [
        MINE_MASTER_GPKG,
        RWR_INPUT,
        PRIORITY_FILE,
        SPECIES_OVERLAP_GPKG,
    ]:
        require_file(path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mines = gpd.read_file(MINE_MASTER_GPKG, layer=MINE_MASTER_LAYER)
    rwr = pd.read_csv(RWR_INPUT, low_memory=False)
    priority = pd.read_csv(PRIORITY_FILE, low_memory=False)
    species = gpd.read_file(SPECIES_OVERLAP_GPKG, layer=SPECIES_OVERLAP_LAYER)

    mines["property_id"] = normalise_id(mines["property_id"])
    rwr["mine_id"] = normalise_id(rwr["mine_id"])
    priority["species_id"] = normalise_id(priority["species_id"])
    species["species_id"] = normalise_id(species["species_id"])

    rwr_fields = [
        "mine_id",
        "reference_lon",
        "reference_lat",
        "reference_inside_attributed_country",
        "grid_id",
        "weighted_richness",
        "rwr_bin",
        "is_top10",
        "rwr_percentile",
        "grid_match_method",
        "mine_to_grid_centre_km",
    ]
    missing_rwr_fields = [field for field in rwr_fields if field not in rwr.columns]
    if missing_rwr_fields:
        raise ValueError(f"RWR input is missing columns: {missing_rwr_fields}")

    priority_fields = ["species_id", "scientific_name", "taxon", "status"]
    missing_priority_fields = [
        field for field in priority_fields if field not in priority.columns
    ]
    if missing_priority_fields:
        raise ValueError(f"Priority list is missing columns: {missing_priority_fields}")

    priority = priority[priority_fields].drop_duplicates("species_id").copy()
    species = species[species["species_id"].isin(set(priority["species_id"]))].copy()
    # The spatial table also contains taxon/status fields. Keep only its stable
    # identifier and geometry so the screened list is the single authoritative
    # source for all species descriptors and pandas cannot create _x/_y suffixes.
    species = species[["species_id", "geometry"]].merge(
        priority, on="species_id", how="left", validate="one_to_one"
    )

    mines_eq = mines.to_crs(EQUAL_AREA_CRS).copy()
    species_eq = species.to_crs(EQUAL_AREA_CRS).copy()
    mines_eq["geometry"] = mines_eq.geometry.map(clean_polygon)
    species_eq["geometry"] = species_eq.geometry.map(clean_polygon)
    mines_eq = mines_eq[mines_eq.geometry.notna() & ~mines_eq.geometry.is_empty].copy()
    species_eq = species_eq[
        species_eq.geometry.notna() & ~species_eq.geometry.is_empty
    ].copy()

    species_sindex = species_eq.sindex
    link_rows = []
    for mine in mines_eq.itertuples(index=False):
        candidate_positions = list(
            species_sindex.query(mine.geometry, predicate="intersects")
        )
        for position in candidate_positions:
            record = species_eq.iloc[int(position)]
            intersection = mine.geometry.intersection(record.geometry)
            if intersection is None or intersection.is_empty:
                continue
            area_km2 = float(intersection.area / 1_000_000.0)
            if not np.isfinite(area_km2) or area_km2 <= MIN_INTERSECTION_AREA_KM2:
                continue
            link_rows.append(
                {
                    "property_id": str(mine.property_id),
                    "property_name": mine.property_name,
                    "resource_country": mine.resource_country,
                    "resource_iso3": mine.resource_iso3,
                    "species_id": str(record["species_id"]),
                    "scientific_name": record["scientific_name"],
                    "taxon": record["taxon"],
                    "status": record["status"],
                    "overlap_area_km2": area_km2,
                }
            )

    links = pd.DataFrame(link_rows)
    if links.empty:
        raise RuntimeError("No positive-area mine/species intersections were found.")
    links = links.sort_values(
        ["property_id", "taxon", "scientific_name", "species_id"],
        kind="stable",
    ).reset_index(drop=True)

    mine_species = (
        links.groupby("property_id", as_index=False)
        .agg(
            num_priority_screened_species=("species_id", "nunique"),
            priority_screened_species_ids=(
                "species_id",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            priority_screened_species_names=(
                "scientific_name",
                lambda x: ";".join(sorted(set(map(str, x)))),
            ),
            priority_screened_species_sum_overlap_area_km2=(
                "overlap_area_km2",
                "sum",
            ),
        )
    )

    bridge = mines.merge(
        rwr[rwr_fields].rename(columns={"mine_id": "property_id"}),
        on="property_id",
        how="left",
        validate="one_to_one",
    )
    bridge = bridge.merge(
        mine_species,
        on="property_id",
        how="left",
        validate="one_to_one",
    )
    bridge["num_priority_screened_species"] = (
        bridge["num_priority_screened_species"].fillna(0).astype(int)
    )
    bridge["priority_screened_species_ids"] = bridge[
        "priority_screened_species_ids"
    ].fillna("")
    bridge["priority_screened_species_names"] = bridge[
        "priority_screened_species_names"
    ].fillna("")
    bridge["priority_screened_species_sum_overlap_area_km2"] = bridge[
        "priority_screened_species_sum_overlap_area_km2"
    ].fillna(0.0)
    bridge["is_high_rwr"] = pd.to_numeric(bridge["is_top10"], errors="coerce").eq(1)
    bridge["has_priority_screened_species_overlap"] = bridge[
        "num_priority_screened_species"
    ].gt(0)

    bridge["ecological_attention_category"] = np.select(
        [
            bridge["is_high_rwr"] & bridge["has_priority_screened_species_overlap"],
            bridge["is_high_rwr"] & ~bridge["has_priority_screened_species_overlap"],
            ~bridge["is_high_rwr"] & bridge["has_priority_screened_species_overlap"],
        ],
        [
            "Dual ecological attention",
            "High RWR only",
            "Priority-screened species only",
        ],
        default="Neither",
    )
    bridge["ecological_attention_category"] = pd.Categorical(
        bridge["ecological_attention_category"],
        categories=CATEGORY_ORDER,
        ordered=True,
    )

    mine_ids = set(mines["property_id"])
    rwr_ids = set(rwr["mine_id"])
    priority_ids = set(priority["species_id"])
    spatial_priority_ids = set(species_eq["species_id"])
    linked_priority_ids = set(links["species_id"])
    tongshanling_link = (
        links["property_id"].eq("54179") & links["species_id"].eq("145401912")
    ).any()

    qa = [
        qa_row(
            "mine_count",
            pass_fail(len(mines) == EXPECTED_MINES),
            len(mines),
            EXPECTED_MINES,
            "All mines from the ownership master must be retained.",
        ),
        qa_row(
            "mine_id_unique",
            pass_fail(mines["property_id"].is_unique),
            int(mines["property_id"].nunique()),
            len(mines),
            "property_id must be a one-to-one mine key.",
        ),
        qa_row(
            "rwr_id_exact_match",
            pass_fail(mine_ids == rwr_ids),
            len(mine_ids & rwr_ids),
            len(mine_ids),
            f"missing_in_rwr={len(mine_ids-rwr_ids)}; extra_in_rwr={len(rwr_ids-mine_ids)}",
        ),
        qa_row(
            "priority_species_count",
            pass_fail(len(priority_ids) == EXPECTED_PRIORITY_SPECIES),
            len(priority_ids),
            EXPECTED_PRIORITY_SPECIES,
            "The retained priority-screened list is fixed by the prior analysis.",
        ),
        qa_row(
            "priority_species_spatial_coverage",
            pass_fail(priority_ids == spatial_priority_ids),
            len(priority_ids & spatial_priority_ids),
            len(priority_ids),
            f"missing_geometry={len(priority_ids-spatial_priority_ids)}",
        ),
        qa_row(
            "tongshanling_missing_link_recovered",
            pass_fail(tongshanling_link),
            int(tongshanling_link),
            1,
            "Mine 54179 and species 145401912 must be present in the all-mine bridge.",
        ),
        qa_row(
            "category_complete",
            pass_fail(
                bridge["ecological_attention_category"].notna().all()
                and len(bridge) == EXPECTED_MINES
            ),
            int(bridge["ecological_attention_category"].notna().sum()),
            EXPECTED_MINES,
            "Every mine must receive exactly one category.",
        ),
        qa_row(
            "linked_priority_species",
            "INFO",
            len(linked_priority_ids),
            EXPECTED_PRIORITY_SPECIES,
            "Number of retained species intersecting at least one of all 1,278 mine zones.",
        ),
    ]
    qa_df = pd.DataFrame(qa)
    qa_df.to_csv(OUT_QA, index=False, encoding="utf-8-sig")

    failures = qa_df[qa_df["status"].eq("FAIL")]
    if not failures.empty:
        raise RuntimeError(
            "Mine-bridge QA failed: " + ", ".join(failures["check"].astype(str))
        )

    summary = (
        bridge.groupby("ecological_attention_category", observed=False)
        .agg(
            mine_count=("property_id", "size"),
            mine_share_percent=("property_id", lambda x: len(x) / len(bridge) * 100.0),
            priority_screened_species=(
                "priority_screened_species_ids",
                lambda x: len(
                    {
                        item
                        for value in x
                        for item in str(value).split(";")
                        if item
                    }
                ),
            ),
        )
        .reset_index()
    )
    summary["ecological_attention_category"] = summary[
        "ecological_attention_category"
    ].astype("string")

    bridge = bridge.sort_values("source_order", kind="stable").reset_index(drop=True)
    summary.to_csv(OUT_GLOBAL, index=False, encoding="utf-8-sig")
    links.to_csv(OUT_LINKS, index=False, encoding="utf-8-sig")
    bridge.drop(columns="geometry").to_csv(
        OUT_MINE_CSV, index=False, encoding="utf-8-sig"
    )

    if OUT_MINE_GPKG.exists():
        OUT_MINE_GPKG.unlink()
    bridge.to_crs(OUTPUT_CRS).to_file(
        OUT_MINE_GPKG,
        layer="mine_ecological_attention",
        driver="GPKG",
        index=False,
    )

    print("13_mine_species_links completed.")
    print(f"Mines: {len(bridge):,}")
    print(f"Positive mine-species links: {len(links):,}")
    print(f"Linked priority-screened species: {len(linked_priority_ids):,}")
    print(summary.to_string(index=False))
    print(f"QA: {OUT_QA}")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
