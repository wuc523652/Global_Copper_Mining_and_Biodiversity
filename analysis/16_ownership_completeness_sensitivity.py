# -*- coding: utf-8 -*-
"""Ownership-missingness sensitivity analysis.

The main analysis retains all mines, does not normalize known equity shares,
and assigns unknown or unresolved shares to an unattributed category. This
script compares that specification with a complete-case restriction containing
only mines for which 100% of project equity and all corresponding owner-country
information are resolved. Spatial unions are recomputed for the restricted mine
set; no missing ownership is imputed.
"""

from __future__ import annotations

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


BASE_DIR = Path(_nc_project(""))
OUTPUT_ROOT = BASE_DIR / "16_ownership_sensitivity"
COMPLETE_DIR = OUTPUT_ROOT / "complete_case"
COMPARISON_DIR = OUTPUT_ROOT / "comparison"
LOG_DIR = OUTPUT_ROOT / "logs"

COMMON_MODULE_PATH = Path(_nc_script("15_buffer_distance_sensitivity.py"))
ATTRIBUTION_MODULE_PATH = (
    Path(_nc_script("12_attribute_exposure_by_ownership.py"))
)
MINE_MASTER_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Master_30km.csv"
OWNER_LEDGER_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Owner_Share_Ledger.csv"
OWNERSHIP_COVERAGE_CSV = BASE_DIR / "10_mine_ownership" / "Mine_Ownership_Coverage.csv"
MINE_EXPOSURE_CSV = BASE_DIR / "11_mine_exposure" / "Mine_Exposure_30km.csv"
EXPOSURE_GPKG = BASE_DIR / "11_mine_exposure" / "Exposure_Union_30km.gpkg"

MAIN_ATTRIBUTION_DIR = BASE_DIR / "12_ownership_attribution"
MAIN_GLOBAL = MAIN_ATTRIBUTION_DIR / "Global_Responsibility_Summary.csv"
MAIN_INVESTOR = MAIN_ATTRIBUTION_DIR / "Investor_Country_Summary_All.csv"
MAIN_RESOURCE = MAIN_ATTRIBUTION_DIR / "Resource_Country_Summary_All.csv"
MAIN_BILATERAL = MAIN_ATTRIBUTION_DIR / "Bilateral_Responsibility_Flow_All.csv"

EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_CRS = "EPSG:4326"
COMPLETE_SHARE_TOLERANCE_PCT = 1e-6
AREA_TOLERANCE_KM2 = 1e-5
DISPLAY_OMITTED_ISO3 = set()


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ensure_inputs() -> None:
    required = [
        COMMON_MODULE_PATH, ATTRIBUTION_MODULE_PATH, MINE_MASTER_CSV, OWNER_LEDGER_CSV,
        OWNERSHIP_COVERAGE_CSV, MINE_EXPOSURE_CSV, EXPOSURE_GPKG,
        MAIN_GLOBAL, MAIN_INVESTOR, MAIN_RESOURCE, MAIN_BILATERAL,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    for directory in [OUTPUT_ROOT, COMPLETE_DIR, COMPARISON_DIR, LOG_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def build_complete_case_exposure(common, complete_ids: set[str], mine_exposure: pd.DataFrame):
    print("[1/5] Recomputing spatial unions for complete-case mines", flush=True)
    geometry = gpd.read_file(
        EXPOSURE_GPKG, layer="mine_exposure_union", engine="pyogrio"
    ).to_crs(EQUAL_AREA_CRS)
    geometry["property_id"] = common.normalise_id(geometry["property_id"])
    geometry = geometry[geometry["property_id"].isin(complete_ids)].copy()
    geometry["geometry"] = geometry.geometry.map(common.clean_polygon)
    geometry = geometry[geometry.geometry.notna() & ~geometry.geometry.is_empty].copy()

    complete_mine = mine_exposure[mine_exposure["property_id"].isin(complete_ids)].copy()
    if set(geometry["property_id"]) != set(complete_mine.loc[
        complete_mine["mine_exposure_area_km2"] > 0, "property_id"
    ]):
        raise ValueError("Complete-case exposure geometry and table have different positive-mine sets.")

    country_rows = []
    country_geometries = []
    meta = complete_mine.groupby(
        ["resource_iso3", "resource_country", "resource_continent"], as_index=False
    ).agg(
        num_mines=("property_id", "nunique"),
        num_exposed_mines=("mine_exposure_area_km2", lambda values: int((values > 0).sum())),
        mine_exposure_area_sum_km2=("mine_exposure_area_km2", "sum"),
    )
    geometry_groups = {key: group for key, group in geometry.groupby("resource_iso3", sort=False)}
    for record in meta.itertuples(index=False):
        group = geometry_groups.get(record.resource_iso3)
        union_geometry = None if group is None else common.union_in_batches(
            group.geometry.tolist(), label=f"complete-case {record.resource_iso3}"
        )
        union_area = 0.0 if union_geometry is None else float(union_geometry.area / 1_000_000.0)
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
        country_geometries.append(union_geometry)
    country = pd.DataFrame(country_rows)
    country["within_country_cross_mine_duplicate_pct_of_union"] = np.where(
        country["physical_exposure_union_area_km2"] > 0,
        country["within_country_cross_mine_duplicate_area_km2"]
        / country["physical_exposure_union_area_km2"] * 100.0,
        0.0,
    )
    country_gdf = gpd.GeoDataFrame(country, geometry=country_geometries, crs=EQUAL_AREA_CRS)

    strict_global_geometry = common.union_in_batches(
        geometry.geometry.tolist(), label="complete-case global exposure"
    )
    if strict_global_geometry is None:
        raise RuntimeError("Complete-case global exposure union is empty.")
    strict_global_area = float(strict_global_geometry.area / 1_000_000.0)
    country_sum = float(country["physical_exposure_union_area_km2"].sum())
    mine_sum = float(complete_mine["mine_exposure_area_km2"].sum())
    global_table = pd.DataFrame([{
        "global_strict_physical_exposure_area_km2": strict_global_area,
        "country_attributed_union_sum_km2": country_sum,
        "cross_resource_country_duplicate_area_km2": max(0.0, country_sum - strict_global_area),
        "mine_exposure_area_sum_km2": mine_sum,
        "within_country_cross_mine_duplicate_area_km2": max(0.0, mine_sum - country_sum),
        "num_mines": complete_mine["property_id"].nunique(),
        "num_exposed_mines": int((complete_mine["mine_exposure_area_km2"] > 0).sum()),
        "definition": "complete-case strict union across species and selected mines",
    }])

    complete_mine.to_csv(COMPLETE_DIR / "Mine_Exposure_30km_CompleteCase.csv", index=False, encoding="utf-8-sig")
    country.to_csv(COMPLETE_DIR / "Country_Exposure_Union_30km_CompleteCase.csv", index=False, encoding="utf-8-sig")
    global_table.to_csv(COMPLETE_DIR / "Global_Exposure_Union_30km_CompleteCase.csv", index=False, encoding="utf-8-sig")

    gpkg = COMPLETE_DIR / "Exposure_Union_30km_CompleteCase.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    gpd.GeoDataFrame(
        [{"scope": "Complete case", "area_km2": strict_global_area}],
        geometry=[strict_global_geometry], crs=EQUAL_AREA_CRS,
    ).to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="global_exposure_union", driver="GPKG", engine="pyogrio", index=False
    )
    geometry.to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="mine_exposure_union", driver="GPKG", engine="pyogrio", index=False, mode="a"
    )
    country_gdf[country_gdf["physical_exposure_union_area_km2"] > 0].to_crs(OUTPUT_CRS).to_file(
        gpkg, layer="country_exposure_union", driver="GPKG", engine="pyogrio", index=False, mode="a"
    )
    return complete_mine, country, global_table


def flow_table_with_key(table: pd.DataFrame) -> pd.DataFrame:
    output = table[table["relationship"].eq("Foreign")].copy()
    output["flow_key"] = output["investor_iso3"].astype(str) + "->" + output["resource_iso3"].astype(str)
    return output






def main() -> None:
    started = time.time()
    ensure_inputs()
    common = import_module(COMMON_MODULE_PATH, "buffer_sensitivity_common")
    attribution = import_module(ATTRIBUTION_MODULE_PATH, "ownership_attribution")

    master = common.read_table(MINE_MASTER_CSV)
    ledger = common.read_table(OWNER_LEDGER_CSV)
    coverage = common.read_table(OWNERSHIP_COVERAGE_CSV)
    mine_exposure = common.read_table(MINE_EXPOSURE_CSV)
    for table in [master, ledger, coverage, mine_exposure]:
        table["property_id"] = common.normalise_id(table["property_id"])

    complete_mask = (
        pd.to_numeric(coverage["explicit_country_resolved_share_pct"], errors="coerce")
        .sub(100.0).abs().le(COMPLETE_SHARE_TOLERANCE_PCT)
        & pd.to_numeric(coverage["total_unattributed_share_pct"], errors="coerce")
        .abs().le(COMPLETE_SHARE_TOLERANCE_PCT)
    )
    complete_ids = set(coverage.loc[complete_mask, "property_id"].astype(str))
    if not complete_ids:
        raise RuntimeError("No complete-case mines were identified.")
    complete_selection = coverage.assign(is_complete_case=complete_mask.astype(int))[[
        "property_id", "property_name", "resource_country", "resource_iso3",
        "num_owner_records", "num_explicit_owner_records", "explicit_reported_share_pct",
        "explicit_country_resolved_share_pct", "explicit_country_unresolved_share_pct",
        "unreported_share_pct", "total_unattributed_share_pct", "ownership_source",
        "is_complete_case",
    ]]
    complete_selection.to_csv(
        COMPARISON_DIR / "Ownership_CompleteCase_Selection.csv", index=False, encoding="utf-8-sig"
    )

    complete_mine, complete_country, complete_global_exposure = build_complete_case_exposure(
        common, complete_ids, mine_exposure
    )
    complete_master = master[master["property_id"].isin(complete_ids)].copy()
    complete_coverage = coverage[coverage["property_id"].isin(complete_ids)].copy()
    complete_ledger = ledger[ledger["property_id"].isin(complete_ids)].copy()

    print("[2/5] Re-running country-balanced ownership attribution", flush=True)
    complete_attribution_dir = COMPLETE_DIR / "ownership_attribution"
    complete_global, complete_investor, complete_resource, complete_bilateral, complete_balance = common.run_attribution(
        attribution,
        complete_attribution_dir,
        complete_master,
        complete_ledger,
        complete_coverage,
        complete_mine,
        complete_country,
        complete_global_exposure,
    )

    print("[3/5] Comparing main and complete-case results", flush=True)
    main_global = common.read_table(MAIN_GLOBAL)
    main_investor = common.read_table(MAIN_INVESTOR)
    main_resource = common.read_table(MAIN_RESOURCE)
    main_bilateral = common.read_table(MAIN_BILATERAL)
    main_record = main_global.iloc[0]
    complete_record = complete_global.iloc[0]
    scenario = pd.DataFrame([
        {
            "scenario": "All mines",
            **{column: main_record[column] for column in main_global.columns},
        },
        {
            "scenario": "Complete case",
            **{column: complete_record[column] for column in complete_global.columns},
        },
    ])
    scenario["mine_count_coverage_vs_main_pct"] = scenario["num_mines"] / float(main_record["num_mines"]) * 100.0
    scenario["strict_physical_area_coverage_vs_main_pct"] = (
        scenario["global_strict_physical_exposure_area_km2"]
        / float(main_record["global_strict_physical_exposure_area_km2"]) * 100.0
    )
    scenario["country_attributed_area_coverage_vs_main_pct"] = (
        scenario["resource_country_attributed_physical_area_km2"]
        / float(main_record["resource_country_attributed_physical_area_km2"]) * 100.0
    )

    investor_compare = main_investor[[
        "investor_iso3", "investor_country", "foreign_project_attributed_area_km2"
    ]].rename(columns={"foreign_project_attributed_area_km2": "main_foreign_area_km2"}).merge(
        complete_investor[["investor_iso3", "foreign_project_attributed_area_km2"]].rename(
            columns={"foreign_project_attributed_area_km2": "complete_case_foreign_area_km2"}
        ), on="investor_iso3", how="outer"
    )
    investor_compare[["main_foreign_area_km2", "complete_case_foreign_area_km2"]] = investor_compare[[
        "main_foreign_area_km2", "complete_case_foreign_area_km2"
    ]].fillna(0.0)
    investor_compare["main_rank"] = investor_compare["main_foreign_area_km2"].rank(
        method="min", ascending=False
    )
    investor_compare["complete_case_rank"] = investor_compare["complete_case_foreign_area_km2"].rank(
        method="min", ascending=False
    )

    main_flows = flow_table_with_key(main_bilateral)
    complete_flows = flow_table_with_key(complete_bilateral)
    bilateral_compare = main_flows[[
        "flow_key", "investor_iso3", "resource_iso3", "attributed_exposure_area_km2"
    ]].rename(columns={"attributed_exposure_area_km2": "main_foreign_area_km2"}).merge(
        complete_flows[["flow_key", "attributed_exposure_area_km2"]].rename(
            columns={"attributed_exposure_area_km2": "complete_case_foreign_area_km2"}
        ), on="flow_key", how="outer"
    )
    bilateral_compare[["main_foreign_area_km2", "complete_case_foreign_area_km2"]] = bilateral_compare[[
        "main_foreign_area_km2", "complete_case_foreign_area_km2"
    ]].fillna(0.0)
    bilateral_compare["main_rank"] = bilateral_compare["main_foreign_area_km2"].rank(
        method="min", ascending=False
    )
    bilateral_compare["complete_case_rank"] = bilateral_compare["complete_case_foreign_area_km2"].rank(
        method="min", ascending=False
    )

    investor_rank_correlation = float(investor_compare["main_foreign_area_km2"].corr(
        investor_compare["complete_case_foreign_area_km2"], method="spearman"
    ))
    bilateral_rank_correlation = float(bilateral_compare["main_foreign_area_km2"].corr(
        bilateral_compare["complete_case_foreign_area_km2"], method="spearman"
    ))
    investor_top10_main = set(investor_compare.nsmallest(10, "main_rank")["investor_iso3"].astype(str))
    investor_top10_complete = set(investor_compare.nsmallest(10, "complete_case_rank")["investor_iso3"].astype(str))
    flow_top10_main = set(bilateral_compare.nsmallest(10, "main_rank")["flow_key"].astype(str))
    flow_top10_complete = set(bilateral_compare.nsmallest(10, "complete_case_rank")["flow_key"].astype(str))
    metrics = pd.DataFrame([{
        "complete_case_mines": len(complete_ids),
        "complete_case_mine_share_pct": len(complete_ids) / len(master) * 100.0,
        "complete_case_strict_physical_area_coverage_pct": float(
            scenario.loc[scenario["scenario"].eq("Complete case"), "strict_physical_area_coverage_vs_main_pct"].iloc[0]
        ),
        "investor_rank_correlation": investor_rank_correlation,
        "investor_top10_shared": len(investor_top10_main & investor_top10_complete),
        "investor_top10_jaccard": len(investor_top10_main & investor_top10_complete)
        / len(investor_top10_main | investor_top10_complete),
        "bilateral_rank_correlation": bilateral_rank_correlation,
        "bilateral_top10_shared": len(flow_top10_main & flow_top10_complete),
        "bilateral_top10_jaccard": len(flow_top10_main & flow_top10_complete)
        / len(flow_top10_main | flow_top10_complete),
    }])

    scenario.to_csv(COMPARISON_DIR / "Table_S_OwnershipMissingness_Scenarios.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(COMPARISON_DIR / "Table_S_OwnershipMissingness_Stability.csv", index=False, encoding="utf-8-sig")
    investor_display = investor_compare[
        ~investor_compare["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    bilateral_display = bilateral_compare[
        ~bilateral_compare["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
        & ~bilateral_compare["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    investor_compare.sort_values("main_rank").to_csv(
        COMPARISON_DIR / "OwnershipMissingness_InvestorRankComparison_All.csv",
        index=False, encoding="utf-8-sig",
    )
    bilateral_compare.sort_values("main_rank").to_csv(
        COMPARISON_DIR / "OwnershipMissingness_BilateralRankComparison_All.csv",
        index=False, encoding="utf-8-sig",
    )
    investor_display.sort_values("main_rank").to_csv(
        COMPARISON_DIR / "OwnershipMissingness_InvestorRankComparison_NamedDisplay.csv",
        index=False, encoding="utf-8-sig",
    )
    bilateral_display.sort_values("main_rank").to_csv(
        COMPARISON_DIR / "OwnershipMissingness_BilateralRankComparison_NamedDisplay.csv",
        index=False, encoding="utf-8-sig",
    )

    complete_unattributed = float(complete_record["unattributed_area_km2"])
    complete_named = common.read_table(
        complete_attribution_dir / "Investor_Country_Summary_NamedDisplay.csv"
    )
    qa = pd.DataFrame([
        {
            "check": "complete_case_selection_nonempty",
            "value": len(complete_ids),
            "expected": ">0",
            "status": "PASS" if complete_ids else "FAIL",
        },
        {
            "check": "complete_case_resolved_share_exact_100",
            "value": float(complete_coverage["explicit_country_resolved_share_pct"].sub(100).abs().max()),
            "expected": f"<={COMPLETE_SHARE_TOLERANCE_PCT}",
            "status": "PASS" if complete_coverage["explicit_country_resolved_share_pct"].sub(100).abs().max() <= COMPLETE_SHARE_TOLERANCE_PCT else "FAIL",
        },
        {
            "check": "complete_case_unattributed_area_zero",
            "value": complete_unattributed,
            "expected": f"<={AREA_TOLERANCE_KM2}",
            "status": "PASS" if complete_unattributed <= AREA_TOLERANCE_KM2 else "FAIL",
        },
        {
            "check": "complete_case_global_union_not_greater_than_main",
            "value": float(complete_global_exposure.iloc[0]["global_strict_physical_exposure_area_km2"]),
            "expected": f"<={float(main_record['global_strict_physical_exposure_area_km2'])}",
            "status": "PASS" if float(complete_global_exposure.iloc[0]["global_strict_physical_exposure_area_km2"]) <= float(main_record["global_strict_physical_exposure_area_km2"]) + AREA_TOLERANCE_KM2 else "FAIL",
        },
        {
            "check": "complete_case_named_display_retains_all_registrations",
            "value": complete_named["investor_iso3"].nunique(),
            "expected": complete_investor["investor_iso3"].nunique(),
            "status": "PASS" if set(complete_named["investor_iso3"]) == set(complete_investor["investor_iso3"]) else "FAIL",
        },
    ])
    qa.to_csv(COMPARISON_DIR / "OwnershipMissingness_QA.csv", index=False, encoding="utf-8-sig")
    if qa["status"].eq("FAIL").any():
        raise RuntimeError(
            "Ownership-missingness QA failed: "
            + ", ".join(qa.loc[qa["status"].eq("FAIL"), "check"])
        )

    with pd.ExcelWriter(COMPARISON_DIR / "Table_S_OwnershipMissingness_Sensitivity.xlsx") as writer:
        scenario.to_excel(writer, sheet_name="Scenarios", index=False)
        metrics.to_excel(writer, sheet_name="Stability", index=False)
        investor_display.sort_values("main_rank").to_excel(writer, sheet_name="Investor ranks", index=False)
        bilateral_display.sort_values("main_rank").to_excel(writer, sheet_name="Bilateral ranks", index=False)
        qa.to_excel(writer, sheet_name="QA", index=False)


    print("[5/5] Writing reproducibility record", flush=True)
    record = {
        "analysis_id": "ownership_missingness_sensitivity",
        "status": "completed",
        "main_rule": "all mines; known shares not normalized; unknown share unattributed",
        "complete_case_rule": "100% equity share and owner country resolved",
        "complete_case_mines": len(complete_ids),
        "imputation": "none",
        "spatial_union": "recomputed within complete-case sample",
        "registration_display_rule": "all countries and regions retained in named display tables",
        "python": sys.version,
        "platform": platform.platform(),
        "geopandas": gpd.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "elapsed_seconds": time.time() - started,
        "qa_failures": int(qa["status"].eq("FAIL").sum()),
        "output_root": str(OUTPUT_ROOT),
    }
    (LOG_DIR / "OwnershipMissingness_Sensitivity_RunRecord.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Completed. Output: {OUTPUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
