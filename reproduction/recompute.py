"""Recompute manuscript statistics from archived spatial checkpoints.

No expected-result table is read until comparison. The companion raw-input
pipeline is required to regenerate the spatial checkpoints themselves.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import traceback

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STEPS = ("spatial", "randomization", "priority", "ownership", "portfolios", "distance", "complete-case")
DEPENDENCIES = {"portfolios": {"priority"}, "distance": {"priority"}, "complete-case": {"ownership"}}


def module(prefix: str):
    path, = (ROOT / "analysis").glob(prefix + "_*.py")
    name = "recompute_" + path.stem
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def read(path: Path) -> pd.DataFrame:
    from _ownership_labels import normalize_registration_names
    return normalize_registration_names(pd.read_csv(path, dtype={"property_id": str, "mine_id": str, "species_id": str}, low_memory=False))


def save(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def stage(data: Path, work: Path) -> None:
    checkpoint = data / "11Reanalysis_Checkpoints"
    for name in ("Mine_Master_30km.csv", "Mine_Owner_Share_Ledger.csv", "Mine_Ownership_Coverage.csv"):
        save(read(checkpoint / "ownership" / name), work / "10_mine_ownership" / name)
    for kind in ("Mine_Exposure", "Country_Exposure_Union", "Global_Exposure_Union"):
        save(read(data / "5Spatial_Overlap" / f"{kind}_30km_Primary.csv"), work / "11_mine_exposure" / f"{kind}_30km.csv")
    shutil.copy2(checkpoint / "spatial/Exposure_Union_30km.gpkg", work / "11_mine_exposure/Exposure_Union_30km.gpkg")


def randomization(data: Path, work: Path) -> None:
    analysis = module("06")
    observations = read(data / "6Randomization/Mine_Grid_Observations.csv")
    background = read(data / "6Randomization/Background_Grid_Values.csv")
    # Preserve original grid and within-country pool order for seeded draws.
    background = background.sort_values("grid_order", kind="stable").rename(columns={"grid_order": "_grid_pos"})
    countries = analysis.build_country_descriptive(observations, background)
    analysis.build_observed_table(observations, countries)
    pools = analysis.build_country_random_pools(background, countries)
    iterations = analysis.run_randomizations(countries, pools)
    analysis.build_test_summary(observations, iterations)


def spatial(data: Path, work: Path) -> None:
    """Aggregate spatial area checkpoints; intersection geometries are inputs."""
    analysis = module("05")
    output = work / "05_rwr_overlap"
    output.mkdir(parents=True, exist_ok=True)
    grid = pd.read_csv(data / "3Grid_and_Biodiversity/Global_10km_Grid_Biodiversity.csv", usecols=["weighted_richness"])
    thresholds = analysis.compute_thresholds(grid)
    country_links = read(data / "11Reanalysis_Checkpoints/spatial/Country_CopperOverlap_BySpecies.csv")
    counts = country_links.groupby("Country", as_index=False).agg(num_species=("species_id", "nunique"))
    country = read(data / "5Spatial_Overlap/Country_Exposure_Union_30km_Primary.csv").rename(columns={
        "resource_country": "Country", "resource_continent": "Continent", "num_mines": "total_mine_count", "physical_exposure_union_area_km2": "overlap_area_km2"})
    country = country.merge(counts, on="Country", validate="one_to_one")
    global_area = read(data / "5Spatial_Overlap/Global_Exposure_Union_30km_Primary.csv").iloc[0]["global_strict_physical_exposure_area_km2"]
    inventory = read(data / "2Species_Inventory/species_inventory.csv")
    exposure = read(data / "5Spatial_Overlap/Species_Exposure_30km_Primary.csv")
    global_meta = pd.DataFrame([{"Country": "Global", "Continent": "Global", "overlap_area_km2": global_area,
        "total_mine_count": int(country["total_mine_count"].sum()), "num_species": int(exposure["total_overlap_area_km2"].gt(0).sum()),
        "input_species_count": len(inventory), "mine_host_country_count": len(country)}])
    for label, meta in (("Country", country), ("Global", global_meta)):
        intersections = read(data / f"11Reanalysis_Checkpoints/spatial/{label}_Key_RWR_Intersected_Area_Long.csv")
        intersections["rwr_bin"] = intersections["weighted_richness"].map(lambda value: analysis.classify_rwr(value, thresholds))
        summary = analysis.summarize_area(intersections, meta)
        save(summary, output / f"{label}_RWR_Spatial_Summary.csv")


def priority(data: Path, work: Path) -> None:
    analysis = module("07")
    habitat = read(data / "2Species_Inventory/species_inventory.csv").rename(columns={
        "iucn_species_id": "species_id", "taxonomic_group": "taxon", "iucn_threat_category": "status",
        "global_range_area_km2": "total_habitat_area_km2"})
    habitat["species_id"] = habitat["species_id"].astype(str)
    habitat = habitat[["species_id", "scientific_name", "taxon", "status", "total_habitat_area_km2"]]
    exposure = read(data / "5Spatial_Overlap/Species_Exposure_30km_Primary.csv")
    # Deliberately exclude archived ratios and classifications: recompute them.
    exposure = exposure[["species_id", "total_overlap_area_km2", "country_accumulated_overlap_area_km2", "overlap_duplicate_area_km2"]]
    exposure = exposure.merge(habitat[["species_id", "taxon", "status"]], on="species_id", validate="one_to_one")
    exposure = exposure[exposure["total_overlap_area_km2"] > 0].copy()
    country_links = read(data / "11Reanalysis_Checkpoints/spatial/Country_CopperOverlap_BySpecies.csv")
    country_counts = country_links.groupby(["species_id", "taxon", "status"], as_index=False).agg(
        num_countries=("Country", "nunique"), country_occurrences=("Country", "size"),
        countries=("Country", lambda values: ";".join(sorted(set(values.dropna().astype(str))))))
    exposure = exposure.merge(country_counts, on=["species_id", "taxon", "status"], validate="one_to_one")
    # Supply the spatial integration endpoints to the unchanged scientific routine.
    analysis.calculate_all_habitat_areas = lambda: habitat.copy()
    analysis.load_species_overlap_data = lambda: exposure.copy()
    analysis.generate_priority_species_analysis()
    module("08").main()
    module("09").main()


def ownership(data: Path, work: Path) -> None:
    module("12").main()


def portfolios(data: Path, work: Path) -> None:
    master = read(work / "10_mine_ownership/Mine_Master_30km.csv")
    observations = read(data / "6Randomization/Mine_Grid_Observations.csv")
    links = read(data / "11Reanalysis_Checkpoints/spatial/Mine_PriorityScreenedSpecies_Link_AllMines.csv")
    priorities = read(work / "07_priority_species/Species_PrioritySpecies_Priority.csv")
    if set(links["species_id"]) != set(priorities["species_id"]):
        raise ValueError("Spatial links do not cover exactly the recomputed priority species set")
    groups = links.groupby("property_id", as_index=False).agg(
        num_priority_screened_species=("species_id", "nunique"),
        priority_screened_species_sum_overlap_area_km2=("overlap_area_km2", "sum"))
    bridge = master.merge(observations[["mine_id", "is_top10"]].rename(columns={"mine_id": "property_id"}),
                          on="property_id", validate="one_to_one", how="left")
    if bridge["is_top10"].isna().any():
        raise ValueError("Missing mine RWR observations")
    bridge = bridge.merge(groups, on="property_id", how="left", validate="one_to_one")
    bridge["num_priority_screened_species"] = bridge["num_priority_screened_species"].fillna(0).astype(int)
    bridge["priority_screened_species_sum_overlap_area_km2"] = bridge["priority_screened_species_sum_overlap_area_km2"].fillna(0)
    bridge["is_high_rwr"] = bridge["is_top10"].eq(1)
    bridge["has_priority_screened_species_overlap"] = bridge["num_priority_screened_species"].gt(0)
    high, species = bridge["is_high_rwr"], bridge["has_priority_screened_species_overlap"]
    bridge["ecological_attention_category"] = np.select(
        [high & species, high & ~species, ~high & species],
        ["Dual ecological attention", "High RWR only", "Priority-screened species only"], default="Neither")
    save(bridge, work / "13_mine_species_links/Mine_EcologicalAttention_30km.csv")
    save(links, work / "13_mine_species_links/Mine_PriorityScreenedSpecies_Link_AllMines.csv")
    # Recompute the aggregate tables used in Supplementary Data 5.
    save(links.groupby(["species_id", "scientific_name", "taxon", "status"], as_index=False).agg(
        associated_mine_count=("property_id", "nunique"), mine_host_country_count=("resource_iso3", "nunique")),
        work / "13_mine_species_links/Species_Summary.csv")
    country = links.groupby(["resource_iso3", "resource_country"], as_index=False).agg(
        associated_mine_count=("property_id", "nunique"), priority_species_count=("species_id", "nunique"),
        mine_species_association_count=("species_id", "size"))
    maximum = bridge.groupby("resource_iso3")["num_priority_screened_species"].max()
    country["maximum_priority_species_per_mine"] = country["resource_iso3"].map(maximum)
    country["associations_per_mine"] = country["mine_species_association_count"] / country["associated_mine_count"]
    save(country, work / "13_mine_species_links/Country_Summary.csv")
    module("14").main()


def distance(data: Path, work: Path) -> None:
    common, attribution = module("15"), module("12")
    master = read(work / "10_mine_ownership/Mine_Master_30km.csv")
    ledger = read(work / "10_mine_ownership/Mine_Owner_Share_Ledger.csv")
    coverage = read(work / "10_mine_ownership/Mine_Ownership_Coverage.csv")
    baseline = read(work / "07_priority_species/Species_PrioritySpecies_Full.csv")
    meta = baseline[["species_id", "scientific_name", "taxon", "status", "Range", "Conflict_Target"]]
    inputs = {name: read(data / "5Spatial_Overlap" / f"{name}_Sensitivity_10_50km.csv")
              for name in ("Species_Exposure", "Mine_Exposure", "Country_Exposure", "Global_Exposure")}
    results, core = {}, []
    for radius in (10, 20, 30, 40, 50):
        out = work / f"15_buffer_sensitivity/distance_{radius}km"
        species = inputs["Species_Exposure"].loc[lambda x: x["radius_km"].eq(radius),
                    ["species_id", "total_habitat_area_km2", "total_overlap_area_km2"]].merge(meta, on="species_id", validate="one_to_one")
        species = common.calculate_priority(species)
        species["buffer_distance_km"] = radius
        priority_rows = species[species["Priority"].eq("Yes")]
        save(species, out / f"Species_Overlap_and_Priority_{radius}km.csv")
        save(priority_rows, out / f"Priority_Species_{radius}km.csv")
        mine, country, global_exp = [inputs[name].loc[lambda x: x["radius_km"].eq(radius)].drop(columns="radius_km").copy()
                                    for name in ("Mine_Exposure", "Country_Exposure", "Global_Exposure")]
        global_owner, investor, resource, bilateral, _ = common.run_attribution(
            attribution, out / "ownership_attribution", master, ledger, coverage, mine, country, global_exp)
        record = global_owner.iloc[0]
        row: dict[str, int | float] = {"buffer_distance_km": radius, "affected_species_count": int(species["Is_Copper_Affected"].sum()),
               "priority_species_count": len(priority_rows)}
        row.update({f"priority_{category.lower()}_count": int(priority_rows["Range"].eq(category).sum())
                    for category in ("Narrow", "Intermediate", "Widespread")})
        for key in ("global_strict_physical_exposure_area_km2", "resource_country_attributed_physical_area_km2",
                    "traceable_equity_associated_area_km2", "unattributed_area_km2", "traceable_rate_pct",
                    "domestic_equity_associated_area_km2", "foreign_equity_associated_area_km2",
                    "domestic_share_of_traceable_pct", "foreign_share_of_traceable_pct", "num_mines", "num_exposed_mines"):
            row[key] = record[key]
        core.append(row)
        bilateral = bilateral[bilateral["relationship"].eq("Foreign")].copy()
        bilateral["flow_key"] = bilateral["investor_iso3"] + "->" + bilateral["resource_iso3"]
        results[radius] = {"ids": set(priority_rows["species_id"]), "investor": investor, "resource": resource, "bilateral": bilateral}
    stability = []
    reference = results[30]
    for radius, current in results.items():
        shared, union = reference["ids"] & current["ids"], reference["ids"] | current["ids"]
        row = {"buffer_distance_km": radius, "priority_species_shared_with_30km": len(shared),
               "priority_species_gained_vs_30km": len(current["ids"] - reference["ids"]),
               "priority_species_lost_vs_30km": len(reference["ids"] - current["ids"]),
               "priority_set_jaccard_vs_30km": len(shared) / len(union) if union else 1.0}
        for name, key, value in (("investor", "investor_iso3", "foreign_project_attributed_area_km2"),
                                 ("resource", "resource_iso3", "foreign_attributed_area_km2"),
                                 ("bilateral", "flow_key", "attributed_exposure_area_km2")):
            count, jaccard = common.top_set_overlap(reference[name], current[name], key, value, 10)
            row[f"{name}_rank_correlation_vs_30km"] = common.rank_correlation(reference[name], current[name], key, value)
            row[f"{name}_top10_shared_vs_30km"] = count
            row[f"{name}_top10_jaccard_vs_30km"] = jaccard
        stability.append(row)
    save(pd.DataFrame(core), work / "15_buffer_sensitivity/comparison/Table_S_BufferDistance_CoreResults.csv")
    save(pd.DataFrame(stability), work / "15_buffer_sensitivity/comparison/Table_S_BufferDistance_SetAndRankStability.csv")
    # Recompute the full buffer-union areas, including marine areas, for the
    # marginal screening-yield diagnostics in Supplementary Data 9.
    import geopandas as gpd
    areas = {}
    for radius in (10, 20, 30, 40, 50):
        buffers = gpd.read_file(data / f"4Mine_Buffers/Mine_Buffers_{radius}km.gpkg", engine="pyogrio").to_crs("ESRI:54009")
        areas[radius] = float(buffers.geometry.union_all().area / 1e6)
    core_by_radius = {row["buffer_distance_km"]: row for row in core}
    diagnostics = []
    for start, end in zip((10, 20, 30, 40), (20, 30, 40, 50)):
        before, after = core_by_radius[start], core_by_radius[end]
        species_gain = after["priority_species_count"] - before["priority_species_count"]
        area_gain = areas[end] - areas[start]
        share_start = 100 * before["priority_species_count"] / before["affected_species_count"]
        share_end = 100 * after["priority_species_count"] / after["affected_species_count"]
        diagnostics.append({"start_distance_km": start, "end_distance_km": end,
            "priority_species_start": before["priority_species_count"], "priority_species_end": after["priority_species_count"],
            "affected_species_start": before["affected_species_count"], "affected_species_end": after["affected_species_count"],
            "full_buffer_union_area_start_km2": areas[start], "full_buffer_union_area_end_km2": areas[end],
            "additional_priority_species": species_gain, "additional_full_buffer_union_area_km2": area_gain,
            "additional_priority_species_per_100000_km2": species_gain / area_gain * 100000,
            "priority_species_share_start_pct": share_start, "priority_species_share_end_pct": share_end,
            "priority_species_share_change_percentage_points": share_end - share_start})
    save(pd.DataFrame(diagnostics), work / "15_buffer_sensitivity/comparison/Distance_Diagnostics.csv")


def complete_case(data: Path, work: Path) -> None:
    # This executes a new union of 959 selected mine-exposure geometries,
    # followed by attribution and rank analyses. It does not sum mine areas.
    module("16").main()


def compare_frames(actual, expected, label, keys=None, columns=None, atol=1e-7, rtol=1e-10):
    columns = columns if columns is not None else [c for c in expected.select_dtypes(include="number")
        if c not in (keys or []) and c != "elapsed_seconds"]
    required = set(keys or []) | set(columns)
    missing = sorted(required - set(actual.columns))
    missing_reference = sorted(required - set(expected.columns))
    if missing or missing_reference:
        return {"check": label, "status": "FAIL", "missing_columns": missing,
                "missing_reference_columns": missing_reference, "detail": "required comparison fields absent"}
    if len(actual) != len(expected):
        return {"check": label, "status": "FAIL", "detail": f"row count {len(actual)} != {len(expected)}"}
    if keys:
        for frame in (actual, expected):
            for key in keys:
                frame[key] = frame[key].astype("string").str.replace(r"\.0$", "", regex=True)
        if actual.duplicated(keys).any() or expected.duplicated(keys).any():
            return {"check": label, "status": "FAIL", "detail": "duplicate comparison keys"}
        actual = actual.sort_values(keys).reset_index(drop=True)
        expected = expected.sort_values(keys).reset_index(drop=True)
        if not actual[keys].equals(expected[keys]):
            return {"check": label, "status": "FAIL", "detail": "identifier sets differ"}
    else:
        actual, expected = actual.reset_index(drop=True), expected.reset_index(drop=True)
    mismatches, maximum = [], 0.0
    for column in columns:
        left = pd.to_numeric(actual[column], errors="raise").to_numpy(dtype=float)
        right = pd.to_numeric(expected[column], errors="raise").to_numpy(dtype=float)
        if not np.allclose(left, right, atol=atol, rtol=rtol, equal_nan=True):
            mismatches.append(column)
        finite = np.isfinite(left) & np.isfinite(right)
        if finite.any():
            maximum = max(maximum, float(np.max(abs(left[finite] - right[finite]))))
    return {"check": label, "status": "FAIL" if mismatches or not columns else "PASS", "rows": len(actual),
            "numeric_columns": columns, "mismatching_columns": mismatches, "max_absolute_difference": maximum,
            "atol": atol, "rtol": rtol}


def comparisons(data: Path, work: Path, steps: set[str]) -> list[dict]:
    checks = []
    ref = data / "11Reanalysis_Checkpoints/expected"
    choices: list[tuple[str, list[str] | None]] = []
    if "randomization" in steps:
        for name, key in (("Iterations", "iteration"), ("Summary", "metric")):
            checks.append(compare_frames(read(work / f"06_randomization/Mine_RWR_Randomization_{name}.csv"),
                          read(data / f"6Randomization/Randomization_{name}.csv"), f"randomization_{name}", [key], atol=1e-10))
    if "priority" in steps:
        choices += [("07_priority_species/Species_PrioritySpecies_Full.csv", ["species_id"]),
                    ("07_priority_species/Species_PrioritySpecies_Priority.csv", ["species_id"]),
                    ("08_priority_robustness/Robustness_SpeciesResults.csv", ["Species_ID"]),
                    ("08_priority_robustness/Robustness_ByTaxon.csv", ["Group"]),
                    ("08_priority_robustness/Robustness_ByCurrentRange.csv", ["Group"]),
                    ("09_priority_diagnostics/Species_PriorityScreening_Diagnostics_203.csv", ["species_id"])]
    for step, directory in (("ownership", "12_ownership_attribution"), ("portfolios", "14_conservation_portfolios")):
        if step in steps:
            suffix = "_EcologicalPortfolio" if step == "portfolios" else ""
            for filename, keys in ((f"Investor_Country{suffix}_Summary_All.csv", ["investor_iso3"]),
                                   (f"Resource_Country{suffix}_Summary_All.csv", ["resource_iso3"]),
                                   ("Bilateral_EcologicalPortfolio_Summary_All.csv" if step == "portfolios" else "Bilateral_Responsibility_Flow_All.csv", ["investor_iso3", "resource_iso3"]),
                                   ("Global_EcologicalPortfolio_Summary.csv" if step == "portfolios" else "Global_Responsibility_Summary.csv", None)):
                choices.append((f"{directory}/{filename}", keys))
    if "distance" in steps:
        choices += [(f"15_buffer_sensitivity/comparison/Table_S_BufferDistance_{kind}.csv", ["buffer_distance_km"])
                    for kind in ("CoreResults", "SetAndRankStability")]
        expected = pd.read_excel(data / "8Final_Results/Supplementary_Data_9.xlsx", sheet_name="Priority_Lists")
        for radius in (10, 20, 30, 40, 50):
            actual = read(work / f"15_buffer_sensitivity/distance_{radius}km/Priority_Species_{radius}km.csv")
            archived = expected[expected["Buffer distance (km)"].eq(radius)]["IUCN species ID"].astype(str)
            equal = set(actual["species_id"]) == set(archived)
            checks.append({"check": f"priority_species_set_{radius}km", "status": "PASS" if equal else "FAIL", "species_count": len(actual)})
    if "complete-case" in steps:
        choices += [("16_ownership_sensitivity/comparison/Table_S_OwnershipMissingness_Scenarios.csv", ["scenario"]),
                    ("16_ownership_sensitivity/comparison/OwnershipMissingness_InvestorRankComparison_All.csv", ["investor_iso3"]),
                    ("16_ownership_sensitivity/comparison/OwnershipMissingness_BilateralRankComparison_All.csv", ["flow_key"])]
    for relative, keys in choices:
        checks.append(compare_frames(read(work / relative), read(ref / relative), relative, keys))
    from publication_checks import check_publication_tables
    checks.extend(check_publication_tables(data, work, steps, read, compare_frames, module))
    return checks


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "data",
                        help="Complete dataset directory (default: data in the code package)")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--steps", nargs="+", choices=STEPS, default=list(STEPS))
    args = parser.parse_args(argv)
    sys.path.insert(0, str(ROOT / "analysis"))
    from _paths import validate_roots
    data, work = validate_roots(args.data_root, args.output_root)
    if work.is_relative_to(ROOT) or ROOT.is_relative_to(work):
        parser.error("Output must be outside the entire code package and must not contain it")
    if work.exists() and any(work.iterdir()):
        parser.error("Use a new or empty output directory; existing results are not overwritten")
    work.mkdir(parents=True, exist_ok=True)
    os.environ.update(NC_DATA_ROOT=str(data / "10Raw_Inputs"), NC_WORK_ROOT=str(work),
                      NC_CAP_CHECK="1", MPLBACKEND="Agg", MPLCONFIGDIR=str(work / ".matplotlib"))
    chosen = set(args.steps)
    while True:
        expanded = chosen | set().union(*(DEPENDENCIES.get(step, set()) for step in chosen))
        if expanded == chosen:
            break
        chosen = expanded
    report = {"mode": "statistics_from_archived_spatial_checkpoints", "steps_requested": args.steps,
              "steps_executed": [], "scope": "Recalculates numerical inference, classifications and attribution; spatial checkpoints are inputs except complete-case and buffer-distance diagnostic unions. Plotting uses archived publication tables checked here, not the generated CSV files.",
              "checks": [], "status": "RUNNING"}
    try:
        stage(data, work)
        for name in STEPS:
            if name in chosen:
                print(f"RECOMPUTE {name}", flush=True)
                globals()[name.replace("-", "_")](data, work)
                report["steps_executed"].append(name)
        report["checks"] = comparisons(data, work, chosen)
        report["status"] = "PASS" if report["checks"] and all(c["status"] == "PASS" for c in report["checks"]) else "FAIL"
    except Exception as error:
        report["status"] = "ERROR"
        report["error"] = f"{type(error).__name__}: {error}"
        traceback.print_exc()
    (work / "comparison_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(f"{report['status']}: {work / 'comparison_report.json'}", flush=True)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(main())
