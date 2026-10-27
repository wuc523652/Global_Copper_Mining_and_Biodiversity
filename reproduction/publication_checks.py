"""Compare newly computed values against the submitted supplementary workbooks.

Mappings are explicit so the report states exactly which published columns
were checked. Workbooks are comparison targets, never computation inputs.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def check_publication_tables(data, work, steps, read, compare, module):
    checks = []

    def add(actual, number, sheet, mapping, keys=None, exact=()):
        missing = set(mapping) - set(actual.columns)
        if missing:
            raise ValueError(f"Missing recomputed fields in Data {number}/{sheet}: {missing}")
        actual = actual[list(mapping)].rename(columns=mapping).copy()
        destination = work / "publication_numeric_tables" / f"Data_{number}_{sheet}.csv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        actual.to_csv(destination, index=False, encoding="utf-8-sig")
        expected = pd.read_excel(data / f"8Final_Results/Supplementary_Data_{number}.xlsx", sheet_name=sheet)
        keys = [mapping[key] for key in keys] if keys else None
        columns = [mapping[c] for c in mapping if mapping[c] not in (keys or []) and mapping[c] not in exact]
        checks.append(compare(actual, expected, f"Supplementary_Data_{number}/{sheet}", keys, columns))
        if exact:
            left, right = actual, expected
            if keys:
                left, right = left.sort_values(keys).reset_index(drop=True), right.sort_values(keys).reset_index(drop=True)
            for column in exact:
                equal = left[column].astype(str).str.lower().equals(right[column].astype(str).str.lower())
                checks.append({"check": f"Supplementary_Data_{number}/{sheet}/{column}", "status": "PASS" if equal else "FAIL", "comparison": "exact labels"})

    if "spatial" in steps:
        common = {"total_mine_count": "Mine count", "num_species": "Exposed species count",
            "area_total_rwr_km2": "Influence area in non-zero RWR cells (km²)",
            "pct_ge90": "Classified influence area in top-10% RWR cells (%)", "pct_ge99": "Classified influence area in top-1% RWR cells (%)"}
        country_map = {"resource_iso3": "Mine host iso3", **common,
            "overlap_area_km2": "Unioned exposed range area (km²)",
            "area_le50_km2": "Rwr 0 50 area (km²)", "area_50_75_km2": "Rwr 50 75 area (km²)",
            "area_75_90_km2": "Rwr 75 90 area (km²)", "area_90_99_km2": "Rwr 90 99 area (km²)",
            "area_99_100_km2": "Rwr 99 100 area (km²)"}
        add(read(work / "05_rwr_overlap/Country_RWR_Spatial_Summary.csv"), 3, "Country_Spatial", country_map, ["resource_iso3"])
        add(read(work / "05_rwr_overlap/Global_RWR_Spatial_Summary.csv"), 3, "Global_Spatial", {**common,
            "overlap_area_km2": "Global union exposed range area (km²)", "input_species_count": "Input species count",
            "mine_host_country_count": "Mine host country count", "area_ge90_km2": "Top10 influence area (km²)", "area_ge99_km2": "Top1 influence area (km²)"})
    if "randomization" in steps:
        actual = read(work / "06_randomization/Mine_RWR_Randomization_Iterations.csv")
        add(actual, 3, "Randomization_Iterations", {"iteration": "Iteration", "top10_mine_count": "Top10 mine count",
            "pct_top10_mines": "Reference points in top-10% RWR cells (%)", "mean_rwr_percentile": "Mean global RWR percentile"}, ["iteration"])

    if "priority" in steps:
        mapping = {"species_id": "IUCN species ID", "total_habitat_area_km2": "Mapped global range area (km²)",
            "total_overlap_area_km2": "Unioned overlap area (km²)", "Impact_Percent": "Overlapping range (%)",
            "Protection_Target_Percent": "Representation target (%)", "Conflict_Target_Percent": "Main-rule overlap threshold (%)",
            "Priority": "Main-rule priority species"}
        for filename, sheet in (("Species_CopperAffected_Full.csv", "Exposed_Species"), ("Species_PrioritySpecies_Priority.csv", "Priority_Species")):
            actual = read(work / "07_priority_species" / filename)
            if sheet == "Priority_Species":
                actual["Priority"] = "Yes"
            add(actual, 4, sheet, mapping, ["species_id"], ["Main-rule priority species"])
        # The published main and alternative rules are uncapped. The cap is a
        # separate sensitivity; derive uncapped thresholds with the original
        # rule function instead of copying the already archived comparisons.
        robust = module("08")
        original = read(work / "07_priority_species/Species_CopperAffected_Full.csv")
        rule = original[["species_id", "total_habitat_area_km2", "Impact_Percent", "Conflict_Target_Percent", "Priority"]].copy()
        protection = robust.representation_target_percent(rule["total_habitat_area_km2"].to_numpy(), robust.SCHEMES[robust.TRAD], apply_cap=False)
        rule["alternative_threshold"] = 100 - protection
        rule["alternative_priority"] = np.where(rule["Impact_Percent"] > rule["alternative_threshold"] + robust.EPS_PERCENT, "Yes", "No")
        rule["main_margin"] = rule["Impact_Percent"] - rule["Conflict_Target_Percent"]
        rule["alternative_margin"] = rule["Impact_Percent"] - rule["alternative_threshold"]
        add(rule, 4, "Rule_Comparison", {"species_id": "IUCN species ID", "total_habitat_area_km2": "Mapped global range area (km²)",
            "Impact_Percent": "Overlapping range (%)", "Conflict_Target_Percent": "Main-rule overlap threshold (%)",
            "alternative_threshold": "Alternative-rule overlap threshold (%)", "main_margin": "Margin above main threshold (percentage points)",
            "alternative_margin": "Margin above alternative threshold (percentage points)", "Priority": "Main-rule priority species",
            "alternative_priority": "Alternative-rule priority species"}, ["species_id"], ["Main-rule priority species", "Alternative-rule priority species"])

    owner_map = {"investor_iso3": "Owner registration iso3", "total_attributed_area_km2": "Total attributed area (km²)",
        "domestic_project_attributed_area_km2": "Domestic attributed area (km²)", "foreign_project_attributed_area_km2": "Cross border attributed area (km²)",
        "foreign_project_share_pct": "Cross border share (%)", "num_mines": "Num mines", "num_resource_countries": "Num mine host countries",
        "num_foreign_resource_countries": "Num cross border mine host countries", "foreign_area_rank": "Cross border area rank",
        "share_of_global_foreign_area_pct": "Share of global cross border area (%)"}
    host_map = {"resource_iso3": "Mine host iso3", "num_mines": "Num mines", "num_exposed_mines": "Num exposed mines",
        "physical_exposure_union_area_km2": "Physical exposure union area (km²)", "traceable_area_km2": "Traceable area (km²)",
        "unattributed_area_km2": "Unattributed area (km²)", "domestic_attributed_area_km2": "Domestic attributed area (km²)",
        "foreign_attributed_area_km2": "Cross border attributed area (km²)", "traceable_rate_pct": "Traceable rate (%)",
        "foreign_share_of_traceable_pct": "Cross border share of traceable (%)", "foreign_area_rank": "Cross border area rank"}
    global_map = {"num_mines": "Num mines", "global_strict_physical_exposure_area_km2": "Global strict physical exposure area (km²)",
        "resource_country_attributed_physical_area_km2": "Mine host country attributed physical area (km²)",
        "traceable_equity_associated_area_km2": "Traceable equity associated area (km²)", "unattributed_area_km2": "Unattributed area (km²)",
        "traceable_rate_pct": "Traceable rate (%)", "domestic_equity_associated_area_km2": "Domestic equity associated area (km²)",
        "foreign_equity_associated_area_km2": "Cross border equity associated area (km²)", "domestic_share_of_traceable_pct": "Domestic share of traceable (%)",
        "foreign_share_of_traceable_pct": "Cross border share of traceable (%)"}
    for step, directory, prefix in (("ownership", "12_ownership_attribution", ""),
                                  ("complete-case", "16_ownership_sensitivity/complete_case/ownership_attribution", "CompleteCase_")):
        if step not in steps:
            continue
        add(read(work / directory / "Investor_Country_Summary_All.csv"), 6, prefix + "Owners" if prefix else "Owner_Countries", owner_map, ["investor_iso3"])
        add(read(work / directory / "Resource_Country_Summary_All.csv"), 6, prefix + "Hosts" if prefix else "Mine_Host_Countries", host_map, ["resource_iso3"])
        if not prefix:
            add(read(work / directory / "Global_Responsibility_Summary.csv"), 6, "Global_Attribution", global_map)
        flows = read(work / directory / "Bilateral_Responsibility_Flow_All.csv")
        flows = flows[flows["relationship"].eq("Foreign")].copy()
        flows["share"] = 100 * flows["attributed_exposure_area_km2"] / flows["attributed_exposure_area_km2"].sum()
        flows["rank"] = flows["attributed_exposure_area_km2"].rank(ascending=False, method="average")
        add(flows, 7, "CompleteCase_Flows" if prefix else "CrossBorder_Flows", {
            "investor_iso3": "Owner registration iso3", "resource_iso3": "Mine host iso3",
            "attributed_exposure_area_km2": "Attributed exposure area (km²)", "num_mines": "Num mines",
            "num_owner_records": "Num owner records", "equity_share_sum_across_mines_pct": "Equity share sum across mines (%)",
            "share": "Share of global cross border area (%)", "rank": "Cross border area rank"}, ["investor_iso3", "resource_iso3"])
    if "complete-case" in steps:
        actual = read(work / "16_ownership_sensitivity/comparison/Table_S_OwnershipMissingness_Stability.csv")
        add(actual, 6, "Completeness_Stability", {"complete_case_mines": "Complete case mines", "complete_case_mine_share_pct": "Complete case mine share (%)",
            "complete_case_strict_physical_area_coverage_pct": "Complete case strict physical area coverage (%)",
            "investor_rank_correlation": "Owner registration country or region rank correlation",
            "investor_top10_shared": "Owner registration country or region top10 shared", "investor_top10_jaccard": "Owner registration country or region top10 jaccard",
            "bilateral_rank_correlation": "Bilateral rank correlation", "bilateral_top10_shared": "Bilateral top10 shared", "bilateral_top10_jaccard": "Bilateral top10 jaccard"})
    if "portfolios" in steps:
        bridge = read(work / "13_mine_species_links/Mine_EcologicalAttention_30km.csv")
        mine = bridge[bridge["num_priority_screened_species"] > 0].copy()
        add(mine, 5, "Mine_Summary", {"property_id": "Mine ID", "num_priority_screened_species": "Priority species count"}, ["property_id"])
        add(read(work / "13_mine_species_links/Species_Summary.csv"), 5, "Species_Summary", {
            "species_id": "IUCN species ID", "associated_mine_count": "Associated mine count", "mine_host_country_count": "Mine host country count"}, ["species_id"])
        add(read(work / "13_mine_species_links/Country_Summary.csv"), 5, "Country_Summary", {
            "resource_iso3": "Mine host iso3", "associated_mine_count": "Associated mine count", "priority_species_count": "Priority species count",
            "mine_species_association_count": "Mine species association count", "maximum_priority_species_per_mine": "Maximum priority species per mine",
            "associations_per_mine": "Associations per mine"}, ["resource_iso3"])
        add(bridge, 8, "Mine_Priority_Categories", {"property_id": "Mine ID", "num_priority_screened_species": "Priority species count",
            "is_high_rwr": "Regional priority", "has_priority_screened_species_overlap": "Species priority"}, ["property_id"])
        portfolio = read(work / "14_conservation_portfolios/Investor_Country_EcologicalPortfolio_Summary_All.csv")
        portfolio["any"] = portfolio["high_rwr_only_mine_equivalent"] + portfolio["priority_screened_species_only_mine_equivalent"] + portfolio["dual_ecological_attention_mine_equivalent"]
        portfolio["any_share"] = portfolio["any"] / portfolio["total_traceable_mine_equivalent"] * 100
        mapping = {"investor_iso3": "Owner registration iso3", "num_unique_mines": "Num unique mines", "total_traceable_mine_equivalent": "Total traceable mine equivalent",
            "ecological_attention_unique_mines": "Conservation priority unique mines", "any": "Mine equivalents meeting either priority criterion",
            "any_share": "Portfolio meeting either priority criterion (%)"}
        for source, target in (("high_rwr_only", "High rwr only"), ("priority_screened_species_only", "Priority species only"),
                               ("dual_ecological_attention", "Dual conservation priority"), ("neither", "Neither")):
            mapping[f"{source}_mine_equivalent"] = f"{target} mine equivalent"
            mapping[f"{source}_portfolio_share_percent"] = f"{target} portfolio share (%)"
        for source, target in (("domestic", "Domestic"), ("foreign", "Cross border")):
            mapping[f"{source}_unique_mines"] = f"{target} unique mines"
            mapping[f"{source}_ecological_attention_unique_mines"] = f"{target} conservation priority unique mines"
            mapping[f"{source}_traceable_mine_equivalent"] = f"{target} traceable mine equivalent"
        add(portfolio, 8, "Owner_Portfolios", mapping, ["investor_iso3"])
    if "distance" in steps:
        actual = read(work / "15_buffer_sensitivity/comparison/Distance_Diagnostics.csv")
        add(actual, 9, "Distance_Diagnostics", {
            "start_distance_km": "Start distance (km)", "end_distance_km": "End distance (km)",
            "priority_species_start": "Priority species start", "priority_species_end": "Priority species end",
            "affected_species_start": "Affected species start", "affected_species_end": "Affected species end",
            "full_buffer_union_area_start_km2": "Full buffer union area start (km²)", "full_buffer_union_area_end_km2": "Full buffer union area end (km²)",
            "additional_priority_species": "Additional priority species", "additional_full_buffer_union_area_km2": "Additional full buffer union area (km²)",
            "additional_priority_species_per_100000_km2": "Additional priority species per 100000 (km²)",
            "priority_species_share_start_pct": "Priority species share start (%)", "priority_species_share_end_pct": "Priority species share end (%)",
            "priority_species_share_change_percentage_points": "Priority species share change percentage points"}, ["start_distance_km"])
    return checks
