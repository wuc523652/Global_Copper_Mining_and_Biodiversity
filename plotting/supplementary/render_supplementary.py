"""Reproduce Supplementary Figures 1–9 from the accompanying source tables.

Python/Matplotlib only. These are deterministic summaries; no sampling,
uncertainty estimates, or new statistical analyses are performed here.
"""
from __future__ import annotations

import argparse
import sys
sys.dont_write_bytecode = True
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import supplementary_1
import supplementary_2
import supplementary_3
import supplementary_4
import supplementary_5
import supplementary_9
import supplementary_6
import supplementary_7
import supplementary_8

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parent))
from data_paths import external_output_root, resolve_data_root

FIGURE_INPUTS = None
ZENODO = None

COMMON_FIELDS = {
    "Buffer distance (km)": "buffer_distance_km",
    "Affected species count": "affected_species_count",
    "Priority species count": "priority_species_count",
    "Num mines": "num_mines",
    "Global strict physical exposure area (km²)": "global_strict_physical_exposure_area_km2",
    "Domestic equity associated area (km²)": "domestic_equity_associated_area_km2",
    "Cross border equity associated area (km²)": "foreign_equity_associated_area_km2",
    "Unattributed area (km²)": "unattributed_area_km2",
    "Owner registration iso3": "investor_iso3",
    "Owner registration country": "investor_country",
    "Owner registration country or region": "investor_country",
    "Mine host iso3": "resource_iso3",
    "Mine host country": "resource_country",
    "Main rank": "main_rank",
    "Complete case rank": "complete_case_rank",
}


def read_table(tables, number, sheet, mapping=None):
    frame = pd.read_excel(tables / f"Supplementary_Data_{number}.xlsx", sheet_name=sheet)
    frame = frame.rename(columns={**COMMON_FIELDS, **(mapping or {})})
    if frame.columns.duplicated().any():
        raise ValueError(f"Ambiguous duplicate fields in Data {number}/{sheet}")
    return frame


def require_finite(frame, columns):
    values = frame[columns].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Missing or non-finite values in {columns}")


def prepare_1(tables):
    import geopandas as gpd
    grid = gpd.read_file(ZENODO / "3Grid_and_Biodiversity" / "Global_10km_Grid_Biodiversity.gpkg",
                         columns=["grid_id", "total_richness"])
    if grid.crs is None:
        raise ValueError("The biodiversity grid must have a defined CRS")
    grid = grid.to_crs(epsg=4326)
    if len(grid) != 1404410 or not grid["grid_id"].is_unique:
        raise ValueError("Unexpected global biodiversity grid")
    require_finite(grid, ["total_richness"])
    world = gpd.read_file(FIGURE_INPUTS / "world_s1.gpkg")
    return grid, world


def prepare_2(tables):
    data = read_table(tables, 3, "Country_Spatial", {
        "Mine host country": "Country", "Mine count": "conflict_mine_count",
        "Unioned exposed range area (km²)": "overlap_area_km2",
        "Rwr 0 50 area (km²)": "area_le50_km2",
        "Rwr 50 75 area (km²)": "area_50_75_km2",
        "Rwr 75 90 area (km²)": "area_75_90_km2",
        "Rwr 90 99 area (km²)": "area_90_99_km2",
        "Rwr 99 100 area (km²)": "area_99_100_km2",
        "Classified influence area in top-10% RWR cells (%)": "pct_ge90",
    })
    require_finite(data, supplementary_2.BIN_FIELDS + ["overlap_area_km2", "pct_ge90", "conflict_mine_count"])
    return (data,)


def prepare_3(tables):
    country = read_table(tables, 3, "Randomization_Countries", {
        "Mine host country": "Country", "Mine count": "mine_count",
        "Top10 mine (%)": "observed_top10_mine_pct",
        "Background top10 grid (%)": "background_top10_grid_pct",
        "Top10 observed minus background (percentage points)": "top10_pct_difference_observed_minus_background",
    })
    iterations = read_table(tables, 3, "Randomization_Iterations", {
        "Mean global RWR percentile": "mean_rwr_percentile",
    })
    summary = read_table(tables, 3, "Randomization_Global", {
        "Metric": "metric", "Observed value": "observed",
        "Randomization 2.5th percentile": "null_q025",
        "Randomization 97.5th percentile": "null_q975",
        "Standardized effect size": "ses", "Upper-tail empirical probability": "p_upper_one_sided",
    })
    require_finite(country, ["mine_count", "observed_top10_mine_pct", "background_top10_grid_pct"])
    require_finite(iterations, ["mean_rwr_percentile"])
    return country, iterations, summary


def prepare_4(tables):
    data = read_table(tables, 4, "Exposed_Species", {
        "Mapped global range area (km²)": "total_habitat_area_km2",
        "Overlapping range (%)": "Impact_Percent",
        "Main-rule priority species": "Priority_ForPlot",
    })
    require_finite(data, ["total_habitat_area_km2", "Impact_Percent"])
    if len(data) != 891 or data["IUCN species ID"].duplicated().any():
        raise ValueError("Unexpected affected-species records")
    return (data,)


def prepare_5(tables):
    data = read_table(tables, 4, "Rule_Comparison", {
        "IUCN species ID": "Species_ID", "Taxonomic group": "Taxon",
        "Mapped global range area (km²)": "Range_Area_km2",
        "Overlapping range (%)": "Impact_Percent",
        "Main-rule priority species": f"{supplementary_5.MAIN}_Priority",
        "Alternative-rule priority species": f"{supplementary_5.ALT}_Priority",
    })
    require_finite(data, ["Range_Area_km2", "Impact_Percent"])
    data["Taxon"] = data["Taxon"].replace({"Freshwater fishes": "Fish"})
    area = data["Range_Area_km2"].to_numpy()
    data["Current_Range_Class"] = np.where(area <= 200, "Narrow", np.where(area >= 250000, "Widespread", "Intermediate"))
    for prefix, scheme in supplementary_5.SCHEMES.items():
        column = f"{prefix}_Priority"
        data[column] = data[column].astype(str).str.lower().isin(["true", "1", "yes"])
        expected = data["Impact_Percent"].to_numpy() > (supplementary_5.screening_threshold(area, scheme, apply_cap=False) + 1e-10)
        if not np.array_equal(data[column], expected):
            raise ValueError(f"Stored screening classification differs from {prefix}")
    main = data[f"{supplementary_5.MAIN}_Priority"]
    alternative = data[f"{supplementary_5.ALT}_Priority"]
    comparison = {
        "Stable_Core_N": int((main & alternative).sum()),
        "Current_Only_N": int((main & ~alternative).sum()),
        "Traditional_Only_N": int((~main & alternative).sum()),
        "Jaccard_Index": float((main & alternative).sum() / (main | alternative).sum()),
    }
    if len(data) != 891 or comparison["Stable_Core_N"] != 202:
        raise ValueError("Unexpected rule-comparison data")
    return data, comparison


def prepare_9(tables):
    return ()


def prepare_6(tables):
    import geopandas as gpd
    world = gpd.read_file(FIGURE_INPUTS / "world_s6.gpkg")
    if world.crs.to_epsg() != 4326 or world["iso3"].duplicated().any():
        raise ValueError("Expected original unique EPSG:4326 map boundaries")
    hosts = read_table(tables, 6, "Mine_Host_Countries", {
        "Physical exposure union area (km²)": "physical_exposure_union_area_km2",
        "Cross border share of traceable (%)": "foreign_share_of_traceable_pct",
        "Cross border attributed area (km²)": "foreign_attributed_area_km2",
    })
    owners = read_table(tables, 6, "Owner_Countries", {
        "Cross border attributed area (km²)": "foreign_project_attributed_area_km2",
    })
    require_finite(hosts, ["physical_exposure_union_area_km2", "foreign_attributed_area_km2"])
    require_finite(owners, ["foreign_project_attributed_area_km2"])
    if owners["investor_iso3"].duplicated().any() or hosts["resource_iso3"].duplicated().any():
        raise ValueError("Duplicated registration or host identifiers")
    ordered = owners.sort_values("foreign_project_attributed_area_km2", ascending=False)
    hkg = ordered.loc[ordered["investor_iso3"].eq("HKG")].iloc[0]
    assert ordered.iloc[7]["investor_iso3"] == "HKG"
    assert hkg["investor_country"] == "Hong Kong, China"
    assert np.isclose(hkg["foreign_project_attributed_area_km2"], 21344.260108647366)
    return world, hosts, owners


def prepare_7(tables):
    scenarios = read_table(tables, 6, "Completeness_Scenarios")
    if scenarios["Scenario"].tolist() != ["All mines", "Complete case"]:
        raise ValueError("Expected all-mines then complete-case scenarios")
    owners = read_table(tables, 6, "Owner_Rank_Comparison")
    pairs = read_table(tables, 7, "Pair_Rank_Comparison")
    pairs["flow_key"] = pairs["investor_iso3"] + "->" + pairs["resource_iso3"]
    for frame in (owners, pairs):
        require_finite(frame, ["main_rank", "complete_case_rank"])
    top_owners = owners.nsmallest(20, "main_rank")
    top_pairs = pairs.nsmallest(20, "main_rank")
    hkg = top_owners.loc[top_owners["investor_iso3"].eq("HKG")].iloc[0]
    hkg_chn = top_pairs.loc[top_pairs["flow_key"].eq("HKG->CHN")].iloc[0]
    assert (hkg["main_rank"], hkg["complete_case_rank"]) == (8, 10)
    assert (hkg_chn["main_rank"], hkg_chn["complete_case_rank"]) == (6, 26)
    assert "BMU" not in top_owners["investor_iso3"].tolist()
    assert "AUS->PNG" not in top_pairs["flow_key"].tolist()
    return scenarios, owners, pairs


def prepare_8(tables):
    core = read_table(tables, 9, "Core_Results")
    stability = read_table(tables, 9, "Set_Rank_Stability", {
        "Priority set jaccard vs 30km": "priority_set_jaccard_vs_30km",
        "Owner registration country rank correlation vs 30km": "investor_rank_correlation_vs_30km",
        "Owner registration country or region rank correlation vs 30km": "investor_rank_correlation_vs_30km",
    })
    diagnostics = read_table(tables, 9, "Distance_Diagnostics", {
        "Additional priority species per 100000 (km²)": "additional_priority_species_per_100000_km2",
        "Priority species share change percentage points": "priority_species_share_change_percentage_points",
    })
    for frame in (core, stability):
        if frame["buffer_distance_km"].tolist() != [10, 20, 30, 40, 50]:
            raise ValueError("Distances must be exactly 10, 20, 30, 40, 50 km in order")
    if diagnostics["Start distance (km)"].tolist() != [10, 20, 30, 40]:
        raise ValueError("Unexpected buffer intervals")
    n = core["affected_species_count"].to_numpy()
    p = core["priority_species_count"].to_numpy()
    assert (n > 0).all() and (p >= 0).all() and (p <= n).all()
    increments = diagnostics["Additional full buffer union area (km²)"].to_numpy()
    assert (increments > 0).all()
    assert np.allclose(np.diff(p) / increments * 100000,
                       diagnostics["additional_priority_species_per_100000_km2"])
    assert np.allclose(np.diff(p / n * 100),
                       diagnostics["priority_species_share_change_percentage_points"])
    require_finite(stability, ["priority_set_jaccard_vs_30km", "investor_rank_correlation_vs_30km"])
    return core, stability, diagnostics


def style_6():
    plt.rcParams.update({
        "font.family": "Times New Roman", "font.size": 6.2,
        "axes.labelsize": 6.5, "axes.titlesize": 7.0,
        "xtick.labelsize": 5.8, "ytick.labelsize": 5.8, "legend.fontsize": 5.8,
        "axes.linewidth": 0.65, "xtick.major.width": 0.55,
        "ytick.major.width": 0.55, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.right": False, "axes.spines.top": False,
        "axes.unicode_minus": False, "legend.frameon": False,
    })


# Published canvas crops in inches. These are layout coordinates, not image inputs.
EXPORT_CROPS = {2: [0.8606958876722428, 0.15769837499999886, 6.66667, 11.873336666666667], 3: [0.2510968678915141, 0.121488472222222, 6.910003333333334, 6.03667], 4: [0.346451608705162, 0.22900833333333362, 6.813336666666666, 4.573336666666666], 5: [0.24096432390395653, 0.2591432624113468, 7.003333333333333, 6.323333333333333]}


def main():
    global FIGURE_INPUTS, ZENODO
    from matplotlib.transforms import Bbox
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", nargs="+", type=int, choices=range(1, 10), default=list(range(1, 10)))
    parser.add_argument("--data-root", type=Path, help="Dataset root; overrides COPPER_DATA_ROOT and bundled data/")
    parser.add_argument("--output", "--output-dir", type=Path, required=True, help="Directory outside both downloaded packages")
    parser.add_argument("--formats", nargs="+", choices=["png", "pdf", "svg"], default=["png"])
    args = parser.parse_args()
    try:
        ZENODO = resolve_data_root(args.data_root)
        if ZENODO is None:
            raise ValueError("Supply --data-root or set COPPER_DATA_ROOT")
        args.output = external_output_root(args.output, ZENODO)
    except ValueError as error:
        parser.error(str(error))
    FIGURE_INPUTS = ZENODO / "9Figure_Inputs" / "Supplementary"
    tables = ZENODO / "8Final_Results"
    args.output.mkdir(parents=True, exist_ok=True)
    for number in dict.fromkeys(args.figures):
        with plt.rc_context():
            plt.rcdefaults()
            if number in [1, 2, 3, 4, 6]:
                style_6()
            plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none"})
            inputs = globals()[f"prepare_{number}"](tables)
            figure = globals()[f"supplementary_{number}"].build_plot(*inputs)
            figure.canvas.draw()
            if number in [1, 6]:
                export_options = {"bbox_inches": "tight", "pad_inches": 0.06}
            elif number in EXPORT_CROPS:
                export_options = {"bbox_inches": Bbox.from_bounds(*EXPORT_CROPS[number]), "pad_inches": 0}
            elif number == 7:
                export_options = {"bbox_inches": "tight"}
            else:
                export_options = {}
            dpi = 600 if number == 9 else 300 if number in [2, 3, 4, 5] else 400
            stem = f"Supplementary_Figure_{number}"
            for extension in dict.fromkeys(args.formats):
                figure.savefig(args.output / f"{stem}.{extension}", dpi=dpi, **export_options)
            plt.close(figure)
            print(f"Rendered {stem}: {', '.join(dict.fromkeys(args.formats))}")


if __name__ == "__main__":
    main()
