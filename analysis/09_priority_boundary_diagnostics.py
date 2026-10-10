# -*- coding: utf-8 -*-
"""
Create transparent diagnostics for the study-defined priority screening rule.

This script does not change the 203-species selection. It makes the algebra and
decision margin explicit without treating a 30 km analysis-zone intersection as
observed habitat loss or as an established mining-impact threshold.
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(_nc_project(""))
INPUT_FILE = (
    BASE_DIR
    / "07_priority_species"
    / "Species_PrioritySpecies_Full.csv"
)

OUTPUT_DIR = BASE_DIR / "09_priority_diagnostics"

OUTPUT_FULL = OUTPUT_DIR / "Species_PriorityScreening_Diagnostics_Full.csv"
OUTPUT_PRIORITY = OUTPUT_DIR / "Species_PriorityScreening_Diagnostics_203.csv"
OUTPUT_SUMMARY = OUTPUT_DIR / "PriorityScreening_Diagnostics_Summary.csv"
OUTPUT_QA = OUTPUT_DIR / "PriorityScreening_Diagnostics_QA.csv"

REQUIRED_COLUMNS = [
    "species_id",
    "scientific_name",
    "taxon",
    "status",
    "Range",
    "total_habitat_area_km2",
    "total_overlap_area_km2",
    "Impact_Percent",
    "Protection_Target_Percent",
    "Conflict_Target_Percent",
    "Priority",
]

NUMERIC_COLUMNS = [
    "total_habitat_area_km2",
    "total_overlap_area_km2",
    "Impact_Percent",
    "Protection_Target_Percent",
    "Conflict_Target_Percent",
]

# Source script uses EPS=1e-12 on proportions; convert it to percentage points.
SCREENING_EPS_PP = 1e-10


def normalize_priority(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("yes")


def require_columns(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(
            f"Input is missing required columns {missing}: {INPUT_FILE}"
        )


def median_or_nan(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.median()) if not values.empty else np.nan


def build_summary(frame: pd.DataFrame) -> pd.DataFrame:
    group_specs = [
        ("All", []),
        ("Range", ["Range"]),
        ("Taxon", ["taxon"]),
        ("Taxon x range", ["taxon", "Range"]),
    ]
    rows = []

    for scope, columns in group_specs:
        groups = [((), frame)] if not columns else frame.groupby(columns, dropna=False)
        for keys, group in groups:
            if not isinstance(keys, tuple):
                keys = (keys,)
            priority = group[group["Priority_Recomputed"]]
            row = {
                "scope": scope,
                "taxon": "All",
                "Range": "All",
                "candidate_species_n": int(group["species_id"].nunique()),
                "analysis_zone_intersecting_species_n": int(
                    (group["Analysis_Zone_Coverage_Percent"] > 0).sum()
                ),
                "priority_screened_species_n": int(
                    priority["species_id"].nunique()
                ),
                "priority_share_percent": (
                    float(priority["species_id"].nunique())
                    / float(group["species_id"].nunique())
                    * 100.0
                ),
                "median_analysis_zone_coverage_percent_priority": median_or_nan(
                    priority["Analysis_Zone_Coverage_Percent"]
                ),
                "median_threshold_exceedance_pp_priority": median_or_nan(
                    priority["Threshold_Exceedance_pp"]
                ),
                "median_representation_margin_pp_priority": median_or_nan(
                    priority["Representation_Margin_pp"]
                ),
            }
            for column, key in zip(columns, keys):
                row[column] = key
            rows.append(row)

    return pd.DataFrame(rows)


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input does not exist: {INPUT_FILE}")

    frame = pd.read_csv(INPUT_FILE)
    require_columns(frame)
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    if frame[NUMERIC_COLUMNS].isna().any().any():
        bad_columns = frame[NUMERIC_COLUMNS].columns[
            frame[NUMERIC_COLUMNS].isna().any()
        ].tolist()
        raise ValueError(f"Non-numeric or missing diagnostic values in {bad_columns}")

    frame["Analysis_Zone_Coverage_Percent"] = frame["Impact_Percent"]
    frame["Range_Outside_Zone_Percent"] = (
        100.0 - frame["Analysis_Zone_Coverage_Percent"]
    )
    frame["Threshold_Exceedance_pp"] = (
        frame["Analysis_Zone_Coverage_Percent"]
        - frame["Conflict_Target_Percent"]
    )
    frame["Representation_Margin_pp"] = (
        frame["Range_Outside_Zone_Percent"]
        - frame["Protection_Target_Percent"]
    )
    frame["Priority_Recomputed"] = (
        frame["Threshold_Exceedance_pp"] > SCREENING_EPS_PP
    )
    frame["Priority_Recorded"] = normalize_priority(frame["Priority"])
    frame["Priority_Rule_Equivalent"] = (
        "analysis-zone coverage > (100 - representation target)"
    )
    frame["Interpretation_Boundary"] = (
        "screening under a conservative non-coexistence assumption; "
        "not observed habitat loss or an established mining-impact threshold"
    )

    identity_error = (
        frame["Representation_Margin_pp"] + frame["Threshold_Exceedance_pp"]
    ).abs()
    mismatch = frame["Priority_Recomputed"] != frame["Priority_Recorded"]

    qa = pd.DataFrame(
        [
            {"check": "input_species_rows", "value": int(len(frame)), "status": "INFO"},
            {
                "check": "input_unique_species",
                "value": int(frame["species_id"].nunique()),
                "status": "INFO",
            },
            {
                "check": "recorded_priority_species",
                "value": int(frame.loc[frame["Priority_Recorded"], "species_id"].nunique()),
                "status": "INFO",
            },
            {
                "check": "recomputed_priority_species",
                "value": int(frame.loc[frame["Priority_Recomputed"], "species_id"].nunique()),
                "status": "PASS" if not mismatch.any() else "FAIL",
            },
            {
                "check": "priority_classification_mismatches",
                "value": int(mismatch.sum()),
                "status": "PASS" if not mismatch.any() else "FAIL",
            },
            {
                "check": "max_algebraic_identity_error_pp",
                "value": float(identity_error.max()),
                "status": "PASS" if float(identity_error.max()) <= 1e-8 else "FAIL",
            },
        ]
    )

    if (qa["status"] == "FAIL").any():
        raise RuntimeError(
            "Priority screening diagnostic QA failed; outputs were not written."
        )

    output_columns = [
        "species_id",
        "scientific_name",
        "taxon",
        "status",
        "Range",
        "total_habitat_area_km2",
        "total_overlap_area_km2",
        "Analysis_Zone_Coverage_Percent",
        "Range_Outside_Zone_Percent",
        "Protection_Target_Percent",
        "Conflict_Target_Percent",
        "Threshold_Exceedance_pp",
        "Representation_Margin_pp",
        "Priority_Recorded",
        "Priority_Recomputed",
        "Priority_Rule_Equivalent",
        "Interpretation_Boundary",
    ]
    output_columns = [column for column in output_columns if column in frame.columns]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame[output_columns].to_csv(OUTPUT_FULL, index=False, encoding="utf-8-sig")
    frame.loc[frame["Priority_Recomputed"], output_columns].to_csv(
        OUTPUT_PRIORITY, index=False, encoding="utf-8-sig"
    )
    build_summary(frame).to_csv(OUTPUT_SUMMARY, index=False, encoding="utf-8-sig")
    qa.to_csv(OUTPUT_QA, index=False, encoding="utf-8-sig")

    print(f"Input species: {frame['species_id'].nunique():,}")
    print(f"Priority-screened species: {frame['Priority_Recomputed'].sum():,}")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
