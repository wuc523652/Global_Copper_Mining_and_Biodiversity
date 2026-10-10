"""Fast checks using real species-derived data and a synthetic ownership example.

This demonstrates selected scientific calculations, not a full spatial reanalysis.
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name, alias):
    spec = importlib.util.spec_from_file_location(alias, ROOT / "analysis" / name)
    if spec is None or spec.loader is None:
        raise ImportError(name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def run(output):
    started = time.perf_counter()
    os.environ["NC_WORK_ROOT"] = str(output)
    os.environ["NC_DATA_ROOT"] = str(ROOT / "inputs")
    os.environ["NC_CAP_CHECK"] = "0"
    os.environ.pop("COPPER_BIODIV_STAGING_ROOT", None)
    output.mkdir(parents=True, exist_ok=True)

    rules = load("08_compare_priority_rules.py", "demo_rules")
    rules.INPUT_FILE = ROOT / "demo/data/exposed_species.csv"
    rules.IMPACT_UNIT = "percent"
    rules.STRICT_EXPECTED_COUNT = True
    species, _ = rules.read_affected_species()
    results = rules.calculate_results(species)
    current = set(results.loc[results[f"{rules.MAIN}_Priority"], "Species_ID"])
    alternative = set(results.loc[results[f"{rules.TRAD}_Priority"], "Species_ID"])
    expected_ids = set(pd.read_csv(ROOT / "demo/data/expected_priority_ids.csv", dtype=str)["species_id"])
    check(len(species) == 891, "Expected 891 exposed species")
    check(current == expected_ids, "Main priority set differs from the 203-species manuscript set")
    check(len(alternative) == 247 and len(current & alternative) == 202, "Alternative rule comparison differs")
    results.to_csv(output / "species_rule_results.csv", index=False)
    rules.comparison_summary(results).to_csv(output / "species_rule_comparison.csv", index=False)

    # Check the same main rule implemented in the original species-analysis source.
    priority = load("07_screen_priority_species.py", "demo_priority")
    thresholds = np.array([priority.calculate_conflict_target(a) for a in species["Range_Area_km2"]])
    decisions = np.array([priority.determine_priority({"Cal_IAP": impact / 100, "Conflict_Target": limit})
                          for impact, limit in zip(species["Impact_Percent"], thresholds)])
    check(set(species.loc[decisions == "Yes", "Species_ID"]) == current, "Primary and comparison implementations disagree")
    check(priority.determine_priority({"Cal_IAP": 0.7, "Conflict_Target": 0.7}) == "No", "Equality is not priority")
    np.testing.assert_allclose([priority.calculate_protection_target(a) for a in [200, 250000]], [1, .3])

    attribution = load("12_attribute_exposure_by_ownership.py", "demo_attribution")
    attribution.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # Illustrative inputs only. No commercial owner identities or raw shares are included.
    master = pd.DataFrame({"property_id": ["demo_1", "demo_2"], "property_name": ["Synthetic 1", "Synthetic 2"],
                           "resource_country": ["Host A"] * 2, "resource_iso3": ["AAA"] * 2,
                           "resource_continent": ["Example"] * 2})
    coverage = pd.DataFrame({"property_id": ["demo_1", "demo_2"], "explicit_reported_share_pct": [60., 100.],
        "explicit_country_resolved_share_pct": [60., 100.], "explicit_country_unresolved_share_pct": [0., 0.],
        "unreported_share_pct": [40., 0.], "total_unattributed_share_pct": [40., 0.], "ownership_source": ["synthetic"] * 2})
    exposure = pd.DataFrame({"property_id": ["demo_1", "demo_2"], "mine_exposure_area_km2": [100., 100.],
                             "num_exposed_species": [1, 1], "buffer_area_km2": [200., 200.]})
    countries = pd.DataFrame({"resource_iso3": ["AAA"], "physical_exposure_union_area_km2": [150.],
                              "mine_exposure_area_sum_km2": [200.]})
    owners = pd.DataFrame({"property_id": ["demo_1", "demo_2"], "reported_share_pct": [60., 100.],
        "investor_iso3": ["HKG", "AAA"], "investor_country": ["Hong Kong", "Host A"],
        "investor_continent": ["Asia", "Example"], "owner_name": ["Synthetic owner 1", "Synthetic owner 2"]})
    with contextlib.redirect_stdout(io.StringIO()):
        base, factors = attribution.build_mine_base(master, coverage, exposure, countries)
        allocated, balance, _ = attribution.build_investor_attribution(base, owners)
    np.testing.assert_allclose(base["attribution_base_area_km2"], [75., 75.])
    np.testing.assert_allclose(balance["traceable_attributed_area_km2"], [45., 75.])
    np.testing.assert_allclose(balance["unattributed_area_km2"], [30., 0.])
    check("HKG" in set(allocated["investor_iso3"]), "HKG must remain in calculations")
    check(np.isclose(factors["mine_overlap_balance_factor"].iloc[0], .75), "Country union balancing failed")

    randomization = load("06_country_stratified_randomization.py", "demo_randomization")
    randomization.ACTIVE_N_ITERATIONS = 99
    country_counts = pd.DataFrame({"Country": ["Host A", "Host B"], "mine_count": [2, 1]})
    pools = {"Host A": {"is_top10": np.array([0, 1]), "rwr_percentile": np.array([20., 95.])},
             "Host B": {"is_top10": np.array([1]), "rwr_percentile": np.array([99.])}}
    with contextlib.redirect_stdout(io.StringIO()):
        first = randomization.run_randomizations(country_counts, pools)
        second = randomization.run_randomizations(country_counts, pools)
    pd.testing.assert_frame_equal(first, second)
    check(first["top10_mine_count"].between(1, 3).all(), "Country-stratified counts changed")
    stats = randomization.empirical_stats(100., np.array([0., 50., 75.]))
    check(stats["p_upper_one_sided"] == .25, "Empirical test must include plus-one correction")

    summary = {"status": "PASS", "scope": "real derived species rule validation; synthetic ownership and randomization checks",
        "full_spatial_reanalysis": False, "exposed_species": len(species), "main_priority_species": len(current),
        "alternative_priority_species": len(alternative), "shared_priority_species": len(current & alternative),
        "synthetic_country_union_area_km2": 150., "synthetic_traceable_area_km2": 120.,
        "synthetic_unattributed_area_km2": 30., "hkg_retained": True, "seed_reproducible": True,
        "elapsed_seconds": round(time.perf_counter() - started, 3)}
    (output / "demo_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "demo_output")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output == ROOT or output.is_relative_to(ROOT / "analysis") or output.is_relative_to(ROOT / "demo"):
        parser.error("Use a separate output directory, not package source or demonstration data")
    if (output / "0code").is_dir() or (output / ".git").is_dir():
        parser.error("Use a separate output directory, not the original project")
    run(output)


if __name__ == "__main__":
    main()
