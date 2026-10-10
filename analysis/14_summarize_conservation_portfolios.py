# -*- coding: utf-8 -*-
"""
14_summarize_conservation_portfolios.py

Aggregate the four mine-level ecological-attention categories from step 13
to investor-country, resource-country and bilateral portfolios.

Two explicitly separated accounting units are retained:

1. equity-weighted mine equivalents for the four-category portfolio
   composition; and
2. unique mines for the ecological-attention counts shown in Fig. 5.

The equity-weighted unit is:

    reported ownership share (%) / 100

Only explicitly reported ownership shares that can be resolved to an investor
country are attributed. Missing, unreported or country-unresolved shares are
not redistributed. Consequently, portfolio shares describe the traceable
country-resolved portion of each investor portfolio, while the global table
also reports unattributed equity. All country-resolved records, including HKG,
enter the calculations and global totals; HKG is omitted only from the
separately generated named-country display selection.

The input tables from steps 10 and 13 are read without modification.
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings

import numpy as np
import pandas as pd


warnings.filterwarnings("ignore")


BASE_DIR = Path(_nc_project(""))
OUTPUT_DIR = BASE_DIR / "14_conservation_portfolios"

MINE_BRIDGE = (
    BASE_DIR
    / "13_mine_species_links"
    / "Mine_EcologicalAttention_30km.csv"
)
OWNER_LEDGER = (
    BASE_DIR
    / "10_mine_ownership"
    / "Mine_Owner_Share_Ledger.csv"
)
OWNERSHIP_COVERAGE = (
    BASE_DIR
    / "10_mine_ownership"
    / "Mine_Ownership_Coverage.csv"
)

OUT_LEDGER = OUTPUT_DIR / "Mine_InvestorCountry_EcologicalLedger_All.csv"
OUT_INVESTOR = OUTPUT_DIR / "Investor_Country_EcologicalPortfolio_Summary_All.csv"
OUT_RESOURCE = OUTPUT_DIR / "Resource_Country_EcologicalPortfolio_Summary_All.csv"
OUT_BILATERAL = OUTPUT_DIR / "Bilateral_EcologicalPortfolio_Summary_All.csv"
OUT_INVESTOR_DISPLAY = OUTPUT_DIR / "Investor_Country_EcologicalPortfolio_NamedDisplay.csv"
OUT_RESOURCE_DISPLAY = OUTPUT_DIR / "Resource_Country_EcologicalPortfolio_NamedDisplay.csv"
OUT_BILATERAL_DISPLAY = OUTPUT_DIR / "Bilateral_EcologicalPortfolio_NamedDisplay.csv"
OUT_GLOBAL = OUTPUT_DIR / "Global_EcologicalPortfolio_Summary.csv"
OUT_SELECTION = OUTPUT_DIR / "Fig5_Display_Selection.csv"
OUT_QA = OUTPUT_DIR / "Fig5_EquityAndMineCount_QA.csv"
OUT_METHOD = OUTPUT_DIR / "Fig5_Method_Definition.md"
OUT_DISPLAY_OMISSION_AUDIT = OUTPUT_DIR / "NamedDisplay_Omission_Audit.csv"

EXPECTED_MINES = 1278
SHARE_TOLERANCE = 1e-6
MIN_PORTFOLIO_MINES_FIG5A = 10
N_FIG5A_INVESTORS = 12
N_FIG5B_INVESTORS = 8
N_FIG5B_RESOURCES = 10
DISPLAY_OMITTED_ISO3 = {"HKG"}

CATEGORY_ORDER = [
    "Neither",
    "High RWR only",
    "Priority-screened species only",
    "Dual ecological attention",
]
CATEGORY_SLUG = {
    "Neither": "neither",
    "High RWR only": "high_rwr_only",
    "Priority-screened species only": "priority_screened_species_only",
    "Dual ecological attention": "dual_ecological_attention",
}


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input not found: {path}")


def normalise_id(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.strip()
    return values.str.replace(r"\.0$", "", regex=True)


def first_notna(series: pd.Series):
    values = series.dropna()
    return values.iloc[0] if len(values) else np.nan


def safe_percent(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    numerator = pd.to_numeric(numerator, errors="coerce").fillna(0.0)
    denominator = pd.to_numeric(denominator, errors="coerce").fillna(0.0)
    return np.where(denominator > 0, numerator / denominator * 100.0, 0.0)


def add_category_columns(
    base: pd.DataFrame,
    data: pd.DataFrame,
    group_fields: list[str],
    value_field: str = "equity_weighted_mine_equivalent",
) -> pd.DataFrame:
    grouped = (
        data.groupby(
            group_fields + ["ecological_attention_category"],
            as_index=False,
            dropna=False,
            observed=False,
        )[value_field]
        .sum()
    )
    wide = grouped.pivot_table(
        index=group_fields,
        columns="ecological_attention_category",
        values=value_field,
        aggfunc="sum",
        fill_value=0.0,
        observed=False,
    ).reset_index()
    wide.columns.name = None
    rename = {
        category: f"{slug}_mine_equivalent"
        for category, slug in CATEGORY_SLUG.items()
    }
    wide = wide.rename(columns=rename)
    for column in rename.values():
        if column not in wide.columns:
            wide[column] = 0.0
    return base.merge(wide, on=group_fields, how="left", validate="one_to_one")


def portfolio_summary(
    data: pd.DataFrame,
    unit_prefix: str,
    name_field: str,
    iso_field: str,
    continent_field: str,
) -> pd.DataFrame:
    group_fields = [iso_field, name_field, continent_field]
    base = (
        data.groupby(group_fields, as_index=False, dropna=False)
        .agg(
            total_traceable_mine_equivalent=(
                "equity_weighted_mine_equivalent",
                "sum",
            ),
            num_unique_mines=("property_id", "nunique"),
        )
    )
    result = add_category_columns(base, data, group_fields)

    category_columns = [
        f"{CATEGORY_SLUG[category]}_mine_equivalent" for category in CATEGORY_ORDER
    ]
    result[category_columns] = result[category_columns].fillna(0.0)
    for category in CATEGORY_ORDER:
        slug = CATEGORY_SLUG[category]
        result[f"{slug}_portfolio_share_percent"] = safe_percent(
            result[f"{slug}_mine_equivalent"],
            result["total_traceable_mine_equivalent"],
        )

    result["high_rwr_related_mine_equivalent"] = (
        result["high_rwr_only_mine_equivalent"]
        + result["dual_ecological_attention_mine_equivalent"]
    )
    result["priority_species_related_mine_equivalent"] = (
        result["priority_screened_species_only_mine_equivalent"]
        + result["dual_ecological_attention_mine_equivalent"]
    )
    result["high_rwr_related_portfolio_share_percent"] = safe_percent(
        result["high_rwr_related_mine_equivalent"],
        result["total_traceable_mine_equivalent"],
    )
    result["priority_species_related_portfolio_share_percent"] = safe_percent(
        result["priority_species_related_mine_equivalent"],
        result["total_traceable_mine_equivalent"],
    )

    attention_counts = (
        data[~data["ecological_attention_category"].eq("Neither")]
        .groupby(group_fields, as_index=False, dropna=False)
        .agg(ecological_attention_unique_mines=("property_id", "nunique"))
    )
    result = result.merge(
        attention_counts, on=group_fields, how="left", validate="one_to_one"
    )

    domestic = data[data["is_cross_border"].eq(False)]
    foreign = data[data["is_cross_border"].eq(True)]
    for label, subset in [("domestic", domestic), ("foreign", foreign)]:
        totals = (
            subset.groupby(group_fields, as_index=False, dropna=False)
            .agg(
                **{
                    f"{label}_traceable_mine_equivalent": (
                        "equity_weighted_mine_equivalent",
                        "sum",
                    ),
                    f"{label}_unique_mines": ("property_id", "nunique"),
                }
            )
        )
        dual = (
            subset[
                subset["ecological_attention_category"].eq(
                    "Dual ecological attention"
                )
            ]
            .groupby(group_fields, as_index=False, dropna=False)
            .agg(
                **{
                    f"{label}_dual_ecological_attention_mine_equivalent": (
                        "equity_weighted_mine_equivalent",
                        "sum",
                    ),
                    f"{label}_dual_unique_mines": ("property_id", "nunique"),
                }
            )
        )
        result = result.merge(totals, on=group_fields, how="left")
        result = result.merge(dual, on=group_fields, how="left")
        attention = (
            subset[~subset["ecological_attention_category"].eq("Neither")]
            .groupby(group_fields, as_index=False, dropna=False)
            .agg(
                **{
                    f"{label}_ecological_attention_unique_mines": (
                        "property_id", "nunique"
                    )
                }
            )
        )
        result = result.merge(attention, on=group_fields, how="left")

    numeric = result.select_dtypes(include=[np.number]).columns
    result[numeric] = result[numeric].fillna(0.0)
    result[f"{unit_prefix}_iso3"] = result[iso_field]
    return result


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
    for path in [MINE_BRIDGE, OWNER_LEDGER, OWNERSHIP_COVERAGE]:
        require_file(path)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    bridge = pd.read_csv(MINE_BRIDGE, low_memory=False)
    owners = pd.read_csv(OWNER_LEDGER, low_memory=False)
    coverage = pd.read_csv(OWNERSHIP_COVERAGE, low_memory=False)

    bridge["property_id"] = normalise_id(bridge["property_id"])
    owners["property_id"] = normalise_id(owners["property_id"])
    coverage["property_id"] = normalise_id(coverage["property_id"])

    bridge_fields = [
        "property_id",
        "property_name",
        "resource_country",
        "resource_iso3",
        "resource_continent",
        "is_high_rwr",
        "has_priority_screened_species_overlap",
        "num_priority_screened_species",
        "ecological_attention_category",
    ]
    missing = [field for field in bridge_fields if field not in bridge.columns]
    if missing:
        raise ValueError(f"Mine bridge is missing columns: {missing}")

    mine_ids = set(bridge["property_id"])
    coverage_ids = set(coverage["property_id"])

    owners["reported_share_pct"] = pd.to_numeric(
        owners["reported_share_pct"], errors="coerce"
    )
    owners["has_reported_share"] = pd.to_numeric(
        owners["has_reported_share"], errors="coerce"
    ).fillna(0)
    owners["investor_country_resolved"] = pd.to_numeric(
        owners["investor_country_resolved"], errors="coerce"
    ).fillna(0)
    resolved = owners[
        owners["has_reported_share"].eq(1)
        & owners["investor_country_resolved"].eq(1)
        & owners["reported_share_pct"].gt(0)
        & owners["investor_iso3"].notna()
    ].copy()
    if resolved.empty:
        raise RuntimeError("No country-resolved reported ownership shares were found.")

    # Several company rows can map to the same investor country in one mine.
    # Combine them before attribution to make the country the statistical unit.
    country_ledger = (
        resolved.groupby(
            ["property_id", "investor_iso3"],
            as_index=False,
            dropna=False,
        )
        .agg(
            property_name=("property_name", first_notna),
            resource_country=("resource_country", first_notna),
            resource_iso3=("resource_iso3", first_notna),
            resource_continent=("resource_continent", first_notna),
            investor_country=("investor_country", first_notna),
            investor_continent=("investor_continent", first_notna),
            reported_share_pct=("reported_share_pct", "sum"),
            num_owner_records_combined=("owner_name", "size"),
            owner_names=(
                "owner_name",
                lambda x: "; ".join(sorted(set(x.dropna().astype(str)))),
            ),
        )
    )

    country_ledger = country_ledger.drop(
        columns=[
            "property_name",
            "resource_country",
            "resource_iso3",
            "resource_continent",
        ]
    ).merge(
        bridge[bridge_fields],
        on="property_id",
        how="left",
        validate="many_to_one",
    )
    country_ledger["equity_weighted_mine_equivalent"] = (
        country_ledger["reported_share_pct"] / 100.0
    )
    country_ledger["is_cross_border"] = ~country_ledger["investor_iso3"].eq(
        country_ledger["resource_iso3"]
    )
    country_ledger = country_ledger.sort_values(
        ["property_id", "investor_iso3"], kind="stable"
    ).reset_index(drop=True)

    display_omitted_country_ledger = country_ledger[
        country_ledger["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
        | country_ledger["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    display_omitted_country_ledger["display_rule"] = (
        "included in calculations and global totals; omitted only from named-country displays"
    )

    investor = portfolio_summary(
        country_ledger,
        unit_prefix="investor",
        name_field="investor_country",
        iso_field="investor_iso3",
        continent_field="investor_continent",
    )
    investor["eligible_fig5a_min_10_mines"] = investor["num_unique_mines"].ge(
        MIN_PORTFOLIO_MINES_FIG5A
    )
    investor = investor.sort_values(
        [
            "ecological_attention_unique_mines",
            "total_traceable_mine_equivalent",
            "investor_iso3",
        ],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    investor_display = investor[
        ~investor["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()

    resource = portfolio_summary(
        country_ledger,
        unit_prefix="resource",
        name_field="resource_country",
        iso_field="resource_iso3",
        continent_field="resource_continent",
    )
    resource = resource.sort_values(
        [
            "ecological_attention_unique_mines",
            "total_traceable_mine_equivalent",
            "resource_iso3",
        ],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    resource_display = resource[
        ~resource["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()

    bilateral_base = (
        country_ledger.groupby(
            [
                "investor_iso3",
                "investor_country",
                "investor_continent",
                "resource_iso3",
                "resource_country",
                "resource_continent",
                "is_cross_border",
            ],
            as_index=False,
            dropna=False,
        )
        .agg(
            total_traceable_mine_equivalent=(
                "equity_weighted_mine_equivalent",
                "sum",
            ),
            num_unique_mines=("property_id", "nunique"),
            ecological_attention_unique_mines=(
                "property_id",
                lambda x: x[
                    ~country_ledger.loc[
                        x.index, "ecological_attention_category"
                    ].eq("Neither")
                ].nunique(),
            ),
        )
    )
    bilateral_group_fields = [
        "investor_iso3",
        "investor_country",
        "investor_continent",
        "resource_iso3",
        "resource_country",
        "resource_continent",
        "is_cross_border",
    ]
    bilateral = add_category_columns(
        bilateral_base, country_ledger, bilateral_group_fields
    )
    bilateral["dual_portfolio_share_percent"] = safe_percent(
        bilateral["dual_ecological_attention_mine_equivalent"],
        bilateral["total_traceable_mine_equivalent"],
    )
    bilateral = bilateral.sort_values(
        [
            "is_cross_border",
            "dual_ecological_attention_mine_equivalent",
            "investor_iso3",
            "resource_iso3",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    bilateral_display = bilateral[
        ~bilateral["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
        & ~bilateral["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()

    coverage_numeric = coverage[
        [
            "property_id",
            "explicit_country_resolved_share_pct",
            "total_unattributed_share_pct",
        ]
    ].copy()
    coverage_numeric["explicit_country_resolved_share_pct"] = pd.to_numeric(
        coverage_numeric["explicit_country_resolved_share_pct"], errors="coerce"
    ).fillna(0.0)
    coverage_numeric["total_unattributed_share_pct"] = pd.to_numeric(
        coverage_numeric["total_unattributed_share_pct"], errors="coerce"
    ).fillna(0.0)
    mine_equity = (
        country_ledger.groupby("property_id", as_index=False)
        .agg(
            ledger_country_resolved_share_pct=("reported_share_pct", "sum"),
            traceable_country_resolved_mine_equivalent=(
                "equity_weighted_mine_equivalent",
                "sum",
            ),
        )
    )
    mine_balance = (
        bridge[["property_id", "ecological_attention_category"]]
        .merge(coverage_numeric, on="property_id", how="left", validate="one_to_one")
        .merge(mine_equity, on="property_id", how="left", validate="one_to_one")
    )
    mine_balance[
        [
            "ledger_country_resolved_share_pct",
            "traceable_country_resolved_mine_equivalent",
        ]
    ] = mine_balance[
        [
            "ledger_country_resolved_share_pct",
            "traceable_country_resolved_mine_equivalent",
        ]
    ].fillna(0.0)
    mine_balance["unattributed_mine_equivalent"] = (
        mine_balance["total_unattributed_share_pct"] / 100.0
    )
    mine_balance["equity_balance"] = (
        mine_balance["traceable_country_resolved_mine_equivalent"]
        + mine_balance["unattributed_mine_equivalent"]
    )

    global_rows = []
    for category in CATEGORY_ORDER:
        subset = mine_balance[
            mine_balance["ecological_attention_category"].eq(category)
        ]
        global_rows.append(
            {
                "ecological_attention_category": category,
                "mine_count": len(subset),
                "mine_share_percent": len(subset) / len(mine_balance) * 100.0,
                "traceable_country_resolved_mine_equivalent": subset[
                    "traceable_country_resolved_mine_equivalent"
                ].sum(),
                "unattributed_mine_equivalent": subset[
                    "unattributed_mine_equivalent"
                ].sum(),
                "traceable_equity_coverage_percent": subset[
                    "traceable_country_resolved_mine_equivalent"
                ].sum()
                / len(subset)
                * 100.0
                if len(subset)
                else 0.0,
            }
        )
    global_rows.append(
        {
            "ecological_attention_category": "Total",
            "mine_count": len(mine_balance),
            "mine_share_percent": 100.0,
            "traceable_country_resolved_mine_equivalent": mine_balance[
                "traceable_country_resolved_mine_equivalent"
            ].sum(),
            "unattributed_mine_equivalent": mine_balance[
                "unattributed_mine_equivalent"
            ].sum(),
            "traceable_equity_coverage_percent": mine_balance[
                "traceable_country_resolved_mine_equivalent"
            ].sum()
            / len(mine_balance)
            * 100.0,
        }
    )
    global_summary = pd.DataFrame(global_rows)

    fig5a = (
        investor_display[investor_display["eligible_fig5a_min_10_mines"]]
        .sort_values(
            [
                "ecological_attention_unique_mines",
                "num_unique_mines",
                "investor_iso3",
            ],
            ascending=[False, False, True],
            kind="stable",
        )
        .head(N_FIG5A_INVESTORS)
        .copy()
    )
    fig5a["selection_rank"] = np.arange(1, len(fig5a) + 1)

    fig5b_investor = (
        investor_display[
            investor_display["foreign_ecological_attention_unique_mines"].gt(0)
        ]
        .sort_values(
            [
                "foreign_ecological_attention_unique_mines",
                "investor_iso3",
            ],
            ascending=[False, True],
            kind="stable",
        )
        .head(N_FIG5B_INVESTORS)
        .copy()
    )
    fig5b_investor["selection_rank"] = np.arange(1, len(fig5b_investor) + 1)

    # Resource columns must be informative for the investor rows actually
    # displayed in Fig. 5b. Rank resources within those eight selected
    # investor portfolios; otherwise a globally top resource can form an
    # entirely empty matrix column when its investors are not displayed.
    fig5b_resource_population = (
        country_ledger[
            country_ledger["is_cross_border"].eq(True)
            & ~country_ledger["ecological_attention_category"].eq("Neither")
            & country_ledger["investor_iso3"].isin(fig5b_investor["investor_iso3"])
            & ~country_ledger["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
        ]
        .groupby(
            ["resource_iso3", "resource_country", "resource_continent"],
            as_index=False,
            dropna=False,
        )
        .agg(
            selected_investor_foreign_ecological_attention_unique_mines=(
                "property_id",
                "nunique",
            )
        )
    )
    fig5b_resource = (
        fig5b_resource_population[
            fig5b_resource_population[
                "selected_investor_foreign_ecological_attention_unique_mines"
            ].gt(0)
        ]
        .sort_values(
            [
                "selected_investor_foreign_ecological_attention_unique_mines",
                "resource_iso3",
            ],
            ascending=[False, True],
            kind="stable",
        )
        .head(N_FIG5B_RESOURCES)
        .copy()
    )
    fig5b_resource["selection_rank"] = np.arange(1, len(fig5b_resource) + 1)

    selection_rows = []
    for row in fig5a.itertuples(index=False):
        selection_rows.append(
            {
                "selection_group": "Fig5a investor portfolio",
                "selection_rank": row.selection_rank,
                "country_iso3": row.investor_iso3,
                "country_name": row.investor_country,
                "selection_metric": "ecological_attention_unique_mines",
                "selection_metric_value": row.ecological_attention_unique_mines,
                "eligibility_rule": "at least 10 unique mines with traceable country-resolved equity",
            }
        )
    for row in fig5b_investor.itertuples(index=False):
        selection_rows.append(
            {
                "selection_group": "Fig5b investor country",
                "selection_rank": row.selection_rank,
                "country_iso3": row.investor_iso3,
                "country_name": row.investor_country,
                "selection_metric": "foreign_ecological_attention_unique_mines",
                "selection_metric_value": row.foreign_ecological_attention_unique_mines,
                "eligibility_rule": "positive count of cross-border mines meeting either ecological-attention condition",
            }
        )
    for row in fig5b_resource.itertuples(index=False):
        selection_rows.append(
            {
                "selection_group": "Fig5b resource country",
                "selection_rank": row.selection_rank,
                "country_iso3": row.resource_iso3,
                "country_name": row.resource_country,
                "selection_metric": "selected_investor_foreign_ecological_attention_unique_mines",
                "selection_metric_value": row.selected_investor_foreign_ecological_attention_unique_mines,
                "eligibility_rule": "positive count of cross-border mines meeting either condition within the eight selected investor countries",
            }
        )
    selection = pd.DataFrame(selection_rows)

    category_columns = [
        f"{CATEGORY_SLUG[category]}_mine_equivalent" for category in CATEGORY_ORDER
    ]
    investor_category_error = (
        investor[category_columns].sum(axis=1)
        - investor["total_traceable_mine_equivalent"]
    ).abs()
    resource_category_error = (
        resource[category_columns].sum(axis=1)
        - resource["total_traceable_mine_equivalent"]
    ).abs()
    domestic_foreign_error = (
        investor["domestic_traceable_mine_equivalent"]
        + investor["foreign_traceable_mine_equivalent"]
        - investor["total_traceable_mine_equivalent"]
    ).abs()
    mine_share_error = (
        mine_balance["ledger_country_resolved_share_pct"]
        - mine_balance["explicit_country_resolved_share_pct"]
    ).abs()
    mine_equity_error = (mine_balance["equity_balance"] - 1.0).abs()
    ledger_total = float(country_ledger["equity_weighted_mine_equivalent"].sum())
    display_omitted_total = float(
        display_omitted_country_ledger["equity_weighted_mine_equivalent"].sum()
    )
    unattributed_total = float(mine_balance["unattributed_mine_equivalent"].sum())
    eligible_count = int(investor_display["eligible_fig5a_min_10_mines"].sum())

    qa = pd.DataFrame(
        [
            qa_row(
                "mine_bridge_count",
                pass_fail(len(bridge) == EXPECTED_MINES),
                len(bridge),
                EXPECTED_MINES,
                "All mine-level ecological labels must be retained.",
            ),
            qa_row(
                "coverage_id_exact_match",
                pass_fail(mine_ids == coverage_ids),
                len(mine_ids & coverage_ids),
                EXPECTED_MINES,
                f"missing={len(mine_ids-coverage_ids)}; extra={len(coverage_ids-mine_ids)}",
            ),
            qa_row(
                "ecological_label_complete",
                pass_fail(
                    country_ledger["ecological_attention_category"].notna().all()
                ),
                int(country_ledger["ecological_attention_category"].notna().sum()),
                len(country_ledger),
                "Every attributed mine-country row must inherit one mine label.",
            ),
            qa_row(
                "resolved_share_matches_coverage",
                pass_fail(mine_share_error.max() <= SHARE_TOLERANCE),
                float(mine_share_error.max()),
                f"<={SHARE_TOLERANCE}",
                "Country-ledger reported shares must reproduce step 10 coverage.",
            ),
            qa_row(
                "mine_equity_conservation",
                pass_fail(mine_equity_error.max() <= SHARE_TOLERANCE),
                float(mine_equity_error.max()),
                f"<={SHARE_TOLERANCE}",
                "Traceable country equity plus unattributed equity must equal one mine.",
            ),
            qa_row(
                "global_equity_conservation",
                pass_fail(
                    abs(
                        ledger_total
                        + unattributed_total
                        - EXPECTED_MINES
                    )
                    <= SHARE_TOLERANCE
                ),
                ledger_total + unattributed_total,
                EXPECTED_MINES,
                "Global mine-equivalent balance.",
            ),
            qa_row(
                "investor_category_sum",
                pass_fail(investor_category_error.max() <= SHARE_TOLERANCE),
                float(investor_category_error.max()),
                f"<={SHARE_TOLERANCE}",
                "Four categories must sum to each investor's traceable portfolio.",
            ),
            qa_row(
                "resource_category_sum",
                pass_fail(resource_category_error.max() <= SHARE_TOLERANCE),
                float(resource_category_error.max()),
                f"<={SHARE_TOLERANCE}",
                "Four categories must sum to each resource country's traceable portfolio.",
            ),
            qa_row(
                "investor_total_matches_ledger",
                pass_fail(
                    abs(investor["total_traceable_mine_equivalent"].sum() - ledger_total)
                    <= SHARE_TOLERANCE
                ),
                float(investor["total_traceable_mine_equivalent"].sum()),
                ledger_total,
                "Investor aggregation must preserve the country ledger total.",
            ),
            qa_row(
                "resource_total_matches_ledger",
                pass_fail(
                    abs(resource["total_traceable_mine_equivalent"].sum() - ledger_total)
                    <= SHARE_TOLERANCE
                ),
                float(resource["total_traceable_mine_equivalent"].sum()),
                ledger_total,
                "Resource aggregation must preserve the country ledger total.",
            ),
            qa_row(
                "bilateral_total_matches_ledger",
                pass_fail(
                    abs(
                        bilateral["total_traceable_mine_equivalent"].sum()
                        - ledger_total
                    )
                    <= SHARE_TOLERANCE
                ),
                float(bilateral["total_traceable_mine_equivalent"].sum()),
                ledger_total,
                "Bilateral aggregation must preserve the country ledger total.",
            ),
            qa_row(
                "domestic_foreign_sum",
                pass_fail(domestic_foreign_error.max() <= SHARE_TOLERANCE),
                float(domestic_foreign_error.max()),
                f"<={SHARE_TOLERANCE}",
                "Domestic plus cross-border equity must equal each investor portfolio.",
            ),
            qa_row(
                "fig5a_selection_count",
                pass_fail(len(fig5a) == min(N_FIG5A_INVESTORS, eligible_count)),
                len(fig5a),
                min(N_FIG5A_INVESTORS, eligible_count),
                f"Eligible investor countries={eligible_count}.",
            ),
            qa_row(
                "fig5b_investor_selection_count",
                pass_fail(len(fig5b_investor) == N_FIG5B_INVESTORS),
                len(fig5b_investor),
                N_FIG5B_INVESTORS,
                "Top investor countries by cross-border ecological-attention mine count.",
            ),
            qa_row(
                "fig5b_resource_selection_count",
                pass_fail(len(fig5b_resource) == N_FIG5B_RESOURCES),
                len(fig5b_resource),
                N_FIG5B_RESOURCES,
                "Top resource countries by displayed cross-border ecological-attention mine count.",
            ),
            qa_row(
                "hkg_included_in_calculations",
                pass_fail(
                    country_ledger["investor_iso3"].isin(DISPLAY_OMITTED_ISO3).any()
                    and investor["investor_iso3"].isin(DISPLAY_OMITTED_ISO3).any()
                    and bilateral["investor_iso3"].isin(DISPLAY_OMITTED_ISO3).any()
                ),
                int(country_ledger["investor_iso3"].isin(DISPLAY_OMITTED_ISO3).sum()),
                ">0",
                "HKG records must enter the analytical ledger and aggregate totals.",
            ),
            qa_row(
                "hkg_omitted_from_named_display_selection",
                pass_fail(not selection["country_iso3"].isin(DISPLAY_OMITTED_ISO3).any()),
                int(selection["country_iso3"].isin(DISPLAY_OMITTED_ISO3).sum()),
                0,
                "HKG must be absent only from the named-country display selection.",
            ),
            qa_row(
                "traceable_country_equity",
                "INFO",
                ledger_total,
                EXPECTED_MINES,
                "Country-resolved reported equity in mine-equivalent units; not a completeness requirement.",
            ),
            qa_row(
                "named_display_omitted_equity",
                "INFO",
                display_omitted_total,
                EXPECTED_MINES,
                "Included in calculations; omitted only from named-country display outputs.",
            ),
            qa_row(
                "unattributed_equity",
                "INFO",
                unattributed_total,
                EXPECTED_MINES,
                "Unreported or country-unresolved equity in mine-equivalent units.",
            ),
        ]
    )
    qa.to_csv(OUT_QA, index=False, encoding="utf-8-sig")

    failures = qa[qa["status"].eq("FAIL")]
    if not failures.empty:
        raise RuntimeError(
            "Ownership ecological-portfolio QA failed: "
            + ", ".join(failures["check"].astype(str))
        )

    country_ledger.to_csv(OUT_LEDGER, index=False, encoding="utf-8-sig")
    display_omitted_country_ledger.to_csv(
        OUT_DISPLAY_OMISSION_AUDIT, index=False, encoding="utf-8-sig"
    )
    investor.to_csv(OUT_INVESTOR, index=False, encoding="utf-8-sig")
    resource.to_csv(OUT_RESOURCE, index=False, encoding="utf-8-sig")
    bilateral.to_csv(OUT_BILATERAL, index=False, encoding="utf-8-sig")
    investor_display.to_csv(OUT_INVESTOR_DISPLAY, index=False, encoding="utf-8-sig")
    resource_display.to_csv(OUT_RESOURCE_DISPLAY, index=False, encoding="utf-8-sig")
    bilateral_display.to_csv(OUT_BILATERAL_DISPLAY, index=False, encoding="utf-8-sig")
    global_summary.to_csv(OUT_GLOBAL, index=False, encoding="utf-8-sig")
    selection.to_csv(OUT_SELECTION, index=False, encoding="utf-8-sig")

    OUT_METHOD.write_text(
        "# Fig. 5 ecological-attention portfolio definition\n\n"
        "- Mine unit: all 1,278 copper-mine 30-km potential spatial interaction zones.\n"
        "- High RWR: mine reference point in the global top 10% of RWR grid cells.\n"
        "- Priority-screened species: positive-area range intersection with at least one of the 203 screened species.\n"
        "- Attribution unit: reported ownership share divided by 100 (equity-weighted mine equivalent).\n"
        "- Investor-country denominator: all explicitly reported shares resolved to that investor country.\n"
        "- Missing, unreported and country-unresolved shares are not redistributed.\n"
        "- All country-resolved records enter calculations and global totals; the named-country display selection is produced separately.\n"
        "- Cross-border link: investor ISO3 differs from resource-country ISO3.\n"
        "- Ecological-attention mine: High RWR only, priority-screened species only, or dual ecological attention; equivalently, either ecological condition is met.\n"
        "- Fig. 5a eligibility: at least 10 unique mines with traceable country-resolved equity; display is selected by the number of unique ecological-attention mines.\n"
        "- Fig. 5b selection: leading investor countries by cross-border ecological-attention mine count, followed by leading resource countries within those selected portfolios; each mine is counted once per investor-resource pair regardless of ownership percentage.\n"
        "- Cross-portfolio counting: if the same mine has reported ownership resolved to more than one investor country, it is counted once in each relevant investor-country portfolio; counts therefore describe country portfolios and are not summed as a global number of distinct mines.\n",
        encoding="utf-8",
    )

    print("19Ownership_EcologicalPortfolio_Statistics completed.")
    print(f"Traceable country-resolved mine equivalents: {ledger_total:.3f}")
    print(f"Named-display omitted mine equivalents: {display_omitted_total:.3f}")
    print(f"Unattributed mine equivalents: {unattributed_total:.3f}")
    print(f"Investor countries: {len(investor)}")
    print(f"Fig. 5a eligible countries: {eligible_count}; displayed: {len(fig5a)}")
    print("Fig. 5a selection:")
    print(
        fig5a[
            [
                "selection_rank",
                "investor_country",
                "num_unique_mines",
                "ecological_attention_unique_mines",
            ]
        ].to_string(index=False)
    )
    print("Fig. 5b investor countries:")
    print(
        fig5b_investor[
            [
                "selection_rank",
                "investor_country",
                "foreign_ecological_attention_unique_mines",
            ]
        ].to_string(index=False)
    )
    print("Fig. 5b resource countries:")
    print(
        fig5b_resource[
            [
                "selection_rank",
                "resource_country",
                "selected_investor_foreign_ecological_attention_unique_mines",
            ]
        ].to_string(index=False)
    )
    print(f"QA: {OUT_QA}")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
