# -*- coding: utf-8 -*-
"""
12_attribute_exposure_by_ownership.py

研究内容四——股权归因、全球/国家统计和双边流向基础表。

核心公式：
    attributed_area(m, j) = attribution_base_area(m) * share(m, j) / 100

面积定义：步骤 11 已经在单矿区内对全部物种栖息地做了空间union，不按物种叠加。

跨矿区重叠的默认处理：
    ATTRIBUTION_AREA_MODE = "country_union_balanced"
同一资源拥有国内，先以该国严格union实际面积为总量，再按照各矿区单独union面积
的相对权重分配到矿区。这样股权归因总量与该国实际物理暴露面积严格平衡，避免
相邻矿区30 km范围重叠导致面积重复。公式为：
    E*_m = E_m * UnionArea_country / sum(E_m in country)

若需要完全复现传统“每个项目单独计算后直接相加”的口径，可将参数改为
"mine_union_unbalanced"；该模式仍不按物种叠加，但会重复计算不同矿区共享的空间。

所有未知股权和投资国未解析股权均保留为 unattributed，不会把已知股权归一化到100%。
全部可解析注册地记录进入计算、全球汇总和具名展示表。
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script
from _ownership_labels import normalize_registration_names

from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd


warnings.filterwarnings("ignore")


# =============================================================================
# 路径与参数
# =============================================================================

BASE_DIR = Path(_nc_project(""))
PREP_DIR = BASE_DIR / "10_mine_ownership"
EXPOSURE_DIR = BASE_DIR / "11_mine_exposure"
OUTPUT_DIR = BASE_DIR / "12_ownership_attribution"

MINE_MASTER_CSV = PREP_DIR / "Mine_Master_30km.csv"
OWNER_LEDGER_CSV = PREP_DIR / "Mine_Owner_Share_Ledger.csv"
OWNERSHIP_COVERAGE_CSV = PREP_DIR / "Mine_Ownership_Coverage.csv"
MINE_EXPOSURE_CSV = EXPOSURE_DIR / "Mine_Exposure_30km.csv"
COUNTRY_EXPOSURE_CSV = EXPOSURE_DIR / "Country_Exposure_Union_30km.csv"
GLOBAL_EXPOSURE_CSV = EXPOSURE_DIR / "Global_Exposure_Union_30km.csv"

OUT_MINE_ATTRIBUTION = OUTPUT_DIR / "Mine_InvestorCountry_Attribution.csv"
OUT_MINE_BALANCE = OUTPUT_DIR / "Mine_Attribution_Balance.csv"
OUT_GLOBAL = OUTPUT_DIR / "Global_Responsibility_Summary.csv"
OUT_INVESTOR = OUTPUT_DIR / "Investor_Country_Summary_All.csv"
OUT_RESOURCE = OUTPUT_DIR / "Resource_Country_Summary_All.csv"
OUT_BILATERAL = OUTPUT_DIR / "Bilateral_Responsibility_Flow_All.csv"
OUT_INVESTOR_DISPLAY = OUTPUT_DIR / "Investor_Country_Summary_NamedDisplay.csv"
OUT_RESOURCE_DISPLAY = OUTPUT_DIR / "Resource_Country_Summary_NamedDisplay.csv"
OUT_BILATERAL_DISPLAY = OUTPUT_DIR / "Bilateral_Responsibility_Flow_NamedDisplay.csv"
OUT_DISPLAY_OMISSION_AUDIT = OUTPUT_DIR / "NamedDisplay_Omission_Audit.csv"
OUT_QA = OUTPUT_DIR / "Attribution_Balance_QA.csv"
OUT_METHOD = OUTPUT_DIR / "Attribution_Method_Metadata.txt"

# 默认确保资源拥有国层面的归因总面积等于该国严格union实际面积。
ATTRIBUTION_AREA_MODE = "country_union_balanced"
VALID_AREA_MODES = {"country_union_balanced", "mine_union_unbalanced"}
BALANCE_TOLERANCE_KM2 = 1e-5
SHARE_TOLERANCE_PCT = 1e-6
DISPLAY_OMITTED_ISO3 = set()


# =============================================================================
# 读取与检查
# =============================================================================

def read_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"缺少{label}：{path}")
    return normalize_registration_names(
        pd.read_csv(path, encoding="utf-8-sig", dtype={"property_id": str})
    )


def require_columns(df, columns, label):
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f"{label}缺少字段：{sorted(missing)}")


def safe_pct(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    return np.divide(
        numerator * 100.0, denominator,
        out=np.zeros_like(numerator, dtype=float), where=denominator > 0,
    )


def load_inputs():
    print("\n[1/6] 读取所有权与实际暴露结果")
    master = read_csv(MINE_MASTER_CSV, "矿区主表")
    ledger = read_csv(OWNER_LEDGER_CSV, "所有者股权账本")
    coverage = read_csv(OWNERSHIP_COVERAGE_CSV, "股权覆盖率表")
    mine_exp = read_csv(MINE_EXPOSURE_CSV, "矿区实际暴露表")
    country_exp = read_csv(COUNTRY_EXPOSURE_CSV, "国家union暴露表")
    global_exp = read_csv(GLOBAL_EXPOSURE_CSV, "全球union暴露表")

    require_columns(master, [
        "property_id", "property_name", "resource_country", "resource_iso3", "resource_continent"
    ], "矿区主表")
    require_columns(ledger, [
        "property_id", "reported_share_pct", "investor_country", "investor_iso3",
        "investor_continent", "resource_iso3"
    ], "所有者账本")
    require_columns(coverage, [
        "property_id", "explicit_country_resolved_share_pct", "total_unattributed_share_pct"
    ], "股权覆盖率表")
    require_columns(mine_exp, [
        "property_id", "mine_exposure_area_km2", "resource_iso3", "num_exposed_species"
    ], "矿区实际暴露表")
    require_columns(country_exp, [
        "resource_iso3", "physical_exposure_union_area_km2", "mine_exposure_area_sum_km2"
    ], "国家实际暴露表")
    require_columns(global_exp, ["global_strict_physical_exposure_area_km2"], "全球实际暴露表")

    for df, label in [(master, "master"), (coverage, "coverage"), (mine_exp, "mine exposure")]:
        if df["property_id"].duplicated().any():
            raise ValueError(f"{label} 的 property_id 重复。")
    return master, ledger, coverage, mine_exp, country_exp, global_exp


# =============================================================================
# 矿区归因底数与股权归因
# =============================================================================

def build_mine_base(master, coverage, mine_exp, country_exp):
    print("\n[2/6] 建立矿区归因面积底数")
    if ATTRIBUTION_AREA_MODE not in VALID_AREA_MODES:
        raise ValueError(f"ATTRIBUTION_AREA_MODE必须是：{sorted(VALID_AREA_MODES)}")

    base = master[[
        "property_id", "property_name", "resource_country", "resource_iso3", "resource_continent"
    ]].merge(
        mine_exp[[
            "property_id", "mine_exposure_area_km2", "num_exposed_species", "buffer_area_km2"
        ]], on="property_id", how="left", validate="one_to_one"
    ).merge(
        coverage[[
            "property_id", "explicit_reported_share_pct", "explicit_country_resolved_share_pct",
            "explicit_country_unresolved_share_pct", "unreported_share_pct",
            "total_unattributed_share_pct", "ownership_source"
        ]], on="property_id", how="left", validate="one_to_one"
    )
    if base["mine_exposure_area_km2"].isna().any():
        raise ValueError("部分矿区缺少步骤 11 的实际暴露面积。")

    country_sum = base.groupby("resource_iso3", as_index=False).agg(
        mine_area_sum_recalculated_km2=("mine_exposure_area_km2", "sum")
    )
    factors = country_exp[[
        "resource_iso3", "physical_exposure_union_area_km2", "mine_exposure_area_sum_km2"
    ]].merge(country_sum, on="resource_iso3", how="outer", validate="one_to_one")
    factors[[
        "physical_exposure_union_area_km2", "mine_exposure_area_sum_km2",
        "mine_area_sum_recalculated_km2"
    ]] = factors[[
        "physical_exposure_union_area_km2", "mine_exposure_area_sum_km2",
        "mine_area_sum_recalculated_km2"
    ]].fillna(0.0)
    diff = (
        factors["mine_exposure_area_sum_km2"] - factors["mine_area_sum_recalculated_km2"]
    ).abs()
    if (diff > np.maximum(BALANCE_TOLERANCE_KM2, factors["mine_area_sum_recalculated_km2"] * 1e-8)).any():
        raise ValueError("步骤 11 的矿区面积合计与国家表不一致。")

    if ATTRIBUTION_AREA_MODE == "country_union_balanced":
        impossible = factors[
            factors["physical_exposure_union_area_km2"]
            > factors["mine_area_sum_recalculated_km2"]
            + np.maximum(BALANCE_TOLERANCE_KM2, factors["mine_area_sum_recalculated_km2"] * 1e-8)
        ]
        if not impossible.empty:
            raise ValueError(
                "国家union面积大于该国矿区单独union面积之和；请检查步骤 11 几何结果。"
            )
        factors["mine_overlap_balance_factor"] = np.where(
            factors["mine_area_sum_recalculated_km2"] > 0,
            factors["physical_exposure_union_area_km2"]
            / factors["mine_area_sum_recalculated_km2"],
            1.0,
        )
    else:
        factors["mine_overlap_balance_factor"] = 1.0

    base = base.merge(
        factors[["resource_iso3", "physical_exposure_union_area_km2", "mine_overlap_balance_factor"]],
        on="resource_iso3", how="left", validate="many_to_one"
    )
    base["attribution_base_area_km2"] = (
        base["mine_exposure_area_km2"] * base["mine_overlap_balance_factor"]
    )
    base["attribution_area_mode"] = ATTRIBUTION_AREA_MODE

    if ATTRIBUTION_AREA_MODE == "country_union_balanced":
        check = base.groupby("resource_iso3")["attribution_base_area_km2"].sum()
        target = factors.set_index("resource_iso3")["physical_exposure_union_area_km2"]
        aligned = pd.concat([check.rename("actual"), target.rename("target")], axis=1).fillna(0)
        if ((aligned["actual"] - aligned["target"]).abs() > BALANCE_TOLERANCE_KM2).any():
            raise RuntimeError("国家union平衡后，矿区归因底数仍未与国家实际面积平衡。")
    print(f"  面积模式：{ATTRIBUTION_AREA_MODE}")
    print(f"  归因面积底数：{base['attribution_base_area_km2'].sum():,.2f} km²")
    return base, factors


def build_investor_attribution(base, ledger):
    print("\n[3/6] 按投资国合并股权并分配面积")
    ledger = normalize_registration_names(ledger)
    ledger["reported_share_pct"] = pd.to_numeric(ledger["reported_share_pct"], errors="coerce")
    resolved = ledger[
        ledger["reported_share_pct"].notna()
        & ledger["reported_share_pct"].gt(0)
        & ledger["investor_iso3"].fillna("").astype(str).str.strip().ne("")
    ].copy()
    country_share_all = resolved.groupby([
        "property_id", "investor_iso3", "investor_country", "investor_continent"
    ], as_index=False).agg(
        investor_country_share_pct=("reported_share_pct", "sum"),
        num_owner_records=("owner_name", "size"),
        owners=("owner_name", lambda x: "; ".join(sorted(set(map(str, x))))),
    )
    excessive = country_share_all.groupby("property_id")["investor_country_share_pct"].sum()
    if (excessive > 100.0 + SHARE_TOLERANCE_PCT).any():
        bad = excessive[excessive > 100.0 + SHARE_TOLERANCE_PCT]
        raise ValueError(f"投资国合并后的矿区股权超过100%：{bad.to_dict()}")

    attribution = country_share_all.merge(
        base, on="property_id", how="left", validate="many_to_one"
    )
    if attribution["attribution_base_area_km2"].isna().any():
        raise ValueError("部分股权记录无法匹配矿区归因底数。")
    attribution["attributed_exposure_area_km2"] = (
        attribution["attribution_base_area_km2"]
        * attribution["investor_country_share_pct"] / 100.0
    )
    attribution["is_domestic"] = (
        attribution["investor_iso3"] == attribution["resource_iso3"]
    ).astype(int)
    attribution["relationship"] = np.where(
        attribution["is_domestic"].eq(1), "Domestic", "Foreign"
    )
    attribution = attribution.sort_values([
        "resource_country", "property_name", "attributed_exposure_area_km2"
    ], ascending=[True, True, False])
    attribution.to_csv(OUT_MINE_ATTRIBUTION, index=False, encoding="utf-8-sig")

    display_omitted_attribution = attribution[
        attribution["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    display_omitted_attribution["display_rule"] = (
        "all countries and regions retained in calculations and named displays"
    )
    display_omitted_attribution.to_csv(
        OUT_DISPLAY_OMISSION_AUDIT, index=False, encoding="utf-8-sig"
    )

    resolved_by_mine = attribution.groupby("property_id", as_index=False).agg(
        resolved_share_recalculated_pct=("investor_country_share_pct", "sum"),
        traceable_attributed_area_km2=("attributed_exposure_area_km2", "sum"),
        domestic_attributed_area_km2=(
            "attributed_exposure_area_km2",
            lambda x: x[attribution.loc[x.index, "is_domestic"].eq(1)].sum(),
        ),
        foreign_attributed_area_km2=(
            "attributed_exposure_area_km2",
            lambda x: x[attribution.loc[x.index, "is_domestic"].eq(0)].sum(),
        ),
    )
    balance = base.merge(
        resolved_by_mine, on="property_id", how="left", validate="one_to_one"
    )
    for col in [
        "resolved_share_recalculated_pct", "traceable_attributed_area_km2",
        "domestic_attributed_area_km2", "foreign_attributed_area_km2",
    ]:
        balance[col] = balance[col].fillna(0.0)
    balance["unattributed_area_km2"] = (
        balance["attribution_base_area_km2"]
        - balance["traceable_attributed_area_km2"]
    ).clip(lower=0.0)
    balance["balance_total_km2"] = (
        balance["traceable_attributed_area_km2"]
        + balance["unattributed_area_km2"]
    )
    balance["balance_difference_km2"] = (
        balance["balance_total_km2"] - balance["attribution_base_area_km2"]
    )
    balance["share_crosscheck_difference_pct"] = (
        balance["resolved_share_recalculated_pct"]
        - balance["explicit_country_resolved_share_pct"]
    )
    if (balance["balance_difference_km2"].abs() > BALANCE_TOLERANCE_KM2).any():
        raise RuntimeError("矿区股权归因面积未能平衡。")
    if (balance["share_crosscheck_difference_pct"].abs() > SHARE_TOLERANCE_PCT).any():
        raise RuntimeError("所有者账本与股权覆盖率表的已解析股权不一致。")
    balance.to_csv(OUT_MINE_BALANCE, index=False, encoding="utf-8-sig")
    print(f"  已解析投资国归因记录：{len(attribution):,}")
    print(
        f"  Named-display omission records: {len(display_omitted_attribution):,}; "
        f"area={display_omitted_attribution['attributed_exposure_area_km2'].sum():,.2f} km2"
    )
    return attribution, balance, display_omitted_attribution


# =============================================================================
# 汇总表
# =============================================================================

def build_global_summary(balance, global_exp):
    traceable = balance["traceable_attributed_area_km2"].sum()
    domestic = balance["domestic_attributed_area_km2"].sum()
    foreign = balance["foreign_attributed_area_km2"].sum()
    unattributed = balance["unattributed_area_km2"].sum()
    base_total = balance["attribution_base_area_km2"].sum()
    raw_mine_total = balance["mine_exposure_area_km2"].sum()
    strict_global = float(global_exp.iloc[0]["global_strict_physical_exposure_area_km2"])
    summary = pd.DataFrame([{
        "attribution_area_mode": ATTRIBUTION_AREA_MODE,
        "global_strict_physical_exposure_area_km2": strict_global,
        "resource_country_attributed_physical_area_km2": base_total,
        "raw_sum_of_mine_union_areas_km2_not_used_in_default": raw_mine_total,
        "traceable_equity_associated_area_km2": traceable,
        "unattributed_area_km2": unattributed,
        "traceable_rate_pct": traceable / base_total * 100 if base_total > 0 else 0.0,
        "domestic_equity_associated_area_km2": domestic,
        "foreign_equity_associated_area_km2": foreign,
        "domestic_share_of_traceable_pct": domestic / traceable * 100 if traceable > 0 else 0.0,
        "foreign_share_of_traceable_pct": foreign / traceable * 100 if traceable > 0 else 0.0,
        "foreign_share_of_total_area_pct": foreign / base_total * 100 if base_total > 0 else 0.0,
        "num_mines": balance["property_id"].nunique(),
        "num_exposed_mines": int((balance["attribution_base_area_km2"] > 0).sum()),
        "interpretation_boundary": (
            "equity-associated shared governance relationship; not legal liability or observed habitat loss"
        ),
    }])
    summary.to_csv(OUT_GLOBAL, index=False, encoding="utf-8-sig")
    return summary


def build_investor_summary(attribution):
    rows = []
    for keys, group in attribution.groupby([
        "investor_iso3", "investor_country", "investor_continent"
    ]):
        iso3, country, continent = keys
        domestic = group.loc[group["is_domestic"].eq(1), "attributed_exposure_area_km2"].sum()
        foreign = group.loc[group["is_domestic"].eq(0), "attributed_exposure_area_km2"].sum()
        rows.append({
            "investor_iso3": iso3,
            "investor_country": country,
            "investor_continent": continent,
            "total_attributed_area_km2": domestic + foreign,
            "domestic_project_attributed_area_km2": domestic,
            "foreign_project_attributed_area_km2": foreign,
            "num_mines": group["property_id"].nunique(),
            "num_resource_countries": group["resource_iso3"].nunique(),
            "num_foreign_resource_countries": group.loc[
                group["is_domestic"].eq(0), "resource_iso3"
            ].nunique(),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError("没有可汇总的投资国归因结果。")
    out["foreign_project_share_pct"] = safe_pct(
        out["foreign_project_attributed_area_km2"], out["total_attributed_area_km2"]
    )
    out = out.sort_values(
        ["foreign_project_attributed_area_km2", "total_attributed_area_km2", "investor_country"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    global_foreign = out["foreign_project_attributed_area_km2"].sum()
    out["foreign_area_rank"] = np.arange(1, len(out) + 1)
    out["share_of_global_foreign_area_pct"] = safe_pct(
        out["foreign_project_attributed_area_km2"], global_foreign
    )
    out["cumulative_global_foreign_area_pct"] = out[
        "share_of_global_foreign_area_pct"
    ].cumsum().clip(upper=100.0)
    out.to_csv(OUT_INVESTOR, index=False, encoding="utf-8-sig")
    return out


def build_resource_summary(balance, country_exp):
    grouped = balance.groupby([
        "resource_iso3", "resource_country", "resource_continent"
    ], as_index=False).agg(
        attribution_base_area_km2=("attribution_base_area_km2", "sum"),
        raw_sum_of_mine_union_areas_km2=("mine_exposure_area_km2", "sum"),
        traceable_area_km2=("traceable_attributed_area_km2", "sum"),
        domestic_attributed_area_km2=("domestic_attributed_area_km2", "sum"),
        foreign_attributed_area_km2=("foreign_attributed_area_km2", "sum"),
        unattributed_area_km2=("unattributed_area_km2", "sum"),
        num_mines_in_attribution=("property_id", "nunique"),
        num_exposed_mines_in_attribution=("attribution_base_area_km2", lambda x: int((x > 0).sum())),
    )
    out = country_exp.merge(
        grouped, on=["resource_iso3", "resource_country", "resource_continent"],
        how="left", validate="one_to_one"
    )
    numeric = [
        "attribution_base_area_km2", "raw_sum_of_mine_union_areas_km2",
        "traceable_area_km2", "domestic_attributed_area_km2",
        "foreign_attributed_area_km2", "unattributed_area_km2",
    ]
    out[numeric] = out[numeric].fillna(0.0)
    out["traceable_rate_pct"] = safe_pct(
        out["traceable_area_km2"], out["attribution_base_area_km2"]
    )
    out["foreign_share_of_traceable_pct"] = np.where(
        out["traceable_area_km2"] > 0,
        out["foreign_attributed_area_km2"] / out["traceable_area_km2"] * 100.0,
        np.nan,
    )
    out["foreign_share_of_total_area_pct"] = safe_pct(
        out["foreign_attributed_area_km2"], out["attribution_base_area_km2"]
    )
    out = out.sort_values(
        ["foreign_attributed_area_km2", "physical_exposure_union_area_km2", "resource_country"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    global_foreign = out["foreign_attributed_area_km2"].sum()
    out["foreign_area_rank"] = np.arange(1, len(out) + 1)
    out["share_of_global_foreign_area_pct"] = safe_pct(
        out["foreign_attributed_area_km2"], global_foreign
    )
    out["cumulative_global_foreign_area_pct"] = out[
        "share_of_global_foreign_area_pct"
    ].cumsum().clip(upper=100.0)
    out.to_csv(OUT_RESOURCE, index=False, encoding="utf-8-sig")
    return out


def build_bilateral_flows(attribution):
    bilateral = attribution.groupby([
        "resource_iso3", "resource_country", "resource_continent",
        "investor_iso3", "investor_country", "investor_continent", "is_domestic", "relationship"
    ], as_index=False).agg(
        attributed_exposure_area_km2=("attributed_exposure_area_km2", "sum"),
        num_mines=("property_id", "nunique"),
        num_owner_records=("num_owner_records", "sum"),
        equity_share_sum_across_mines_pct=("investor_country_share_pct", "sum"),
    )
    bilateral = bilateral.sort_values(
        ["is_domestic", "attributed_exposure_area_km2"], ascending=[True, False]
    ).reset_index(drop=True)
    bilateral.to_csv(OUT_BILATERAL, index=False, encoding="utf-8-sig")
    return bilateral


def save_qa(base, factors, balance, global_summary, display_omitted_attribution):
    qa_rows = []
    qa_rows.append({
        "check": "max_absolute_mine_balance_difference_km2",
        "value": balance["balance_difference_km2"].abs().max(),
        "status": "PASS" if balance["balance_difference_km2"].abs().max() <= BALANCE_TOLERANCE_KM2 else "FAIL",
    })
    qa_rows.append({
        "check": "max_absolute_share_crosscheck_difference_pct",
        "value": balance["share_crosscheck_difference_pct"].abs().max(),
        "status": "PASS" if balance["share_crosscheck_difference_pct"].abs().max() <= SHARE_TOLERANCE_PCT else "FAIL",
    })
    qa_rows.append({
        "check": "num_mines_without_traceable_investor_country",
        "value": int((balance["resolved_share_recalculated_pct"] <= 0).sum()),
        "status": "INFO",
    })
    qa_rows.append({
        "check": "num_mines_with_partial_traceable_share",
        "value": int(balance["resolved_share_recalculated_pct"].between(0, 100, inclusive="neither").sum()),
        "status": "INFO",
    })
    qa_rows.append({
        "check": "num_resource_countries_with_cross_mine_overlap",
        "value": int((factors["mine_overlap_balance_factor"] < 1 - 1e-10).sum()),
        "status": "INFO",
    })
    qa_rows.append({
        "check": "named_display_omitted_attribution_records",
        "value": len(display_omitted_attribution),
        "status": "INFO",
    })
    qa_rows.append({
        "check": "named_display_omitted_attributed_area_km2",
        "value": display_omitted_attribution["attributed_exposure_area_km2"].sum(),
        "status": "INFO",
    })
    qa_rows.append({
        "check": "attribution_area_mode",
        "value": ATTRIBUTION_AREA_MODE,
        "status": "INFO",
    })
    pd.DataFrame(qa_rows).to_csv(OUT_QA, index=False, encoding="utf-8-sig")

    metadata = {
        "area_mode": ATTRIBUTION_AREA_MODE,
        "species_overlap_rule": "union across species; each physical location counted once",
        "cross_mine_rule": (
            "country physical union balanced proportionally across mine-level union areas"
            if ATTRIBUTION_AREA_MODE == "country_union_balanced"
            else "mine-level union areas summed; cross-mine overlap may repeat"
        ),
        "unknown_share_rule": "not normalized; retained as unattributed",
        "named_country_display_rule": (
            "All countries and regions are retained in calculations, global totals "
            "and named display tables"
        ),
        "ownership_interpretation": "equity-associated shared governance, not legal liability",
        "global_summary": global_summary.iloc[0].to_dict(),
    }
    OUT_METHOD.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    master, ledger, coverage, mine_exp, country_exp, global_exp = load_inputs()
    base, factors = build_mine_base(master, coverage, mine_exp, country_exp)
    attribution, balance, display_omitted_attribution = build_investor_attribution(base, ledger)

    print("\n[4/6] 汇总全球、投资国和资源拥有国")
    global_summary = build_global_summary(balance, global_exp)
    investor = build_investor_summary(attribution)
    resource = build_resource_summary(balance, country_exp)

    print("\n[5/6] 生成完整双边流向表")
    bilateral = build_bilateral_flows(attribution)

    investor_display = investor[
        ~investor["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    resource_display = resource[
        ~resource["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    bilateral_display = bilateral[
        ~bilateral["investor_iso3"].isin(DISPLAY_OMITTED_ISO3)
        & ~bilateral["resource_iso3"].isin(DISPLAY_OMITTED_ISO3)
    ].copy()
    investor_display.to_csv(OUT_INVESTOR_DISPLAY, index=False, encoding="utf-8-sig")
    resource_display.to_csv(OUT_RESOURCE_DISPLAY, index=False, encoding="utf-8-sig")
    bilateral_display.to_csv(OUT_BILATERAL_DISPLAY, index=False, encoding="utf-8-sig")

    print("\n[6/6] 保存平衡QA与方法元数据")
    save_qa(base, factors, balance, global_summary, display_omitted_attribution)
    print(f"  投资国：{len(investor):,}")
    print(f"  资源拥有国：{len(resource):,}")
    print(f"  双边关系：{len(bilateral):,}")
    print(f"  输出目录：{OUTPUT_DIR}")
    print("\n完成。下一步可依次运行15–17制图脚本。")


if __name__ == "__main__":
    main()
