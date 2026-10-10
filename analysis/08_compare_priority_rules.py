# -*- coding: utf-8 -*-
"""
08_compare_priority_rules.py

研究内容二：空间冲突限制目标稳健性分析（轻量版）

本脚本只读取步骤 07 已生成的物种统计表，不进行任何空间相交、栅格计算或
随机化，因此通常可在数秒至数分钟内完成。

比较两套有文献依据的物种代表性目标，并将保护目标转换为最大允许空间冲突比例：

1. 主分析：Sreekar–Watson (2026)
   - 分布面积 <= 200 km²：保护目标 100%
   - 分布面积 >= 250,000 km²：保护目标 30%
   - 中间范围：对数线性插值
   - 默认不设保护目标面积上限；可选 NC_CAP_CHECK=1 启用 1,000,000 km² 检验

2. 敏感性方案：Rodrigues et al. (2004)，亦见 Lewin et al. (2024)
   - 分布面积 <= 1,000 km²：保护目标 100%
   - 分布面积 >= 250,000 km²：保护目标 10%
   - 中间范围：对数线性插值
   - 默认不设保护目标面积上限；可选 NC_CAP_CHECK=1 启用 1,000,000 km² 检验

对任一方案：
    最大允许冲突比例 = 100% - 保护目标比例
    优先干预物种 = Impact_Percent > Conflict_Limit_Percent + EPS

只分析已发生铜矿空间暴露的物种（预期 891 种）。

输入
----
<work-root>/07_priority_species/
    Species_PrioritySpecies_Full.csv

输出
----
<work-root>/08_priority_robustness/
    Robustness_SpeciesResults.csv
    Robustness_SchemeSummary.csv
    Robustness_ComparisonSummary.csv
    Robustness_ByTaxon.csv
    Robustness_ByCurrentRange.csv
    Robustness_BoundarySpecies.csv
    Robustness_RunSummary.txt

依赖
----
numpy, pandas
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings

import numpy as np
import pandas as pd


# =============================================================================
# 1. 路径与可调参数
# =============================================================================
BASE_DIR = Path(_nc_project(""))
INPUT_FILE = (
    BASE_DIR
    / "07_priority_species"
    / "Species_PrioritySpecies_Full.csv"
)
OUTPUT_DIR = BASE_DIR / "08_priority_robustness"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_COPPER_AFFECTED_SPECIES = 891
EXPECTED_MAIN_PRIORITY_SPECIES = 203

# False：数量不一致时发出警告但继续运行；True：数量不一致时停止。
STRICT_EXPECTED_COUNT = False

# "auto" 会自动判断冲突字段是 0–1 比例还是 0–100 百分比。
# 如需手动指定，可改为 "proportion" 或 "percent"。
IMPACT_UNIT = "auto"

# 只有严格超过目标才判为优先物种。
EPS_PERCENT = 1e-10

# 距离主分析阈值不超过 1 或 5 个百分点的物种定义为边界物种。
BOUNDARY_LEVELS_PP = (1.0, 5.0)

FIGURE_PREFIX = "FigS_LimitTarget_Robustness"
FIGURE_DPI = 600
SAVE_TIF = True


# =============================================================================
# 2. 两套限制目标
# =============================================================================
SCHEMES = {
    "Main_SreekarWatson2026": {
        "label": "Sreekar–Watson (2026)",
        "short_label": "Current",
        "rare_area_km2": 200.0,
        "rare_target": 1.0,
        "common_area_km2": 250000.0,
        "common_target": 0.30,
        "cap_area_km2": 1000000.0,
        "color": "#D55E00",
    },
    "Traditional_RodriguesLewin": {
        "label": "Rodrigues–Lewin",
        "short_label": "Traditional",
        "rare_area_km2": 1000.0,
        "rare_target": 1.0,
        "common_area_km2": 250000.0,
        "common_target": 0.10,
        "cap_area_km2": 1000000.0,
        "color": "#0072B2",
    },
}

# The manuscript compares uncapped rules; the reported cap check is optional.
from os import environ as _nc_environ
if _nc_environ.get("NC_CAP_CHECK", "0") != "1":
    for _nc_scheme in SCHEMES.values():
        _nc_scheme["cap_area_km2"] = float("inf")


MAIN = "Main_SreekarWatson2026"
TRAD = "Traditional_RodriguesLewin"

TAXON_ORDER = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
TAXON_ALIASES = {
    "amphibian": "Amphibians",
    "amphibians": "Amphibians",
    "amphibia": "Amphibians",
    "bird": "Birds",
    "birds": "Birds",
    "aves": "Birds",
    "fish": "Fish",
    "fishes": "Fish",
    "actinopterygii": "Fish",
    "mammal": "Mammals",
    "mammals": "Mammals",
    "mammalia": "Mammals",
    "reptile": "Reptiles",
    "reptiles": "Reptiles",
    "reptilia": "Reptiles",
}

TAXON_COLORS = {
    "Amphibians": "#D55E00",
    "Birds": "#E69F00",
    "Fish": "#56B4E9",
    "Mammals": "#009E73",
    "Reptiles": "#7A5195",
}

RANGE_ORDER = ["Narrow", "Intermediate", "Widespread"]
RANGE_COLORS = {
    "Narrow": "#e28740",
    "Intermediate": "#84cfbb",
    "Widespread": "#7fc7ef",
}



# =============================================================================
# 3. 数据读取与字段整理
# =============================================================================
def find_column(df, aliases, label, required=True):
    """不区分大小写查找字段。"""
    lookup = {str(c).strip().lower(): c for c in df.columns}
    for alias in aliases:
        if alias.lower() in lookup:
            return lookup[alias.lower()]
    if required:
        raise ValueError(
            f"无法识别{label}字段。\n"
            f"可接受字段名：{aliases}\n"
            f"实际字段：{list(df.columns)}"
        )
    return None


def clean_id(value):
    """统一物种 ID，避免 123 与 123.0 被视为不同物种。"""
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    if text.endswith(".0"):
        text = text[:-2]
    return text


def numeric_series(series):
    """兼容数字中的逗号与百分号。"""
    return pd.to_numeric(
        series.astype(str)
        .str.replace(",", "", regex=False)
        .str.replace("%", "", regex=False)
        .str.strip(),
        errors="coerce",
    )


def parse_bool(series):
    """兼容 True/False、1/0、Yes/No 等常见布尔写法。"""
    true_values = {"true", "1", "1.0", "yes", "y", "t", "affected", "是"}
    return series.astype(str).str.strip().str.lower().isin(true_values)


def normalize_taxon(series):
    cleaned = series.astype(str).str.strip()
    mapped = cleaned.str.lower().map(TAXON_ALIASES)
    return mapped.fillna(cleaned)


def convert_impact_to_percent(raw_series, column_name):
    """将冲突比例统一为 0–100 百分比。"""
    values = numeric_series(raw_series)
    has_percent_sign = raw_series.astype(str).str.contains("%", regex=False).any()

    if IMPACT_UNIT == "proportion":
        values = values * 100.0
        unit_used = "proportion (0–1)"
    elif IMPACT_UNIT == "percent":
        unit_used = "percent (0–100)"
    elif IMPACT_UNIT == "auto":
        finite = values.dropna()
        if has_percent_sign:
            unit_used = "percent sign detected"
        elif not finite.empty and finite.max() <= 1.0 + 1e-10:
            values = values * 100.0
            unit_used = "auto: proportion (0–1)"
        else:
            unit_used = "auto: percent (0–100)"
    else:
        raise ValueError(
            'IMPACT_UNIT 只能设为 "auto"、"proportion" 或 "percent"。'
        )

    if values.dropna().lt(-EPS_PERCENT).any() or values.dropna().gt(100 + 1e-6).any():
        raise ValueError(
            f"{column_name} 转为百分比后存在小于 0 或大于 100 的值，"
            "请检查 IMPACT_UNIT 设置及原始字段。"
        )
    return values.clip(0.0, 100.0), unit_used


def read_affected_species():
    """读取完整表，并仅保留已发生铜矿空间暴露的物种。"""
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"未找到输入文件：{INPUT_FILE}")

    raw = pd.read_csv(INPUT_FILE, low_memory=False)

    id_col = find_column(
        raw,
        ["species_id", "speciesid", "sisid", "sis_id", "id_no", "id"],
        "物种 ID",
        required=False,
    )
    name_col = find_column(
        raw,
        [
            "scientific_name",
            "species_name",
            "binomial",
            "sci_name",
            "species",
        ],
        "物种学名",
        required=False,
    )
    if id_col is None and name_col is None:
        raise ValueError("输入表至少需要包含物种 ID 或物种学名字段。")

    taxon_col = find_column(
        raw,
        ["taxon", "taxa", "taxonomic_group", "group", "class_name", "class"],
        "类群",
    )
    area_col = find_column(
        raw,
        [
            "total_habitat_area_km2",
            "habitat_area_km2",
            "range_area_km2",
            "global_range_area_km2",
            "area_km2",
        ],
        "全球分布面积",
    )
    impact_col = find_column(
        raw,
        [
            "impact_percent",
            "cal_iap_percent",
            "cal_iap",
            "conflict_ratio",
            "conflict_percent",
            "overlap_percent",
            "impacted_percent",
        ],
        "铜矿空间暴露比例",
    )
    affected_col = find_column(
        raw,
        [
            "is_copper_affected",
            "copper_affected",
            "is_affected",
            "affected",
            "has_copper_overlap",
            "has_overlap",
        ],
        "铜矿空间暴露状态",
        required=False,
    )

    impact_percent, unit_used = convert_impact_to_percent(
        raw[impact_col], impact_col
    )
    if affected_col is not None:
        affected = parse_bool(raw[affected_col])
        affected_source = affected_col
    else:
        affected = impact_percent > EPS_PERCENT
        affected_source = f"{impact_col} > 0"

    species_id = (
        raw[id_col].apply(clean_id)
        if id_col is not None
        else raw[name_col].apply(clean_id)
    )
    scientific_name = (
        raw[name_col].astype(str).str.strip()
        if name_col is not None
        else species_id
    )

    df = pd.DataFrame(
        {
            "Species_ID": species_id,
            "Scientific_Name": scientific_name,
            "Taxon": normalize_taxon(raw[taxon_col]),
            "Range_Area_km2": numeric_series(raw[area_col]),
            "Impact_Percent": impact_percent,
            "Copper_Affected": affected.astype(bool),
        }
    )

    invalid_taxa = sorted(
        set(df.loc[~df["Taxon"].isin(TAXON_ORDER), "Taxon"].dropna())
    )
    if invalid_taxa:
        raise ValueError(
            "存在无法归入五个类群的值："
            f"{invalid_taxa[:20]}。请在 TAXON_ALIASES 中补充映射。"
        )

    invalid = (
        df["Species_ID"].isna()
        | df["Range_Area_km2"].isna()
        | (df["Range_Area_km2"] <= 0)
        | df["Impact_Percent"].isna()
    )
    if invalid.any():
        warnings.warn(
            f"已排除 {int(invalid.sum())} 条缺少物种 ID、有效面积或暴露比例的记录。"
        )
        df = df.loc[~invalid].copy()

    # 同一物种如存在完全一致的重复记录，仅保留一条；不一致则停止。
    duplicated = df[df.duplicated("Species_ID", keep=False)]
    if not duplicated.empty:
        check = duplicated.groupby("Species_ID").agg(
            n_taxon=("Taxon", "nunique"),
            n_area=("Range_Area_km2", "nunique"),
            n_impact=("Impact_Percent", "nunique"),
            n_affected=("Copper_Affected", "nunique"),
        )
        bad = check[
            (check["n_taxon"] > 1)
            | (check["n_area"] > 1)
            | (check["n_impact"] > 1)
            | (check["n_affected"] > 1)
        ]
        if not bad.empty:
            raise ValueError(
                "同一物种存在不一致的类群、面积、暴露比例或暴露状态。示例："
                f"{bad.index[:20].tolist()}"
            )
        df = df.drop_duplicates("Species_ID", keep="first")

    df = df[df["Copper_Affected"]].copy()
    df = df.sort_values(["Taxon", "Scientific_Name", "Species_ID"]).reset_index(
        drop=True
    )

    n_affected = len(df)
    if n_affected != EXPECTED_COPPER_AFFECTED_SPECIES:
        message = (
            f"识别到 {n_affected:,} 个铜矿空间暴露物种，"
            f"与预期 {EXPECTED_COPPER_AFFECTED_SPECIES:,} 不一致；"
            f"暴露状态来源：{affected_source}。"
        )
        if STRICT_EXPECTED_COUNT:
            raise ValueError(message)
        warnings.warn(message)

    metadata = {
        "input_rows": len(raw),
        "affected_rows": n_affected,
        "impact_column": impact_col,
        "impact_unit": unit_used,
        "affected_source": affected_source,
    }
    return df, metadata


# =============================================================================
# 4. 限制目标计算与优先物种判别
# =============================================================================
def representation_target_percent(area_km2, scheme, apply_cap=True):
    """
    计算保护代表性目标百分比。

    中间范围在 log(area) 尺度上线性插值；面积上限施加于保护目标面积：
        capped_target_area = min(relative_target * range_area, cap_area)
    """
    area = np.asarray(area_km2, dtype=float)
    rare_area = scheme["rare_area_km2"]
    common_area = scheme["common_area_km2"]
    rare_target = scheme["rare_target"]
    common_target = scheme["common_target"]

    relative = np.empty_like(area, dtype=float)
    rare = area <= rare_area
    common = area >= common_area
    middle = ~(rare | common)

    relative[rare] = rare_target
    relative[common] = common_target
    if middle.any():
        position = (
            np.log(area[middle]) - np.log(rare_area)
        ) / (np.log(common_area) - np.log(rare_area))
        relative[middle] = rare_target + position * (
            common_target - rare_target
        )

    target_area = relative * area
    if apply_cap:
        target_area = np.minimum(target_area, scheme["cap_area_km2"])
    return np.clip(target_area / area * 100.0, 0.0, 100.0)


def classify_range(area_km2, rare_area=200.0, common_area=250000.0):
    """按当前主分析锚点划分 Narrow、Intermediate 与 Widespread。"""
    area = np.asarray(area_km2, dtype=float)
    return np.select(
        [area <= rare_area, area < common_area],
        ["Narrow", "Intermediate"],
        default="Widespread",
    )


def calculate_results(df):
    result = df.copy()
    result["Current_Range_Class"] = classify_range(result["Range_Area_km2"])

    for prefix, scheme in SCHEMES.items():
        target = representation_target_percent(
            result["Range_Area_km2"], scheme, apply_cap=True
        )
        target_no_cap = representation_target_percent(
            result["Range_Area_km2"], scheme, apply_cap=False
        )
        limit = 100.0 - target
        limit_no_cap = 100.0 - target_no_cap

        result[f"{prefix}_Protection_Target_Percent"] = target
        result[f"{prefix}_Conflict_Limit_Percent"] = limit
        result[f"{prefix}_Exceedance_Margin_pp"] = (
            result["Impact_Percent"] - limit
        )
        result[f"{prefix}_Priority"] = (
            result["Impact_Percent"] > limit + EPS_PERCENT
        )
        result[f"{prefix}_Cap_Applied"] = target < target_no_cap - EPS_PERCENT
        result[f"{prefix}_Priority_WithoutCap"] = (
            result["Impact_Percent"] > limit_no_cap + EPS_PERCENT
        )
        result[f"{prefix}_Cap_Changes_Priority"] = (
            result[f"{prefix}_Priority"]
            != result[f"{prefix}_Priority_WithoutCap"]
        )

    main_priority = result[f"{MAIN}_Priority"]
    trad_priority = result[f"{TRAD}_Priority"]
    result["Robustness_Status"] = np.select(
        [
            main_priority & trad_priority,
            main_priority & ~trad_priority,
            ~main_priority & trad_priority,
        ],
        ["Stable core", "Current only", "Traditional only"],
        default="Neither",
    )

    for level in BOUNDARY_LEVELS_PP:
        suffix = str(level).replace(".", "p")
        result[f"Near_Current_{suffix}pp"] = (
            result[f"{MAIN}_Exceedance_Margin_pp"].abs() <= level
        )
        result[f"Near_Traditional_{suffix}pp"] = (
            result[f"{TRAD}_Exceedance_Margin_pp"].abs() <= level
        )
    return result


# =============================================================================
# 5. 汇总统计
# =============================================================================
def safe_percent(numerator, denominator):
    return numerator / denominator * 100.0 if denominator else np.nan


def scheme_summary(df):
    rows = []
    for prefix, scheme in SCHEMES.items():
        priority = df[f"{prefix}_Priority"]
        n_priority = int(priority.sum())
        ar_n = int(
            (priority & df["Taxon"].isin(["Amphibians", "Reptiles"])).sum()
        )
        rows.append(
            {
                "Scheme": prefix,
                "Scheme_Label": scheme["label"],
                "Copper_Affected_N": len(df),
                "Priority_N": n_priority,
                "Priority_Rate_Percent": safe_percent(n_priority, len(df)),
                "Amphibians_Reptiles_N": ar_n,
                "Amphibians_Reptiles_Share_Percent": safe_percent(
                    ar_n, n_priority
                ),
                "Narrow_Priority_N": int(
                    (priority & (df["Current_Range_Class"] == "Narrow")).sum()
                ),
                "Intermediate_Priority_N": int(
                    (
                        priority
                        & (df["Current_Range_Class"] == "Intermediate")
                    ).sum()
                ),
                "Widespread_Priority_N": int(
                    (
                        priority
                        & (df["Current_Range_Class"] == "Widespread")
                    ).sum()
                ),
                "Within_1pp_N": int(
                    (
                        df[f"{prefix}_Exceedance_Margin_pp"].abs()
                        <= BOUNDARY_LEVELS_PP[0]
                    ).sum()
                ),
                "Within_5pp_N": int(
                    (
                        df[f"{prefix}_Exceedance_Margin_pp"].abs()
                        <= BOUNDARY_LEVELS_PP[1]
                    ).sum()
                ),
                "Cap_Applied_N": int(df[f"{prefix}_Cap_Applied"].sum()),
                "Cap_Changes_Priority_N": int(
                    df[f"{prefix}_Cap_Changes_Priority"].sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def comparison_summary(df):
    main = df[f"{MAIN}_Priority"]
    trad = df[f"{TRAD}_Priority"]
    both = int((main & trad).sum())
    main_only = int((main & ~trad).sum())
    trad_only = int((~main & trad).sum())
    neither = int((~main & ~trad).sum())
    union = both + main_only + trad_only

    return pd.DataFrame(
        [
            {
                "Copper_Affected_N": len(df),
                "Current_Priority_N": int(main.sum()),
                "Traditional_Priority_N": int(trad.sum()),
                "Stable_Core_N": both,
                "Current_Only_N": main_only,
                "Traditional_Only_N": trad_only,
                "Neither_N": neither,
                "Priority_Union_N": union,
                "Jaccard_Index": both / union if union else np.nan,
                "Jaccard_Percent": safe_percent(both, union),
                "Current_Set_Retained_Percent": safe_percent(
                    both, int(main.sum())
                ),
                "Traditional_Set_Within_Current_Percent": safe_percent(
                    both, int(trad.sum())
                ),
                "Overall_Binary_Agreement_Percent": safe_percent(
                    both + neither, len(df)
                ),
            }
        ]
    )


def grouped_summary(df, group_col, group_order):
    main_total = int(df[f"{MAIN}_Priority"].sum())
    trad_total = int(df[f"{TRAD}_Priority"].sum())
    rows = []
    for group in group_order:
        sub = df[df[group_col] == group]
        main = sub[f"{MAIN}_Priority"]
        trad = sub[f"{TRAD}_Priority"]
        both = int((main & trad).sum())
        union = int((main | trad).sum())
        rows.append(
            {
                "Group": group,
                "Copper_Affected_N": len(sub),
                "Current_Priority_N": int(main.sum()),
                "Current_Priority_Rate_Percent": safe_percent(
                    int(main.sum()), len(sub)
                ),
                "Current_Composition_Percent": safe_percent(
                    int(main.sum()), main_total
                ),
                "Traditional_Priority_N": int(trad.sum()),
                "Traditional_Priority_Rate_Percent": safe_percent(
                    int(trad.sum()), len(sub)
                ),
                "Traditional_Composition_Percent": safe_percent(
                    int(trad.sum()), trad_total
                ),
                "Stable_Core_N": both,
                "Current_Only_N": int((main & ~trad).sum()),
                "Traditional_Only_N": int((~main & trad).sum()),
                "Group_Jaccard_Percent": safe_percent(both, union),
            }
        )
    return pd.DataFrame(rows)


def boundary_species(df):
    main_abs = df[f"{MAIN}_Exceedance_Margin_pp"].abs()
    trad_abs = df[f"{TRAD}_Exceedance_Margin_pp"].abs()
    selected = df[
        (main_abs <= BOUNDARY_LEVELS_PP[1])
        | (trad_abs <= BOUNDARY_LEVELS_PP[1])
    ].copy()
    selected["Minimum_Absolute_Margin_pp"] = np.minimum(
        selected[f"{MAIN}_Exceedance_Margin_pp"].abs(),
        selected[f"{TRAD}_Exceedance_Margin_pp"].abs(),
    )
    return selected.sort_values(
        ["Minimum_Absolute_Margin_pp", "Scientific_Name", "Species_ID"]
    )


# =============================================================================
# 6. 补充材料图
# =============================================================================






# =============================================================================
# 7. 输出
# =============================================================================
def write_run_summary(metadata, schemes, comparison):
    main_row = schemes.loc[schemes["Scheme"] == MAIN].iloc[0]
    trad_row = schemes.loc[schemes["Scheme"] == TRAD].iloc[0]
    comp = comparison.iloc[0]
    lines = [
        "限制目标稳健性分析运行摘要",
        "=" * 40,
        f"输入文件：{INPUT_FILE}",
        f"输出目录：{OUTPUT_DIR}",
        f"输入记录数：{metadata['input_rows']}",
        f"铜矿空间暴露物种数：{metadata['affected_rows']}",
        f"暴露比例字段：{metadata['impact_column']}",
        f"暴露比例单位判断：{metadata['impact_unit']}",
        f"铜矿暴露状态来源：{metadata['affected_source']}",
        "",
        f"当前目标优先物种数：{int(main_row['Priority_N'])}",
        f"传统目标优先物种数：{int(trad_row['Priority_N'])}",
        f"稳定核心集合：{int(comp['Stable_Core_N'])}",
        f"仅当前目标识别：{int(comp['Current_Only_N'])}",
        f"仅传统目标识别：{int(comp['Traditional_Only_N'])}",
        f"Jaccard 指数：{comp['Jaccard_Index']:.4f}",
        f"当前集合在传统目标下保留比例："
        f"{comp['Current_Set_Retained_Percent']:.2f}%",
        "",
        f"当前目标 1 个百分点内边界物种：{int(main_row['Within_1pp_N'])}",
        f"当前目标 5 个百分点内边界物种：{int(main_row['Within_5pp_N'])}",
        f"传统目标 1 个百分点内边界物种：{int(trad_row['Within_1pp_N'])}",
        f"传统目标 5 个百分点内边界物种：{int(trad_row['Within_5pp_N'])}",
        "",
        f"当前目标受 100 万 km² 面积上限影响的物种："
        f"{int(main_row['Cap_Applied_N'])}",
        f"面积上限改变当前目标优先判定的物种："
        f"{int(main_row['Cap_Changes_Priority_N'])}",
        f"传统目标受 100 万 km² 面积上限影响的物种："
        f"{int(trad_row['Cap_Applied_N'])}",
        f"面积上限改变传统目标优先判定的物种："
        f"{int(trad_row['Cap_Changes_Priority_N'])}",
    ]
    (OUTPUT_DIR / "Robustness_RunSummary.txt").write_text(
        "\n".join(lines), encoding="utf-8-sig"
    )


def main():
    print("=" * 72)
    print("研究内容二：限制目标稳健性分析")
    print("=" * 72)
    print(f"输入：{INPUT_FILE}")
    print(f"输出：{OUTPUT_DIR}")
    print("说明：本脚本不进行空间运算，只处理已有物种统计表。")

    species, metadata = read_affected_species()
    results = calculate_results(species)

    schemes = scheme_summary(results)
    comparison = comparison_summary(results)
    by_taxon = grouped_summary(results, "Taxon", TAXON_ORDER)
    by_range = grouped_summary(
        results, "Current_Range_Class", RANGE_ORDER
    )
    boundary = boundary_species(results)

    main_n = int(results[f"{MAIN}_Priority"].sum())
    if main_n != EXPECTED_MAIN_PRIORITY_SPECIES:
        message = (
            f"当前目标识别到 {main_n} 个优先物种，"
            f"与预期 {EXPECTED_MAIN_PRIORITY_SPECIES} 不一致。"
        )
        if STRICT_EXPECTED_COUNT:
            raise ValueError(message)
        warnings.warn(message)

    results.to_csv(
        OUTPUT_DIR / "Robustness_SpeciesResults.csv",
        index=False,
        encoding="utf-8-sig",
    )
    schemes.to_csv(
        OUTPUT_DIR / "Robustness_SchemeSummary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    comparison.to_csv(
        OUTPUT_DIR / "Robustness_ComparisonSummary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    by_taxon.to_csv(
        OUTPUT_DIR / "Robustness_ByTaxon.csv",
        index=False,
        encoding="utf-8-sig",
    )
    by_range.to_csv(
        OUTPUT_DIR / "Robustness_ByCurrentRange.csv",
        index=False,
        encoding="utf-8-sig",
    )
    boundary.to_csv(
        OUTPUT_DIR / "Robustness_BoundarySpecies.csv",
        index=False,
        encoding="utf-8-sig",
    )

    write_run_summary(metadata, schemes, comparison)

    comp = comparison.iloc[0]
    print("\n分析完成")
    print(f"  铜矿空间暴露物种：{len(results)}")
    print(f"  当前目标优先物种：{int(comp['Current_Priority_N'])}")
    print(f"  传统目标优先物种：{int(comp['Traditional_Priority_N'])}")
    print(f"  稳定核心集合：{int(comp['Stable_Core_N'])}")
    print(f"  Jaccard 指数：{comp['Jaccard_Index']:.4f}")
    print(f"  结果目录：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
