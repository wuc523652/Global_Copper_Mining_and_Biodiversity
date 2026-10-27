# -*- coding: utf-8 -*-
"""
07_screen_priority_species.py

功能：
基于候选物种全集，筛选栖息地受铜矿开采威胁且超过冲突限制目标的优先救助物种。

关键逻辑：
1. 从 IUCN 原始栖息地数据读取所有候选 CR/EN/VU 脊椎动物物种；
2. 使用 Mollweide 等面积投影 ESRI:54009 计算每个物种完整全球栖息地面积；
3. 从 04_calculate_mine_species_overlap.py 输出的 Species_CopperOverlap_GlobalUnique.csv
   读取全球物种尺度严格空间去重铜矿重叠面积；
4. 以候选物种全集 LEFT JOIN 全球去重重叠面积；
5. 未发生铜矿重叠的物种 total_overlap_area_km2 = 0，因此 Cal_IAP = 0；
6. 计算：
   Cal_IAP = total_overlap_area_km2 / total_habitat_area_km2
7. 严格按照 Sreekar & Watson (2026) 更新的物种代表性目标，将保护目标转换为
   铜矿空间冲突限制目标：Conflict_Target = 1 - Protection_Target：
   - Narrow:       A_i <= 200 km²，保护目标 = 100%，冲突上限 = 0%；
   - Intermediate: 200 < A_i < 250,000 km²，保护目标在 log(A_i) 空间由
                   100% 线性降至 30%，冲突上限相应由 0% 升至 70%；
   - Widespread:   A_i >= 250,000 km²，保护目标 = 30%，冲突上限 = 70%。
8. 只有 Cal_IAP > Conflict_Target 才识别为优先干预物种；等于阈值不算超限。

重要说明：
- Species_PrioritySpecies_Full.csv 是候选物种全集表；
- Species_CopperAffected_Full.csv 是发生铜矿重叠的物种表；
- Species_PrioritySpecies_Priority.csv 是优先救助物种表；
- 本版暂不计算“超限幅度”，也不再生成按国家制图所需的旧版结果。

方法来源：
Sreekar, R. & Watson, J. E. M. Updating species representation targets for
protected and conserved area planning. Nature Reviews Biodiversity (2026).
https://doi.org/10.1038/s44358-026-00169-7
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import unary_union

try:
    from shapely.validation import make_valid
except Exception:
    def make_valid(geom):
        if geom is None:
            return None
        try:
            return geom.buffer(0)
        except Exception:
            return geom


warnings.filterwarnings("ignore")


# ============================================================
# 配置参数
# ============================================================
BASE_DIR = Path(_nc_project(""))

INPUT_DIR = BASE_DIR / "04_mine_species_overlap"
IUCN_DIR = Path(_nc_data("IUCNredlist/202502"))

OUTPUT_DIR = BASE_DIR / "07_priority_species"
CACHE_DIR = OUTPUT_DIR / "cache_global_habitat_mollweide_v3_sw2026"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_CRS = "ESRI:54009"

TAXON_LIST = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
STATUS_LIST = ["CR", "EN", "VU"]

SPECIES_ID_FIELD_CANDIDATES = ["id_no", "ID_NO", "Id_no", "ID_No"]
SPECIES_NAME_FIELD_CANDIDATES = [
    "binomial", "BINOMIAL", "Binomial",
    "sci_name", "SCI_NAME", "scientific", "SCIENTIFIC",
]

# Sreekar & Watson (2026) 的更新锚点。
RANGE_NARROW_THRESHOLD = 200.0
RANGE_INTERMEDIATE_THRESHOLD = 250000.0

NARROW_CONFLICT_TARGET = 0.0
WIDESPREAD_CONFLICT_TARGET = 0.70
NARROW_PROTECTION_TARGET = 1.00
WIDESPREAD_PROTECTION_TARGET = 0.30

TARGET_RULE_VERSION = "Sreekar_Watson_2026"
TARGET_RULE_DOI = "10.1038/s44358-026-00169-7"

# 阈值比较容差：仅用于“是否超过冲突限制目标”的严格边界判断。
EPS = 1e-12

# 比例边界容差：Cal_IAP 理论上应位于 [0, 1]。
# 1e-8 的比例容差等价于 1e-6 个百分点，可容纳空间面积运算产生的
# 100.0000000000001% 一类浮点误差，但不会掩盖具有实际意义的越界。
CAL_IAP_BOUND_TOL = 1e-8

FORCE_RECALCULATE_HABITAT_AREA = True

IMPACT_PERCENT_BINS = [-np.inf, 0, 0.1, 1, 5, 10, 20, 50, 90, 100, np.inf]
IMPACT_PERCENT_LABELS = [
    "0",
    "(0,0.1]",
    "(0.1,1]",
    "(1,5]",
    "(5,10]",
    "(10,20]",
    "(20,50]",
    "(50,90]",
    "(90,100]",
    ">100"
]

RANGE_ORDER = ["Narrow", "Intermediate", "Widespread", "Unknown"]
STATUS_ORDER = ["CR", "EN", "VU"]
TAXON_ORDER = TAXON_LIST


# ============================================================
# 基础工具函数
# ============================================================
def clean_species_id(species_id):
    """清理物种 ID，移除浮点数末尾的 .0。"""
    if pd.isna(species_id):
        return None

    species_id_str = str(species_id).strip()

    if species_id_str == "" or species_id_str.lower() == "nan":
        return None

    if species_id_str.endswith(".0"):
        species_id_str = species_id_str[:-2]

    return species_id_str


def read_vector(path):
    """读取矢量文件，优先使用 pyogrio。"""
    try:
        return gpd.read_file(path, engine="pyogrio")
    except Exception:
        return gpd.read_file(path)


def detect_species_id_field(columns):
    """自动识别 IUCN 数据中的物种 ID 字段。"""
    for field in SPECIES_ID_FIELD_CANDIDATES:
        if field in columns:
            return field

    raise ValueError(
        f"未找到物种 ID 字段。候选字段为 {SPECIES_ID_FIELD_CANDIDATES}，"
        f"当前数据字段为：{list(columns)}"
    )


def detect_species_name_field(columns):
    """自动识别 IUCN 数据中的学名字段；找不到时返回 None。"""
    lower_map = {str(col).lower(): col for col in columns}
    for field in SPECIES_NAME_FIELD_CANDIDATES:
        if field in columns:
            return field
        if field.lower() in lower_map:
            return lower_map[field.lower()]
    return None


def ensure_projected(gdf, target_crs, name="vector"):
    """确保数据具有 CRS，并投影到目标 CRS。"""
    if gdf.crs is None:
        warnings.warn(
            f"[CRS] {name} 没有 CRS，已默认设置为 EPSG:4326。"
            f"请确认原始数据是否确实为经纬度坐标。"
        )
        gdf = gdf.set_crs(epsg=4326)

    if str(gdf.crs) != str(target_crs):
        print(f"  投影转换: {gdf.crs} -> {target_crs}")
        gdf = gdf.to_crs(target_crs)

    return gdf


def fix_geometries(gdf, name="vector"):
    """修复几何，并删除空几何。"""
    before = len(gdf)

    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[~gdf.geometry.is_empty].copy()

    def _fix_one(geom):
        if geom is None or geom.is_empty:
            return None
        try:
            return make_valid(geom)
        except Exception:
            try:
                return geom.buffer(0)
            except Exception:
                return None

    gdf["geometry"] = gdf.geometry.apply(_fix_one)

    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[~gdf.geometry.is_empty].copy()

    after = len(gdf)
    if after < before:
        print(f"  [GEOM] {name}: 删除无效或空几何 {before - after} 条")

    return gdf


def safe_union(geometries):
    """对同一物种的多个栖息地斑块进行空间合并，避免物种内部重叠重复计算。"""
    valid_geoms = [
        geom for geom in geometries
        if geom is not None and not geom.is_empty
    ]

    if len(valid_geoms) == 0:
        return None

    try:
        union_geom = unary_union(valid_geoms)
    except Exception:
        fixed_geoms = []
        for geom in valid_geoms:
            try:
                fixed_geoms.append(make_valid(geom))
            except Exception:
                try:
                    fixed_geoms.append(geom.buffer(0))
                except Exception:
                    pass

        if len(fixed_geoms) == 0:
            return None

        union_geom = unary_union(fixed_geoms)

    if union_geom is None or union_geom.is_empty:
        return None

    return union_geom


def count_and_percent(df, group_cols, within_cols=None, count_name="num_species"):
    """
    按 group_cols 统计物种数量，并计算百分比。
    """
    if df.empty:
        return pd.DataFrame(columns=group_cols + [count_name, "percent"])

    out = (
        df.groupby(group_cols, observed=False)["species_id"]
        .nunique()
        .reset_index(name=count_name)
    )

    if within_cols is None:
        total = out[count_name].sum()
        out["percent"] = np.where(total > 0, out[count_name] / total * 100, 0)
    else:
        denom = out.groupby(within_cols, observed=False)[count_name].transform("sum")
        out["percent"] = np.where(denom > 0, out[count_name] / denom * 100, 0)

    return out


def aggregate_impact_stats(df, group_cols):
    """统计重叠/冲突面积占比的均值、中位数、分位数等信息。"""
    if df.empty:
        return pd.DataFrame(columns=group_cols)

    def q25(x):
        return x.quantile(0.25)

    def q75(x):
        return x.quantile(0.75)

    out = (
        df.groupby(group_cols, observed=False)
        .agg(
            num_species=("species_id", "nunique"),
            num_copper_affected_species=("Is_Copper_Affected", "sum"),
            total_overlap_area_km2=("total_overlap_area_km2", "sum"),
            total_habitat_area_km2=("total_habitat_area_km2", "sum"),
            mean_impact_percent=("Impact_Percent", "mean"),
            median_impact_percent=("Impact_Percent", "median"),
            q25_impact_percent=("Impact_Percent", q25),
            q75_impact_percent=("Impact_Percent", q75),
            min_impact_percent=("Impact_Percent", "min"),
            max_impact_percent=("Impact_Percent", "max"),
            mean_conflict_target_percent=("Conflict_Target_Percent", "mean"),
            median_conflict_target_percent=("Conflict_Target_Percent", "median")
        )
        .reset_index()
    )

    out["aggregate_impact_percent"] = np.where(
        out["total_habitat_area_km2"] > 0,
        out["total_overlap_area_km2"] / out["total_habitat_area_km2"] * 100,
        np.nan
    )

    out["copper_affected_percent"] = np.where(
        out["num_species"] > 0,
        out["num_copper_affected_species"] / out["num_species"] * 100,
        np.nan
    )

    return out


# ============================================================
# 栖息地面积计算
# ============================================================
def calculate_species_habitat_area(taxon, status):
    """
    计算指定 taxon-status 组合中每个物种的完整全球栖息地面积。
    """
    print(f"\n处理完整全球栖息地面积: {taxon} - {status}")

    cache_file = CACHE_DIR / f"global_habitat_area_{taxon}_{status}.pkl"

    if cache_file.exists() and not FORCE_RECALCULATE_HABITAT_AREA:
        print(f"  从缓存读取: {cache_file.name}")
        return pd.read_pickle(cache_file)

    habitat_file = IUCN_DIR / f"{taxon}-{status}" / "data_0.shp"

    if not habitat_file.exists():
        print(f"  警告：IUCN 原始栖息地文件不存在: {habitat_file}")
        return pd.DataFrame(columns=["species_id", "scientific_name", "total_habitat_area_km2"])

    print(f"  读取 IUCN 原始栖息地数据: {habitat_file}")

    gdf = read_vector(habitat_file)

    if gdf.empty:
        print(f"  警告：{habitat_file} 为空")
        return pd.DataFrame(columns=["species_id", "scientific_name", "total_habitat_area_km2"])

    id_field = detect_species_id_field(gdf.columns)
    name_field = detect_species_name_field(gdf.columns)

    gdf["species_id"] = gdf[id_field].apply(clean_species_id)
    if name_field is None:
        gdf["scientific_name"] = pd.NA
        print("  提示：未识别到学名字段，scientific_name 将留空。")
    else:
        gdf["scientific_name"] = gdf[name_field].astype("string").str.strip()
        gdf.loc[
            gdf["scientific_name"].isin(["", "nan", "None", "<NA>"]),
            "scientific_name",
        ] = pd.NA
    gdf = gdf.dropna(subset=["species_id"]).copy()

    print(f"  物种数量: {gdf['species_id'].nunique():,}")

    gdf = ensure_projected(gdf, TARGET_CRS, name=f"{taxon}-{status}")
    gdf = fix_geometries(gdf, name=f"{taxon}-{status}")

    gdf = gdf[["species_id", "scientific_name", "geometry"]].copy()

    print("  按物种 dissolve/union 计算完整全球栖息地面积...")

    species_areas = []
    n_species = gdf["species_id"].nunique()

    for i, (species_id, sub) in enumerate(gdf.groupby("species_id"), 1):
        if i % 200 == 0 or i == 1:
            print(f"    {taxon}-{status}: {i}/{n_species}")

        union_geom = safe_union(sub.geometry)

        if union_geom is None:
            area_km2 = np.nan
        else:
            area_km2 = union_geom.area / 1e6

        names = sub["scientific_name"].dropna().astype(str)
        scientific_name = names.iloc[0] if not names.empty else pd.NA

        species_areas.append({
            "species_id": species_id,
            "scientific_name": scientific_name,
            "total_habitat_area_km2": area_km2
        })

    result_df = pd.DataFrame(species_areas)

    result_df.to_pickle(cache_file)

    print(f"  完成，共 {len(result_df):,} 个物种")
    print(f"  缓存已保存: {cache_file.name}")

    return result_df


def calculate_all_habitat_areas():
    """计算所有 taxon-status 组合中每个物种的完整全球栖息地面积。"""
    print("\n" + "=" * 80)
    print("Step 1: 计算所有候选物种的完整全球栖息地面积")
    print("=" * 80)

    all_habitat_areas = []

    for taxon in TAXON_LIST:
        for status in STATUS_LIST:
            df = calculate_species_habitat_area(taxon, status)

            if not df.empty:
                df["taxon"] = taxon
                df["status"] = status
                all_habitat_areas.append(df)

    if len(all_habitat_areas) == 0:
        raise ValueError("没有成功计算任何物种的栖息地面积，请检查 IUCN_DIR 路径。")

    result_df = pd.concat(all_habitat_areas, ignore_index=True)

    result_df["species_id"] = result_df["species_id"].apply(clean_species_id)
    result_df = result_df.dropna(subset=["species_id"]).copy()

    result_df = (
        result_df
        .sort_values(["taxon", "status", "species_id"])
        .drop_duplicates(subset=["species_id", "taxon", "status"], keep="first")
        .copy()
    )

    print(f"\n候选物种全集记录数: {len(result_df):,}")
    print(f"候选物种唯一 ID 数: {result_df['species_id'].nunique():,}")

    return result_df


# ============================================================
# 铜矿重叠面积数据读取
# ============================================================
def load_species_overlap_data():
    """
    读取全球物种尺度严格空间去重铜矿重叠面积：
    Species_CopperOverlap_GlobalUnique.csv

    新大纲不再绘制旧版国家优先物种地图，因此本步骤不再依赖
    Country_CopperOverlap_BySpecies.csv。global-unique 表中已有的国家计数字段
    仍会原样保留，供物种名单表参考。
    """
    print("\n" + "=" * 80)
    print("Step 2: 加载全球物种尺度严格空间去重铜矿重叠面积")
    print("=" * 80)

    global_overlap_file = INPUT_DIR / "Species_CopperOverlap_GlobalUnique.csv"

    if not global_overlap_file.exists():
        raise FileNotFoundError(
            f"输入文件不存在: {global_overlap_file}\n"
            "请先运行 04_calculate_mine_species_overlap.py，"
            "生成 Species_CopperOverlap_GlobalUnique.csv。"
        )

    df = pd.read_csv(global_overlap_file)

    required_cols = [
        "species_id",
        "taxon",
        "status",
        "global_unique_overlap_area_km2"
    ]

    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"{global_overlap_file.name} 缺少必要字段: {col}")

    df["species_id"] = df["species_id"].apply(clean_species_id)

    df["global_unique_overlap_area_km2"] = pd.to_numeric(
        df["global_unique_overlap_area_km2"], errors="coerce"
    ).fillna(0)

    species_overlap = df.rename(
        columns={
            "global_unique_overlap_area_km2": "total_overlap_area_km2"
        }
    )

    if "num_countries" not in species_overlap.columns:
        species_overlap["num_countries"] = 0

    if "country_occurrences" not in species_overlap.columns:
        species_overlap["country_occurrences"] = 0

    keep_cols = [
        "species_id",
        "taxon",
        "status",
        "total_overlap_area_km2",
        "num_countries",
        "country_occurrences"
    ]

    if "country_accumulated_overlap_area_km2" in species_overlap.columns:
        keep_cols.append("country_accumulated_overlap_area_km2")

    if "overlap_duplicate_area_km2" in species_overlap.columns:
        keep_cols.append("overlap_duplicate_area_km2")

    if "countries" in species_overlap.columns:
        keep_cols.append("countries")

    species_overlap = species_overlap[keep_cols].copy()

    print(f"  读取全球去重重叠物种数: {len(species_overlap):,}")
    print(f"  类群分布:\n{species_overlap['taxon'].value_counts()}")
    print(f"  濒危等级分布:\n{species_overlap['status'].value_counts()}")

    return species_overlap


# ============================================================
# 优先救助筛选规则
# ============================================================
def classify_range(area_km2):
    if pd.isna(area_km2) or area_km2 <= 0:
        return "Unknown"

    if area_km2 <= RANGE_NARROW_THRESHOLD:
        return "Narrow"
    elif area_km2 < RANGE_INTERMEDIATE_THRESHOLD:
        return "Intermediate"
    else:
        return "Widespread"


def calculate_protection_target(area_km2):
    """
    根据物种完整全球栖息地面积计算代表性保护目标 T_i（比例，0--1）。

    严格采用 Sreekar & Watson (2026) 更新的两个锚点：
    - A_i <= 200 km²: T_i = 1.00；
    - 200 < A_i < 250,000 km²: T_i 在 log(A_i) 空间由 1.00 降至 0.30；
    - A_i >= 250,000 km²: T_i = 0.30。

    中间区间公式：
    T_i = 1 - 0.70 * [ln(A_i)-ln(200)] / [ln(250000)-ln(200)]
    """
    if pd.isna(area_km2) or area_km2 <= 0:
        return np.nan

    if area_km2 <= RANGE_NARROW_THRESHOLD:
        return NARROW_PROTECTION_TARGET
    if area_km2 >= RANGE_INTERMEDIATE_THRESHOLD:
        return WIDESPREAD_PROTECTION_TARGET

    log_fraction = (
        (np.log(area_km2) - np.log(RANGE_NARROW_THRESHOLD))
        / (np.log(RANGE_INTERMEDIATE_THRESHOLD) - np.log(RANGE_NARROW_THRESHOLD))
    )
    return NARROW_PROTECTION_TARGET - (
        NARROW_PROTECTION_TARGET - WIDESPREAD_PROTECTION_TARGET
    ) * log_fraction


def calculate_conflict_target(area_km2):
    """
    将代表性保护目标严格转换为最大允许空间冲突比例 L_i。

    L_i = 1 - T_i，因此：
    - A_i <= 200 km²: L_i = 0；
    - 200 < A_i < 250,000 km²:
      L_i = 0.70 * [ln(A_i)-ln(200)]/[ln(250000)-ln(200)]；
    - A_i >= 250,000 km²: L_i = 0.70。

    使用 ln、log10 或任意同底对数会得到相同插值比例。
    """
    protection_target = calculate_protection_target(area_km2)
    if pd.isna(protection_target):
        return np.nan
    return 1.0 - protection_target


def determine_priority(row):
    """判断是否超过冲突限制目标：只有 C_i > L_i 时才识别为优先救助物种。"""
    cal_iap = row["Cal_IAP"]
    target = row["Conflict_Target"]

    if pd.isna(cal_iap) or pd.isna(target):
        return "Unknown"

    if cal_iap > target + EPS:
        return "Yes"

    return "No"


def determine_meet_target(row):
    """
    判断是否满足/未超过冲突限制目标。

    注意：本版本中只有 C_i > L_i 才为 Priority = Yes，
    因此 C_i <= L_i 时视为 Meet_Target = Yes。
    """
    cal_iap = row["Cal_IAP"]
    target = row["Conflict_Target"]

    if pd.isna(cal_iap) or pd.isna(target):
        return "Unknown"

    if cal_iap <= target + EPS:
        return "Yes"

    return "No"


def cleanup_stale_boundary_qa_files():
    """
    删除仅与 Cal_IAP 边界检查有关的旧版/上次运行 QA 文件。

    这样，当本次运行没有真正越界记录时，输出目录中不会残留旧的警告表。
    这里只处理明确列出的派生 QA 文件，不删除任何输入或主要结果。
    """
    stale_names = [
        "Warning_Cal_IAP_GreaterThan1.csv",  # 旧版文件名
        "Warning_Cal_IAP_Raw_GreaterThan1.csv",
        "Warning_Cal_IAP_Raw_LessThan0.csv",
        "Cal_IAP_Boundary_QA.csv",
    ]
    for name in stale_names:
        path = OUTPUT_DIR / name
        if path.exists() and path.is_file():
            path.unlink()


# ============================================================
# 结果统计表生成
# ============================================================
def generate_all_candidate_range_impact_excel(merged_df):
    print("\n生成 Excel: 候选物种全集范围分布与冲突比例统计...")

    output_xlsx = OUTPUT_DIR / "AllCandidateSpecies_Range_Impact_ByTaxon.xlsx"

    df = merged_df.copy()

    df["Range"] = pd.Categorical(df["Range"], categories=RANGE_ORDER, ordered=True)
    df["taxon"] = pd.Categorical(df["taxon"], categories=TAXON_ORDER, ordered=True)

    range_by_taxon = count_and_percent(
        df,
        group_cols=["taxon", "Range"],
        within_cols=["taxon"],
        count_name="num_species"
    )

    range_overall = count_and_percent(
        df,
        group_cols=["Range"],
        within_cols=None,
        count_name="num_species"
    )

    impact_bin_by_taxon = count_and_percent(
        df,
        group_cols=["taxon", "Impact_Bin"],
        within_cols=["taxon"],
        count_name="num_species"
    )

    range_impact_by_taxon = count_and_percent(
        df,
        group_cols=["taxon", "Range", "Impact_Bin"],
        within_cols=["taxon", "Range"],
        count_name="num_species"
    )

    impact_stats_by_taxon_range = aggregate_impact_stats(
        df,
        group_cols=["taxon", "Range"]
    )

    impact_stats_by_taxon = aggregate_impact_stats(
        df,
        group_cols=["taxon"]
    )

    affected_df = df[df["Is_Copper_Affected"]].copy()

    affected_range_by_taxon = count_and_percent(
        affected_df,
        group_cols=["taxon", "Range"],
        within_cols=["taxon"],
        count_name="num_copper_affected_species"
    )

    affected_impact_stats_by_taxon_range = aggregate_impact_stats(
        affected_df,
        group_cols=["taxon", "Range"]
    )

    detail_cols = [
        "species_id", "scientific_name", "taxon", "status",
        "total_habitat_area_km2",
        "total_overlap_area_km2",
        "Cal_IAP_raw", "Cal_IAP",
        "Impact_Percent_raw", "Impact_Percent",
        "Range",
        "Protection_Target", "Protection_Target_Percent",
        "Conflict_Target", "Conflict_Target_Percent",
        "Priority", "Meet_Target",
        "Is_Copper_Affected",
        "num_countries",
        "country_occurrences",
        "Impact_Bin",
        "Warning_Flag", "Target_Rule_Version", "Target_Rule_DOI"
    ]

    extra_cols = [
        "country_accumulated_overlap_area_km2",
        "overlap_duplicate_area_km2",
        "countries"
    ]

    for col in extra_cols:
        if col in df.columns and col not in detail_cols:
            detail_cols.append(col)

    detail_cols = [c for c in detail_cols if c in df.columns]
    detail = df[detail_cols].copy()

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        range_by_taxon.to_excel(writer, sheet_name="Range_by_Taxon", index=False)
        range_overall.to_excel(writer, sheet_name="Range_Overall", index=False)
        impact_bin_by_taxon.to_excel(writer, sheet_name="ImpactBin_by_Taxon", index=False)
        range_impact_by_taxon.to_excel(writer, sheet_name="Range_Impact_Taxon", index=False)
        impact_stats_by_taxon_range.to_excel(writer, sheet_name="ImpactStats_TaxonRange", index=False)
        impact_stats_by_taxon.to_excel(writer, sheet_name="ImpactStats_Taxon", index=False)
        affected_range_by_taxon.to_excel(writer, sheet_name="Affected_Range_Taxon", index=False)
        affected_impact_stats_by_taxon_range.to_excel(writer, sheet_name="Affected_ImpactStats", index=False)
        detail.to_excel(writer, sheet_name="AllCandidate_Detail", index=False)

    print(f"  已保存: {output_xlsx}")

    return {
        "range_by_taxon": range_by_taxon,
        "range_overall": range_overall,
        "impact_bin_by_taxon": impact_bin_by_taxon,
        "range_impact_by_taxon": range_impact_by_taxon,
        "impact_stats_by_taxon_range": impact_stats_by_taxon_range,
        "impact_stats_by_taxon": impact_stats_by_taxon,
        "affected_range_by_taxon": affected_range_by_taxon,
        "affected_impact_stats_by_taxon_range": affected_impact_stats_by_taxon_range
    }


def generate_priority_range_status_excel(priority_species):
    print("\n生成 Excel: 优先救助物种范围分布与濒危等级统计...")

    output_xlsx = OUTPUT_DIR / "PrioritySpecies_Range_Status_ByAllAndTaxon.xlsx"

    df = priority_species.copy()

    df["Range"] = pd.Categorical(df["Range"], categories=RANGE_ORDER, ordered=True)
    df["status"] = pd.Categorical(df["status"], categories=STATUS_ORDER, ordered=True)
    df["taxon"] = pd.Categorical(df["taxon"], categories=TAXON_ORDER, ordered=True)

    overall_range = count_and_percent(
        df,
        group_cols=["Range"],
        within_cols=None,
        count_name="num_priority_species"
    )

    overall_status = count_and_percent(
        df,
        group_cols=["status"],
        within_cols=None,
        count_name="num_priority_species"
    )

    overall_range_status = count_and_percent(
        df,
        group_cols=["Range", "status"],
        within_cols=None,
        count_name="num_priority_species"
    )

    taxon_total = count_and_percent(
        df,
        group_cols=["taxon"],
        within_cols=None,
        count_name="num_priority_species"
    )

    taxon_range = count_and_percent(
        df,
        group_cols=["taxon", "Range"],
        within_cols=["taxon"],
        count_name="num_priority_species"
    )

    taxon_status = count_and_percent(
        df,
        group_cols=["taxon", "status"],
        within_cols=["taxon"],
        count_name="num_priority_species"
    )

    taxon_range_status = count_and_percent(
        df,
        group_cols=["taxon", "Range", "status"],
        within_cols=["taxon"],
        count_name="num_priority_species"
    )

    detail_cols = [
        "species_id", "scientific_name", "taxon", "status",
        "total_habitat_area_km2",
        "total_overlap_area_km2",
        "Cal_IAP_raw", "Cal_IAP",
        "Impact_Percent_raw", "Impact_Percent",
        "Range",
        "Protection_Target", "Protection_Target_Percent",
        "Conflict_Target", "Conflict_Target_Percent",
        "Is_Copper_Affected",
        "num_countries",
        "country_occurrences",
        "Impact_Bin",
        "Warning_Flag", "Target_Rule_Version", "Target_Rule_DOI"
    ]

    detail_cols = [c for c in detail_cols if c in df.columns]
    detail = df[detail_cols].copy()

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        overall_range.to_excel(writer, sheet_name="Overall_Range", index=False)
        overall_status.to_excel(writer, sheet_name="Overall_Status", index=False)
        overall_range_status.to_excel(writer, sheet_name="Overall_Range_Status", index=False)
        taxon_total.to_excel(writer, sheet_name="Taxon_Total", index=False)
        taxon_range.to_excel(writer, sheet_name="Taxon_Range", index=False)
        taxon_status.to_excel(writer, sheet_name="Taxon_Status", index=False)
        taxon_range_status.to_excel(writer, sheet_name="Taxon_Range_Status", index=False)
        detail.to_excel(writer, sheet_name="Priority_Detail", index=False)

    print(f"  已保存: {output_xlsx}")

    return {
        "overall_range": overall_range,
        "overall_status": overall_status,
        "overall_range_status": overall_range_status,
        "taxon_total": taxon_total,
        "taxon_range": taxon_range,
        "taxon_status": taxon_status,
        "taxon_range_status": taxon_range_status
    }


# ============================================================
# 主分析流程
# ============================================================
def generate_priority_species_analysis():
    print("\n" + "=" * 80)
    print("开始优先救助物种分析：候选物种全集 + 全球物种去重重叠面积")
    print("=" * 80)

    cleanup_stale_boundary_qa_files()

    habitat_areas = calculate_all_habitat_areas()

    overlap_data = load_species_overlap_data()

    print("\n" + "=" * 80)
    print("Step 3: 候选物种全集合并全球去重铜矿重叠面积并计算冲突比例")
    print("=" * 80)

    merged_df = pd.merge(
        habitat_areas,
        overlap_data,
        on=["species_id", "taxon", "status"],
        how="left"
    )

    check_unmatched = overlap_data.merge(
        habitat_areas[["species_id", "taxon", "status"]],
        on=["species_id", "taxon", "status"],
        how="left",
        indicator=True
    )

    unmatched_overlap = check_unmatched[check_unmatched["_merge"] == "left_only"].copy()

    if not unmatched_overlap.empty:
        unmatched_file = OUTPUT_DIR / "OverlapSpecies_NotInCandidateHabitatAreas.csv"
        unmatched_overlap.to_csv(unmatched_file, index=False, encoding="utf-8-sig")
        print(
            f"  警告：有 {len(unmatched_overlap):,} 条 overlap 物种记录未在候选物种面积表中找到，"
            f"已输出: {unmatched_file}"
        )

    fill_zero_cols = ["total_overlap_area_km2", "num_countries", "country_occurrences"]
    for col in fill_zero_cols:
        if col not in merged_df.columns:
            merged_df[col] = 0
        merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce").fillna(0)

    optional_numeric_cols = [
        "country_accumulated_overlap_area_km2",
        "overlap_duplicate_area_km2"
    ]

    for col in optional_numeric_cols:
        if col in merged_df.columns:
            merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce").fillna(0)

    merged_df["total_habitat_area_km2"] = pd.to_numeric(
        merged_df["total_habitat_area_km2"], errors="coerce"
    )

    merged_df["Is_Copper_Affected"] = merged_df["total_overlap_area_km2"] > 0

    # 先保存未经截断的原始比例，供 QA 和复现检查使用。
    merged_df["Cal_IAP_raw"] = np.where(
        merged_df["total_habitat_area_km2"] > 0,
        merged_df["total_overlap_area_km2"] / merged_df["total_habitat_area_km2"],
        np.nan
    )

    # 理论上 Cal_IAP 必须位于 [0, 1]。正式分析字段采用截断后的值，
    # 从而避免极小的浮点误差产生 >100% 的 Impact_Percent；原始值仍保留。
    merged_df["Cal_IAP"] = merged_df["Cal_IAP_raw"].clip(lower=0.0, upper=1.0)

    merged_df["Range"] = merged_df["total_habitat_area_km2"].apply(classify_range)

    merged_df["Protection_Target"] = merged_df["total_habitat_area_km2"].apply(
        calculate_protection_target
    )

    merged_df["Conflict_Target"] = merged_df["total_habitat_area_km2"].apply(
        calculate_conflict_target
    )

    merged_df["Impact_Percent_raw"] = merged_df["Cal_IAP_raw"] * 100
    merged_df["Impact_Percent"] = merged_df["Cal_IAP"] * 100
    merged_df["Protection_Target_Percent"] = merged_df["Protection_Target"] * 100
    merged_df["Conflict_Target_Percent"] = merged_df["Conflict_Target"] * 100

    # 将方法版本写入每一行，避免后续表格与旧规则结果混淆。
    merged_df["Target_Rule_Version"] = TARGET_RULE_VERSION
    merged_df["Target_Rule_DOI"] = TARGET_RULE_DOI

    merged_df["Priority"] = merged_df.apply(determine_priority, axis=1)
    merged_df["Meet_Target"] = merged_df.apply(determine_meet_target, axis=1)

    merged_df["Impact_Bin"] = pd.cut(
        merged_df["Impact_Percent"],
        bins=IMPACT_PERCENT_BINS,
        labels=IMPACT_PERCENT_LABELS,
        right=True,
        include_lowest=True
    )

    merged_df["Warning_Flag"] = ""

    merged_df.loc[
        merged_df["total_habitat_area_km2"].isna(),
        "Warning_Flag"
    ] = merged_df["Warning_Flag"] + "Missing_habitat_area;"

    merged_df.loc[
        merged_df["total_habitat_area_km2"] <= 0,
        "Warning_Flag"
    ] = merged_df["Warning_Flag"] + "Invalid_habitat_area;"

    true_over_one = merged_df["Cal_IAP_raw"] > (1.0 + CAL_IAP_BOUND_TOL)
    true_below_zero = merged_df["Cal_IAP_raw"] < (0.0 - CAL_IAP_BOUND_TOL)
    floating_point_clipped = (
        merged_df["Cal_IAP_raw"].notna()
        & ~true_over_one
        & ~true_below_zero
        & (
            (merged_df["Cal_IAP_raw"] > 1.0)
            | (merged_df["Cal_IAP_raw"] < 0.0)
        )
    )

    merged_df.loc[
        true_over_one,
        "Warning_Flag"
    ] = merged_df["Warning_Flag"] + "Cal_IAP_raw_greater_than_1_beyond_tolerance;"

    merged_df.loc[
        true_below_zero,
        "Warning_Flag"
    ] = merged_df["Warning_Flag"] + "Cal_IAP_raw_less_than_0_beyond_tolerance;"

    print(f"  候选物种记录数: {len(merged_df):,}")
    print(f"  已发生铜矿重叠物种数: {merged_df['Is_Copper_Affected'].sum():,}")
    print(f"  未发生铜矿重叠物种数: {(~merged_df['Is_Copper_Affected']).sum():,}")

    print("\nRange 分布：")
    print(merged_df["Range"].value_counts(dropna=False))

    print("\nPriority 分布：")
    print(merged_df["Priority"].value_counts(dropna=False))

    print("\nMeet_Target 分布：")
    print(merged_df["Meet_Target"].value_counts(dropna=False))

    n_true_over_1 = int(true_over_one.sum())
    n_true_below_0 = int(true_below_zero.sum())
    n_floating_point_clipped = int(floating_point_clipped.sum())

    qa_bounds = pd.DataFrame([{
        "n_candidate_species": int(len(merged_df)),
        "cal_iap_bound_tolerance": CAL_IAP_BOUND_TOL,
        "n_raw_above_1_within_tolerance": int(
            ((merged_df["Cal_IAP_raw"] > 1.0) & ~true_over_one).sum()
        ),
        "n_raw_above_1_beyond_tolerance": n_true_over_1,
        "n_raw_below_0_within_tolerance": int(
            ((merged_df["Cal_IAP_raw"] < 0.0) & ~true_below_zero).sum()
        ),
        "n_raw_below_0_beyond_tolerance": n_true_below_0,
        "n_rows_clipped_within_tolerance": n_floating_point_clipped,
        "maximum_cal_iap_raw": float(merged_df["Cal_IAP_raw"].max(skipna=True)),
        "minimum_cal_iap_raw": float(merged_df["Cal_IAP_raw"].min(skipna=True)),
        "maximum_impact_percent_raw": float(
            merged_df["Impact_Percent_raw"].max(skipna=True)
        ),
        "formal_cal_iap_range": "[0,1] after clip",
        "formal_impact_percent_range": "[0,100] after clip",
    }])
    qa_bounds_file = OUTPUT_DIR / "Cal_IAP_Boundary_QA.csv"
    qa_bounds.to_csv(qa_bounds_file, index=False, encoding="utf-8-sig")
    print(f"\nCal_IAP 边界 QA 已保存: {qa_bounds_file}")
    print(f"  容差内浮点截断记录数: {n_floating_point_clipped}")
    print(f"  真正 > 1 的记录数: {n_true_over_1}")
    print(f"  真正 < 0 的记录数: {n_true_below_0}")

    if n_true_over_1 > 0:
        over1_file = OUTPUT_DIR / "Warning_Cal_IAP_Raw_GreaterThan1.csv"
        merged_df[true_over_one].to_csv(
            over1_file,
            index=False,
            encoding="utf-8-sig"
        )
        print(
            f"\n警告：有 {n_true_over_1} 个物种的 Cal_IAP_raw "
            f"> 1 + tolerance，已输出: {over1_file}"
        )
    else:
        print("\n检查通过：没有物种的 Cal_IAP_raw 真正超过 1。")

    if n_true_below_0 > 0:
        below0_file = OUTPUT_DIR / "Warning_Cal_IAP_Raw_LessThan0.csv"
        merged_df[true_below_zero].to_csv(
            below0_file,
            index=False,
            encoding="utf-8-sig"
        )
        print(
            f"警告：有 {n_true_below_0} 个物种的 Cal_IAP_raw "
            f"< 0 - tolerance，已输出: {below0_file}"
        )

    print("\n" + "=" * 80)
    print("Step 4: 保存物种完整结果、受铜矿影响物种结果与优先救助物种清单")
    print("=" * 80)

    full_cols = [
        "species_id", "scientific_name", "taxon", "status",
        "total_habitat_area_km2",
        "total_overlap_area_km2",
        "Cal_IAP_raw", "Cal_IAP",
        "Impact_Percent_raw", "Impact_Percent",
        "Range",
        "Protection_Target", "Protection_Target_Percent",
        "Conflict_Target", "Conflict_Target_Percent",
        "Priority", "Meet_Target",
        "Is_Copper_Affected",
        "num_countries",
        "country_occurrences",
        "Impact_Bin",
        "Warning_Flag",
        "Target_Rule_Version", "Target_Rule_DOI"
    ]

    for col in [
        "country_accumulated_overlap_area_km2",
        "overlap_duplicate_area_km2",
        "countries"
    ]:
        if col in merged_df.columns:
            full_cols.append(col)

    full_cols = [c for c in full_cols if c in merged_df.columns]

    full_output = (
        merged_df[full_cols]
        .sort_values(
            ["Is_Copper_Affected", "Priority", "Impact_Percent"],
            ascending=[False, True, False]
        )
    )

    full_output_file = OUTPUT_DIR / "Species_PrioritySpecies_Full.csv"
    full_output.to_csv(full_output_file, index=False, encoding="utf-8-sig")
    print(f"  已保存候选全集完整表: {full_output_file}")
    print(f"  包含候选物种数: {len(full_output):,}")

    affected_species = merged_df[merged_df["Is_Copper_Affected"]].copy()
    affected_output_file = OUTPUT_DIR / "Species_CopperAffected_Full.csv"
    affected_species[full_cols].to_csv(affected_output_file, index=False, encoding="utf-8-sig")
    print(f"  已保存已发生铜矿重叠物种表: {affected_output_file}")
    print(f"  受铜矿影响物种数: {len(affected_species):,}")

    priority_species = merged_df[merged_df["Priority"] == "Yes"].copy()

    priority_cols = [
        "species_id", "scientific_name", "taxon", "status",
        "total_habitat_area_km2",
        "total_overlap_area_km2",
        "Cal_IAP_raw", "Cal_IAP",
        "Impact_Percent_raw", "Impact_Percent",
        "Range",
        "Protection_Target", "Protection_Target_Percent",
        "Conflict_Target", "Conflict_Target_Percent",
        "Is_Copper_Affected",
        "num_countries",
        "country_occurrences",
        "Impact_Bin",
        "Warning_Flag",
        "Target_Rule_Version", "Target_Rule_DOI"
    ]

    priority_cols = [c for c in priority_cols if c in priority_species.columns]

    priority_output = (
        priority_species[priority_cols]
        .sort_values(["Impact_Percent", "total_overlap_area_km2"], ascending=[False, False])
    )

    priority_output_file = OUTPUT_DIR / "Species_PrioritySpecies_Priority.csv"
    priority_output.to_csv(priority_output_file, index=False, encoding="utf-8-sig")
    print(f"  已保存优先救助物种清单: {priority_output_file}")
    print(f"  优先救助物种数: {len(priority_output):,}")

    print("\n" + "=" * 80)
    print("Step 5: 生成用于绘图的 Excel 统计表")
    print("=" * 80)

    all_candidate_stats = generate_all_candidate_range_impact_excel(merged_df)

    priority_stats = generate_priority_range_status_excel(priority_species)

    print("\n" + "=" * 80)
    print("Step 6: 保存 CSV 汇总统计")
    print("=" * 80)

    summary_taxon_status = (
        priority_species.groupby(["taxon", "status"], observed=False)
        .agg(
            num_priority_species=("species_id", "nunique"),
            total_overlap_area_km2=("total_overlap_area_km2", "sum"),
            total_habitat_area_km2=("total_habitat_area_km2", "sum"),
            total_country_occurrences=("num_countries", "sum"),
            mean_impact_percent=("Impact_Percent", "mean"),
            median_impact_percent=("Impact_Percent", "median"),
            mean_conflict_target_percent=("Conflict_Target_Percent", "mean"),
            median_conflict_target_percent=("Conflict_Target_Percent", "median")
        )
        .reset_index()
    )

    summary_taxon_status["aggregate_impact_percent"] = np.where(
        summary_taxon_status["total_habitat_area_km2"] > 0,
        summary_taxon_status["total_overlap_area_km2"] / summary_taxon_status["total_habitat_area_km2"] * 100,
        np.nan
    )

    summary_taxon_status_file = OUTPUT_DIR / "Priority_Summary_By_Taxon_Status.csv"
    summary_taxon_status.to_csv(summary_taxon_status_file, index=False, encoding="utf-8-sig")
    print(f"  已保存: {summary_taxon_status_file}")

    summary_range = (
        priority_species.groupby("Range", observed=False)
        .agg(
            num_priority_species=("species_id", "nunique"),
            total_overlap_area_km2=("total_overlap_area_km2", "sum"),
            total_habitat_area_km2=("total_habitat_area_km2", "sum"),
            mean_impact_percent=("Impact_Percent", "mean"),
            median_impact_percent=("Impact_Percent", "median"),
            mean_conflict_target_percent=("Conflict_Target_Percent", "mean"),
            median_conflict_target_percent=("Conflict_Target_Percent", "median")
        )
        .reset_index()
    )

    summary_range["aggregate_impact_percent"] = np.where(
        summary_range["total_habitat_area_km2"] > 0,
        summary_range["total_overlap_area_km2"] / summary_range["total_habitat_area_km2"] * 100,
        np.nan
    )

    summary_range_file = OUTPUT_DIR / "Priority_Summary_By_Range.csv"
    summary_range.to_csv(summary_range_file, index=False, encoding="utf-8-sig")
    print(f"  已保存: {summary_range_file}")

    print("\n" + "=" * 80)
    print("分析完成 - 关键统计")
    print("=" * 80)

    total_candidate = len(merged_df)
    total_affected = int(merged_df["Is_Copper_Affected"].sum())
    total_priority = len(priority_species)

    print(f"候选物种全集数量: {total_candidate:,}")
    print(f"已发生铜矿重叠物种数: {total_affected:,}")
    print(f"优先救助物种数: {total_priority:,}")

    if total_candidate > 0:
        print(f"铜矿重叠物种占候选全集比例: {total_affected / total_candidate * 100:.2f}%")
        print(f"优先救助物种占候选全集比例: {total_priority / total_candidate * 100:.2f}%")

    if total_affected > 0:
        print(f"优先救助物种占已发生铜矿重叠物种比例: {total_priority / total_affected * 100:.2f}%")

    print("\n候选全集按类群统计:")
    print(merged_df["taxon"].value_counts())

    print("\n已发生铜矿重叠物种按类群统计:")
    print(affected_species["taxon"].value_counts())

    print("\n优先救助物种按类群统计:")
    print(priority_species["taxon"].value_counts())

    print("\n优先救助物种按濒危等级统计:")
    print(priority_species["status"].value_counts())

    print("\n优先救助物种按分布范围统计:")
    print(priority_species["Range"].value_counts())

    print("\n输出目录:")
    print(OUTPUT_DIR)

    print("\n关键输出文件:")
    print(f"1. {OUTPUT_DIR / 'Species_PrioritySpecies_Full.csv'}")
    print(f"2. {OUTPUT_DIR / 'Species_CopperAffected_Full.csv'}")
    print(f"3. {OUTPUT_DIR / 'Species_PrioritySpecies_Priority.csv'}")
    print(f"4. {OUTPUT_DIR / 'AllCandidateSpecies_Range_Impact_ByTaxon.xlsx'}")
    print(f"5. {OUTPUT_DIR / 'PrioritySpecies_Range_Status_ByAllAndTaxon.xlsx'}")
    print(f"6. {OUTPUT_DIR / 'Priority_Summary_By_Taxon_Status.csv'}")
    print(f"7. {OUTPUT_DIR / 'Priority_Summary_By_Range.csv'}")
    print(f"8. {OUTPUT_DIR / 'Cal_IAP_Boundary_QA.csv'}")

    print("=" * 80)

    return merged_df, priority_species


# ============================================================
# 程序入口
# ============================================================
if __name__ == "__main__":
    try:
        full_data, priority_data = generate_priority_species_analysis()

        print("\n✓ 所有分析完成！")
        print(f"✓ 结果已保存至: {OUTPUT_DIR}")

    except Exception as e:
        print(f"\n✗ 错误: {e}")
        import traceback
        traceback.print_exc()
        raise
