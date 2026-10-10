# -*- coding: utf-8 -*-
"""
05_summarize_rwr_overlap.py

研究内容一：国家与全球铜矿—物种栖息地冲突区的 RWR 构成。

本版不再用“相交网格数量”近似面积。所有空间相交和面积计算均在
Mollweide 等面积投影（ESRI:54009）中完成，并计算每个国家冲突区实际落入
五个 RWR 百分位等级的面积（km²）及面积占比。

国家口径保持 04_calculate_mine_species_overlap.py 的定义：矿区按矿山属性中的所属国
归属，直接使用其 country_overlap_union 图层，不按缓冲区跨境部分重新分配。

Global 行采用所有国家冲突区几何的严格空间 union 后重新计算，不对国家面积
直接求和，因此不会重复计算跨国或重叠部分。

主要输出：
1. Country_Key_RWR_Area_Proportion.csv
2. Country_Key_RWR_Area_Proportion.gpkg
3. Country_Key_RWR_Intersected_Area_Long.csv
4. Global_Key_RWR_Intersected_Area_Long.csv
5. Country_Key_RWR_BinThresholds.csv
6. Country_Key_RWR_BinArea_Long.csv
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from pyproj import CRS
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.ops import transform, unary_union
from shapely.validation import make_valid


# ============================================================
# 路径与参数
# ============================================================

BASE_DIR = Path(_nc_project(""))

COUNTRY_GPKG = (
    BASE_DIR
    / "04_mine_species_overlap"
    / "Country_CopperOverlap_Summary.gpkg"
)
COUNTRY_LAYER = "country_overlap_union"

RWR_GPKG = (
    BASE_DIR
    / "02_rwr"
    / "Global_RWR_10km_grid.gpkg"
)
RWR_LAYER = "Global_RWR_10km_grid"

OUT_DIR = BASE_DIR / "05_rwr_overlap"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_SUMMARY_CSV = OUT_DIR / "Country_Key_RWR_Area_Proportion.csv"
OUT_SUMMARY_GPKG = OUT_DIR / "Country_Key_RWR_Area_Proportion.gpkg"
OUT_COUNTRY_INTERSECTION_CSV = OUT_DIR / "Country_Key_RWR_Intersected_Area_Long.csv"
OUT_GLOBAL_INTERSECTION_CSV = OUT_DIR / "Global_Key_RWR_Intersected_Area_Long.csv"
OUT_THRESHOLDS_CSV = OUT_DIR / "Country_Key_RWR_BinThresholds.csv"
OUT_BIN_AREA_LONG_CSV = OUT_DIR / "Country_Key_RWR_BinArea_Long.csv"

COUNTRY_FIELD = "Country"
CONTINENT_FIELD = "Continent"
GRID_ID_FIELD = "grid_id"
RWR_FIELD = "weighted_richness"

EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_GPKG_CRS = "EPSG:4326"
PERCENTILE_BASE = "all"  # 可改为 "nonzero"

# 小于该阈值的相交碎片视为浮点误差，不纳入统计。
MIN_INTERSECTION_AREA_KM2 = 1e-8
# 当 RWR 网格覆盖面积与冲突区几何面积偏差超过以下任一阈值时发出警告。
QA_ABS_TOL_KM2 = 1.0
QA_REL_TOL_PCT = 0.5

warnings.filterwarnings("ignore")


# ============================================================
# RWR 分级字段
# ============================================================

BIN_ORDER = ["le50", "p50_75", "p75_90", "p90_99", "p99_100"]
BIN_DISPLAY = {
    "le50": "≤50%",
    "p50_75": "50–75%",
    "p75_90": "75–90%",
    "p90_99": "90–99%",
    "p99_100": "99–100%",
}
AREA_COL_MAP = {
    "le50": "area_le50_km2",
    "p50_75": "area_50_75_km2",
    "p75_90": "area_75_90_km2",
    "p90_99": "area_90_99_km2",
    "p99_100": "area_99_100_km2",
}
PCT_COL_MAP = {
    "le50": "pct_le50",
    "p50_75": "pct_50_75",
    "p75_90": "pct_75_90",
    "p90_99": "pct_90_99",
    "p99_100": "pct_99_100",
}


# ============================================================
# 基础工具
# ============================================================

def read_vector_layer(path: Path, layer: str) -> gpd.GeoDataFrame:
    """读取 GPKG 图层，兼容 ArcGIS 显示的 main.layer 写法。"""
    if not path.exists():
        raise FileNotFoundError(f"文件不存在：{path}")
    layer = str(layer).split(".", 1)[-1] if str(layer).lower().startswith("main.") else str(layer)
    try:
        return gpd.read_file(path, layer=layer, engine="pyogrio")
    except Exception:
        return gpd.read_file(path, layer=layer)


def detect_field(columns, candidates, label):
    for field in candidates:
        if field in columns:
            return field
    raise ValueError(f"未找到{label}字段。候选：{candidates}；实际：{list(columns)}")


def drop_zm(geom):
    if geom is None or geom.is_empty:
        return None
    try:
        return transform(lambda x, y, *args: (x, y), geom)
    except Exception:
        return geom


def extract_polygonal(geom):
    """仅保留 Polygon/MultiPolygon。"""
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "Polygon":
        return geom if geom.area > 0 else None
    if geom.geom_type == "MultiPolygon":
        parts = [g for g in geom.geoms if not g.is_empty and g.area > 0]
    elif geom.geom_type == "GeometryCollection":
        parts = []
        for item in geom.geoms:
            part = extract_polygonal(item)
            if part is None:
                continue
            parts.extend(list(part.geoms) if part.geom_type == "MultiPolygon" else [part])
    else:
        return None
    if not parts:
        return None
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def repair_polygonal(geom):
    if geom is None or geom.is_empty:
        return None
    geom = drop_zm(geom)
    try:
        geom = make_valid(geom)
    except Exception:
        try:
            geom = geom.buffer(0)
        except Exception:
            return None
    return extract_polygonal(geom)


def clean_polygon_gdf(gdf: gpd.GeoDataFrame, name: str, repair_all=False) -> gpd.GeoDataFrame:
    """删除空几何；国家面全部修复，规则网格仅修复无效记录。"""
    before = len(gdf)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if repair_all:
        gdf["geometry"] = gdf.geometry.apply(repair_polygonal)
    else:
        invalid = ~gdf.geometry.is_valid
        if invalid.any():
            print(f"[GEOM] {name}: 修复无效几何 {int(invalid.sum()):,} 条")
            gdf.loc[invalid, "geometry"] = gdf.loc[invalid].geometry.apply(repair_polygonal)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if len(gdf) < before:
        print(f"[GEOM] {name}: 删除无效/空几何 {before - len(gdf):,} 条")
    return gdf


def safe_union(geometries):
    valid = [g for g in geometries if g is not None and not g.is_empty]
    if not valid:
        return None
    try:
        result = unary_union(valid)
    except Exception:
        result = unary_union([repair_polygonal(g) for g in valid if repair_polygonal(g) is not None])
    return repair_polygonal(result)


def project_equal_area(gdf: gpd.GeoDataFrame, name: str) -> gpd.GeoDataFrame:
    if gdf.crs is None:
        raise ValueError(f"{name} 缺少 CRS，不能进行实际面积计算。")
    if not CRS.from_user_input(gdf.crs).equals(CRS.from_user_input(EQUAL_AREA_CRS)):
        print(f"[CRS] {name}: {gdf.crs} -> {EQUAL_AREA_CRS}")
        gdf = gdf.to_crs(EQUAL_AREA_CRS)
    return gdf


# ============================================================
# 数据准备与百分位阈值
# ============================================================

def prepare_country_overlap():
    print("\n[1/8] 读取并准备 country_overlap_union")
    gdf = read_vector_layer(COUNTRY_GPKG, COUNTRY_LAYER)
    country_field = detect_field(
        gdf.columns, [COUNTRY_FIELD, "COUNTRY", "country", "NAME_0", "NAME"], "国家名称"
    )
    if country_field != COUNTRY_FIELD:
        gdf = gdf.rename(columns={country_field: COUNTRY_FIELD})
    if CONTINENT_FIELD not in gdf.columns:
        gdf[CONTINENT_FIELD] = "Unknown"
    gdf[COUNTRY_FIELD] = gdf[COUNTRY_FIELD].astype(str).str.strip()
    gdf[CONTINENT_FIELD] = gdf[CONTINENT_FIELD].fillna("Unknown").astype(str)
    gdf = clean_polygon_gdf(gdf, "country_overlap_union", repair_all=True)

    # 保存原脚本报告的面积用于 QA；本脚本的 overlap_area_km2 由几何重新计算。
    if "overlap_area_km2" in gdf.columns:
        gdf = gdf.rename(columns={"overlap_area_km2": "reported_overlap_area_km2"})

    if gdf[COUNTRY_FIELD].duplicated().any():
        print("  同一国家存在多条记录：按 Country dissolve 后计算实际面积。")
        attrs = gdf.drop(columns="geometry").drop_duplicates(COUNTRY_FIELD)
        geoms = gdf[[COUNTRY_FIELD, "geometry"]].dissolve(by=COUNTRY_FIELD, as_index=False)
        gdf = geoms.merge(attrs, on=COUNTRY_FIELD, how="left")
        gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs=geoms.crs)

    gdf = project_equal_area(gdf, "country_overlap_union")
    gdf = clean_polygon_gdf(gdf, "country_overlap_union projected", repair_all=True)
    gdf["overlap_area_km2"] = gdf.geometry.area / 1_000_000.0
    if (gdf["overlap_area_km2"] <= 0).any():
        raise RuntimeError("存在实际冲突面积不大于 0 的国家几何。")
    print(f"  国家数: {len(gdf):,}")
    return gdf


def prepare_rwr_grid():
    print("\n[2/8] 读取并准备 Global_RWR_10km_grid")
    grid = read_vector_layer(RWR_GPKG, RWR_LAYER)
    for field in [GRID_ID_FIELD, RWR_FIELD]:
        if field not in grid.columns:
            raise ValueError(f"RWR 网格缺少字段 {field}；实际字段：{list(grid.columns)}")
    grid = grid[[GRID_ID_FIELD, RWR_FIELD, "geometry"]].copy()
    grid[GRID_ID_FIELD] = grid[GRID_ID_FIELD].astype(str)
    if grid[GRID_ID_FIELD].duplicated().any():
        raise ValueError(f"{GRID_ID_FIELD} 必须唯一；发现重复 ID。")
    grid[RWR_FIELD] = pd.to_numeric(grid[RWR_FIELD], errors="coerce").fillna(0.0)
    grid = clean_polygon_gdf(grid, "Global_RWR_10km_grid", repair_all=False)
    grid = project_equal_area(grid, "Global_RWR_10km_grid")
    print(f"  RWR 网格数: {len(grid):,}")
    return grid.reset_index(drop=True)


def compute_thresholds(grid: gpd.GeoDataFrame):
    print("\n[3/8] 计算 RWR 百分位阈值")
    values = grid[RWR_FIELD].dropna()
    if PERCENTILE_BASE == "nonzero":
        values_for_q = values[values > 0]
    elif PERCENTILE_BASE == "all":
        values_for_q = values
    else:
        raise ValueError("PERCENTILE_BASE 只能为 'all' 或 'nonzero'。")
    if values_for_q.empty:
        raise RuntimeError("用于计算 RWR 百分位的数据为空。")
    q = values_for_q.quantile([0.50, 0.75, 0.90, 0.99, 1.00])
    thresholds = {
        "p50": float(q.loc[0.50]), "p75": float(q.loc[0.75]),
        "p90": float(q.loc[0.90]), "p99": float(q.loc[0.99]),
        "p100": float(q.loc[1.00]), "min": float(values.min()),
        "percentile_base": PERCENTILE_BASE,
        "n_grid_used": int(len(values_for_q)),
    }
    rows = []
    bounds = [
        ("le50", "≤50%", thresholds["min"], thresholds["p50"]),
        ("p50_75", "50–75%", thresholds["p50"], thresholds["p75"]),
        ("p75_90", "75–90%", thresholds["p75"], thresholds["p90"]),
        ("p90_99", "90–99%", thresholds["p90"], thresholds["p99"]),
        ("p99_100", "99–100%", thresholds["p99"], thresholds["p100"]),
    ]
    for order, (code, label, low, high) in enumerate(bounds, 1):
        rows.append({
            "bin_order": order, "bin_code": code, "bin_label": label,
            "lower_value": low, "upper_value": high,
            "lower_inclusive": code == "le50", "upper_inclusive": True,
            "percentile_base": PERCENTILE_BASE,
            "n_grid_used_for_percentile": len(values_for_q),
        })
    pd.DataFrame(rows).to_csv(OUT_THRESHOLDS_CSV, index=False, encoding="utf-8-sig")
    print(pd.DataFrame(rows)[["bin_label", "lower_value", "upper_value"]].to_string(index=False))
    return thresholds


def classify_rwr(value, t):
    value = float(value)
    if value <= t["p50"]:
        return "le50"
    if value <= t["p75"]:
        return "p50_75"
    if value <= t["p90"]:
        return "p75_90"
    if value <= t["p99"]:
        return "p90_99"
    return "p99_100"


# ============================================================
# 实际相交面积计算
# ============================================================

def intersect_units_with_grid(
    units: gpd.GeoDataFrame,
    grid: gpd.GeoDataFrame,
    thresholds: dict,
    progress_label: str,
) -> pd.DataFrame:
    """逐统计单元使用空间索引筛选候选网格，再计算精确相交面积。"""

    if units.crs is None or grid.crs is None:
        raise ValueError(
            f"统计单元或 RWR 网格缺少 CRS："
            f"统计单元={units.crs}；RWR 网格={grid.crs}"
        )

    units_crs = CRS.from_user_input(units.crs)
    grid_crs = CRS.from_user_input(grid.crs)

    if not units_crs.equals(grid_crs):
        raise ValueError(
            f"统计单元与 RWR 网格 CRS 确实不一致："
            f"统计单元={units.crs}；RWR 网格={grid.crs}"
        )

    sindex = grid.sindex
    records = []
    n_units = len(units)

    for i, row in enumerate(units.itertuples(index=False), 1):
        country = getattr(row, COUNTRY_FIELD)
        continent = getattr(row, CONTINENT_FIELD)
        geom = row.geometry
        if i == 1 or i % 20 == 0 or i == n_units:
            print(f"  {progress_label}: {i:,}/{n_units:,} - {country}")

        candidate_pos = np.asarray(sindex.query(geom, predicate="intersects"), dtype=int)
        if candidate_pos.size == 0:
            warnings.warn(f"{country} 未与任何 RWR 网格相交。")
            continue
        candidate = grid.iloc[candidate_pos][[GRID_ID_FIELD, RWR_FIELD, "rwr_bin", "geometry"]].copy()
        try:
            intersection = candidate.geometry.intersection(geom)
        except Exception:
            intersection = candidate.geometry.apply(
                lambda grid_geom: repair_polygonal(grid_geom).intersection(geom)
            )
        area = intersection.area.to_numpy(dtype=float) / 1_000_000.0
        keep = np.isfinite(area) & (area > MIN_INTERSECTION_AREA_KM2)
        if not keep.any():
            warnings.warn(f"{country} 的候选网格精确相交面积均为 0。")
            continue

        selected = candidate.loc[keep, [GRID_ID_FIELD, RWR_FIELD, "rwr_bin"]].copy()
        selected[COUNTRY_FIELD] = country
        selected[CONTINENT_FIELD] = continent
        selected["rwr_bin_label"] = selected["rwr_bin"].map(BIN_DISPLAY)
        selected["intersect_area_km2"] = area[keep]
        records.append(selected[[
            COUNTRY_FIELD, CONTINENT_FIELD, GRID_ID_FIELD, RWR_FIELD,
            "rwr_bin", "rwr_bin_label", "intersect_area_km2",
        ]])

    if not records:
        raise RuntimeError(f"{progress_label} 没有产生有效的精确相交面积。")
    result = pd.concat(records, ignore_index=True)
    return result.sort_values([COUNTRY_FIELD, "rwr_bin", GRID_ID_FIELD]).reset_index(drop=True)


def summarize_area(long_df: pd.DataFrame, unit_meta: pd.DataFrame) -> pd.DataFrame:
    """将逐网格精确相交面积汇总为国家/Global 的宽表。"""
    pivot = long_df.pivot_table(
        index=[COUNTRY_FIELD, CONTINENT_FIELD], columns="rwr_bin",
        values="intersect_area_km2", aggfunc="sum", fill_value=0.0,
    )
    for code in BIN_ORDER:
        if code not in pivot.columns:
            pivot[code] = 0.0
    pivot = pivot[BIN_ORDER].rename(columns=AREA_COL_MAP).reset_index()

    area_cols = [AREA_COL_MAP[c] for c in BIN_ORDER]
    pivot["area_total_rwr_km2"] = pivot[area_cols].sum(axis=1)
    for code in BIN_ORDER:
        pivot[PCT_COL_MAP[code]] = np.where(
            pivot["area_total_rwr_km2"] > 0,
            pivot[AREA_COL_MAP[code]] / pivot["area_total_rwr_km2"] * 100.0,
            0.0,
        )
    pivot["area_ge90_km2"] = pivot["area_90_99_km2"] + pivot["area_99_100_km2"]
    pivot["pct_ge90"] = pivot["pct_90_99"] + pivot["pct_99_100"]
    pivot["area_ge99_km2"] = pivot["area_99_100_km2"]
    pivot["pct_ge99"] = pivot["pct_99_100"]

    counts = long_df.groupby(COUNTRY_FIELD)[GRID_ID_FIELD].nunique().rename("n_intersected_grids")
    weighted = long_df.assign(
        _weighted=long_df[RWR_FIELD] * long_df["intersect_area_km2"]
    ).groupby(COUNTRY_FIELD).agg(
        rwr_min=(RWR_FIELD, "min"), rwr_max=(RWR_FIELD, "max"),
        weighted_sum=("_weighted", "sum"),
    )
    pivot = pivot.merge(counts, on=COUNTRY_FIELD, how="left").merge(
        weighted.reset_index(), on=COUNTRY_FIELD, how="left"
    )
    pivot["rwr_area_weighted_mean"] = (
        pivot["weighted_sum"] / pivot["area_total_rwr_km2"].replace(0, np.nan)
    )
    pivot = pivot.drop(columns="weighted_sum")

    meta = unit_meta.copy()
    if CONTINENT_FIELD in meta.columns:
        meta = meta.drop(columns=CONTINENT_FIELD)
    pivot = pivot.merge(meta, on=COUNTRY_FIELD, how="left")
    pivot["area_balance_difference_km2"] = (
        pivot["area_total_rwr_km2"] - pivot["overlap_area_km2"]
    )
    pivot["area_balance_difference_pct"] = np.where(
        pivot["overlap_area_km2"] > 0,
        pivot["area_balance_difference_km2"] / pivot["overlap_area_km2"] * 100.0,
        np.nan,
    )
    pivot["rwr_grid_coverage_pct"] = np.where(
        pivot["overlap_area_km2"] > 0,
        pivot["area_total_rwr_km2"] / pivot["overlap_area_km2"] * 100.0,
        np.nan,
    )
    pivot["uncovered_area_km2"] = (
        pivot["overlap_area_km2"] - pivot["area_total_rwr_km2"]
    ).clip(lower=0.0)
    return pivot


def build_bin_area_long(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in summary.iterrows():
        for order, code in enumerate(BIN_ORDER, 1):
            rows.append({
                COUNTRY_FIELD: row[COUNTRY_FIELD],
                CONTINENT_FIELD: row[CONTINENT_FIELD],
                "rwr_bin_order": order,
                "rwr_bin": code,
                "rwr_bin_label": BIN_DISPLAY[code],
                "area_km2": row[AREA_COL_MAP[code]],
                "area_pct": row[PCT_COL_MAP[code]],
                "area_total_rwr_km2": row["area_total_rwr_km2"],
                "overlap_area_km2": row["overlap_area_km2"],
            })
    return pd.DataFrame(rows)


def report_area_qa(summary: pd.DataFrame):
    problematic = summary[
        (summary["area_balance_difference_km2"].abs() > QA_ABS_TOL_KM2)
        & (summary["area_balance_difference_pct"].abs() > QA_REL_TOL_PCT)
    ]
    if problematic.empty:
        print("[QA] 所有统计单元的 RWR 分级面积和与冲突区几何面积在容差内一致。")
    else:
        warnings.warn(
            f"{len(problematic)} 个统计单元的面积平衡偏差超过 "
            f"{QA_ABS_TOL_KM2} km² 且 {QA_REL_TOL_PCT}%。"
            "这通常表示冲突区部分超出 RWR 陆地网格覆盖范围；详见输出的 QA 字段。"
        )
        print(problematic[[
            COUNTRY_FIELD, "overlap_area_km2", "area_total_rwr_km2",
            "area_balance_difference_km2", "area_balance_difference_pct",
        ]].head(20).to_string(index=False))


def write_summary_gpkg(
    country: gpd.GeoDataFrame,
    global_geom,
    summary: pd.DataFrame,
):
    country_geom = country[[COUNTRY_FIELD, "geometry"]].merge(
        summary[summary[COUNTRY_FIELD] != "Global"], on=COUNTRY_FIELD, how="inner"
    )
    global_gdf = gpd.GeoDataFrame(
        [{COUNTRY_FIELD: "Global", "geometry": global_geom}],
        geometry="geometry", crs=country.crs,
    ).merge(summary[summary[COUNTRY_FIELD] == "Global"], on=COUNTRY_FIELD, how="inner")
    out = pd.concat([global_gdf, country_geom], ignore_index=True, sort=False)
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=country.crs).to_crs(OUTPUT_GPKG_CRS)
    for col in out.select_dtypes(include="object").columns:
        out[col] = out[col].fillna("").astype(str)
    if OUT_SUMMARY_GPKG.exists():
        OUT_SUMMARY_GPKG.unlink()
    out.to_file(
        OUT_SUMMARY_GPKG,
        layer="country_key_rwr_area_proportion",
        driver="GPKG",
        engine="pyogrio",
    )


# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 90)
    print("国家与全球铜矿冲突区的 RWR 等级实际面积及面积占比")
    print("=" * 90)

    country = prepare_country_overlap()
    grid = prepare_rwr_grid()
    thresholds = compute_thresholds(grid)
    grid["rwr_bin"] = grid[RWR_FIELD].apply(lambda value: classify_rwr(value, thresholds))

    print("\n[4/8] 计算各国家冲突区 × RWR 网格的精确相交面积")
    country_units = country[[COUNTRY_FIELD, CONTINENT_FIELD, "geometry"]].copy()
    country_long = intersect_units_with_grid(country_units, grid, thresholds, "国家")
    country_long.to_csv(OUT_COUNTRY_INTERSECTION_CSV, index=False, encoding="utf-8-sig")
    print(f"  国家逐网格实际面积表: {OUT_COUNTRY_INTERSECTION_CSV}")

    print("\n[5/8] 严格 union 后计算 Global 精确相交面积")
    global_geom = safe_union(country.geometry)
    if global_geom is None:
        raise RuntimeError("所有国家冲突区的 Global union 为空。")
    global_units = gpd.GeoDataFrame(
        [{COUNTRY_FIELD: "Global", CONTINENT_FIELD: "Global", "geometry": global_geom}],
        geometry="geometry", crs=country.crs,
    )
    global_long = intersect_units_with_grid(global_units, grid, thresholds, "Global")
    global_long.to_csv(OUT_GLOBAL_INTERSECTION_CSV, index=False, encoding="utf-8-sig")
    print(f"  Global 逐网格实际面积表: {OUT_GLOBAL_INTERSECTION_CSV}")

    print("\n[6/8] 汇总国家与 Global 五级 RWR 的实际面积及占比")
    country_meta = country.drop(columns="geometry").copy()
    country_summary = summarize_area(country_long, country_meta)
    global_meta = pd.DataFrame([{
        COUNTRY_FIELD: "Global",
        "overlap_area_km2": global_geom.area / 1_000_000.0,
    }])
    global_summary = summarize_area(global_long, global_meta)
    summary = pd.concat([global_summary, country_summary], ignore_index=True, sort=False)

    preferred = [
        COUNTRY_FIELD, CONTINENT_FIELD, "overlap_area_km2", "area_total_rwr_km2",
        "area_le50_km2", "area_50_75_km2", "area_75_90_km2",
        "area_90_99_km2", "area_99_100_km2", "area_ge90_km2", "area_ge99_km2",
        "pct_le50", "pct_50_75", "pct_75_90", "pct_90_99", "pct_99_100",
        "pct_ge90", "pct_ge99", "n_intersected_grids", "rwr_min", "rwr_max",
        "rwr_area_weighted_mean", "rwr_grid_coverage_pct",
        "area_balance_difference_km2", "area_balance_difference_pct", "uncovered_area_km2",
    ]
    other = [c for c in summary.columns if c not in preferred]
    summary = summary[preferred + other]
    summary["_global_order"] = np.where(summary[COUNTRY_FIELD] == "Global", 0, 1)
    summary = summary.sort_values(
        ["_global_order", "overlap_area_km2", "pct_ge90"],
        ascending=[True, False, False],
    ).drop(columns="_global_order").reset_index(drop=True)
    report_area_qa(summary)
    summary.to_csv(OUT_SUMMARY_CSV, index=False, encoding="utf-8-sig")
    build_bin_area_long(summary).to_csv(OUT_BIN_AREA_LONG_CSV, index=False, encoding="utf-8-sig")

    print("\n[7/8] 输出带几何的国家与 Global 统计结果")
    write_summary_gpkg(country, global_geom, summary)

    print("\n[8/8] 完成")
    print(f"1. {OUT_SUMMARY_CSV}")
    print(f"2. {OUT_SUMMARY_GPKG}")
    print(f"3. {OUT_COUNTRY_INTERSECTION_CSV}")
    print(f"4. {OUT_GLOBAL_INTERSECTION_CSV}")
    print(f"5. {OUT_THRESHOLDS_CSV}")
    print(f"6. {OUT_BIN_AREA_LONG_CSV}")
    print("\n说明：所有 area_*_km2 均为等面积投影下的精确几何相交面积；")
    print("pct_ge90 = pct_90_99 + pct_99_100，即冲突面积落入全球 RWR 前 10% 的比例。")


if __name__ == "__main__":
    main()
