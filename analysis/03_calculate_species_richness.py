# -*- coding: utf-8 -*-
"""
03_calculate_species_richness.py

统计全球 10 km x 10 km 陆地网格中受威胁(CR/EN/VU)脊椎动物物种丰富度。

对每个网格输出:
  - total_richness        : 涉及的所有物种数 (去重)
  - n_Amphibians / n_Birds / n_Fish / n_Mammals / n_Reptiles : 各类群物种数
  - n_CR / n_EN / n_VU    : 各濒危等级物种数

输入:
  网格 : <work-root>/01_grid\\Global_10km_grid.gpkg (layer=grid)
  IUCN : <data-root>/IUCNredlist\\202502\\{Taxon}-{Status}\\data_0.shp

输出:
  <work-root>/03_species_richness\\Global_10km_grid.gpkg
  <work-root>/03_species_richness\\Global_10km_grid.csv

说明:
  - 同一物种与同一网格无论相交多少个 polygon, 只计 1 次 (按 ID_NO 去重)
  - PRESENCE / ORIGIN / SEASONAL 不做筛选, 全部纳入
  - 仅判断相交关系 (predicate=intersects), 不计算相交面积
  - 中间结果按 (taxon, status) 写出 parquet, 支持断点续算
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings
import gc
import time

import pandas as pd
import geopandas as gpd


# ============================================================
# 路径与参数
# ============================================================
GRID_PATH = Path(_nc_project("01_grid/Global_10km_grid.gpkg"))
GRID_LAYER = "grid"

IUCN_ROOT = Path(_nc_data("IUCNredlist/202502"))

OUT_DIR = Path(_nc_project("03_species_richness"))
OUT_GPKG = OUT_DIR / "Global_10km_grid.gpkg"
OUT_CSV = OUT_DIR / "Global_10km_grid.csv"
INTERIM_DIR = OUT_DIR / "interim"

TAXA = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
STATUSES = ["CR", "EN", "VU"]

# 物种识别字段优先级: 优先 ID_NO, 不存在再退回 SCI_NAME
ID_FIELD_CANDIDATES = ["ID_NO", "SCI_NAME"]


# ============================================================
# 工具函数
# ============================================================
def detect_id_field(columns):
    """从一个 IUCN shapefile 的字段列表中选择物种识别字段."""
    for f in ID_FIELD_CANDIDATES:
        if f in columns:
            return f
    raise ValueError(
        f"未找到物种识别字段, 期望 {ID_FIELD_CANDIDATES}, "
        f"实际字段: {list(columns)}"
    )


def process_one_group(taxon, status, grid):
    """
    处理一个 taxon-status 组合, 返回 (grid_id, species_id, taxon, status) 长表。
    若该组合数据缺失或为空, 返回 None。
    若已存在中间 parquet, 直接读取返回, 实现断点续算。
    """
    sub = IUCN_ROOT / f"{taxon}-{status}"
    shp = sub / "data_0.shp"
    interim_path = INTERIM_DIR / f"{taxon}_{status}.parquet"

    if interim_path.exists():
        print(f"[CACHE] 已有中间结果, 跳过空间叠加: {interim_path.name}")
        return pd.read_parquet(interim_path)

    if not shp.exists():
        warnings.warn(f"[SKIP] 数据文件不存在: {shp}")
        return None

    t0 = time.time()
    print(f"[READ ] {taxon}-{status}: {shp}")
    sp = gpd.read_file(shp, engine="pyogrio")
    if sp.empty:
        warnings.warn(f"[EMPTY] {shp} 为空")
        return None

    id_field = detect_id_field(sp.columns)

    # 投影到网格 CRS (网格通常为 Mollweide 等面积投影)
    if sp.crs is None:
        sp = sp.set_crs(epsg=4326)
    if sp.crs != grid.crs:
        sp = sp.to_crs(grid.crs)

    # 修复无效几何 (跨日界线/极地投影后可能产生)
    try:
        sp["geometry"] = sp.geometry.make_valid()
    except Exception:
        sp["geometry"] = sp.geometry.buffer(0)
    sp = sp[sp.geometry.notna() & ~sp.geometry.is_empty].copy()

    # 仅保留必要字段
    sp = sp[[id_field, "geometry"]].rename(columns={id_field: "species_id"})
    # species_id 转字符串, 避免 ID_NO 数值/SCI_NAME 字符串混用
    sp["species_id"] = sp["species_id"].astype(str)

    print(f"[SJOIN] {taxon}-{status}: {len(sp)} polygons -> sjoin with {len(grid)} grids")
    joined = gpd.sjoin(
        grid[["grid_id", "geometry"]],
        sp,
        how="inner",
        predicate="intersects",
    )

    df = (
        joined[["grid_id", "species_id"]]
        .drop_duplicates(subset=["grid_id", "species_id"])
        .reset_index(drop=True)
    )
    df["taxon"] = taxon
    df["status"] = status

    df.to_parquet(interim_path, index=False)
    dt = time.time() - t0
    print(
        f"[OK   ] {taxon}-{status}: {len(df):>10,d} (grid_id, species_id) 对, "
        f"unique species = {df['species_id'].nunique():,d}, "
        f"用时 {dt:.1f}s -> {interim_path.name}"
    )

    del sp, joined
    gc.collect()
    return df


# ============================================================
# 主流程
# ============================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[GRID ] 读取网格: {GRID_PATH} | layer={GRID_LAYER}")
    grid = gpd.read_file(GRID_PATH, layer=GRID_LAYER, engine="pyogrio")
    if "grid_id" not in grid.columns:
        raise ValueError("网格图层缺少 grid_id 字段")
    print(f"[GRID ] 网格数量: {len(grid):,d}, CRS: {grid.crs}")

    # ------------------------------------------------------------
    # 1. 对 5 类群 x 3 等级 = 15 组逐一做空间叠加
    # ------------------------------------------------------------
    all_pairs = []
    for taxon in TAXA:
        for status in STATUSES:
            df = process_one_group(taxon, status, grid)
            if df is not None and len(df) > 0:
                all_pairs.append(df)

    if not all_pairs:
        raise RuntimeError("没有成功处理任何 taxon-status 组合, 请检查 IUCN 数据路径")

    pairs = pd.concat(all_pairs, ignore_index=True)

    # 跨组合不应出现重复 (一个物种只属于一个 taxon 和一个 status),
    # 但为稳健起见再去一次重
    pairs = pairs.drop_duplicates(subset=["grid_id", "species_id"]).reset_index(drop=True)
    print(f"[MERGE] 合并后 (grid_id, species_id) 对: {len(pairs):,d}")
    print(f"[MERGE] 涉及 unique species: {pairs['species_id'].nunique():,d}")
    print(f"[MERGE] 涉及 unique grid:    {pairs['grid_id'].nunique():,d}")

    # ------------------------------------------------------------
    # 2. 统计指标
    # ------------------------------------------------------------
    # 物种总数
    total = (
        pairs.groupby("grid_id")["species_id"]
        .nunique()
        .rename("total_richness")
    )

    # 各类群物种数
    per_taxon = (
        pairs.groupby(["grid_id", "taxon"])["species_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    per_taxon.columns = [f"n_{c}" for c in per_taxon.columns]
    for t in TAXA:
        col = f"n_{t}"
        if col not in per_taxon.columns:
            per_taxon[col] = 0
    per_taxon = per_taxon[[f"n_{t}" for t in TAXA]]

    # 各濒危等级物种数
    per_status = (
        pairs.groupby(["grid_id", "status"])["species_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    per_status.columns = [f"n_{c}" for c in per_status.columns]
    for s in STATUSES:
        col = f"n_{s}"
        if col not in per_status.columns:
            per_status[col] = 0
    per_status = per_status[[f"n_{s}" for s in STATUSES]]

    stats = pd.concat([total, per_taxon, per_status], axis=1).reset_index()
    int_cols = [c for c in stats.columns if c != "grid_id"]
    stats[int_cols] = stats[int_cols].fillna(0).astype("int32")

    print(f"[STATS] 含物种的网格数: {len(stats):,d} / {len(grid):,d}")
    print(f"[STATS] 单网格物种数 max: {stats['total_richness'].max()}, "
          f"mean: {stats['total_richness'].mean():.2f}")

    # ------------------------------------------------------------
    # 3. 合并回网格 + 写出
    # ------------------------------------------------------------
    out = grid.merge(stats, on="grid_id", how="left")
    for c in int_cols:
        out[c] = out[c].fillna(0).astype("int32")

    print(f"[WRITE] GeoPackage -> {OUT_GPKG}")
    out.to_file(OUT_GPKG, layer="Global_10km_grid", driver="GPKG")

    print(f"[WRITE] CSV       -> {OUT_CSV}")
    out_csv = pd.DataFrame(out.drop(columns="geometry"))
    out_csv.to_csv(OUT_CSV, index=False)

    print("[DONE ] 全部完成")


if __name__ == "__main__":
    main()
