# -*- coding: utf-8 -*-
"""
02_calculate_rarity_weighted_richness.py

按照 IUCN / IBAT 对 Rarity-Weighted Richness 的核心含义，
对全球 10 km x 10 km 陆地网格计算受威胁(CR/EN/VU)脊椎动物的面积比例型稀有度加权丰富度。

核心公式：
  对物种 s 和网格 g：
      contribution_{s,g} = A_{s,g} / A_s

  其中：
      A_{s,g} = 物种 s 的分布范围与网格 g 的实际相交面积
      A_s     = 物种 s 的全球分布范围总面积

  网格 g 的 rarity-weighted richness：
      RWR_g = sum_s contribution_{s,g}

与旧版 1 / N_s 网格占据数算法的区别：
  - 旧版：只要物种与网格相交，无论相交面积多大，都给该网格 1 / 覆盖网格数。
  - 本版：按真实相交面积分配权重，边界擦碰但面积为 0 的网格不贡献权重。

输入：
  网格 : <work-root>/01_grid\\Global_10km_grid.gpkg (layer=grid)
  IUCN : <data-root>/IUCNredlist\\202502\\{Taxon}-{Status}\\data_0.shp

输出：
  <work-root>/2Global_WeightedRichness_10km_grid\\Global_RWR_10km_grid.gpkg
  <work-root>/2Global_WeightedRichness_10km_grid\\Global_RWR_10km_grid.csv
  <work-root>/2Global_WeightedRichness_10km_grid\\Species_RangeArea_Table.csv
  <work-root>/2Global_WeightedRichness_10km_grid\\interim\\{Taxon}_{Status}_area_rwr.parquet

说明：
  - 物种识别字段优先级：ID_NO > SCI_NAME，与此前丰富度代码一致。
  - PRESENCE / ORIGIN / SEASONAL 不做筛选，与此前丰富度代码保持一致。
  - 同一物种的多个 polygon 会先 dissolve/union，避免同一物种内部重叠造成重复计面积。
  - CRS 统一使用网格 CRS，通常为 Mollweide 等面积投影，因此面积单位为 m²。
  - 本脚本的分母 A_s 使用物种完整 IUCN 分布范围面积，包括可能存在的海洋/水域部分；
    因此对海陆双分布物种，其在所有陆地网格中的贡献总和允许小于 1。
  - grid_coverage_ratio / rwr_contribution_sum 表示该物种完整分布范围中被当前陆地网格覆盖的比例；
    小于 1 通常不是错误，而是说明该物种部分分布范围位于陆地网格之外。
  - 如果需要生态学更严格版本，可在 read_species_data() 中统一加入 PRESENCE/ORIGIN/SEASONAL 筛选，
    但应与 Species Richness 指标同步重算，避免两个指标口径不一致。
"""

from __future__ import annotations

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

import gc
import shutil
import time
import warnings
from pathlib import Path
from typing import Optional, Tuple, List

import geopandas as gpd
import pandas as pd
import numpy as np


# ============================================================
# 路径与参数
# ============================================================
GRID_PATH = Path(_nc_project("01_grid/Global_10km_grid.gpkg"))
GRID_LAYER = "grid"

IUCN_ROOT = Path(_nc_data("IUCNredlist/202502"))

OUT_DIR = Path(_nc_project("02_rwr"))
OUT_GPKG = OUT_DIR / "Global_RWR_10km_grid.gpkg"
OUT_CSV = OUT_DIR / "Global_RWR_10km_grid.csv"
SPECIES_AREA_CSV = OUT_DIR / "Species_RangeArea_Table.csv"
INTERIM_DIR = OUT_DIR / "interim"

TAXA = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
STATUSES = ["CR", "EN", "VU"]

ID_FIELD_CANDIDATES = ["ID_NO", "SCI_NAME"]

# 面积阈值：用于剔除线接触/点接触导致的 0 面积相交。
# 设为 0 表示只剔除严格 0 面积；如果想避免极小碎片，可改为 1.0 或 100.0（单位 m²）。
MIN_OVERLAP_AREA_M2 = 0.0
MIN_SPECIES_AREA_M2 = 0.0

# 每处理多少个物种写一个临时 part parquet，降低内存占用。
BATCH_SPECIES = 50

# 是否删除每组的 part 临时文件夹。
REMOVE_PARTS_AFTER_MERGE = True


# ============================================================
# 工具函数
# ============================================================
def detect_id_field(columns) -> str:
    """从 IUCN shapefile 字段列表中选择物种识别字段。"""
    for f in ID_FIELD_CANDIDATES:
        if f in columns:
            return f
    raise ValueError(
        f"未找到物种识别字段，期望 {ID_FIELD_CANDIDATES}，"
        f"实际字段: {list(columns)}"
    )


def repair_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """修复无效几何并剔除空几何。"""
    try:
        gdf["geometry"] = gdf.geometry.make_valid()
    except Exception:
        gdf["geometry"] = gdf.geometry.buffer(0)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    return gdf


def read_species_data(taxon: str, status: str, target_crs) -> Optional[gpd.GeoDataFrame]:
    """读取一个 taxon-status 组合的 IUCN 数据，统一 CRS 并 dissolve 到物种级。"""
    shp = IUCN_ROOT / f"{taxon}-{status}" / "data_0.shp"
    if not shp.exists():
        warnings.warn(f"[SKIP] 数据文件不存在: {shp}")
        return None

    print(f"[READ ] {taxon}-{status}: {shp}")
    sp = gpd.read_file(shp, engine="pyogrio")
    if sp.empty:
        warnings.warn(f"[EMPTY] {shp} 为空")
        return None

    id_field = detect_id_field(sp.columns)

    # 如需更严格版本，可在此处打开筛选条件；但应同步重算 Species Richness。
    # if "PRESENCE" in sp.columns:
    #     sp = sp[sp["PRESENCE"].isin([1, 2])]
    # if "ORIGIN" in sp.columns:
    #     sp = sp[sp["ORIGIN"].isin([1, 2])]

    sp = sp[[id_field, "geometry"]].rename(columns={id_field: "species_id"})
    sp = sp[sp["species_id"].notna()].copy()
    sp["species_id"] = sp["species_id"].astype(str)

    if sp.crs is None:
        sp = sp.set_crs(epsg=4326)
    if sp.crs != target_crs:
        sp = sp.to_crs(target_crs)

    sp = repair_geometry(sp)
    if sp.empty:
        warnings.warn(f"[EMPTY] {taxon}-{status} 修复几何后为空")
        return None

    # 关键：同一物种多个 polygon 先 union/dissolve，避免季节、多部分、多来源 polygon 内部重叠造成重复面积。
    print(f"[DISS ] {taxon}-{status}: {len(sp):,d} polygons -> dissolve by species_id")
    sp = sp.dissolve(by="species_id", as_index=False)
    sp = repair_geometry(sp)
    sp["species_total_area_m2"] = sp.geometry.area.astype("float64")
    sp = sp[sp["species_total_area_m2"] > MIN_SPECIES_AREA_M2].copy()

    if sp.empty:
        warnings.warn(f"[EMPTY] {taxon}-{status} dissolve 后没有有效面积物种")
        return None

    sp["taxon"] = taxon
    sp["status"] = status
    print(
        f"[DISS ] {taxon}-{status}: unique species = {len(sp):,d}, "
        f"area range km² = [{sp['species_total_area_m2'].min()/1e6:.6f}, "
        f"{sp['species_total_area_m2'].max()/1e6:.2f}]"
    )
    return sp[["species_id", "taxon", "status", "species_total_area_m2", "geometry"]]


def query_candidate_cells(grid: gpd.GeoDataFrame, grid_sindex, geom) -> gpd.GeoDataFrame:
    """利用空间索引提取与物种几何相交的候选网格。"""
    try:
        idx = grid_sindex.query(geom, predicate="intersects")
        if len(idx) == 0:
            return grid.iloc[[]]
        return grid.iloc[idx]
    except TypeError:
        # 兼容较旧 geopandas/shapely：先按 bbox 查，再手动 intersects。
        idx = list(grid_sindex.intersection(geom.bounds))
        if len(idx) == 0:
            return grid.iloc[[]]
        cand = grid.iloc[idx]
        return cand[cand.geometry.intersects(geom)]


def write_part(records: List[pd.DataFrame], part_dir: Path, part_idx: int) -> int:
    """把当前批次记录写成 part parquet。"""
    if not records:
        return part_idx
    part = pd.concat(records, ignore_index=True)
    part_path = part_dir / f"part_{part_idx:04d}.parquet"
    part.to_parquet(part_path, index=False)
    print(f"[PART ] 写出 {part_path.name}: {len(part):,d} rows")
    del part
    records.clear()
    gc.collect()
    return part_idx + 1


def compute_species_grid_contributions(
    species_gdf: gpd.GeoDataFrame,
    grid: gpd.GeoDataFrame,
    taxon: str,
    status: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    对一个 taxon-status 组合计算 species-grid 面积贡献。

    返回：
      contrib: grid_id, species_id, taxon, status, overlap_area_km2,
               species_total_area_km2, rwr_contribution
      species_summary: species_id, taxon, status, species_total_area_km2,
                       n_cells_positive_area, overlap_area_sum_km2, grid_coverage_ratio
    """
    final_path = INTERIM_DIR / f"{taxon}_{status}_area_rwr.parquet"
    species_summary_path = INTERIM_DIR / f"{taxon}_{status}_species_area.csv"

    if final_path.exists() and species_summary_path.exists():
        print(f"[CACHE] 复用面积型中间结果: {final_path.name}")
        return pd.read_parquet(final_path), pd.read_csv(species_summary_path)

    t0 = time.time()
    part_dir = INTERIM_DIR / f"{taxon}_{status}_parts"
    if part_dir.exists():
        shutil.rmtree(part_dir)
    part_dir.mkdir(parents=True, exist_ok=True)

    grid_sindex = grid.sindex
    records: List[pd.DataFrame] = []
    part_idx = 0
    species_rows = []

    print(f"[AREA ] {taxon}-{status}: start species-cell intersection")
    for k, row in enumerate(species_gdf.itertuples(index=False), start=1):
        species_id = row.species_id
        geom = row.geometry
        total_area_m2 = float(row.species_total_area_m2)

        if geom is None or geom.is_empty or total_area_m2 <= MIN_SPECIES_AREA_M2:
            continue

        cand = query_candidate_cells(grid, grid_sindex, geom)
        if cand.empty:
            species_rows.append({
                "species_id": species_id,
                "taxon": taxon,
                "status": status,
                "species_total_area_km2": total_area_m2 / 1e6,
                "n_cells_positive_area": 0,
                "overlap_area_sum_km2": 0.0,
                "grid_coverage_ratio": 0.0,
            })
            continue

        # 精确相交面积。边界线/点接触的面积为 0，会在下面过滤掉。
        inter = cand.geometry.intersection(geom)
        area_m2 = inter.area.to_numpy(dtype="float64")
        mask = area_m2 > MIN_OVERLAP_AREA_M2

        if not np.any(mask):
            species_rows.append({
                "species_id": species_id,
                "taxon": taxon,
                "status": status,
                "species_total_area_km2": total_area_m2 / 1e6,
                "n_cells_positive_area": 0,
                "overlap_area_sum_km2": 0.0,
                "grid_coverage_ratio": 0.0,
            })
            continue

        overlap_area_m2 = area_m2[mask]
        grid_ids = cand["grid_id"].to_numpy()[mask]
        contrib = overlap_area_m2 / total_area_m2

        df = pd.DataFrame({
            "grid_id": grid_ids,
            "species_id": species_id,
            "taxon": taxon,
            "status": status,
            "overlap_area_km2": overlap_area_m2 / 1e6,
            "species_total_area_km2": total_area_m2 / 1e6,
            "rwr_contribution": contrib,
        })
        records.append(df)

        overlap_sum_km2 = float(overlap_area_m2.sum() / 1e6)
        species_rows.append({
            "species_id": species_id,
            "taxon": taxon,
            "status": status,
            "species_total_area_km2": total_area_m2 / 1e6,
            "n_cells_positive_area": int(mask.sum()),
            "overlap_area_sum_km2": overlap_sum_km2,
            "grid_coverage_ratio": overlap_sum_km2 / (total_area_m2 / 1e6),
        })

        if k % BATCH_SPECIES == 0:
            part_idx = write_part(records, part_dir, part_idx)
            print(f"[AREA ] {taxon}-{status}: processed {k:,d}/{len(species_gdf):,d} species")

    part_idx = write_part(records, part_dir, part_idx)

    part_files = sorted(part_dir.glob("part_*.parquet"))
    if part_files:
        contrib = pd.concat([pd.read_parquet(p) for p in part_files], ignore_index=True)
    else:
        contrib = pd.DataFrame(columns=[
            "grid_id", "species_id", "taxon", "status", "overlap_area_km2",
            "species_total_area_km2", "rwr_contribution"
        ])

    # 理论上每个 species-grid 只出现一次；这里再稳健合并，避免异常重复。
    if not contrib.empty:
        dup_n = contrib.duplicated(subset=["grid_id", "species_id"]).sum()
        if dup_n > 0:
            warnings.warn(f"[{taxon}-{status}] 发现重复 species-grid 记录 {dup_n:,d} 条，将按面积和贡献合并")
            contrib = (
                contrib.groupby(["grid_id", "species_id", "taxon", "status"], as_index=False)
                .agg(
                    overlap_area_km2=("overlap_area_km2", "sum"),
                    species_total_area_km2=("species_total_area_km2", "first"),
                    rwr_contribution=("rwr_contribution", "sum"),
                )
            )

    species_summary = pd.DataFrame(species_rows)
    contrib.to_parquet(final_path, index=False)
    species_summary.to_csv(species_summary_path, index=False)

    if REMOVE_PARTS_AFTER_MERGE and part_dir.exists():
        shutil.rmtree(part_dir)

    dt = time.time() - t0
    print(
        f"[OK   ] {taxon}-{status}: contrib rows = {len(contrib):,d}, "
        f"species = {species_summary['species_id'].nunique() if not species_summary.empty else 0:,d}, "
        f"time = {dt:.1f}s -> {final_path.name}"
    )
    return contrib, species_summary


def aggregate_rwr(contrib: pd.DataFrame) -> pd.DataFrame:
    """按 grid_id 汇总 RWR 总量、类群分量和濒危等级分量。"""
    if contrib.empty:
        raise RuntimeError("contrib 为空，无法计算 RWR")

    # 若同一物种异常出现在多个 taxon/status 组合，避免总 RWR 重复计数。
    dup_n = contrib.duplicated(subset=["grid_id", "species_id"]).sum()
    if dup_n > 0:
        warnings.warn(
            f"合并所有组合后发现重复 grid_id-species_id 记录 {dup_n:,d} 条。"
            f"将保留第一条，避免同一物种重复贡献。"
        )
        contrib = contrib.drop_duplicates(subset=["grid_id", "species_id"], keep="first").copy()

    rwr_total = contrib.groupby("grid_id")["rwr_contribution"].sum().rename("rwr_total")
    # 兼容旧脚本命名，下游如果继续用 weighted_richness 也可直接读取。
    weighted_richness = rwr_total.rename("weighted_richness")

    area_based_richness = (
        contrib.groupby("grid_id")["species_id"].nunique().rename("area_based_richness")
    )

    per_taxon = (
        contrib.groupby(["grid_id", "taxon"])["rwr_contribution"]
        .sum()
        .unstack(fill_value=0.0)
    )
    per_taxon.columns = [f"rwr_{c}" for c in per_taxon.columns]
    for t in TAXA:
        col = f"rwr_{t}"
        if col not in per_taxon.columns:
            per_taxon[col] = 0.0
    per_taxon = per_taxon[[f"rwr_{t}" for t in TAXA]]

    per_status = (
        contrib.groupby(["grid_id", "status"])["rwr_contribution"]
        .sum()
        .unstack(fill_value=0.0)
    )
    per_status.columns = [f"rwr_{c}" for c in per_status.columns]
    for s in STATUSES:
        col = f"rwr_{s}"
        if col not in per_status.columns:
            per_status[col] = 0.0
    per_status = per_status[[f"rwr_{s}" for s in STATUSES]]

    stats = pd.concat(
        [area_based_richness, rwr_total, weighted_richness, per_taxon, per_status],
        axis=1,
    ).reset_index()

    stats["area_based_richness"] = stats["area_based_richness"].fillna(0).astype("int32")
    float_cols = [c for c in stats.columns if c not in ["grid_id", "area_based_richness"]]
    stats[float_cols] = stats[float_cols].fillna(0.0).astype("float64")
    return stats


def build_species_area_table(contrib: pd.DataFrame, species_info: pd.DataFrame) -> pd.DataFrame:
    """构建物种尺度汇总表，便于核验每个物种的权重是否被正确分配。"""
    if species_info.empty:
        return species_info

    # species_info 已在单组中记录过 coverage；这里重新聚合，处理跨组合合并后的情况。
    if contrib.empty:
        out = species_info.copy()
        out["rwr_contribution_sum"] = 0.0
        return out

    agg = (
        contrib.groupby(["species_id", "taxon", "status"], as_index=False)
        .agg(
            n_cells_positive_area=("grid_id", "nunique"),
            overlap_area_sum_km2=("overlap_area_km2", "sum"),
            species_total_area_km2=("species_total_area_km2", "first"),
            rwr_contribution_sum=("rwr_contribution", "sum"),
        )
    )
    agg["grid_coverage_ratio"] = agg["overlap_area_sum_km2"] / agg["species_total_area_km2"]

    all_species = species_info[["species_id", "taxon", "status", "species_total_area_km2"]].drop_duplicates()
    out = all_species.merge(
        agg,
        on=["species_id", "taxon", "status", "species_total_area_km2"],
        how="left",
    )
    out["n_cells_positive_area"] = out["n_cells_positive_area"].fillna(0).astype("int32")
    for c in ["overlap_area_sum_km2", "rwr_contribution_sum", "grid_coverage_ratio"]:
        out[c] = out[c].fillna(0.0).astype("float64")
    return out


# ============================================================
# 主流程
# ============================================================
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[GRID ] 读取网格: {GRID_PATH} | layer={GRID_LAYER}")
    grid = gpd.read_file(GRID_PATH, layer=GRID_LAYER, engine="pyogrio")
    if "grid_id" not in grid.columns:
        raise ValueError("网格图层缺少 grid_id 字段")
    print(f"[GRID ] 网格数量: {len(grid):,d}, CRS: {grid.crs}")

    # 只保留必要字段，降低内存占用。
    grid_small = grid[["grid_id", "geometry"]].copy()
    grid_small = repair_geometry(grid_small)

    all_contrib = []
    all_species_info = []

    for taxon in TAXA:
        for status in STATUSES:
            sp = read_species_data(taxon, status, grid_small.crs)
            if sp is None or sp.empty:
                continue
            contrib, species_info = compute_species_grid_contributions(sp, grid_small, taxon, status)
            if not contrib.empty:
                all_contrib.append(contrib)
            if species_info is not None and not species_info.empty:
                all_species_info.append(species_info)
            del sp, contrib, species_info
            gc.collect()

    if not all_contrib:
        raise RuntimeError("没有成功计算任何 species-grid 面积贡献，请检查 IUCN 数据路径或网格范围")

    contrib_all = pd.concat(all_contrib, ignore_index=True)
    species_info_all = pd.concat(all_species_info, ignore_index=True) if all_species_info else pd.DataFrame()

    print(f"[MERGE] species-grid area contribution rows: {len(contrib_all):,d}")
    print(f"[MERGE] unique species: {contrib_all['species_id'].nunique():,d}")
    print(f"[MERGE] unique grid:    {contrib_all['grid_id'].nunique():,d}")

    # ------------------------------------------------------------
    # 汇总到网格尺度
    # ------------------------------------------------------------
    stats = aggregate_rwr(contrib_all)

    print(
        f"[STATS] 含正面积物种分布的网格数: {len(stats):,d} / {len(grid):,d}; "
        f"rwr_total max = {stats['rwr_total'].max():.6f}, "
        f"mean(nonzero) = {stats.loc[stats['rwr_total'] > 0, 'rwr_total'].mean():.6f}"
    )

    # ------------------------------------------------------------
    # 物种尺度核验表
    # ------------------------------------------------------------
    species_area = build_species_area_table(contrib_all, species_info_all)
    if not species_area.empty:
        print(f"[WRITE] Species area table -> {SPECIES_AREA_CSV}")
        species_area.to_csv(SPECIES_AREA_CSV, index=False)
        print(
            f"[CHECK] species contribution sum total = {species_area['rwr_contribution_sum'].sum():.6f}; "
            f"species n = {species_area['species_id'].nunique():,d}; "
            f"mean grid coverage ratio = {species_area['grid_coverage_ratio'].mean():.6f}"
        )
        low_cov = (species_area["grid_coverage_ratio"] < 0.99).sum()
        if low_cov > 0:
            print(
                f"[NOTE ] {low_cov:,d} species have grid_coverage_ratio < 0.99. "
                f"This is expected when complete IUCN ranges include marine/water areas "
                f"or other areas outside the terrestrial grid domain. "
                f"These species are not forced to sum to 1 over land grids."
            )

    # ------------------------------------------------------------
    # 合并回网格 + 写出
    # ------------------------------------------------------------
    out = grid.merge(stats, on="grid_id", how="left")
    out["area_based_richness"] = out["area_based_richness"].fillna(0).astype("int32")
    float_cols = [
        "rwr_total", "weighted_richness",
        *[f"rwr_{t}" for t in TAXA],
        *[f"rwr_{s}" for s in STATUSES],
    ]
    for c in float_cols:
        if c not in out.columns:
            out[c] = 0.0
        out[c] = out[c].fillna(0.0).astype("float64")

    print(f"[WRITE] GeoPackage -> {OUT_GPKG}")
    out.to_file(OUT_GPKG, layer="Global_RWR_10km_grid", driver="GPKG")

    print(f"[WRITE] CSV       -> {OUT_CSV}")
    out.drop(columns="geometry").to_csv(OUT_CSV, index=False)

    print("[DONE ] IUCN-style area-based RWR calculation completed")


if __name__ == "__main__":
    main()
