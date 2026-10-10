# -*- coding: utf-8 -*-
"""
11_calculate_mine_exposure.py

研究内容四——30 km范围内的矿区尺度“实际潜在暴露面积”。

关键定义（已按最新要求修正）：
    先将全部发生空间暴露的受威胁物种几何做空间并集，再与矿区30 km范围相交。
    同一地表无论覆盖多少个物种、多少个栖息地斑块，都只计一次。

本脚本直接复用研究内容一输出的：
    04_mine_species_overlap/Species_CopperOverlap_GlobalUnique.gpkg
其中每个物种的几何已经是该物种在全部铜矿30 km范围内的严格去重暴露几何。
将这些物种几何再次 union，正好得到“至少涉及一种目标物种”的实际暴露范围，
无须重新读取全部5,890种原始IUCN面。

面积层级：
1. mine_exposure_area_km2：单座矿区内部严格去重；
2. physical_exposure_union_area_km2：同一资源拥有国内不同矿区再次union；
3. global_strict_physical_exposure_area_km2：全球所有国家再次union。

注意：按矿区所属国归属的国家面积允许跨越国界；不同资源拥有国的缓冲区若共享
同一物理空间，国家归属面积之和可大于全球严格去重面积。该差额会单列QA。
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import time
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import unary_union

try:
    from shapely.validation import make_valid
except Exception:
    def make_valid(geom):
        return geom.buffer(0)


warnings.filterwarnings("ignore")


# =============================================================================
# 路径与参数
# =============================================================================

BASE_DIR = Path(_nc_project(""))

MINE_MASTER_GPKG = (
    BASE_DIR / "10_mine_ownership" / "Mine_Master_30km.gpkg"
)
MINE_MASTER_LAYER = "mine_master"

SPECIES_EXPOSURE_GPKG = (
    BASE_DIR / "04_mine_species_overlap" / "Species_CopperOverlap_GlobalUnique.gpkg"
)
SPECIES_EXPOSURE_LAYER = "global_unique_species_overlap"

OUTPUT_DIR = BASE_DIR / "11_mine_exposure"
OUT_MINE_CSV = OUTPUT_DIR / "Mine_Exposure_30km.csv"
OUT_COUNTRY_CSV = OUTPUT_DIR / "Country_Exposure_Union_30km.csv"
OUT_GLOBAL_CSV = OUTPUT_DIR / "Global_Exposure_Union_30km.csv"
OUT_SPECIES_COUNT = OUTPUT_DIR / "Mine_ExposedSpecies_Count.csv"
OUT_GEOMETRY_QA = OUTPUT_DIR / "Exposure_Geometry_QA.csv"
OUT_SPECIES_QA = OUTPUT_DIR / "Species_Geometry_Loading_QA.csv"
OUT_GPKG = OUTPUT_DIR / "Exposure_Union_30km.gpkg"

EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_CRS = "EPSG:4326"
UNION_BATCH_SIZE = 50
AREA_EPSILON_M2 = 1.0
CALCULATE_SPECIES_COUNTS = True

SPECIES_ID_CANDIDATES = ["species_id", "ID_NO", "id_no", "SCI_NAME", "sci_name"]


# =============================================================================
# 几何与读取工具
# =============================================================================

def read_vector(path: Path, layer=None) -> gpd.GeoDataFrame:
    if not path.exists():
        raise FileNotFoundError(f"输入空间文件不存在：{path}")
    kwargs = {"layer": layer} if layer else {}
    try:
        return gpd.read_file(path, engine="pyogrio", **kwargs)
    except Exception:
        return gpd.read_file(path, **kwargs)


def clean_polygon(geom):
    if geom is None or geom.is_empty:
        return None
    try:
        geom = make_valid(geom)
    except Exception:
        try:
            geom = geom.buffer(0)
        except Exception:
            return None
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type in {"Polygon", "MultiPolygon"}:
        return geom
    if geom.geom_type == "GeometryCollection":
        parts = [g for g in geom.geoms if g.geom_type in {"Polygon", "MultiPolygon"}]
        return unary_union(parts) if parts else None
    return None


def safe_intersection(left, right):
    if left is None or right is None or left.is_empty or right.is_empty:
        return None
    try:
        geom = left.intersection(right)
    except Exception:
        try:
            geom = make_valid(left).intersection(make_valid(right))
        except Exception:
            return None
    return clean_polygon(geom)


def union_in_batches(geometries, batch_size=UNION_BATCH_SIZE, label="geometry"):
    """分批union，降低一次性处理大量复杂MultiPolygon的内存峰值。"""
    valid = [clean_polygon(g) for g in geometries if g is not None and not g.is_empty]
    valid = [g for g in valid if g is not None and not g.is_empty]
    if not valid:
        return None
    level = valid
    round_no = 0
    while len(level) > 1:
        round_no += 1
        next_level = []
        for start in range(0, len(level), batch_size):
            merged = clean_polygon(unary_union(level[start:start + batch_size]))
            if merged is not None:
                next_level.append(merged)
        print(f"    {label} union 第{round_no}轮：{len(level):,} -> {len(next_level):,}")
        level = next_level
        if not level:
            return None
    return level[0]


def detect_species_id(columns):
    for field in SPECIES_ID_CANDIDATES:
        if field in columns:
            return field
    raise ValueError(f"物种暴露GPKG无法识别物种ID字段：{list(columns)}")


def clean_species_id(value):
    text = "" if value is None or pd.isna(value) else str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def write_layer(gdf, path, layer, mode):
    out = gdf.to_crs(OUTPUT_CRS).copy()
    out["geometry"] = out.geometry.map(clean_polygon)
    out = out[out.geometry.notna() & ~out.geometry.is_empty].copy()
    for col in out.select_dtypes(include="object").columns:
        out[col] = out[col].fillna("").astype(str)
    try:
        out.to_file(path, layer=layer, driver="GPKG", mode=mode, engine="pyogrio")
    except Exception:
        out.to_file(path, layer=layer, driver="GPKG", mode=mode)


# =============================================================================
# 核心计算
# =============================================================================

def load_inputs():
    print("\n[1/6] 读取矿区主表与研究内容一物种暴露几何")
    mines = read_vector(MINE_MASTER_GPKG, MINE_MASTER_LAYER)
    required_mine = {
        "property_id", "property_name", "resource_country", "resource_iso3",
        "resource_continent", "buffer_area_km2",
    }
    if not required_mine.issubset(mines.columns):
        raise ValueError(f"矿区主表缺少字段：{sorted(required_mine - set(mines.columns))}")
    if mines.crs is None:
        raise ValueError("矿区主表缺少CRS。")
    mines = mines.to_crs(EQUAL_AREA_CRS)
    mines["property_id"] = mines["property_id"].astype(str).str.strip()
    mines["geometry"] = mines.geometry.map(clean_polygon)
    mines = mines[mines.geometry.notna() & ~mines.geometry.is_empty].copy()
    if mines["property_id"].duplicated().any():
        raise ValueError("Mine_Master_30km.gpkg 的 property_id 重复。")

    species = read_vector(SPECIES_EXPOSURE_GPKG, SPECIES_EXPOSURE_LAYER)
    if species.crs is None:
        raise ValueError("物种暴露几何缺少CRS。")
    species = species.to_crs(EQUAL_AREA_CRS)
    id_field = detect_species_id(species.columns)
    species = species[[id_field, "geometry"]].rename(columns={id_field: "species_id"})
    species["species_id"] = species["species_id"].map(clean_species_id)
    species["geometry"] = species.geometry.map(clean_polygon)
    species = species[
        species["species_id"].ne("") & species.geometry.notna() & ~species.geometry.is_empty
    ].copy()
    # 理论上每物种一条；若有重复，先在物种内union。
    if species["species_id"].duplicated().any():
        species = species.dissolve(by="species_id", as_index=False)
        species["geometry"] = species.geometry.map(clean_polygon)

    print(f"  矿区：{len(mines):,}")
    print(f"  发生暴露的物种：{species['species_id'].nunique():,}")
    return mines.reset_index(drop=True), species.reset_index(drop=True)


def build_species_union(species):
    print("\n[2/6] 全部物种暴露几何union（物种间重叠只计一次）")
    t0 = time.time()
    exposure_union = union_in_batches(species.geometry.tolist(), label="species exposure")
    if exposure_union is None or exposure_union.is_empty:
        raise RuntimeError("全部物种暴露几何union结果为空。")
    union_area = exposure_union.area / 1e6
    species_area_sum = species.geometry.area.sum() / 1e6
    qa = pd.DataFrame([{
        "num_species_geometries": len(species),
        "num_unique_species": species["species_id"].nunique(),
        "species_geometry_area_sum_km2_not_used": species_area_sum,
        "all_species_union_area_km2_used": union_area,
        "cross_species_duplicate_area_km2": species_area_sum - union_area,
        "cross_species_duplicate_pct_of_union": (
            (species_area_sum - union_area) / union_area * 100 if union_area > 0 else 0.0
        ),
        "elapsed_seconds": time.time() - t0,
        "definition": "union across all exposed species; each physical location counted once",
    }])
    qa.to_csv(OUT_SPECIES_QA, index=False, encoding="utf-8-sig")
    print(f"  实际物理暴露union面积：{union_area:,.2f} km²")
    print(f"  已消除跨物种重复面积：{species_area_sum - union_area:,.2f} km²")
    return exposure_union


def calculate_species_counts(mines, species):
    if not CALCULATE_SPECIES_COUNTS:
        return pd.DataFrame({"property_id": mines["property_id"], "num_exposed_species": np.nan})
    print("\n[3/6] 统计每座矿区涉及的不同物种数（仅计数，不叠加物种面积）")
    left = mines[["property_id", "geometry"]].copy()
    right = species[["species_id", "geometry"]].copy()
    candidates = gpd.sjoin(left, right, how="inner", predicate="intersects")
    if candidates.empty:
        counts = pd.DataFrame(columns=["property_id", "num_exposed_species"])
    else:
        species_geom = species.geometry.to_dict()
        mine_geom = mines.set_index("property_id").geometry.to_dict()
        valid_pairs = []
        for row in candidates[["property_id", "species_id", "index_right"]].itertuples(index=False):
            inter = safe_intersection(mine_geom[row.property_id], species_geom[row.index_right])
            if inter is not None and inter.area > AREA_EPSILON_M2:
                valid_pairs.append((row.property_id, row.species_id))
        pairs = pd.DataFrame(valid_pairs, columns=["property_id", "species_id"]).drop_duplicates()
        counts = pairs.groupby("property_id", as_index=False).agg(
            num_exposed_species=("species_id", "nunique")
        )
    counts = mines[["property_id"]].merge(counts, on="property_id", how="left")
    counts["num_exposed_species"] = counts["num_exposed_species"].fillna(0).astype(int)
    counts.to_csv(OUT_SPECIES_COUNT, index=False, encoding="utf-8-sig")
    print(f"  涉及至少一种目标物种的矿区：{(counts['num_exposed_species'] > 0).sum():,}")
    return counts


def calculate_mine_exposure(mines, exposure_union, species_counts):
    print("\n[4/6] 计算单矿区实际暴露面积")
    rows = []
    geoms = []
    total = len(mines)
    for i, mine in enumerate(mines.itertuples(index=False), 1):
        if i == 1 or i % 100 == 0 or i == total:
            print(f"  矿区 {i:,}/{total:,}")
        inter = safe_intersection(mine.geometry, exposure_union)
        area = 0.0 if inter is None else inter.area / 1e6
        rows.append({
            "property_id": mine.property_id,
            "property_name": mine.property_name,
            "resource_country": mine.resource_country,
            "resource_iso3": mine.resource_iso3,
            "resource_continent": mine.resource_continent,
            "buffer_area_km2": float(mine.buffer_area_km2),
            "mine_exposure_area_km2": area,
            "exposure_share_of_buffer_pct": (
                area / float(mine.buffer_area_km2) * 100 if float(mine.buffer_area_km2) > 0 else 0.0
            ),
        })
        geoms.append(inter)

    table = pd.DataFrame(rows).merge(species_counts, on="property_id", how="left", validate="one_to_one")
    table["num_exposed_species"] = table["num_exposed_species"].fillna(0).astype(int)
    table.to_csv(OUT_MINE_CSV, index=False, encoding="utf-8-sig")

    gdf = gpd.GeoDataFrame(table, geometry=geoms, crs=EQUAL_AREA_CRS)
    positive = gdf[gdf["mine_exposure_area_km2"] > 0].copy()
    print(f"  正面积暴露矿区：{len(positive):,}/{len(table):,}")
    return table, positive


def calculate_country_and_global(mine_table, mine_positive, exposure_union):
    print("\n[5/6] 国家与全球再次union")
    country_rows = []
    country_geoms = []
    meta = mine_table.groupby(["resource_iso3", "resource_country", "resource_continent"], as_index=False).agg(
        num_mines=("property_id", "nunique"),
        num_exposed_mines=("mine_exposure_area_km2", lambda x: int((x > 0).sum())),
        mine_exposure_area_sum_km2=("mine_exposure_area_km2", "sum"),
    )
    positive_groups = {key: group for key, group in mine_positive.groupby("resource_iso3")}
    for row in meta.itertuples(index=False):
        group = positive_groups.get(row.resource_iso3)
        union_geom = None if group is None else union_in_batches(
            group.geometry.tolist(), label=f"{row.resource_iso3} mine exposure"
        )
        union_area = 0.0 if union_geom is None else union_geom.area / 1e6
        country_rows.append({
            "resource_iso3": row.resource_iso3,
            "resource_country": row.resource_country,
            "resource_continent": row.resource_continent,
            "num_mines": int(row.num_mines),
            "num_exposed_mines": int(row.num_exposed_mines),
            "mine_exposure_area_sum_km2": float(row.mine_exposure_area_sum_km2),
            "physical_exposure_union_area_km2": union_area,
            "within_country_cross_mine_duplicate_area_km2": max(
                0.0, float(row.mine_exposure_area_sum_km2) - union_area
            ),
        })
        country_geoms.append(union_geom)

    country = pd.DataFrame(country_rows)
    country["within_country_cross_mine_duplicate_pct_of_union"] = np.where(
        country["physical_exposure_union_area_km2"] > 0,
        country["within_country_cross_mine_duplicate_area_km2"]
        / country["physical_exposure_union_area_km2"] * 100,
        0.0,
    )
    country.to_csv(OUT_COUNTRY_CSV, index=False, encoding="utf-8-sig")
    country_gdf = gpd.GeoDataFrame(country, geometry=country_geoms, crs=EQUAL_AREA_CRS)
    country_positive = country_gdf[country_gdf["physical_exposure_union_area_km2"] > 0].copy()

    strict_global_area = exposure_union.area / 1e6
    country_sum = country["physical_exposure_union_area_km2"].sum()
    mine_sum = mine_table["mine_exposure_area_km2"].sum()
    global_table = pd.DataFrame([{
        "global_strict_physical_exposure_area_km2": strict_global_area,
        "country_attributed_union_sum_km2": country_sum,
        "cross_resource_country_duplicate_area_km2": max(0.0, country_sum - strict_global_area),
        "mine_exposure_area_sum_km2": mine_sum,
        "within_country_cross_mine_duplicate_area_km2": max(0.0, mine_sum - country_sum),
        "num_mines": mine_table["property_id"].nunique(),
        "num_exposed_mines": int((mine_table["mine_exposure_area_km2"] > 0).sum()),
        "definition": "strict union across species; country attribution follows mine host-country attribute",
    }])
    global_table.to_csv(OUT_GLOBAL_CSV, index=False, encoding="utf-8-sig")
    return country, country_positive, global_table


def save_geometries(exposure_union, mine_positive, country_positive):
    print("\n[6/6] 保存实际暴露几何与QA")
    if OUT_GPKG.exists():
        OUT_GPKG.unlink()
    global_gdf = gpd.GeoDataFrame(
        [{"scope": "Global", "area_km2": exposure_union.area / 1e6}],
        geometry=[exposure_union], crs=EQUAL_AREA_CRS,
    )
    write_layer(global_gdf, OUT_GPKG, "global_exposure_union", "w")
    write_layer(mine_positive, OUT_GPKG, "mine_exposure_union", "a")
    write_layer(country_positive, OUT_GPKG, "country_exposure_union", "a")

    qa = pd.DataFrame([
        {"check": "mine_positive_geometry_count", "value": len(mine_positive)},
        {"check": "country_positive_geometry_count", "value": len(country_positive)},
        {"check": "all_species_union_is_valid", "value": int(exposure_union.is_valid)},
        {"check": "all_species_union_area_km2", "value": exposure_union.area / 1e6},
        {"check": "area_definition", "value": "physical union; no summation across species"},
    ])
    qa.to_csv(OUT_GEOMETRY_QA, index=False, encoding="utf-8-sig")
    print(f"  {OUT_GPKG}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mines, species = load_inputs()
    exposure_union = build_species_union(species)
    species_counts = calculate_species_counts(mines, species)
    mine_table, mine_positive = calculate_mine_exposure(mines, exposure_union, species_counts)
    _, country_positive, _ = calculate_country_and_global(
        mine_table, mine_positive, exposure_union
    )
    save_geometries(exposure_union, mine_positive, country_positive)
    print("\n完成。面积已按物理空间union，同一位置涉及多个物种时只计一次。")
    print("下一步运行 12_attribute_exposure_by_ownership.py。")


if __name__ == "__main__":
    main()
