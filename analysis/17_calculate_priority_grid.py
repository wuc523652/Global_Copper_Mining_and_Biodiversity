"""Recalculate priority-species grid associations from the original range geometry.

The scientific calculations and geometry filtering are preserved from the
study's priority-distribution analysis. This stage does not render figures.
The 347 original grid cells discarded after WGS84 conversion are not zeros.
Run through run_analysis.py after stages 01 and 07.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.validation import make_valid
from _paths import data_path, project_path

GRID_FILE = Path(project_path("01_grid/Global_10km_grid.gpkg"))
PRIORITY_FILE = Path(project_path("07_priority_species/Species_PrioritySpecies_Priority.csv"))
IUCN_DIR = Path(data_path("IUCNredlist/202502"))
OUTPUT_DIR = Path(project_path("17_priority_grid"))
TAXON_ORDER = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
STATUS_ORDER = ["CR", "EN", "VU"]
GROUP_ORDER = ["All"] + TAXON_ORDER
SPECIES_ID_FIELD_CANDIDATES = ["id_no", "ID_NO", "Id_no", "ID_No"]
TARGET_CRS = "EPSG:4326"
REMOVE_ANTARCTICA_BY_CENTROID_LAT = True
ANTARCTICA_LAT_THRESHOLD = -60.0
DOMINANT_COUNT_BINS = [
    {"lo": 1, "hi": 1, "label": "1"},
    {"lo": 2, "hi": 2, "label": "2"},
    {"lo": 3, "hi": 5, "label": "3–5"},
    {"lo": 6, "hi": 10, "label": "6–10"},
    {"lo": 11, "hi": None, "label": ">10"},
]


def find_grid_file():
    if not GRID_FILE.is_file():
        raise FileNotFoundError(f"Run stage 01 first: {GRID_FILE}")
    return GRID_FILE


def clean_species_id(x):
    """统一清理 species_id，避免 123.0 与 123 不一致。"""
    if pd.isna(x):
        return None

    s = str(x).strip()

    if s == "" or s.lower() == "nan":
        return None

    if s.endswith(".0"):
        s = s[:-2]

    return s

def read_vector(path):
    """读取矢量文件，优先使用 pyogrio。"""
    try:
        return gpd.read_file(path, engine="pyogrio")
    except Exception:
        return gpd.read_file(path)

def _is_nonempty_geometry(geom):
    """安全判断几何是否存在且非空，避免坏几何触发 GEOS 异常。"""
    try:
        return geom is not None and (not geom.is_empty)
    except Exception:
        return False

def _has_finite_bounds(geom):
    """检查几何 bounds 是否全部为有限数，用于剔除 NaN/Inf 坐标几何。"""
    if not _is_nonempty_geometry(geom):
        return False
    try:
        b = geom.bounds
        if b is None or len(b) != 4:
            return False
        return bool(np.all(np.isfinite(np.asarray(b, dtype="float64"))))
    except Exception:
        return False

def _safe_make_valid(geom):
    """安全修复单个几何；遇到 NaN/Inf 或无法修复的几何返回 None。"""
    if not _is_nonempty_geometry(geom):
        return None
    if not _has_finite_bounds(geom):
        return None
    try:
        fixed = make_valid(geom)
    except Exception:
        try:
            fixed = geom.buffer(0)
        except Exception:
            return None

    if not _is_nonempty_geometry(fixed):
        return None
    if not _has_finite_bounds(fixed):
        return None

    return fixed

def fix_geometries(gdf, name="vector"):
    """修复几何并删除空几何、NaN/Inf 坐标几何。"""
    if gdf.empty:
        return gdf

    before = len(gdf)

    gdf = gdf[gdf.geometry.apply(_is_nonempty_geometry)].copy()

    finite_mask = gdf.geometry.apply(_has_finite_bounds)
    finite_removed = int((~finite_mask).sum())
    if finite_removed > 0:
        print(f"  [GEOM] {name}: 删除含 NaN/Inf bounds 的几何 {finite_removed:,} 条")
    gdf = gdf[finite_mask].copy()

    gdf["geometry"] = gdf.geometry.apply(_safe_make_valid)
    gdf = gdf[gdf.geometry.apply(_is_nonempty_geometry)].copy()

    finite_mask = gdf.geometry.apply(_has_finite_bounds)
    finite_removed_2 = int((~finite_mask).sum())
    if finite_removed_2 > 0:
        print(f"  [GEOM] {name}: 修复后仍含 NaN/Inf bounds，删除 {finite_removed_2:,} 条")
    gdf = gdf[finite_mask].copy()

    after = len(gdf)
    if after < before:
        print(f"  [GEOM] {name}: 共删除无效、空或不可修复几何 {before - after:,} 条")

    return gdf

def clean_for_spatial_join(gdf, name="vector"):
    """空间连接前再次清理几何，专门防止 GEOS 的 NaN/Inf 报错。"""
    cleaned = fix_geometries(gdf, name=name)
    if cleaned.empty:
        raise ValueError(f"{name} 清理后为空，请检查输入几何。")
    return cleaned

def detect_species_id_field(columns):
    """自动识别 IUCN 数据中的物种 ID 字段。"""
    for field in SPECIES_ID_FIELD_CANDIDATES:
        if field in columns:
            return field

    raise ValueError(
        f"未找到物种 ID 字段。候选字段为 {SPECIES_ID_FIELD_CANDIDATES}，"
        f"当前字段为：{list(columns)}"
    )

def ensure_crs(gdf, default_epsg=4326, name="vector"):
    """若数据缺 CRS，则默认设置为 EPSG:4326。"""
    if gdf.crs is None:
        print(f"  [CRS] {name} 缺少 CRS，默认设置为 EPSG:{default_epsg}")
        gdf = gdf.set_crs(epsg=default_epsg)

    return gdf

def detect_grid_id_column(gdf):
    """识别或创建网格 ID 字段。"""
    candidate_cols = [
        "grid_id",
        "GRID_ID",
        "Grid_ID",
        "cell_id",
        "CELL_ID",
        "id",
        "ID",
        "FID",
        "OBJECTID",
    ]

    for col in candidate_cols:
        if col in gdf.columns:
            return col

    gdf["grid_id"] = np.arange(len(gdf), dtype=np.int64)
    return "grid_id"

def load_priority_species():
    """读取优先救助物种清单。"""
    if not PRIORITY_FILE.exists():
        raise FileNotFoundError(
            f"未找到优先救助物种清单：{PRIORITY_FILE}\n"
            "请先运行 07_screen_priority_species.py，生成 Species_PrioritySpecies_Priority.csv。"
        )

    df = pd.read_csv(PRIORITY_FILE)

    required_cols = ["species_id", "taxon", "status"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{PRIORITY_FILE.name} 缺少必要字段：{missing}")

    df["species_id"] = df["species_id"].apply(clean_species_id)
    df["taxon"] = df["taxon"].astype(str).str.strip()
    df["status"] = df["status"].astype(str).str.strip()

    df = df.dropna(subset=["species_id"]).copy()
    df = df[df["taxon"].isin(TAXON_ORDER)].copy()
    df = df[df["status"].isin(STATUS_ORDER)].copy()

    df = (
        df.sort_values(["taxon", "status", "species_id"])
        .drop_duplicates(subset=["species_id", "taxon", "status"], keep="first")
        .copy()
    )

    if df.empty:
        raise ValueError("优先救助物种清单为空，请检查输入文件。")

    print(f"优先救助物种数: {df['species_id'].nunique():,}")
    print("按类群统计：")
    print(df["taxon"].value_counts().reindex(TAXON_ORDER).fillna(0).astype(int))

    return df

def load_grid():
    """读取全球 10 km 陆地网格。"""
    grid_path = find_grid_file()

    print("\n读取全球 10 km 陆地网格：")
    print(f"  {grid_path}")

    grid = read_vector(grid_path)
    if grid.empty:
        raise ValueError(f"网格文件为空：{grid_path}")

    grid = ensure_crs(grid, default_epsg=4326, name="grid")
    grid = fix_geometries(grid, name="grid")

    grid_id_col = detect_grid_id_column(grid)

    if grid.crs.to_string() != TARGET_CRS:
        print(f"  网格投影转换到 {TARGET_CRS}: {grid.crs} -> {TARGET_CRS}")
        grid_4326 = grid[[grid_id_col, "geometry"]].to_crs(TARGET_CRS)
    else:
        grid_4326 = grid[[grid_id_col, "geometry"]].copy()

    # 投影转换后再次清理。极区或投影边界附近的少量网格可能出现 NaN/Inf 坐标，
    # 如果不提前删除，会在 gpd.sjoin/intersects 中触发 GEOSException。
    grid_4326 = fix_geometries(grid_4326, name="grid after to_crs")

    if REMOVE_ANTARCTICA_BY_CENTROID_LAT:
        # 用 bounds 中心纬度判断，避免 geographic CRS 下 centroid 的警告，
        # 也避免坏几何在 centroid 计算时再次触发异常。
        b = grid_4326.geometry.bounds
        center_lat = (b["miny"] + b["maxy"]) / 2.0
        before = len(grid_4326)
        grid_4326 = grid_4326[center_lat >= ANTARCTICA_LAT_THRESHOLD].copy()
        after = len(grid_4326)
        if after < before:
            print(f"  删除参与空间连接的疑似南极洲网格: {before - after:,} 个")
            print("  注意：南极洲底图面仍由国家边界矢量绘制，不影响地图显示南极洲。")

    grid_4326 = grid_4326.rename(columns={grid_id_col: "grid_id"})
    grid_4326["grid_id"] = grid_4326["grid_id"].astype(str)

    print(f"  网格数量: {len(grid_4326):,}")
    print(f"  网格 CRS: {grid_4326.crs}")

    return grid_4326

def load_priority_habitats(priority_df):
    """
    从 IUCN 原始栖息地数据中读取优先救助物种的分布范围。
    """
    print("\n读取优先救助物种栖息地矢量：")

    all_parts = []

    for taxon in TAXON_ORDER:
        for status in STATUS_ORDER:
            sub_ids = set(
                priority_df.loc[
                    (priority_df["taxon"] == taxon) & (priority_df["status"] == status),
                    "species_id"
                ].astype(str)
            )

            if len(sub_ids) == 0:
                continue

            shp = IUCN_DIR / f"{taxon}-{status}" / "data_0.shp"

            if not shp.exists():
                print(f"  警告：文件不存在，跳过：{shp}")
                continue

            print(f"  {taxon}-{status}: 读取 {len(sub_ids):,} 个优先救助物种")

            gdf = read_vector(shp)
            if gdf.empty:
                print(f"    警告：空文件，跳过：{shp}")
                continue

            id_field = detect_species_id_field(gdf.columns)
            gdf["species_id"] = gdf[id_field].apply(clean_species_id)
            gdf = gdf[gdf["species_id"].isin(sub_ids)].copy()

            if gdf.empty:
                print(f"    警告：未匹配到优先救助物种，跳过")
                continue

            gdf["taxon"] = taxon
            gdf["status"] = status

            gdf = ensure_crs(gdf, default_epsg=4326, name=f"{taxon}-{status}")
            gdf = fix_geometries(gdf, name=f"{taxon}-{status}")

            if gdf.crs.to_string() != TARGET_CRS:
                gdf = gdf.to_crs(TARGET_CRS)

            gdf = fix_geometries(gdf, name=f"{taxon}-{status} after to_crs")

            gdf = gdf[["species_id", "taxon", "status", "geometry"]].copy()
            all_parts.append(gdf)

            print(f"    读取斑块数: {len(gdf):,}")
            print(f"    物种数: {gdf['species_id'].nunique():,}")

    if len(all_parts) == 0:
        raise ValueError("没有读取到任何优先救助物种栖息地数据，请检查 IUCN_DIR 和 species_id。")

    habitat = pd.concat(all_parts, ignore_index=True)
    habitat = gpd.GeoDataFrame(habitat, geometry="geometry", crs=TARGET_CRS)
    habitat = fix_geometries(habitat, name="priority habitat merged")

    print(f"\n优先救助物种栖息地总斑块数: {len(habitat):,}")
    print(f"优先救助物种唯一数量: {habitat['species_id'].nunique():,}")

    return habitat

def dominant_count_bin(value):
    """返回 dominant taxon 丰富度的固定分级序号和标签。"""
    if pd.isna(value) or value <= 0:
        return 0, "No priority species"

    value = int(value)
    for index, item in enumerate(DOMINANT_COUNT_BINS, start=1):
        hi = item["hi"]
        if value >= item["lo"] and (hi is None or value <= hi):
            return index, item["label"]

    raise ValueError(f"dominant taxon 数量 {value} 未落入 DOMINANT_COUNT_BINS。")

def add_dominant_taxon_fields(count_df):
    """
    为每个网格确定物种数最高的优先干预类群。

    多类群并列最大值时，np.argmax 会选择 TAXON_ORDER 中最先出现的类群；
    dominant_taxon_tie_count 会保留并列数量，便于在结果表中识别这些网格。
    """
    out = count_df.copy()
    taxa_cols = [f"richness_{taxon}" for taxon in TAXON_ORDER]
    missing = [col for col in taxa_cols if col not in out.columns]
    if missing:
        raise ValueError(f"计算 dominant taxon 时缺少字段：{missing}")

    matrix = out[taxa_cols].to_numpy(dtype=int)
    maximum = matrix.max(axis=1)
    argmax = matrix.argmax(axis=1)
    tie_count = (matrix == maximum[:, None]).sum(axis=1)

    dominant = np.asarray([TAXON_ORDER[index] for index in argmax], dtype=object)
    dominant[maximum <= 0] = None
    tie_count[maximum <= 0] = 0

    out["dominant_taxon"] = dominant
    out["dominant_taxon_richness"] = maximum.astype(int)
    out["dominant_taxon_tie_count"] = tie_count.astype(int)
    out["dominant_taxon_is_tie"] = (tie_count > 1) & (maximum > 0)

    bin_results = out["dominant_taxon_richness"].apply(dominant_count_bin)
    out["dominant_count_bin"] = bin_results.apply(lambda item: int(item[0]))
    out["dominant_count_label"] = bin_results.apply(lambda item: item[1])
    return out

def calculate_grid_priority_richness(grid_4326, habitat_4326):
    """
    统计每个网格内涉及的优先救助物种数量。
    按 All 和五个 taxa 分别统计。
    """
    print("\n开始空间连接：网格 × 优先救助物种栖息地")

    grid_for_join = grid_4326[["grid_id", "geometry"]].copy()
    habitat_for_join = habitat_4326[["species_id", "taxon", "geometry"]].copy()

    grid_for_join = clean_for_spatial_join(grid_for_join, name="grid for spatial join")
    habitat_for_join = clean_for_spatial_join(habitat_for_join, name="habitat for spatial join")

    try:
        joined = gpd.sjoin(
            grid_for_join,
            habitat_for_join,
            how="inner",
            predicate="intersects"
        )
    except Exception as e:
        print("  空间连接第一次失败，正在强制重新清理几何后重试。")
        print(f"  原始错误: {type(e).__name__}: {e}")
        grid_for_join = clean_for_spatial_join(grid_for_join, name="grid for spatial join retry")
        habitat_for_join = clean_for_spatial_join(habitat_for_join, name="habitat for spatial join retry")
        joined = gpd.sjoin(
            grid_for_join,
            habitat_for_join,
            how="inner",
            predicate="intersects"
        )

    if joined.empty:
        raise ValueError("空间连接结果为空：没有网格与优先救助物种栖息地相交。")

    pairs = (
        joined[["grid_id", "species_id", "taxon"]]
        .drop_duplicates()
        .copy()
    )

    pairs["grid_id"] = pairs["grid_id"].astype(str)

    print(f"网格-物种唯一配对数: {len(pairs):,}")
    print(f"涉及网格数: {pairs['grid_id'].nunique():,}")
    print(f"涉及物种数: {pairs['species_id'].nunique():,}")

    out_pair_file = OUTPUT_DIR / "PrioritySpecies_Distribution_GridSpeciesPairs.csv"
    pairs.to_csv(out_pair_file, index=False, encoding="utf-8-sig")
    print(f"已保存网格-物种配对表: {out_pair_file}")

    count_df = grid_4326[["grid_id"]].copy()

    all_counts = (
        pairs.groupby("grid_id")["species_id"]
        .nunique()
        .reset_index(name="richness_All")
    )

    count_df = count_df.merge(all_counts, on="grid_id", how="left")

    for taxon in TAXON_ORDER:
        taxon_counts = (
            pairs[pairs["taxon"] == taxon]
            .groupby("grid_id")["species_id"]
            .nunique()
            .reset_index(name=f"richness_{taxon}")
        )
        count_df = count_df.merge(taxon_counts, on="grid_id", how="left")

    richness_cols = [c for c in count_df.columns if c.startswith("richness_")]
    for col in richness_cols:
        count_df[col] = count_df[col].fillna(0).astype(int)

    count_df = add_dominant_taxon_fields(count_df)

    out_count_file = OUTPUT_DIR / "PrioritySpecies_Distribution_GridCounts.csv"
    count_df.to_csv(out_count_file, index=False, encoding="utf-8-sig")
    print(f"已保存网格统计表: {out_count_file}")

    grid_counts = grid_4326.merge(count_df, on="grid_id", how="left")
    for col in richness_cols:
        grid_counts[col] = grid_counts[col].fillna(0).astype(int)

    nonzero = grid_counts[grid_counts["richness_All"] > 0].copy()
    out_gpkg = OUTPUT_DIR / "PrioritySpecies_Distribution_NonzeroGrid_4326.gpkg"

    if not nonzero.empty:
        nonzero.to_file(out_gpkg, driver="GPKG")
        print(f"已保存非零网格 GPKG: {out_gpkg}")

    summary_rows = []
    for group in GROUP_ORDER:
        col = f"richness_{group}"
        if col not in grid_counts.columns:
            continue

        vals = grid_counts[col]
        nz = vals[vals > 0]

        summary_rows.append({
            "Group": group,
            "num_nonzero_grids": int((vals > 0).sum()),
            "max_priority_species_richness": int(vals.max()),
            "mean_nonzero_priority_species_richness": float(nz.mean()) if len(nz) > 0 else 0,
            "median_nonzero_priority_species_richness": float(nz.median()) if len(nz) > 0 else 0,
            "total_grid_species_records_sum": int(vals.sum()),
        })

    summary_df = pd.DataFrame(summary_rows)
    out_summary = OUTPUT_DIR / "PrioritySpecies_Distribution_Summary.csv"
    summary_df.to_csv(out_summary, index=False, encoding="utf-8-sig")
    print(f"已保存统计摘要: {out_summary}")

    print("\n统计摘要：")
    print(summary_df.to_string(index=False))

    dominant_summary = (
        count_df[count_df["dominant_taxon_richness"] > 0]
        .groupby("dominant_taxon", observed=False)
        .agg(
            num_dominant_grids=("grid_id", "nunique"),
            max_dominant_species_richness=("dominant_taxon_richness", "max"),
            mean_dominant_species_richness=("dominant_taxon_richness", "mean"),
            num_tie_grids=("dominant_taxon_is_tie", "sum"),
        )
        .reindex(TAXON_ORDER)
        .fillna(0)
        .reset_index()
    )
    dominant_summary.to_csv(
        OUTPUT_DIR / "PrioritySpecies_Distribution_DominantTaxon_Summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    count_df[
        [
            "grid_id",
            "dominant_taxon",
            "dominant_taxon_richness",
            "dominant_count_bin",
            "dominant_count_label",
            "dominant_taxon_tie_count",
            "dominant_taxon_is_tie",
        ]
    ].to_csv(
        OUTPUT_DIR / "PrioritySpecies_Distribution_DominantTaxon_Grid.csv",
        index=False,
        encoding="utf-8-sig",
    )

    return grid_counts, pairs, summary_df

def main():
    required = [GRID_FILE, PRIORITY_FILE]
    for taxon in TAXON_ORDER:
        for status in STATUS_ORDER:
            shape = IUCN_DIR / f"{taxon}-{status}" / "data_0.shp"
            required.extend(shape.with_suffix(ext) for ext in (".shp", ".shx", ".dbf", ".prj"))
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input files: " + ", ".join(missing))
    if OUTPUT_DIR.exists() and any(OUTPUT_DIR.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {OUTPUT_DIR}")
    priority = load_priority_species()
    habitat = load_priority_habitats(priority)
    expected_species = set(priority["species_id"].astype(str))
    actual_species = set(habitat["species_id"].astype(str))
    if actual_species != expected_species:
        raise ValueError("Priority range geometries do not cover every selected species: "
                         + ", ".join(sorted(expected_species - actual_species)))
    grid = load_grid()
    if grid["grid_id"].duplicated().any():
        raise ValueError("The input grid contains duplicate identifiers")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    counts, pairs, _ = calculate_grid_priority_richness(grid, habitat)
    if pairs.duplicated(["grid_id", "species_id"]).any():
        raise ValueError("Duplicate grid-species pairs")
    expected_counts = pairs.groupby("grid_id")["species_id"].nunique()
    actual_counts = counts.set_index("grid_id")["richness_All"].reindex(expected_counts.index)
    if not np.array_equal(actual_counts.to_numpy(), expected_counts.to_numpy()):
        raise ValueError("Grid counts disagree with recomputed associations")
    record = {
        "status": "PASS",
        "algorithm": "WGS84 intersects; unique grid_id/species_id; original geometry filtering",
        "valid_grid_count": len(counts),
        "nonzero_grid_count": int(counts["richness_All"].gt(0).sum()),
        "grid_species_pairs": len(pairs),
        "priority_species": len(expected_species),
        "priority_species_in_grid": int(pairs["species_id"].nunique()),
        "source_grid": str(GRID_FILE),
        "source_priority_list": str(PRIORITY_FILE),
        "note": "Cells excluded during WGS84 geometry filtering are absent, not zero-richness cells.",
    }
    (OUTPUT_DIR / "Priority_Grid_RunRecord.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
