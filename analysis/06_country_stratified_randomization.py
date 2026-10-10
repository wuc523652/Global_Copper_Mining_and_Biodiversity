# -*- coding: utf-8 -*-
"""
06_country_stratified_randomization.py

研究内容一的精简版随机位置对照实验：国家分层的 10 km RWR 网格位置对照。

检验问题
--------
在严格保持各国铜矿数量不变的条件下，真实铜矿参考点是否比所属国内随机
陆地网格更容易位于全球 RWR 前 10% 网格？

与旧版脚本的区别
----------------
旧版在每次模拟中平移 1,278 个矿区面，随后反复执行按国 union、全球 union
以及与约 140 万个 RWR 网格的精确相交。新版只在预处理阶段完成一次：

1. 铜矿参考点（30 km 缓冲区在等面积投影下的质心）与 RWR 网格匹配；
2. RWR 网格中心点与国家匹配；
3. 9,999 次随机化仅在 NumPy 数组中按国家有放回抽取网格。

因此，本脚本检验的是“矿山位置是否非随机集中于高 RWR 景观”，而不是再次
随机化 30 km 潜在影响范围的精确面积。30 km 缓冲区面积叠加仍由研究内容一
的主分析代码负责。

主要指标
--------
T = 位于全球 RWR 前 10% 网格的矿山数 / 全部矿山数 × 100%

同时使用矿山所在网格的平均 RWR 全球百分位作为低成本连续敏感性指标。
报告随机均值、95% 区间、观测－随机差、富集倍数、SES，以及采用“加 1”
修正的单尾经验 P 值。

国家尺度结果仅作描述：默认展示矿山数不少于 10 座的国家，比较真实矿山
比例与该国陆地网格背景比例；不进行逐国显著性检验和 FDR 校正。

注意
----
1. 国家归属取自铜矿属性，与研究内容一既定的“按矿区所属国归属”口径一致；
2. 随机抽样有放回，同一网格可被多座矿山抽中；
3. 全局零分布是各国随机结果按真实矿山数加权后的组合，而不是全球不分层抽样；
4. 显著结果不能解释为生物多样性决定矿山选址，也未控制矿床地质条件、
   交通可达性等选址因素。
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import re
import unicodedata
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS, Geod
from shapely.errors import GEOSException
from shapely.ops import transform
from shapely.validation import make_valid


# ============================================================
# 路径与运行参数：沿用旧版脚本
# ============================================================

MINE_BUFFER_PATH = Path(
    _nc_data("Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp")
)
GADM_PATH = Path(_nc_data("Global_GADM/GC_GADM.shp"))
RWR_GPKG = Path(
    _nc_project("02_rwr/Global_RWR_10km_grid.gpkg")
)
RWR_LAYER = "Global_RWR_10km_grid"

PROJECT_DIR = Path(_nc_project(""))
BASE_OUT_DIR = PROJECT_DIR / "06_randomization"

N_ITERATIONS = 9999
RANDOM_SEED = 20260728
PERCENTILE_BASE = "all"  # 与 RWR 主分析一致；可选 "nonzero"
CREATE_FIGURES = False

# 国家图仅展示矿山数不少于该值的国家；仅作描述，不计算国家 P 值。
COUNTRY_MIN_MINES_FOR_FIG = 10
COUNTRY_MAX_FOR_FIG = 30

# 少量沿海网格中心可能落在国界面之外，可在等面积投影下匹配最近国家。
GRID_COUNTRY_NEAREST_MAX_KM = 25.0

# 仅用于先验证流程。正式分析必须保持 False。
TEST_MODE = False
TEST_N_ITERATIONS = 99
TEST_MAX_MINES = 50

ACTIVE_N_ITERATIONS = TEST_N_ITERATIONS if TEST_MODE else N_ITERATIONS
OUT_DIR = BASE_OUT_DIR / "TEST_OUTPUT_DO_NOT_PUBLISH" if TEST_MODE else BASE_OUT_DIR
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 保留旧版主要输出文件名，方便接入现有分析流程。
OUT_OBSERVED = OUT_DIR / "Mine_RWR_Observed.csv"
OUT_ITERATIONS = OUT_DIR / "Mine_RWR_Randomization_Iterations.csv"
OUT_SUMMARY = OUT_DIR / "Mine_RWR_Randomization_Summary.csv"
OUT_THRESHOLDS = OUT_DIR / "Mine_RWR_Randomization_BinThresholds.csv"
OUT_MINE_QA = OUT_DIR / "Mine_Randomization_Input_QA.csv"
OUT_UNMATCHED_COUNTRIES = OUT_DIR / "Unmatched_Mine_Countries.csv"
OUT_CONFIG = OUT_DIR / "Mine_RWR_Randomization_Config.csv"

# 新版增加的轻量 QA/描述性结果。
OUT_COUNTRY_DESC = OUT_DIR / "Mine_RWR_Country_Descriptive.csv"
OUT_GRID_COUNTRY_QA = OUT_DIR / "RWR_Grid_Country_Background_QA.csv"

OUT_GLOBAL_FIG_PNG = OUT_DIR / "FigS_Randomization_Global.png"
OUT_GLOBAL_FIG_PDF = OUT_DIR / "FigS_Randomization_Global.pdf"
# 为兼容旧路径保留 CountryTop30 文件名；实际只纳入 n >= 10 的国家，最多 30 个。
OUT_COUNTRY_FIG_PNG = OUT_DIR / "FigS_Randomization_CountryTop30.png"
OUT_COUNTRY_FIG_PDF = OUT_DIR / "FigS_Randomization_CountryTop30.pdf"
OUT_COMBINED_FIG_PNG = OUT_DIR / "FigS3_Randomization_GridPointControl.png"
OUT_COMBINED_FIG_PDF = OUT_DIR / "FigS3_Randomization_GridPointControl.pdf"

EQUAL_AREA_CRS = "ESRI:54009"
COUNTRY_FIELD = "Country"
MINE_ID_FIELD = "mine_id"
GRID_ID_FIELD = "grid_id"
RWR_FIELD = "weighted_richness"
GEOD = Geod(ellps="WGS84")

GADM_COUNTRY_FIELD_CANDIDATES = [
    "Country", "COUNTRY", "country", "NAME_0", "NAME", "CNTRY_NAME",
    "ADMIN", "admin", "SOVEREIGNT",
]
MINE_ID_FIELD_CANDIDATES = [
    "Property_I", "Property_ID", "PROPERTY_ID",
    "mine_id", "Mine_ID", "MINE_ID", "id", "ID", "OBJECTID", "FID", "FID_1",
    "SNL_ID", "snl_id", "Project_ID", "PROJECT_ID", "project_id",
    "MineName", "MINE_NAME", "mine_name", "Name", "NAME",
]
MINE_COUNTRY_FIELD_CANDIDATES = [
    "Country", "COUNTRY", "country", "Mine_Country", "MINE_COUNTRY", "mine_country",
    "Host_Country", "HOST_COUNTRY", "host_country", "Project_Country",
    "PROJECT_COUNTRY", "project_country", "SNL_Country", "SNL_COUNTRY", "snl_country",
    "CNTRY_NAME", "ADMIN", "admin", "COUNTRY_NA", "Country_Na", "COUNTRYNAM",
    "ResCountry", "RESCOUNTRY", "Host_Cou", "HOST_COU",
]

BIN_ORDER = ["le50", "p50_75", "p75_90", "p90_99", "p99_100"]
BIN_LABELS = {
    "le50": "≤50%",
    "p50_75": "50–75%",
    "p75_90": "75–90%",
    "p90_99": "90–99%",
    "p99_100": "99–100%",
}

warnings.filterwarnings("ignore")


# ============================================================
# 基础工具与国家名匹配
# ============================================================

def read_vector(path: Path, layer=None):
    if not path.exists():
        raise FileNotFoundError(f"文件不存在：{path}")
    try:
        if layer is None:
            return gpd.read_file(path, engine="pyogrio")
        return gpd.read_file(path, layer=layer, engine="pyogrio")
    except Exception:
        if layer is None:
            return gpd.read_file(path)
        return gpd.read_file(path, layer=layer)


def detect_field(columns, candidates, label):
    for field in candidates:
        if field in columns:
            return field
    raise ValueError(f"未找到{label}字段。候选：{candidates}；实际：{list(columns)}")


def normalize_name(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    value = unicodedata.normalize("NFKD", str(value).strip())
    value = "".join(char for char in value if not unicodedata.combining(char)).lower()
    value = value.replace("&", "and")
    value = re.sub(r"[\(\)\[\]\{\},.;:'\"`´’‘]", " ", value)
    value = re.sub(r"[-_/]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


COUNTRY_ALIASES = {
    "usa": "United States", "u s a": "United States", "us": "United States",
    "u s": "United States", "united states of america": "United States",
    "russian federation": "Russia", "russia federation": "Russia",
    "viet nam": "Vietnam", "lao pdr": "Laos",
    "lao people s democratic republic": "Laos",
    "bolivia plurinational state of": "Bolivia",
    "venezuela bolivarian republic of": "Venezuela",
    "iran islamic republic of": "Iran", "iran islamic rep": "Iran",
    "tanzania united republic of": "Tanzania", "united republic of tanzania": "Tanzania",
    "democratic republic of congo": "Democratic Republic of the Congo",
    "democratic republic of the congo": "Democratic Republic of the Congo",
    "dr congo": "Democratic Republic of the Congo", "drc": "Democratic Republic of the Congo",
    "congo kinshasa": "Democratic Republic of the Congo",
    "congo democratic republic": "Democratic Republic of the Congo",
    "dem rep congo": "Democratic Republic of the Congo",
    "republic of congo": "Republic of Congo", "congo brazzaville": "Republic of Congo",
    "cote d ivoire": "Côte d'Ivoire", "ivory coast": "Côte d'Ivoire",
    "turkiye": "Turkey", "czechia": "Czech Republic", "north macedonia": "Macedonia",
    "eswatini": "Swaziland", "brunei darussalam": "Brunei",
    "korea republic of": "South Korea", "republic of korea": "South Korea",
    "korea dem people s republic": "North Korea",
    "democratic people s republic of korea": "North Korea",
    "moldova republic of": "Moldova", "syrian arab republic": "Syria",
    "cabo verde": "Cape Verde", "slovak republic": "Slovakia",
    "myanmar burma": "Myanmar", "burma": "Myanmar",
    "palestine": "Palestina", "state of palestine": "Palestina",
    "east timor": "Timor-Leste", "micronesia federated states of": "Micronesia",
    "federated states of micronesia": "Micronesia",
    "sao tome and principe": "São Tomé and Príncipe",
    "sao tome principe": "São Tomé and Príncipe",
    "saint barthelemy": "Saint-Barthélemy", "saint martin": "Saint-Martin",
    "curacao": "Curaçao", "aland": "Åland", "aland islands": "Åland",
}


def drop_zm(geom):
    if geom is None or geom.is_empty:
        return None
    try:
        return transform(lambda x, y, *args: (x, y), geom)
    except Exception:
        return geom


def repair_geometry(geom):
    if geom is None or geom.is_empty:
        return None
    geom = drop_zm(geom)
    try:
        repaired = make_valid(geom)
        # 极少数复杂国家边界在 make_valid 后仍可能被 GEOS 判为无效。
        # 再用 buffer(0) 作一次保守兜底，避免后续空间谓词直接崩溃。
        if repaired is not None and not repaired.is_empty and not repaired.is_valid:
            repaired = repaired.buffer(0)
        return repaired
    except Exception:
        try:
            return geom.buffer(0)
        except Exception:
            return None


def clean_gdf(gdf, name, repair_invalid=True):
    before = len(gdf)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if repair_invalid:
        invalid = ~gdf.geometry.is_valid
        if invalid.any():
            print(f"[GEOM] {name}: 修复无效几何 {int(invalid.sum()):,} 条")
            gdf.loc[invalid, "geometry"] = gdf.loc[invalid].geometry.apply(repair_geometry)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if repair_invalid:
        still_invalid = ~gdf.geometry.is_valid
        if still_invalid.any():
            bad_rows = gdf.index[still_invalid].tolist()[:10]
            raise ValueError(
                f"{name} 中仍有 {int(still_invalid.sum())} 条几何无法修复；"
                f"行索引示例：{bad_rows}"
            )
    if len(gdf) < before:
        print(f"[GEOM] {name}: 删除无效/空几何 {before - len(gdf):,} 条")
    return gdf


def safe_covers(country, country_geom, point):
    """稳健执行国家面 covers 点，并在 GEOS 拓扑异常时就地修复一次。"""
    try:
        return bool(country_geom.covers(point))
    except GEOSException as exc:
        repaired = repair_geometry(country_geom)
        if repaired is None or repaired.is_empty or not repaired.is_valid:
            raise ValueError(
                f"国家 {country!r} 的边界几何无效，自动修复失败。"
            ) from exc
        try:
            return bool(repaired.covers(point))
        except GEOSException as second_exc:
            raise ValueError(
                f"国家 {country!r} 的边界几何在修复后仍无法执行 covers。"
            ) from second_exc


def project_to(gdf, crs, name):
    if gdf.crs is None:
        raise ValueError(f"{name} 缺少 CRS。")
    if CRS.from_user_input(gdf.crs).equals(CRS.from_user_input(crs)):
        return gdf.copy()
    return gdf.to_crs(crs)




# ============================================================
# 输入准备
# ============================================================

def prepare_countries():
    print("\n[1/8] 读取并整理 GADM 国家面")
    countries = read_vector(GADM_PATH)
    field = detect_field(countries.columns, GADM_COUNTRY_FIELD_CANDIDATES, "GADM 国家名")
    countries = countries[[field, "geometry"]].rename(columns={field: COUNTRY_FIELD})
    countries[COUNTRY_FIELD] = countries[COUNTRY_FIELD].astype(str).str.strip()
    countries = clean_gdf(countries, "GADM", repair_invalid=True)
    countries = countries.dissolve(by=COUNTRY_FIELD, as_index=False)
    countries = clean_gdf(countries, "GADM dissolved", repair_invalid=True)
    if countries[COUNTRY_FIELD].duplicated().any():
        raise ValueError("GADM 国家面 dissolve 后国家名仍不唯一。")
    print(f"  国家数: {len(countries):,}")
    return countries


def build_country_lookup(countries):
    names = countries[COUNTRY_FIELD].tolist()
    lookup = {normalize_name(name): name for name in names}
    valid = set(names)
    for alias, target in COUNTRY_ALIASES.items():
        if target in valid:
            lookup[normalize_name(alias)] = target
    return lookup


def prepare_rwr_grid():
    print("\n[2/8] 读取 RWR 网格并计算全球百分位")
    grid = read_vector(RWR_GPKG, RWR_LAYER)
    for field in [GRID_ID_FIELD, RWR_FIELD]:
        if field not in grid.columns:
            raise ValueError(f"RWR 网格缺少字段 {field}。")

    grid = grid[[GRID_ID_FIELD, RWR_FIELD, "geometry"]].copy()
    grid[GRID_ID_FIELD] = grid[GRID_ID_FIELD].astype(str)
    if grid[GRID_ID_FIELD].duplicated().any():
        raise ValueError("RWR grid_id 不唯一。")
    grid[RWR_FIELD] = pd.to_numeric(grid[RWR_FIELD], errors="coerce").fillna(0.0)
    grid = clean_gdf(grid, "RWR grid", repair_invalid=False).reset_index(drop=True)
    grid["_grid_pos"] = np.arange(len(grid), dtype=np.int64)

    if PERCENTILE_BASE == "all":
        values_for_q = grid[RWR_FIELD]
    elif PERCENTILE_BASE == "nonzero":
        values_for_q = grid.loc[grid[RWR_FIELD] > 0, RWR_FIELD]
    else:
        raise ValueError("PERCENTILE_BASE 只能为 'all' 或 'nonzero'。")
    if values_for_q.empty:
        raise ValueError("用于计算 RWR 百分位的网格为空。")

    quantiles = values_for_q.quantile([0.50, 0.75, 0.90, 0.99, 1.00])
    thresholds = {
        "p50": float(quantiles.loc[0.50]),
        "p75": float(quantiles.loc[0.75]),
        "p90": float(quantiles.loc[0.90]),
        "p99": float(quantiles.loc[0.99]),
        "p100": float(quantiles.loc[1.00]),
    }

    def classify(values):
        return np.select(
            [
                values <= thresholds["p50"],
                values <= thresholds["p75"],
                values <= thresholds["p90"],
                values <= thresholds["p99"],
            ],
            ["le50", "p50_75", "p75_90", "p90_99"],
            default="p99_100",
        )

    grid["rwr_bin"] = classify(grid[RWR_FIELD].to_numpy(float))
    # 与五级划分严格一致：只有数值大于第 90 百分位阈值才属于前 10%。
    grid["is_top10"] = (grid[RWR_FIELD] > thresholds["p90"]).astype(np.int8)

    # 连续敏感性指标使用经验分布函数百分位；RWR <= 0 明确赋值为 0。
    if PERCENTILE_BASE == "all":
        percentile = (
            grid[RWR_FIELD].rank(method="max", pct=True).to_numpy(float) * 100.0
        )
    else:
        percentile = np.zeros(len(grid), dtype=float)
        positive_mask = grid[RWR_FIELD].to_numpy(float) > 0
        percentile[positive_mask] = (
            grid.loc[positive_mask, RWR_FIELD]
            .rank(method="max", pct=True)
            .to_numpy(float)
            * 100.0
        )
    percentile[grid[RWR_FIELD].to_numpy(float) <= 0] = 0.0
    grid["rwr_percentile"] = percentile

    threshold_rows = [
        {
            "rwr_bin": code,
            "rwr_bin_label": BIN_LABELS[code],
            "lower_value": (
                float(grid[RWR_FIELD].min())
                if code == "le50"
                else thresholds[
                    {"p50_75": "p50", "p75_90": "p75", "p90_99": "p90", "p99_100": "p99"}[code]
                ]
            ),
            "upper_value": thresholds[
                {"le50": "p50", "p50_75": "p75", "p75_90": "p90", "p90_99": "p99", "p99_100": "p100"}[code]
            ],
            "percentile_base": PERCENTILE_BASE,
            "n_grid_used_for_percentile": len(values_for_q),
            "n_grid_in_bin": int((grid["rwr_bin"] == code).sum()),
        }
        for code in BIN_ORDER
    ]
    pd.DataFrame(threshold_rows).to_csv(OUT_THRESHOLDS, index=False, encoding="utf-8-sig")

    print(f"  RWR 网格数: {len(grid):,}")
    print(f"  第 90 百分位阈值: {thresholds['p90']:.8g}")
    print(f"  实际前 10% 网格数: {int(grid['is_top10'].sum()):,} ({grid['is_top10'].mean() * 100:.2f}%)")
    return grid, thresholds


def prepare_mines(countries, grid_crs):
    print("\n[3/8] 读取铜矿缓冲区并生成参考点")
    mine = read_vector(MINE_BUFFER_PATH)
    id_field = next((field for field in MINE_ID_FIELD_CANDIDATES if field in mine.columns), None)
    country_field = detect_field(mine.columns, MINE_COUNTRY_FIELD_CANDIDATES, "矿山所属国")

    if id_field is None:
        mine[MINE_ID_FIELD] = [f"mine_{i + 1}" for i in range(len(mine))]
    else:
        mine[MINE_ID_FIELD] = mine[id_field].astype(str).str.strip()
        missing_id = mine[MINE_ID_FIELD].isin(["", "nan", "None", "NULL"])
        replacement = [f"mine_missing_{i + 1}" for i in range(int(missing_id.sum()))]
        mine.loc[missing_id, MINE_ID_FIELD] = replacement

    country_lookup = build_country_lookup(countries)
    raw_country = mine[country_field].astype(str).str.strip()
    mine[COUNTRY_FIELD] = raw_country.map(
        lambda value: country_lookup.get(normalize_name(value), "")
    )
    unmatched = pd.DataFrame({
        "mine_country_raw": raw_country[mine[COUNTRY_FIELD] == ""]
    })
    if not unmatched.empty:
        unmatched = (
            unmatched.value_counts()
            .rename("mine_count")
            .reset_index()
        )
        unmatched.to_csv(OUT_UNMATCHED_COUNTRIES, index=False, encoding="utf-8-sig")
        raise ValueError(
            f"有 {int(unmatched['mine_count'].sum())} 座矿山的国家名无法匹配 GADM；"
            f"清单已保存：{OUT_UNMATCHED_COUNTRIES}"
        )

    mine = mine[[MINE_ID_FIELD, COUNTRY_FIELD, "geometry"]].copy()
    mine = clean_gdf(mine, "mine buffers", repair_invalid=True)
    mine = mine.dissolve(by=[MINE_ID_FIELD, COUNTRY_FIELD], as_index=False)
    mine = clean_gdf(mine, "mine buffers dissolved", repair_invalid=True)

    multi_country = mine.groupby(MINE_ID_FIELD)[COUNTRY_FIELD].nunique()
    if (multi_country > 1).any():
        raise ValueError("存在同一 mine_id 对应多个所属国，无法保持国家分层。")

    if TEST_MODE:
        mine = (
            mine.sort_values([COUNTRY_FIELD, MINE_ID_FIELD])
            .head(TEST_MAX_MINES)
            .copy()
        )
        print(f"[TEST] 仅使用前 {len(mine)} 座矿山，结果不可用于论文。")

    # 在等面积投影下计算缓冲区质心，再转回 RWR 网格 CRS。
    mine_equal_area = project_to(mine, EQUAL_AREA_CRS, "mine buffers")
    reference_equal_area = gpd.GeoDataFrame(
        mine[[MINE_ID_FIELD, COUNTRY_FIELD]].copy(),
        geometry=mine_equal_area.geometry.centroid,
        crs=EQUAL_AREA_CRS,
    )
    reference = project_to(reference_equal_area, grid_crs, "mine reference points")
    reference = reference.reset_index(drop=True)
    reference["_mine_pos"] = np.arange(len(reference), dtype=np.int64)
    reference_wgs84 = reference.to_crs("EPSG:4326")
    reference["reference_lon"] = reference_wgs84.geometry.x.to_numpy()
    reference["reference_lat"] = reference_wgs84.geometry.y.to_numpy()

    countries_in_grid_crs = project_to(countries, grid_crs, "GADM")
    # dissolve 后的国家面在重投影过程中仍可能产生自相交或环冲突，
    # 因此必须在目标 CRS 下再次检查和修复，而不能只修复原始 GADM。
    countries_in_grid_crs = clean_gdf(
        countries_in_grid_crs,
        "GADM projected to RWR CRS",
        repair_invalid=True,
    )
    country_geometry = countries_in_grid_crs.set_index(COUNTRY_FIELD).geometry.to_dict()
    reference["reference_inside_attributed_country"] = [
        safe_covers(country, country_geometry[country], point)
        for country, point in zip(reference[COUNTRY_FIELD], reference.geometry)
    ]

    print(f"  矿山数: {len(reference):,}；所属国数: {reference[COUNTRY_FIELD].nunique():,}")
    print(
        "  参考点位于所属国外（仅 QA，不改变属性归属）: "
        f"{int((~reference['reference_inside_attributed_country']).sum()):,}"
    )
    return reference


# ============================================================
# 一次性空间匹配
# ============================================================

def assign_grid_centres_to_countries(grid, countries):
    print("\n[4/8] 将 RWR 网格中心点划归国家（一次性预处理）")
    # 10 km 网格中心点用于定义国家随机背景；不做网格 × 国家精确面积叠加。
    centres = gpd.GeoDataFrame(
        grid[["_grid_pos"]].copy(),
        geometry=grid.geometry.centroid,
        crs=grid.crs,
    )
    centres_in_country_crs = project_to(centres, countries.crs, "RWR grid centres")

    joined = gpd.sjoin(
        centres_in_country_crs,
        countries[[COUNTRY_FIELD, "geometry"]],
        how="left",
        predicate="within",
    )
    joined = (
        joined.reset_index(drop=True)
        .sort_values(["_grid_pos", COUNTRY_FIELD], na_position="last")
        .drop_duplicates("_grid_pos", keep="first")
    )
    assignment = joined[["_grid_pos", COUNTRY_FIELD]].copy()
    assignment["country_match_method"] = np.where(
        assignment[COUNTRY_FIELD].notna(), "centre_within_country", "unmatched"
    )
    assignment["country_nearest_distance_km"] = np.nan

    unmatched_pos = assignment.loc[
        assignment[COUNTRY_FIELD].isna(), "_grid_pos"
    ].to_numpy(dtype=np.int64)
    if len(unmatched_pos):
        print(f"  中心点未直接落入国家面的网格: {len(unmatched_pos):,}，尝试最近国家匹配")
        unmatched_centres = centres.loc[
            centres["_grid_pos"].isin(unmatched_pos)
        ].copy()
        unmatched_equal_area = project_to(
            unmatched_centres, EQUAL_AREA_CRS, "unmatched grid centres"
        )
        countries_equal_area = project_to(countries, EQUAL_AREA_CRS, "GADM")
        nearest = gpd.sjoin_nearest(
            unmatched_equal_area,
            countries_equal_area[[COUNTRY_FIELD, "geometry"]],
            how="left",
            max_distance=GRID_COUNTRY_NEAREST_MAX_KM * 1000.0,
            distance_col="_distance_m",
        )
        nearest = (
            nearest.reset_index(drop=True)
            .sort_values(["_grid_pos", "_distance_m", COUNTRY_FIELD], na_position="last")
            .drop_duplicates("_grid_pos", keep="first")
        )
        nearest_map = nearest.set_index("_grid_pos")
        update_mask = assignment["_grid_pos"].isin(nearest_map.index)
        mapped_country = assignment.loc[update_mask, "_grid_pos"].map(
            nearest_map[COUNTRY_FIELD]
        )
        mapped_distance = assignment.loc[update_mask, "_grid_pos"].map(
            nearest_map["_distance_m"]
        )
        assignment.loc[update_mask, COUNTRY_FIELD] = mapped_country.to_numpy()
        assignment.loc[update_mask, "country_nearest_distance_km"] = (
            mapped_distance.to_numpy(float) / 1000.0
        )
        assignment.loc[
            update_mask & assignment[COUNTRY_FIELD].notna(), "country_match_method"
        ] = "nearest_country_within_limit"

    grid_background = grid[
        ["_grid_pos", GRID_ID_FIELD, RWR_FIELD, "rwr_bin", "is_top10", "rwr_percentile"]
    ].merge(
        assignment,
        on="_grid_pos",
        how="left",
        validate="one_to_one",
    )

    qa = (
        grid_background.assign(
            assigned=grid_background[COUNTRY_FIELD].notna()
        )
        .groupby(["country_match_method", "assigned"], dropna=False)
        .agg(
            n_grid=("_grid_pos", "size"),
            mean_nearest_distance_km=("country_nearest_distance_km", "mean"),
            max_nearest_distance_km=("country_nearest_distance_km", "max"),
        )
        .reset_index()
    )
    qa.to_csv(OUT_GRID_COUNTRY_QA, index=False, encoding="utf-8-sig")

    still_unmatched = int(grid_background[COUNTRY_FIELD].isna().sum())
    print(f"  已划归国家的网格: {len(grid_background) - still_unmatched:,}")
    print(f"  未划归网格（不进入国家随机池）: {still_unmatched:,}")
    return grid_background


def nearest_grid_pairs(points, grid):
    """返回未相交矿点与最近 RWR 网格的位置索引。"""
    try:
        pairs = np.asarray(
            grid.sindex.nearest(points.geometry, return_all=False),
            dtype=np.int64,
        )
        if pairs.ndim != 2 or pairs.shape[0] != 2:
            raise RuntimeError("空间索引 nearest 返回格式异常。")
        point_positions = pairs[0]
        grid_positions = pairs[1]
        return pd.DataFrame({
            "_point_subset_pos": point_positions,
            "_grid_index": grid_positions,
        })
    except Exception:
        # 兼容不支持 SpatialIndex.nearest 的旧版 GeoPandas。
        nearest = gpd.sjoin_nearest(
            points,
            grid[["_grid_pos", "geometry"]],
            how="left",
        ).reset_index().rename(columns={"index": "_point_subset_pos"})
        nearest = (
            nearest.sort_values(["_point_subset_pos", "_grid_pos"])
            .drop_duplicates("_point_subset_pos", keep="first")
        )
        nearest["_grid_index"] = nearest["index_right"].astype(np.int64)
        return nearest[["_point_subset_pos", "_grid_index"]]


def match_mines_to_rwr_grid(mine_points, grid):
    print("\n[5/8] 将铜矿参考点匹配到 RWR 网格")
    attributes = [
        "_grid_pos", GRID_ID_FIELD, RWR_FIELD, "rwr_bin",
        "is_top10", "rwr_percentile", "geometry",
    ]
    joined = gpd.sjoin(
        mine_points,
        grid[attributes],
        how="left",
        predicate="intersects",
    ).reset_index(drop=True)

    # 点恰好落在网格边界时可能匹配多个网格；按 grid_id 稳定选取一个，不按 RWR 高低选择。
    joined = (
        joined.sort_values(["_mine_pos", GRID_ID_FIELD], na_position="last")
        .drop_duplicates("_mine_pos", keep="first")
    )
    joined["grid_match_method"] = np.where(
        joined["_grid_pos"].notna(), "point_intersects_grid", "unmatched"
    )

    unmatched_mines = joined.loc[joined["_grid_pos"].isna(), "_mine_pos"].to_numpy(dtype=np.int64)
    if len(unmatched_mines):
        print(f"  未直接匹配 RWR 网格的矿山: {len(unmatched_mines):,}，匹配最近网格")
        subset = (
            mine_points.loc[mine_points["_mine_pos"].isin(unmatched_mines)]
            .sort_values("_mine_pos")
            .reset_index(drop=True)
        )
        nearest_pairs = nearest_grid_pairs(subset, grid)
        nearest_pairs["_mine_pos"] = nearest_pairs["_point_subset_pos"].map(
            subset["_mine_pos"]
        )
        nearest_grid = nearest_pairs.merge(
            grid[attributes[:-1]].reset_index().rename(columns={"index": "_grid_index"}),
            on="_grid_index",
            how="left",
            validate="many_to_one",
        ).set_index("_mine_pos")

        joined = joined.set_index("_mine_pos")
        for column in attributes[:-1]:
            joined.loc[nearest_grid.index, column] = nearest_grid[column]
        joined.loc[nearest_grid.index, "grid_match_method"] = "nearest_rwr_grid"
        joined = joined.reset_index()

    if joined["_grid_pos"].isna().any():
        raise RuntimeError("仍有矿山无法匹配 RWR 网格，请检查矿山坐标和 RWR 网格范围。")

    # 计算矿点到所匹配网格中心的球面距离，作为 QA 而非筛选条件。
    # 只为实际匹配到的矿山提取 1,278 个网格中心，避免再次生成约 140 万个点。
    selected_pos = joined["_grid_pos"].to_numpy(dtype=np.int64)
    selected_grid_geometry = grid.geometry.iloc[selected_pos].reset_index(drop=True)
    grid_centres_wgs84 = gpd.GeoSeries(
        selected_grid_geometry.centroid, crs=grid.crs
    ).to_crs("EPSG:4326")
    centre_lon = grid_centres_wgs84.x.to_numpy(float)
    centre_lat = grid_centres_wgs84.y.to_numpy(float)
    _, _, distance_m = GEOD.inv(
        joined["reference_lon"].to_numpy(float),
        joined["reference_lat"].to_numpy(float),
        centre_lon,
        centre_lat,
    )
    joined["mine_to_grid_centre_km"] = np.asarray(distance_m, dtype=float) / 1000.0

    mine_observed = joined[[
        MINE_ID_FIELD, COUNTRY_FIELD, "reference_lon", "reference_lat",
        "reference_inside_attributed_country",
        GRID_ID_FIELD, RWR_FIELD, "rwr_bin", "is_top10", "rwr_percentile",
        "grid_match_method", "mine_to_grid_centre_km",
    ]].copy()
    mine_observed.to_csv(OUT_MINE_QA, index=False, encoding="utf-8-sig")

    print(
        f"  直接匹配: {(mine_observed['grid_match_method'] == 'point_intersects_grid').sum():,}；"
        f"最近网格匹配: {(mine_observed['grid_match_method'] == 'nearest_rwr_grid').sum():,}"
    )
    return mine_observed


# ============================================================
# 描述统计、随机化与检验
# ============================================================

def build_country_descriptive(mine_observed, grid_background):
    mine_stats = (
        mine_observed.groupby(COUNTRY_FIELD)
        .agg(
            mine_count=(MINE_ID_FIELD, "nunique"),
            observed_top10_mine_count=("is_top10", "sum"),
            observed_mean_rwr_percentile=("rwr_percentile", "mean"),
            observed_median_rwr_percentile=("rwr_percentile", "median"),
        )
        .reset_index()
    )
    mine_stats["observed_top10_mine_pct"] = (
        mine_stats["observed_top10_mine_count"]
        / mine_stats["mine_count"]
        * 100.0
    )

    background_stats = (
        grid_background.dropna(subset=[COUNTRY_FIELD])
        .groupby(COUNTRY_FIELD)
        .agg(
            background_grid_count=("_grid_pos", "size"),
            background_top10_grid_count=("is_top10", "sum"),
            background_mean_rwr_percentile=("rwr_percentile", "mean"),
            background_median_rwr_percentile=("rwr_percentile", "median"),
        )
        .reset_index()
    )
    background_stats["background_top10_grid_pct"] = (
        background_stats["background_top10_grid_count"]
        / background_stats["background_grid_count"]
        * 100.0
    )

    country = mine_stats.merge(
        background_stats,
        on=COUNTRY_FIELD,
        how="left",
        validate="one_to_one",
    )
    missing_background = country[country["background_grid_count"].isna()]
    if not missing_background.empty:
        names = ", ".join(missing_background[COUNTRY_FIELD].tolist())
        raise ValueError(f"以下矿山所属国没有可用于抽样的 RWR 网格：{names}")

    country["top10_pct_difference_observed_minus_background"] = (
        country["observed_top10_mine_pct"] - country["background_top10_grid_pct"]
    )
    country["top10_enrichment_ratio"] = np.where(
        country["background_top10_grid_pct"] > 0,
        country["observed_top10_mine_pct"] / country["background_top10_grid_pct"],
        np.nan,
    )
    country["mean_percentile_difference_observed_minus_background"] = (
        country["observed_mean_rwr_percentile"]
        - country["background_mean_rwr_percentile"]
    )
    country = country.sort_values(
        ["mine_count", "top10_pct_difference_observed_minus_background"],
        ascending=[False, False],
    ).reset_index(drop=True)
    country.to_csv(OUT_COUNTRY_DESC, index=False, encoding="utf-8-sig")
    return country


def build_observed_table(mine_observed, country):
    global_row = {
        COUNTRY_FIELD: "Global",
        "mine_count": int(mine_observed[MINE_ID_FIELD].nunique()),
        "observed_top10_mine_count": int(mine_observed["is_top10"].sum()),
        "observed_top10_mine_pct": float(mine_observed["is_top10"].mean() * 100.0),
        "observed_mean_rwr_percentile": float(mine_observed["rwr_percentile"].mean()),
        "observed_median_rwr_percentile": float(mine_observed["rwr_percentile"].median()),
        "background_grid_count": int(country["background_grid_count"].sum()),
        "background_top10_grid_count": int(country["background_top10_grid_count"].sum()),
        # 国家分层随机背景的理论期望，按各国真实矿山数加权。
        "background_top10_grid_pct": float(
            np.average(
                country["background_top10_grid_pct"],
                weights=country["mine_count"],
            )
        ),
        "background_mean_rwr_percentile": float(
            np.average(
                country["background_mean_rwr_percentile"],
                weights=country["mine_count"],
            )
        ),
        "background_median_rwr_percentile": np.nan,
    }
    global_row["top10_pct_difference_observed_minus_background"] = (
        global_row["observed_top10_mine_pct"] - global_row["background_top10_grid_pct"]
    )
    global_row["top10_enrichment_ratio"] = (
        global_row["observed_top10_mine_pct"] / global_row["background_top10_grid_pct"]
        if global_row["background_top10_grid_pct"] > 0 else np.nan
    )
    global_row["mean_percentile_difference_observed_minus_background"] = (
        global_row["observed_mean_rwr_percentile"]
        - global_row["background_mean_rwr_percentile"]
    )
    observed = pd.concat(
        [pd.DataFrame([global_row]), country],
        ignore_index=True,
        sort=False,
    )
    observed.to_csv(OUT_OBSERVED, index=False, encoding="utf-8-sig")
    return observed


def build_country_random_pools(grid_background, country):
    pools = {}
    mining_countries = set(country[COUNTRY_FIELD])
    for country_name, group in grid_background[
        grid_background[COUNTRY_FIELD].isin(mining_countries)
    ].groupby(COUNTRY_FIELD):
        pools[country_name] = {
            "is_top10": group["is_top10"].to_numpy(dtype=np.int8),
            "rwr_percentile": group["rwr_percentile"].to_numpy(dtype=np.float64),
        }
    return pools


def run_randomizations(country, pools):
    print(f"\n[6/8] 运行 {ACTIVE_N_ITERATIONS} 次国家分层数组抽样")
    country_counts = [
        (row.Country, int(row.mine_count))
        for row in country.itertuples(index=False)
    ]
    total_mines = int(country["mine_count"].sum())
    rows = []

    for iteration in range(1, ACTIVE_N_ITERATIONS + 1):
        rng = np.random.default_rng(
            np.random.SeedSequence([RANDOM_SEED, iteration])
        )
        top10_count = 0
        percentile_sum = 0.0
        for country_name, mine_count in country_counts:
            pool = pools[country_name]
            sampled_pos = rng.integers(
                0,
                len(pool["is_top10"]),
                size=mine_count,
            )
            top10_count += int(pool["is_top10"][sampled_pos].sum())
            percentile_sum += float(pool["rwr_percentile"][sampled_pos].sum())

        rows.append({
            "iteration": iteration,
            "top10_mine_count": top10_count,
            "pct_top10_mines": top10_count / total_mines * 100.0,
            "mean_rwr_percentile": percentile_sum / total_mines,
        })

    iterations = pd.DataFrame(rows)
    iterations.to_csv(OUT_ITERATIONS, index=False, encoding="utf-8-sig")
    print("  随机抽样完成。")
    return iterations


def empirical_stats(observed, null_values):
    null_values = np.asarray(null_values, dtype=float)
    null_values = null_values[np.isfinite(null_values)]
    if len(null_values) == 0:
        raise ValueError("随机分布为空，无法计算统计量。")

    n = len(null_values)
    mean = float(np.mean(null_values))
    sd = float(np.std(null_values, ddof=1)) if n > 1 else np.nan
    p_upper = float((1 + np.sum(null_values >= observed)) / (n + 1))
    p_lower = float((1 + np.sum(null_values <= observed)) / (n + 1))
    enrichment = float(observed / mean) if mean != 0 else np.nan
    return {
        "observed": float(observed),
        "null_n": n,
        "null_mean": mean,
        "null_sd": sd,
        "null_median": float(np.median(null_values)),
        "null_q025": float(np.quantile(null_values, 0.025)),
        "null_q975": float(np.quantile(null_values, 0.975)),
        "observed_minus_null_mean": float(observed - mean),
        "enrichment_ratio": enrichment,
        "ses": (
            float((observed - mean) / sd)
            if np.isfinite(sd) and sd > 0
            else np.nan
        ),
        "p_upper_one_sided": p_upper,
        "p_lower_one_sided": p_lower,
        "p_two_sided": float(min(1.0, 2.0 * min(p_upper, p_lower))),
    }


def build_test_summary(mine_observed, iterations):
    metric_definitions = [
        (
            "pct_top10_mines",
            "primary",
            float(mine_observed["is_top10"].mean() * 100.0),
            "Percentage of mine reference points in global top-10% RWR grids",
        ),
        (
            "mean_rwr_percentile",
            "continuous_sensitivity",
            float(mine_observed["rwr_percentile"].mean()),
            "Mean global RWR percentile of mine reference-point grids",
        ),
    ]
    rows = []
    for metric, role, observed_value, description in metric_definitions:
        row = {
            "metric": metric,
            "metric_role": role,
            "metric_description": description,
        }
        row.update(empirical_stats(observed_value, iterations[metric]))
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, index=False, encoding="utf-8-sig")
    return summary


# ============================================================
# 补充材料图
# ============================================================









# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 96)
    print("国家分层的 10 km RWR 网格位置对照：真实铜矿是否非随机集中于高 RWR 景观")
    print("=" * 96)

    countries = prepare_countries()
    grid, thresholds = prepare_rwr_grid()
    mine_points = prepare_mines(countries, grid.crs)
    grid_background = assign_grid_centres_to_countries(grid, countries)
    mine_observed = match_mines_to_rwr_grid(mine_points, grid)

    country = build_country_descriptive(mine_observed, grid_background)
    observed = build_observed_table(mine_observed, country)
    pools = build_country_random_pools(grid_background, country)
    iterations = run_randomizations(country, pools)
    summary = build_test_summary(mine_observed, iterations)

    config = pd.DataFrame([{
        "analysis_name": "country-stratified 10-km RWR grid-point control",
        "n_iterations": ACTIVE_N_ITERATIONS,
        "random_seed": RANDOM_SEED,
        "percentile_base": PERCENTILE_BASE,
        "rwr_p90_threshold": thresholds["p90"],
        "test_mode": TEST_MODE,
        "n_mines": len(mine_observed),
        "n_mine_countries": mine_observed[COUNTRY_FIELD].nunique(),
        "randomization_unit": "RWR grid sampled within attributed mine country",
        "sampling_with_replacement": True,
        "mine_location_proxy": "centroid of 30-km buffer calculated in ESRI:54009",
        "primary_metric": "pct_top10_mines",
        "continuous_sensitivity_metric": "mean_rwr_percentile",
        "country_inference": "descriptive only; no country P values or FDR",
        "country_min_mines_for_figure": COUNTRY_MIN_MINES_FOR_FIG,
        "grid_country_nearest_max_km": GRID_COUNTRY_NEAREST_MAX_KM,
    }])
    config.to_csv(OUT_CONFIG, index=False, encoding="utf-8-sig")


    global_observed = observed.loc[observed[COUNTRY_FIELD] == "Global"].iloc[0]
    primary = summary.set_index("metric").loc["pct_top10_mines"]
    print("\n" + "-" * 96)
    print(
        f"真实矿山位于前 10% RWR 网格的比例: "
        f"{global_observed['observed_top10_mine_pct']:.3f}%"
    )
    print(
        f"国家分层随机均值及 95% 区间: "
        f"{primary['null_mean']:.3f}% "
        f"[{primary['null_q025']:.3f}%, {primary['null_q975']:.3f}%]"
    )
    print(
        f"富集倍数: {primary['enrichment_ratio']:.3f}；"
        f"SES: {primary['ses']:.3f}；"
        f"单尾经验 P: {primary['p_upper_one_sided']:.4f}"
    )
    print("-" * 96)
    print(f"Observed: {OUT_OBSERVED}")
    print(f"Iterations: {OUT_ITERATIONS}")
    print(f"Summary: {OUT_SUMMARY}")
    print(f"Country descriptive: {OUT_COUNTRY_DESC}")
    print(f"Mine QA: {OUT_MINE_QA}")
    if CREATE_FIGURES:
        print(f"Combined figure: {OUT_COMBINED_FIG_PNG}")
    if TEST_MODE:
        print("\n[WARNING] 当前为 TEST_MODE，结果仅用于验证代码，不可用于论文。")


if __name__ == "__main__":
    main()
