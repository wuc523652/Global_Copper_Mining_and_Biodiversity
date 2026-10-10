# -*- coding: utf-8 -*-
"""
10_prepare_mine_ownership.py

研究内容四——矿区主表与所有权账本准备。

本脚本只整理矿区、所有者、股权和国家，不进行生态面积计算。核心规则：
1. 以 Property_I / PROP_ID 作为稳定矿区主键，禁止用矿名模糊匹配；
2. 保留全部 1,278 座矿区，不按开发阶段筛选；
3. Copper_mine_all_info.xlsx 是主要所有权快照；仅当 PROP_ID 在该表完全不存在时，
   才回退到 Global_Copper_Mine-regionrenew.xlsx 的 List_of_Ow；
4. 没有明确百分比的所有者不平均分配股权；不足 100% 的部分保留为未知；
5. 所有者国家优先使用 Copper_mine_ownerlist-countryrenew.xlsx；
6. 任何明确股权合计超过 100% 的矿区都会停止运行并输出 QA。

主要输出：
- Mine_Master_30km.gpkg / .csv
- Mine_Owner_Share_Ledger.csv
- Mine_Ownership_Coverage.csv
- Ownership_Source_Comparison.csv
- Ownership_Country_Unmatched.csv
- Ownership_Parsing_QA.csv
- Mine_ID_Crosswalk_QA.csv
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import re
import unicodedata
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import unary_union

try:
    from shapely.validation import make_valid
except Exception:  # Shapely < 2
    def make_valid(geom):
        return geom.buffer(0)


warnings.filterwarnings("ignore")


# =============================================================================
# 路径与参数
# =============================================================================

BASE_DIR = Path(_nc_project(""))
EXCEL_DIR = Path(_nc_data("Global_Copper_Mine/Excel"))

MINE_BUFFER_FILE = Path(
    _nc_data("Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp")
)
GADM_FILE = Path(_nc_data("Global_GADM/GC_GADM.shp"))

MASTER_XLSX = EXCEL_DIR / "Global_Copper_Mine-regionrenew.xlsx"
ALL_INFO_XLSX = EXCEL_DIR / "Copper_mine_all_info.xlsx"
OWNER_COUNTRY_XLSX = EXCEL_DIR / "Copper_mine_ownerlist-countryrenew.xlsx"

OUTPUT_DIR = BASE_DIR / "10_mine_ownership"
OUT_MASTER_GPKG = OUTPUT_DIR / "Mine_Master_30km.gpkg"
OUT_MASTER_CSV = OUTPUT_DIR / "Mine_Master_30km.csv"
OUT_LEDGER = OUTPUT_DIR / "Mine_Owner_Share_Ledger.csv"
OUT_COVERAGE = OUTPUT_DIR / "Mine_Ownership_Coverage.csv"
OUT_SOURCE_COMPARISON = OUTPUT_DIR / "Ownership_Source_Comparison.csv"
OUT_COUNTRY_UNMATCHED = OUTPUT_DIR / "Ownership_Country_Unmatched.csv"
OUT_PARSING_QA = OUTPUT_DIR / "Ownership_Parsing_QA.csv"
OUT_ID_QA = OUTPUT_DIR / "Mine_ID_Crosswalk_QA.csv"

EQUAL_AREA_CRS = "ESRI:54009"
OUTPUT_CRS = "EPSG:4326"
SHARE_TOLERANCE = 1e-6

PROPERTY_ID_CANDIDATES = [
    "Property_I", "Property_ID", "PROPERTY_ID", "PROP_ID", "property_id"
]
MINE_NAME_CANDIDATES = [
    "Property_N", "Property_Name", "PROP_NAME", "mine_name", "MineName", "NAME"
]
MINE_COUNTRY_CANDIDATES = [
    "Country", "COUNTRY", "country", "COUNTRY_NAME", "Host_Country", "Mine_Country"
]
GADM_NAME_CANDIDATES = ["NAME_0", "Country", "COUNTRY", "ADMIN", "SOVEREIGNT", "NAME"]
GADM_ISO_CANDIDATES = ["GID_0", "ISO3", "ADM0_A3", "SOV_A3"]

RECOGNIZED_OWNER_ROLES = {
    "owner", "venturer", "venture", "optionee", "optionor", "fractional", "carried"
}


# 常见名称差异直接映射到 ISO3；其余名称从 GADM 读取。
COUNTRY_ISO3_OVERRIDES = {
    "america": "USA",
    "united states of america": "USA",
    "u s a": "USA",
    "england": "GBR",
    "uk": "GBR",
    "russian federation": "RUS",
    "south korea": "KOR",
    "republic of korea": "KOR",
    "north korea": "PRK",
    "democratic peoples republic of korea": "PRK",
    "laos": "LAO",
    "lao peoples democratic republic": "LAO",
    "iran islamic republic of": "IRN",
    "bolivia plurinational state of": "BOL",
    "venezuela bolivarian republic of": "VEN",
    "tanzania united republic of": "TZA",
    "viet nam": "VNM",
    "vietnam": "VNM",
    "macedonia": "MKD",
    "north macedonia": "MKD",
    "democratic republic of congo": "COD",
    "democratic republic of the congo": "COD",
    "dr congo": "COD",
    "republic of congo": "COG",
    "republic of the congo": "COG",
    "ivory coast": "CIV",
    "cote d ivoire": "CIV",
    "hong kong sar china": "HKG",
    "taiwan province of china": "TWN",
}


# 仅用于结果汇总和流向图。跨洲国家沿用研究内容一至三的口径。
CONTINENT_ISO3 = {
    "Africa": set("DZA AGO BEN BWA BFA BDI CPV CMR CAF TCD COM COG COD CIV DJI EGY GNQ ERI SWZ ETH GAB GMB GHA GIN GNB KEN LSO LBR LBY MDG MWI MLI MRT MUS MAR MOZ NAM NER NGA RWA STP SEN SYC SLE SOM ZAF SSD SDN TZA TGO TUN UGA ESH ZMB ZWE REU SHN MYT".split()),
    "Asia": set("AFG ARM AZE BHR BGD BTN BRN KHM CHN CYP GEO HKG IND IDN IRN IRQ ISR JPN JOR KAZ KWT KGZ LAO LBN MYS MDV MNG MMR NPL PRK OMN PAK PSE PHL QAT SAU SGP KOR LKA SYR TWN TJK THA TLS TUR TKM ARE UZB VNM YEM MAC RUS".split()),
    "Europe": set("ALB AND AUT BLR BEL BIH BGR HRV CZE DNK EST FRO FIN FRA DEU GIB GRC GGY VAT HUN ISL IRL IMN ITA JEY XKO LVA LIE LTU LUX MKD MLT MDA MCO MNE NLD NOR POL PRT ROU SMR SRB SVK SVN ESP SWE CHE UKR GBR ALA".split()),
    "North America": set("ATG BHS BRB BLZ CAN CRI CUB DMA DOM SLV GRD GTM HTI HND JAM MEX NIC PAN KNA LCA VCT USA ABW AIA BMU BES VGB CYM CUW GRL GLP MTQ MSR PRI SXM TCA VIR SPM MAF".split()),
    "South America": set("ARG BOL BRA CHL COL ECU GUY PRY PER SUR URY VEN GUF FLK".split()),
    "Oceania": set("AUS FJI KIR MHL FSM NRU NZL PLW PNG WSM SLB TON TUV VUT ASM COK PYF GUM NCL NIU NFK MNP WLF".split()),
    "Antarctica": set("ATA BVT HMD SGS ATF".split()),
}


# =============================================================================
# 通用工具
# =============================================================================

def read_vector(path: Path) -> gpd.GeoDataFrame:
    if not path.exists():
        raise FileNotFoundError(f"输入空间文件不存在：{path}")
    try:
        return gpd.read_file(path, engine="pyogrio")
    except Exception:
        return gpd.read_file(path)


def first_existing(columns, candidates, label, required=True):
    for field in candidates:
        if field in columns:
            return field
    if required:
        raise ValueError(f"无法识别{label}字段。候选={candidates}；实际字段={list(columns)}")
    return None


def normalize_text(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKD", str(value).strip())
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = text.casefold().replace("&", " and ")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_property_id(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return re.sub(r"\.0$", "", text)


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


def continent_from_iso3(iso3: str) -> str:
    for continent, codes in CONTINENT_ISO3.items():
        if iso3 in codes:
            return continent
    return "Unknown"


def build_country_lookup(gadm: gpd.GeoDataFrame):
    name_field = first_existing(gadm.columns, GADM_NAME_CANDIDATES, "GADM国家名")
    iso_field = first_existing(gadm.columns, GADM_ISO_CANDIDATES, "GADM ISO3")
    table = gadm[[name_field, iso_field]].dropna().drop_duplicates().copy()
    table[iso_field] = table[iso_field].astype(str).str.strip().str.upper()
    table[name_field] = table[name_field].astype(str).str.strip()
    name_to_iso = {
        normalize_text(name): iso
        for name, iso in zip(table[name_field], table[iso_field])
        if len(iso) == 3
    }
    iso_to_name = {
        iso: name for name, iso in zip(table[name_field], table[iso_field]) if len(iso) == 3
    }
    name_to_iso.update(COUNTRY_ISO3_OVERRIDES)
    return name_to_iso, iso_to_name


def resolve_country(raw_country, name_to_iso, iso_to_name):
    raw = "" if raw_country is None or pd.isna(raw_country) else str(raw_country).strip()
    iso3 = name_to_iso.get(normalize_text(raw), "")
    canonical = iso_to_name.get(iso3, raw) if iso3 else raw
    return canonical, iso3, continent_from_iso3(iso3)


def parse_owner_entry(entry: str):
    """解析“公司名 (角色) 51%”；公司名内部括号不会被误删。"""
    raw = str(entry).strip()
    share = np.nan
    prefix = raw
    share_match = re.search(r"(?P<share>\d+(?:\.\d+)?)\s*%\s*$", prefix)
    if share_match:
        share = float(share_match.group("share"))
        prefix = prefix[:share_match.start()].strip()

    role = ""
    role_match = re.search(r"\s+\((?P<role>[^()]*)\)\s*$", prefix)
    if role_match and normalize_text(role_match.group("role")) in RECOGNIZED_OWNER_ROLES:
        role = role_match.group("role").strip()
        owner_name = prefix[:role_match.start()].strip()
    else:
        owner_name = prefix.strip()

    return {
        "owner_entry_raw": raw,
        "owner_name": owner_name,
        "owner_role": role,
        "reported_share_pct": share,
        "has_reported_share": int(pd.notna(share)),
    }


def build_owner_country_lookup(mapping: pd.DataFrame):
    required = {"OWNER_NAME", "OWNER_COUNTRY"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"所有者国家表缺少字段：{sorted(required - set(mapping.columns))}")

    exact = {}
    normalized_candidates = {}
    for _, row in mapping.iterrows():
        name = "" if pd.isna(row["OWNER_NAME"]) else str(row["OWNER_NAME"]).strip()
        country = "" if pd.isna(row["OWNER_COUNTRY"]) else str(row["OWNER_COUNTRY"]).strip()
        if not name:
            continue
        exact[name.casefold()] = country
        normalized_candidates.setdefault(normalize_text(name), set()).add(country)
    normalized = {
        key: next(iter(values))
        for key, values in normalized_candidates.items()
        if len(values) == 1
    }
    ambiguous = {key for key, values in normalized_candidates.items() if len(values) > 1}
    return exact, normalized, ambiguous


def lookup_owner_country(owner_name, exact, normalized, ambiguous):
    key_exact = str(owner_name).strip().casefold()
    if key_exact in exact:
        return exact[key_exact], "exact"
    key_norm = normalize_text(owner_name)
    if key_norm in ambiguous:
        return "", "normalized_ambiguous"
    if key_norm in normalized:
        return normalized[key_norm], "normalized_unique"
    return "", "unmatched"


# =============================================================================
# 数据准备
# =============================================================================

def prepare_mine_master(gadm, name_to_iso, iso_to_name):
    print("\n[1/5] 建立 Property ID 矿区主表")
    master = pd.read_excel(MASTER_XLSX, sheet_name=0)
    required = {"FID", "Property_I", "Property_N", "Country", "List_of_Ow"}
    if not required.issubset(master.columns):
        raise ValueError(f"矿区主表缺少字段：{sorted(required - set(master.columns))}")
    master = master.copy()
    master["source_fid"] = pd.to_numeric(master["FID"], errors="raise").astype(int)
    master["property_id"] = master["Property_I"].map(clean_property_id)
    if master["property_id"].eq("").any() or master["property_id"].duplicated().any():
        raise ValueError("Global_Copper_Mine-regionrenew.xlsx 的 Property_I 存在空值或重复。")

    mines = read_vector(MINE_BUFFER_FILE)
    if mines.crs is None:
        raise ValueError("30 km矿区缓冲区缺少CRS。")
    mines = mines.reset_index(drop=True)
    mines["source_order"] = np.arange(len(mines), dtype=int)
    pid_field = first_existing(mines.columns, PROPERTY_ID_CANDIDATES, "矿区Property ID", required=False)
    shape_name_field = first_existing(
        mines.columns, MINE_NAME_CANDIDATES, "缓冲区矿名", required=False
    )
    shape_country_field = first_existing(
        mines.columns, MINE_COUNTRY_CANDIDATES, "缓冲区国家", required=False
    )
    auxiliary_fields = [field for field in [shape_name_field, shape_country_field] if field]
    auxiliary_rename = {}
    if shape_name_field:
        auxiliary_rename[shape_name_field] = "shape_property_name"
    if shape_country_field:
        auxiliary_rename[shape_country_field] = "shape_country_raw"

    if pid_field:
        geom = mines[[pid_field, "source_order", *auxiliary_fields, "geometry"]].copy()
        geom = geom.rename(columns=auxiliary_rename)
        geom["property_id"] = geom[pid_field].map(clean_property_id)
        geom = geom.drop(columns=[pid_field])
        joined = geom.merge(
            master, on="property_id", how="left", validate="many_to_one", suffixes=("_shape", "_excel")
        )
        id_method = f"shapefile:{pid_field}"
    else:
        fid_field = "FID" if "FID" in mines.columns else None
        geom = mines[["source_order", *auxiliary_fields, "geometry"]].copy()
        geom = geom.rename(columns=auxiliary_rename)
        geom["source_fid"] = (
            pd.to_numeric(mines[fid_field], errors="coerce").astype("Int64")
            if fid_field else geom["source_order"].astype("Int64")
        )
        joined = geom.merge(master, on="source_fid", how="left", validate="many_to_one")
        id_method = "shapefile:FID" if fid_field else "validated_source_order_to_excel_FID"

    if joined["property_id"].isna().any() or joined["property_id"].eq("").any():
        raise ValueError("部分30 km缓冲区无法连接到 Property_I；请查看输入顺序或字段。")
    if shape_name_field:
        joined["shape_name_matches_excel"] = (
            joined["shape_property_name"].map(normalize_text)
            == joined["Property_N"].map(normalize_text)
        ).astype(int)
        match_rate = joined["shape_name_matches_excel"].mean()
        if not pid_field and match_rate < 0.98:
            raise ValueError(
                f"依赖FID/源顺序连接时，缓冲区矿名与Excel矿名匹配率仅{match_rate:.1%}；"
                "不能安全建立Property ID。"
            )
    else:
        joined["shape_name_matches_excel"] = np.nan
    if shape_country_field:
        joined["shape_country_matches_excel"] = (
            joined["shape_country_raw"].map(normalize_text)
            == joined["Country"].map(normalize_text)
        ).astype(int)
    else:
        joined["shape_country_matches_excel"] = np.nan
    if joined["property_id"].duplicated().any():
        # 极少数 multipart 记录按 Property ID 做几何合并。
        attrs = joined.drop(columns="geometry").drop_duplicates("property_id")
        geoms = joined[["property_id", "geometry"]].dissolve(by="property_id", as_index=False)
        joined = attrs.merge(geoms, on="property_id", how="inner", validate="one_to_one")

    gdf = gpd.GeoDataFrame(joined, geometry="geometry", crs=mines.crs).to_crs(EQUAL_AREA_CRS)
    gdf["geometry"] = gdf.geometry.map(clean_polygon)
    if gdf.geometry.isna().any():
        raise ValueError("矿区缓冲区存在无法修复的非面几何。")
    gdf["buffer_area_km2"] = gdf.geometry.area / 1e6

    resolved = gdf["Country"].apply(lambda x: resolve_country(x, name_to_iso, iso_to_name))
    gdf[["resource_country", "resource_iso3", "resource_continent"]] = pd.DataFrame(
        resolved.tolist(), index=gdf.index
    )
    if gdf["resource_iso3"].eq("").any():
        bad = sorted(gdf.loc[gdf["resource_iso3"].eq(""), "Country"].astype(str).unique())
        raise ValueError(f"以下资源拥有国无法匹配ISO3：{bad}")

    gdf["property_name"] = gdf["Property_N"].fillna("").astype(str).str.strip()
    id_qa = gdf[[
        "property_id", "source_fid", "source_order", "Property_N", "Country",
        "shape_name_matches_excel", "shape_country_matches_excel"
    ]].copy()
    id_qa = id_qa.rename(columns={
        "Property_N": "property_name", "Country": "resource_country_raw"
    })
    id_qa["id_connection_method"] = id_method
    id_qa["id_unique"] = (~id_qa["property_id"].duplicated()).astype(int)
    id_qa.to_csv(OUT_ID_QA, index=False, encoding="utf-8-sig")

    keep = [
        "property_id", "source_fid", "source_order", "property_name",
        "resource_country", "resource_iso3", "resource_continent",
        "buffer_area_km2", "List_of_Ow", "geometry",
    ]
    gdf = gdf[keep].sort_values("source_fid").reset_index(drop=True)

    print(f"  矿区：{len(gdf):,}；Property ID：{gdf['property_id'].nunique():,}")
    print(f"  连接方式：{id_method}")
    return gdf


def prepare_ownership_source(mine_master):
    print("\n[2/5] 选择所有权快照")
    info = pd.read_excel(ALL_INFO_XLSX, sheet_name="Sheet1")
    required = {"PROP_ID", "OWNER_LIST"}
    if not required.issubset(info.columns):
        raise ValueError(f"Copper_mine_all_info.xlsx 缺少字段：{sorted(required - set(info.columns))}")
    info = info[["PROP_ID", "OWNER_LIST"]].copy()
    info["property_id"] = info["PROP_ID"].map(clean_property_id)
    info = info[info["property_id"].isin(mine_master["property_id"])].copy()
    if info["property_id"].duplicated().any():
        dup = info.loc[info["property_id"].duplicated(False), "property_id"].tolist()
        raise ValueError(f"主要所有权表中研究矿区PROP_ID重复：{dup[:20]}")

    source = mine_master.drop(columns="geometry").merge(
        info[["property_id", "OWNER_LIST"]], on="property_id", how="left", indicator=True,
        validate="one_to_one"
    )
    source["ownership_source"] = np.where(
        source["_merge"].eq("both"), "Copper_mine_all_info.xlsx", "Global_Copper_Mine-regionrenew.xlsx_fallback"
    )
    # 注意：在主要快照中存在但 OWNER_LIST 为空时，不回填旧快照。
    source["ownership_text"] = np.where(
        source["_merge"].eq("both"), source["OWNER_LIST"], source["List_of_Ow"]
    )

    comparison = source[[
        "property_id", "property_name", "resource_country", "List_of_Ow", "OWNER_LIST", "_merge"
    ]].copy()
    comparison["old_owner_text_normalized"] = comparison["List_of_Ow"].map(normalize_text)
    comparison["new_owner_text_normalized"] = comparison["OWNER_LIST"].map(normalize_text)
    comparison["same_owner_text"] = (
        comparison["old_owner_text_normalized"] == comparison["new_owner_text_normalized"]
    ).astype(int)
    comparison.to_csv(OUT_SOURCE_COMPARISON, index=False, encoding="utf-8-sig")
    print(f"  主要快照匹配：{(source['_merge'] == 'both').sum():,}")
    print(f"  完全不存在并回退旧表：{(source['_merge'] != 'both').sum():,}")
    return source.drop(columns="_merge")


def build_owner_ledger(source, name_to_iso, iso_to_name):
    print("\n[3/5] 解析 OWNER_LIST 并连接投资国")
    mapping = pd.read_excel(OWNER_COUNTRY_XLSX, sheet_name="Sheet1")
    exact, normalized, ambiguous = build_owner_country_lookup(mapping)

    rows = []
    qa_rows = []
    for _, mine in source.iterrows():
        raw_text = "" if pd.isna(mine["ownership_text"]) else str(mine["ownership_text"]).strip()
        if not raw_text:
            qa_rows.append({
                "property_id": mine["property_id"], "property_name": mine["property_name"],
                "qa_type": "blank_ownership_text", "qa_detail": mine["ownership_source"],
            })
            continue
        entries = [item.strip() for item in re.split(r"\s*;\s*", raw_text) if item.strip()]
        for order, entry in enumerate(entries, 1):
            parsed = parse_owner_entry(entry)
            owner_country_raw, match_method = lookup_owner_country(
                parsed["owner_name"], exact, normalized, ambiguous
            )
            investor_country, investor_iso3, investor_continent = resolve_country(
                owner_country_raw, name_to_iso, iso_to_name
            )
            row = {
                "property_id": mine["property_id"],
                "property_name": mine["property_name"],
                "resource_country": mine["resource_country"],
                "resource_iso3": mine["resource_iso3"],
                "resource_continent": mine["resource_continent"],
                "ownership_source": mine["ownership_source"],
                "owner_order": order,
                **parsed,
                "owner_country_raw": owner_country_raw,
                "owner_country_match_method": match_method,
                "investor_country": investor_country,
                "investor_iso3": investor_iso3,
                "investor_continent": investor_continent,
                "investor_country_resolved": int(bool(investor_iso3)),
            }
            rows.append(row)
            if not parsed["owner_name"]:
                qa_rows.append({
                    "property_id": mine["property_id"], "property_name": mine["property_name"],
                    "qa_type": "blank_parsed_owner", "qa_detail": entry,
                })
            if pd.isna(parsed["reported_share_pct"]):
                qa_rows.append({
                    "property_id": mine["property_id"], "property_name": mine["property_name"],
                    "qa_type": "owner_without_explicit_share", "qa_detail": entry,
                })
            if pd.notna(parsed["reported_share_pct"]) and not investor_iso3:
                qa_rows.append({
                    "property_id": mine["property_id"], "property_name": mine["property_name"],
                    "qa_type": "explicit_share_without_investor_country", "qa_detail": entry,
                })

    ledger = pd.DataFrame(rows)
    if ledger.empty:
        raise RuntimeError("没有解析出任何所有者记录。")
    ledger["reported_share_pct"] = pd.to_numeric(ledger["reported_share_pct"], errors="coerce")
    ledger.to_csv(OUT_LEDGER, index=False, encoding="utf-8-sig")

    qa = pd.DataFrame(qa_rows, columns=["property_id", "property_name", "qa_type", "qa_detail"])
    qa.to_csv(OUT_PARSING_QA, index=False, encoding="utf-8-sig")

    unmatched = ledger[
        ledger["owner_country_match_method"].isin(["unmatched", "normalized_ambiguous"])
        | ledger["investor_iso3"].eq("")
    ].copy()
    unmatched.to_csv(OUT_COUNTRY_UNMATCHED, index=False, encoding="utf-8-sig")
    print(f"  所有者记录：{len(ledger):,}")
    print(f"  有明确股权记录：{ledger['reported_share_pct'].notna().sum():,}")
    print(f"  所有者/投资国未解析记录：{len(unmatched):,}")
    return ledger


def build_coverage(source, ledger):
    print("\n[4/5] 核算股权覆盖率")
    explicit = ledger[ledger["reported_share_pct"].notna()].copy()
    explicit["resolved_share_pct"] = np.where(
        explicit["investor_iso3"].ne(""), explicit["reported_share_pct"], 0.0
    )
    grouped = explicit.groupby("property_id", as_index=False).agg(
        num_explicit_owner_records=("owner_name", "size"),
        explicit_reported_share_pct=("reported_share_pct", "sum"),
        explicit_country_resolved_share_pct=("resolved_share_pct", "sum"),
    )
    entry_count = ledger.groupby("property_id", as_index=False).agg(
        num_owner_records=("owner_name", "size")
    )
    coverage = source[[
        "property_id", "property_name", "resource_country", "resource_iso3",
        "resource_continent", "ownership_source", "ownership_text"
    ]].merge(entry_count, on="property_id", how="left").merge(
        grouped, on="property_id", how="left"
    )
    number_cols = [
        "num_owner_records", "num_explicit_owner_records", "explicit_reported_share_pct",
        "explicit_country_resolved_share_pct",
    ]
    coverage[number_cols] = coverage[number_cols].fillna(0.0)
    coverage[["num_owner_records", "num_explicit_owner_records"]] = coverage[
        ["num_owner_records", "num_explicit_owner_records"]
    ].astype(int)
    coverage["explicit_country_unresolved_share_pct"] = (
        coverage["explicit_reported_share_pct"] - coverage["explicit_country_resolved_share_pct"]
    ).clip(lower=0.0)
    coverage["unreported_share_pct"] = (
        100.0 - coverage["explicit_reported_share_pct"]
    ).clip(lower=0.0)
    coverage["total_unattributed_share_pct"] = (
        100.0 - coverage["explicit_country_resolved_share_pct"]
    ).clip(lower=0.0)
    coverage["ownership_text_present"] = coverage["ownership_text"].fillna("").astype(str).str.strip().ne("").astype(int)
    coverage["has_explicit_share"] = coverage["explicit_reported_share_pct"].gt(0).astype(int)
    coverage["all_share_country_resolved"] = (
        coverage["explicit_reported_share_pct"].gt(0)
        & coverage["explicit_country_unresolved_share_pct"].le(SHARE_TOLERANCE)
    ).astype(int)

    excessive = coverage[coverage["explicit_reported_share_pct"] > 100.0 + SHARE_TOLERANCE]
    if not excessive.empty:
        excessive.to_csv(OUTPUT_DIR / "Ownership_Share_Over100_ERROR.csv", index=False, encoding="utf-8-sig")
        raise ValueError(
            "存在明确股权合计超过100%的矿区；已输出 Ownership_Share_Over100_ERROR.csv。"
        )
    coverage.to_csv(OUT_COVERAGE, index=False, encoding="utf-8-sig")
    print(f"  有明确股权矿区：{coverage['has_explicit_share'].sum():,}")
    print(f"  明确股权全部可解析投资国：{coverage['all_share_country_resolved'].sum():,}")
    return coverage


def save_mine_master(mine_master, coverage):
    print("\n[5/5] 保存矿区主表")
    attrs = coverage.drop(columns=["ownership_text"]).copy()
    out = mine_master.drop(columns=["List_of_Ow"]).merge(
        attrs, on=[
            "property_id", "property_name", "resource_country", "resource_iso3", "resource_continent"
        ], how="left", validate="one_to_one"
    )
    csv = out.drop(columns="geometry").copy()
    csv.to_csv(OUT_MASTER_CSV, index=False, encoding="utf-8-sig")

    gpkg = out.to_crs(OUTPUT_CRS)
    gpkg["geometry"] = gpkg.geometry.map(clean_polygon)
    for col in gpkg.select_dtypes(include="object").columns:
        gpkg[col] = gpkg[col].fillna("").astype(str)
    if OUT_MASTER_GPKG.exists():
        OUT_MASTER_GPKG.unlink()
    gpkg.to_file(OUT_MASTER_GPKG, layer="mine_master", driver="GPKG")
    print(f"  {OUT_MASTER_GPKG}")
    print(f"  {OUT_MASTER_CSV}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in [MASTER_XLSX, ALL_INFO_XLSX, OWNER_COUNTRY_XLSX, GADM_FILE]:
        if not path.exists():
            raise FileNotFoundError(f"输入文件不存在：{path}")

    gadm = read_vector(GADM_FILE)
    name_to_iso, iso_to_name = build_country_lookup(gadm)
    mine_master = prepare_mine_master(gadm, name_to_iso, iso_to_name)
    source = prepare_ownership_source(mine_master)
    ledger = build_owner_ledger(source, name_to_iso, iso_to_name)
    coverage = build_coverage(source, ledger)
    save_mine_master(mine_master, coverage)

    print("\n完成。下一步运行 11_calculate_mine_exposure.py。")


if __name__ == "__main__":
    main()
