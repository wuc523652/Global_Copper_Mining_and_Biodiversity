# -*- coding: utf-8 -*-
"""
04_calculate_mine_species_overlap.py

目的：
基于“矿区所属国家 / 资源拥有国”口径，统计每个国家名下铜矿影响范围与受威胁物种栖息地的空间冲突规模，
并生成 Country_CopperOverlap_Summary.gpkg。

核心统计口径：
1. 国家归属不再由缓冲区实际落入哪个国家边界决定，而是由铜矿数据属性表中的国家字段决定。
2. 铜矿影响范围不按 GADM 国界裁剪。即使某一国家矿区的缓冲影响范围越过国界，并在其他国家领土内
   与物种栖息地发生空间冲突，该冲突面积、受影响物种和冲突矿区仍全部归属于该矿区属性表记录的所属国家。
3. GADM 国家边界仅用于：
   - 统一国家名称与大洲字段；
   - 作为 country_summary 图层的显示底图；
   - 保留所有国家，并对无矿区或无冲突国家填 0。
4. 输出 Country_CopperOverlap_Summary.gpkg：
   - layer = country_summary：GADM 国家边界 + 按矿区所属国家归属的统计属性。
   - layer = country_overlap_union：按矿区所属国家归属后的冲突区域 union 几何。注意该几何可能跨出
     对应国家边界，因为它表示该国矿区在全球范围内造成的影响区域，而不是该国领土内的冲突区域。
5. 统计每个矿区所属国家：
   - 与物种栖息地存在空间冲突的铜矿区数量 conflict_mine_count
   - 归属到该国的空间冲突总面积 overlap_area_km2，单位 km²，国家归属口径下严格空间去重
   - 涉及物种数量 num_species；同一物种可在不同矿区所属国家重复计数
6. GADM 无大洲字段，因此脚本内置国家-大洲映射，并输出大洲汇总表。

注意：
- 面积计算使用 ESRI:54009 Mollweide 等面积投影。
- 物种数据默认只处理 CR / EN / VU，以及 Amphibians / Birds / Fish / Mammals / Reptiles。
- 当前版本使用 30 km 铜矿影响范围。
- 输入铜矿数据必须包含矿区所属国家字段，默认优先识别 Country / COUNTRY / country 等字段。
  若国家名无法匹配到 GADM，脚本会输出 Unmatched_Mine_Countries.csv 并停止，便于手动补充别名。

ArcGIS Pro 兼容性修正：
- 写出 GPKG 前强制修复 Polygon/MultiPolygon 几何。
- 从 make_valid 可能产生的 GeometryCollection 中只提取面要素。
- 删除 Z/M 维度并再次修复。
- 输出 GPKG 几何默认转为 EPSG:4326；面积字段仍基于 ESRI:54009 计算。
"""

from _paths import project_path as _nc_project, data_path as _nc_data, analysis_script as _nc_script

from pathlib import Path
import warnings
import gc
import time
import unicodedata
import re

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.validation import make_valid
from shapely.ops import unary_union, transform
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection


# ============================================================
# 路径配置
# ============================================================

# 当前版本使用 30 km 铜矿影响范围。
MINE_SHP = Path(
    _nc_data("Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp")
)

# GADM 全球国家边界面数据。该数据只需要包含国家名称字段和 geometry。
GADM_PATH = Path(_nc_data("Global_GADM/GC_GADM.shp"))

# IUCN / BirdLife 物种栖息地数据根目录。
IUCN_ROOT = Path(_nc_data("IUCNredlist/202502"))

# 输出目录
OUT_DIR = Path(_nc_project("04_mine_species_overlap"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 核心输出
OUT_COUNTRY_SUMMARY_GPKG = OUT_DIR / "Country_CopperOverlap_Summary.gpkg"
OUT_COUNTRY_SUMMARY_CSV = OUT_DIR / "Country_CopperOverlap_Summary.csv"
OUT_CONTINENT_SUMMARY = OUT_DIR / "Continent_CopperOverlap_Summary.csv"

# 国家 × 物种、国家 × 矿区辅助表
OUT_BY_SPECIES = OUT_DIR / "Country_CopperOverlap_BySpecies.csv"
OUT_BY_MINE = OUT_DIR / "Country_CopperOverlap_ByMine.csv"

# 全球物种尺度严格空间去重结果，供后续物种优先救助模块使用
OUT_BY_SPECIES_GLOBAL_UNIQUE = OUT_DIR / "Species_CopperOverlap_GlobalUnique.csv"
OUT_BY_SPECIES_GLOBAL_UNIQUE_GEOM = OUT_DIR / "Species_CopperOverlap_GlobalUnique.gpkg"

# 全球总览与类群/等级统计
OUT_GLOBAL_OVERVIEW = OUT_DIR / "Global_CopperOverlap_Overview.csv"
OUT_GLOBAL_BY_TAXON = OUT_DIR / "Global_AffectedSpecies_ByTaxon.csv"
OUT_GLOBAL_BY_STATUS = OUT_DIR / "Global_AffectedSpecies_ByStatus.csv"

# 中间结果文件夹
INTERIM_DIR = OUT_DIR / "interim_mine_country_attribution_mollweide"
INTERIM_DIR.mkdir(parents=True, exist_ok=True)

# 等面积投影：所有面积计算都在该投影中完成
TARGET_CRS = "ESRI:54009"

# GPKG 最终输出几何 CRS：为了提高 ArcGIS Pro 兼容性，最终国家面和冲突面输出为 WGS84。
# 注意：所有面积字段仍然是在 TARGET_CRS 中计算得到，不会因为这里转为 WGS84 而改变。
OUTPUT_GPKG_CRS = "EPSG:4326"

# 物种数据设置
TAXA = ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"]
STATUSES = ["CR", "EN", "VU"]

# 字段候选
ID_FIELD_CANDIDATES = ["ID_NO", "id_no", "Id_no", "ID_No", "SCI_NAME", "sci_name"]
GADM_COUNTRY_FIELD_CANDIDATES = [
    "Country", "COUNTRY", "country", "NAME_0", "NAME", "CNTRY_NAME", "ADMIN", "admin", "SOVEREIGNT"
]
MINE_ID_FIELD_CANDIDATES = [
    "mine_id", "Mine_ID", "MINE_ID", "id", "ID", "OBJECTID", "FID", "FID_1",
    "SNL_ID", "snl_id", "Project_ID", "PROJECT_ID", "project_id",
    "MineName", "MINE_NAME", "mine_name", "Name", "NAME"
]

# 铜矿属性表中的“矿区所属国家 / 资源拥有国 / 项目所在国”字段候选。
# 本脚本的国家归属完全依赖该字段，而不再根据缓冲区实际落入哪个国家边界来分配。
MINE_COUNTRY_FIELD_CANDIDATES = [
    "Country", "COUNTRY", "country",
    "Mine_Country", "MINE_COUNTRY", "mine_country",
    "Host_Country", "HOST_COUNTRY", "host_country",
    "Project_Country", "PROJECT_COUNTRY", "project_country",
    "SNL_Country", "SNL_COUNTRY", "snl_country",
    "CNTRY_NAME", "ADMIN", "admin",
    "COUNTRY_NA", "Country_Na", "COUNTRYNAM",
    "ResCountry", "RESCOUNTRY", "Host_Cou", "HOST_COU"
]

MINE_COUNTRY_RAW_FIELD = "mine_country_raw"

COUNTRY_FIELD = "Country"
CONTINENT_FIELD = "Continent"
GADM_GID_FIELD = "GID_0"
MINE_ID_FIELD = "mine_id"

USE_CACHE = False
warnings.filterwarnings("ignore")


# ============================================================
# 国家-大洲映射
# ============================================================

def normalize_name(value):
    """用于国家名匹配的规范化：小写、去重音、去标点、多空格合并。"""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    s = str(value).strip()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = s.lower()
    s = s.replace("&", "and")
    s = re.sub(r"[\(\)\[\]\{\},.;:'\"`´’‘]", " ", s)
    s = re.sub(r"[-_/]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def build_gid_continent_map():
    """
    基于用户提供的 GC_GADM.xlsx 属性表中 GID_0 对全部国家/地区分配大洲。
    这样即使 GADM.shp 本身没有大洲字段，也不会在大洲统计中出现 Unknown。

    大洲分类口径：
    - Africa / Asia / Europe / North America / South America / Oceania / Antarctica
    - 跨洲国家按本研究矿业统计常用口径处理：Russia、Kazakhstan、Turkey、Azerbaijan、Georgia 归入 Asia。
    - 少量海外领地按其地理位置或通常区域归属处理。
    """
    return {
        "AFG": "Asia", "XAD": "Asia", "ALB": "Europe", "DZA": "Africa", "ASM": "Oceania",
        "AND": "Europe", "AGO": "Africa", "AIA": "North America", "ATA": "Antarctica",
        "ATG": "North America", "ARG": "South America", "ARM": "Asia", "ABW": "North America",
        "AUS": "Oceania", "AUT": "Europe", "AZE": "Asia", "BHS": "North America",
        "BHR": "Asia", "BGD": "Asia", "BRB": "North America", "BLR": "Europe",
        "BEL": "Europe", "BLZ": "North America", "BEN": "Africa", "BMU": "North America",
        "BTN": "Asia", "BOL": "South America", "BES": "North America", "BIH": "Europe",
        "BWA": "Africa", "BVT": "Antarctica", "BRA": "South America", "IOT": "Asia",
        "VGB": "North America", "BRN": "Asia", "BGR": "Europe", "BFA": "Africa",
        "BDI": "Africa", "KHM": "Asia", "CMR": "Africa", "CAN": "North America",
        "CPV": "Africa", "XCA": "Asia", "CYM": "North America", "CAF": "Africa",
        "TCD": "Africa", "CHL": "South America", "CHN": "Asia", "CXR": "Oceania",
        "XCL": "North America", "CCK": "Oceania", "COL": "South America", "COM": "Africa",
        "COK": "Oceania", "CRI": "North America", "HRV": "Europe", "CUB": "North America",
        "CUW": "North America", "CYP": "Asia", "CZE": "Europe", "CIV": "Africa",
        "COD": "Africa", "DNK": "Europe", "DJI": "Africa", "DMA": "North America",
        "DOM": "North America", "ECU": "South America", "EGY": "Africa", "SLV": "North America",
        "GNQ": "Africa", "ERI": "Africa", "EST": "Europe", "ETH": "Africa",
        "FLK": "South America", "FRO": "Europe", "FJI": "Oceania", "FIN": "Europe",
        "FRA": "Europe", "GUF": "South America", "PYF": "Oceania", "ATF": "Antarctica",
        "GAB": "Africa", "GMB": "Africa", "GEO": "Asia", "DEU": "Europe",
        "GHA": "Africa", "GIB": "Europe", "GRC": "Europe", "GRL": "North America",
        "GRD": "North America", "GLP": "North America", "GUM": "Oceania", "GTM": "North America",
        "GGY": "Europe", "GIN": "Africa", "GNB": "Africa", "GUY": "South America",
        "HTI": "North America", "HMD": "Antarctica", "HND": "North America", "HKG": "Asia",
        "HUN": "Europe", "ISL": "Europe", "IND": "Asia", "IDN": "Asia",
        "IRN": "Asia", "IRQ": "Asia", "IRL": "Europe", "IMN": "Europe",
        "ISR": "Asia", "ITA": "Europe", "JAM": "North America", "JPN": "Asia",
        "JEY": "Europe", "JOR": "Asia", "KAZ": "Asia", "KEN": "Africa",
        "KIR": "Oceania", "XKO": "Europe", "KWT": "Asia", "KGZ": "Asia",
        "LAO": "Asia", "LVA": "Europe", "LBN": "Asia", "LSO": "Africa",
        "LBR": "Africa", "LBY": "Africa", "LIE": "Europe", "LTU": "Europe",
        "LUX": "Europe", "MAC": "Asia", "MKD": "Europe", "MDG": "Africa",
        "MWI": "Africa", "MYS": "Asia", "MDV": "Asia", "MLI": "Africa",
        "MLT": "Europe", "MHL": "Oceania", "MTQ": "North America", "MRT": "Africa",
        "MUS": "Africa", "MYT": "Africa", "MEX": "North America", "FSM": "Oceania",
        "MDA": "Europe", "MCO": "Europe", "MNG": "Asia", "MNE": "Europe",
        "MSR": "North America", "MAR": "Africa", "MOZ": "Africa", "MMR": "Asia",
        "NAM": "Africa", "NRU": "Oceania", "NPL": "Asia", "NLD": "Europe",
        "NCL": "Oceania", "NZL": "Oceania", "NIC": "North America", "NER": "Africa",
        "NGA": "Africa", "NIU": "Oceania", "NFK": "Oceania", "PRK": "Asia",
        "XNC": "Asia", "MNP": "Oceania", "NOR": "Europe", "OMN": "Asia",
        "PAK": "Asia", "PLW": "Oceania", "PSE": "Asia", "PAN": "North America",
        "PNG": "Oceania", "XPI": "Asia", "PRY": "South America", "PER": "South America",
        "PHL": "Asia", "PCN": "Oceania", "POL": "Europe", "PRT": "Europe",
        "PRI": "North America", "QAT": "Asia", "COG": "Africa", "REU": "Africa",
        "ROU": "Europe", "RUS": "Asia", "RWA": "Africa", "SHN": "Africa",
        "KNA": "North America", "LCA": "North America", "SPM": "North America", "VCT": "North America",
        "BLM": "North America", "MAF": "North America", "WSM": "Oceania", "SMR": "Europe",
        "SAU": "Asia", "SEN": "Africa", "SRB": "Europe", "SYC": "Africa",
        "SLE": "Africa", "SGP": "Asia", "SXM": "North America", "SVK": "Europe",
        "SVN": "Europe", "SLB": "Oceania", "SOM": "Africa", "ZAF": "Africa",
        "SGS": "Antarctica", "KOR": "Asia", "SSD": "Africa", "ESP": "Europe",
        "XSP": "Asia", "LKA": "Asia", "SDN": "Africa", "SUR": "South America",
        "SJM": "Europe", "SWZ": "Africa", "SWE": "Europe", "CHE": "Europe",
        "SYR": "Asia", "STP": "Africa", "TWN": "Asia", "TJK": "Asia",
        "TZA": "Africa", "THA": "Asia", "TLS": "Asia", "TGO": "Africa",
        "TKL": "Oceania", "TON": "Oceania", "TTO": "North America", "TUN": "Africa",
        "TUR": "Asia", "TKM": "Asia", "TCA": "North America", "TUV": "Oceania",
        "UGA": "Africa", "UKR": "Europe", "ARE": "Asia", "GBR": "Europe",
        "USA": "North America", "UMI": "Oceania", "URY": "South America", "UZB": "Asia",
        "VUT": "Oceania", "VAT": "Europe", "VEN": "South America", "VNM": "Asia",
        "VIR": "North America", "WLF": "Oceania", "ESH": "Africa", "YEM": "Asia",
        "ZMB": "Africa", "ZWE": "Africa", "ALA": "Europe",
    }


def build_country_continent_map():
    """
    国家名-大洲映射，作为 GID_0 不存在时的备用方案。
    这些名称来自用户提供的 GC_GADM.xlsx 属性表中的 NAME_0。
    """
    gid_map = build_gid_continent_map()
    name_to_gid = {
        "Afghanistan": "AFG", "Akrotiri and Dhekelia": "XAD", "Albania": "ALB", "Algeria": "DZA",
        "American Samoa": "ASM", "Andorra": "AND", "Angola": "AGO", "Anguilla": "AIA",
        "Antarctica": "ATA", "Antigua and Barbuda": "ATG", "Argentina": "ARG", "Armenia": "ARM",
        "Aruba": "ABW", "Australia": "AUS", "Austria": "AUT", "Azerbaijan": "AZE",
        "Bahamas": "BHS", "Bahrain": "BHR", "Bangladesh": "BGD", "Barbados": "BRB",
        "Belarus": "BLR", "Belgium": "BEL", "Belize": "BLZ", "Benin": "BEN", "Bermuda": "BMU",
        "Bhutan": "BTN", "Bolivia": "BOL", "Bonaire, Sint Eustatius and Saba": "BES",
        "Bosnia and Herzegovina": "BIH", "Botswana": "BWA", "Bouvet Island": "BVT", "Brazil": "BRA",
        "British Indian Ocean Territory": "IOT", "British Virgin Islands": "VGB", "Brunei": "BRN", "Bulgaria": "BGR",
        "Burkina Faso": "BFA", "Burundi": "BDI", "Cambodia": "KHM", "Cameroon": "CMR", "Canada": "CAN",
        "Cape Verde": "CPV", "Caspian Sea": "XCA", "Cayman Islands": "CYM", "Central African Republic": "CAF",
        "Chad": "TCD", "Chile": "CHL", "China": "CHN", "Christmas Island": "CXR",
        "Clipperton Island": "XCL", "Cocos Islands": "CCK", "Colombia": "COL", "Comoros": "COM",
        "Cook Islands": "COK", "Costa Rica": "CRI", "Croatia": "HRV", "Cuba": "CUB", "Curaçao": "CUW",
        "Cyprus": "CYP", "Czech Republic": "CZE", "Côte d'Ivoire": "CIV", "Democratic Republic of the Congo": "COD",
        "Denmark": "DNK", "Djibouti": "DJI", "Dominica": "DMA", "Dominican Republic": "DOM",
        "Ecuador": "ECU", "Egypt": "EGY", "El Salvador": "SLV", "Equatorial Guinea": "GNQ",
        "Eritrea": "ERI", "Estonia": "EST", "Ethiopia": "ETH", "Falkland Islands": "FLK", "Faroe Islands": "FRO",
        "Fiji": "FJI", "Finland": "FIN", "France": "FRA", "French Guiana": "GUF", "French Polynesia": "PYF",
        "French Southern Territories": "ATF", "Gabon": "GAB", "Gambia": "GMB", "Georgia": "GEO", "Germany": "DEU",
        "Ghana": "GHA", "Gibraltar": "GIB", "Greece": "GRC", "Greenland": "GRL", "Grenada": "GRD",
        "Guadeloupe": "GLP", "Guam": "GUM", "Guatemala": "GTM", "Guernsey": "GGY", "Guinea": "GIN",
        "Guinea-Bissau": "GNB", "Guyana": "GUY", "Haiti": "HTI", "Heard Island and McDonald Islands": "HMD",
        "Honduras": "HND", "Hong Kong": "HKG", "Hungary": "HUN", "Iceland": "ISL", "India": "IND",
        "Indonesia": "IDN", "Iran": "IRN", "Iraq": "IRQ", "Ireland": "IRL", "Isle of Man": "IMN",
        "Israel": "ISR", "Italy": "ITA", "Jamaica": "JAM", "Japan": "JPN", "Jersey": "JEY",
        "Jordan": "JOR", "Kazakhstan": "KAZ", "Kenya": "KEN", "Kiribati": "KIR", "Kosovo": "XKO",
        "Kuwait": "KWT", "Kyrgyzstan": "KGZ", "Laos": "LAO", "Latvia": "LVA", "Lebanon": "LBN",
        "Lesotho": "LSO", "Liberia": "LBR", "Libya": "LBY", "Liechtenstein": "LIE", "Lithuania": "LTU",
        "Luxembourg": "LUX", "Macao": "MAC", "Macedonia": "MKD", "Madagascar": "MDG", "Malawi": "MWI",
        "Malaysia": "MYS", "Maldives": "MDV", "Mali": "MLI", "Malta": "MLT", "Marshall Islands": "MHL",
        "Martinique": "MTQ", "Mauritania": "MRT", "Mauritius": "MUS", "Mayotte": "MYT", "Mexico": "MEX",
        "Micronesia": "FSM", "Moldova": "MDA", "Monaco": "MCO", "Mongolia": "MNG", "Montenegro": "MNE",
        "Montserrat": "MSR", "Morocco": "MAR", "Mozambique": "MOZ", "Myanmar": "MMR", "Namibia": "NAM",
        "Nauru": "NRU", "Nepal": "NPL", "Netherlands": "NLD", "New Caledonia": "NCL", "New Zealand": "NZL",
        "Nicaragua": "NIC", "Niger": "NER", "Nigeria": "NGA", "Niue": "NIU", "Norfolk Island": "NFK",
        "North Korea": "PRK", "Northern Cyprus": "XNC", "Northern Mariana Islands": "MNP", "Norway": "NOR",
        "Oman": "OMN", "Pakistan": "PAK", "Palau": "PLW", "Palestina": "PSE", "Panama": "PAN",
        "Papua New Guinea": "PNG", "Paracel Islands": "XPI", "Paraguay": "PRY", "Peru": "PER", "Philippines": "PHL",
        "Pitcairn Islands": "PCN", "Poland": "POL", "Portugal": "PRT", "Puerto Rico": "PRI", "Qatar": "QAT",
        "Republic of Congo": "COG", "Reunion": "REU", "Romania": "ROU", "Russia": "RUS", "Rwanda": "RWA",
        "Saint Helena": "SHN", "Saint Kitts and Nevis": "KNA", "Saint Lucia": "LCA", "Saint Pierre and Miquelon": "SPM",
        "Saint Vincent and the Grenadines": "VCT", "Saint-Barthélemy": "BLM", "Saint-Martin": "MAF",
        "Samoa": "WSM", "San Marino": "SMR", "Saudi Arabia": "SAU", "Senegal": "SEN", "Serbia": "SRB",
        "Seychelles": "SYC", "Sierra Leone": "SLE", "Singapore": "SGP", "Sint Maarten": "SXM", "Slovakia": "SVK",
        "Slovenia": "SVN", "Solomon Islands": "SLB", "Somalia": "SOM", "South Africa": "ZAF",
        "South Georgia and the South Sandwich Islands": "SGS", "South Korea": "KOR", "South Sudan": "SSD",
        "Spain": "ESP", "Spratly Islands": "XSP", "Sri Lanka": "LKA", "Sudan": "SDN", "Suriname": "SUR",
        "Svalbard and Jan Mayen": "SJM", "Swaziland": "SWZ", "Sweden": "SWE", "Switzerland": "CHE",
        "Syria": "SYR", "São Tomé and Príncipe": "STP", "Taiwan": "TWN", "Tajikistan": "TJK",
        "Tanzania": "TZA", "Thailand": "THA", "Timor-Leste": "TLS", "Togo": "TGO", "Tokelau": "TKL",
        "Tonga": "TON", "Trinidad and Tobago": "TTO", "Tunisia": "TUN", "Turkey": "TUR", "Turkmenistan": "TKM",
        "Turks and Caicos Islands": "TCA", "Tuvalu": "TUV", "Uganda": "UGA", "Ukraine": "UKR",
        "United Arab Emirates": "ARE", "United Kingdom": "GBR", "United States": "USA", "United States Minor Outlying Islands": "UMI",
        "Uruguay": "URY", "Uzbekistan": "UZB", "Vanuatu": "VUT", "Vatican City": "VAT",
        "Venezuela": "VEN", "Vietnam": "VNM", "Virgin Islands, U.S.": "VIR", "Wallis and Futuna": "WLF",
        "Western Sahara": "ESH", "Yemen": "YEM", "Zambia": "ZMB", "Zimbabwe": "ZWE", "Åland": "ALA",
    }
    mapping = {normalize_name(name): gid_map[gid] for name, gid in name_to_gid.items()}
    aliases = {
        "cote d ivoire": "Africa", "ivory coast": "Africa", "dr congo": "Africa", "drc": "Africa",
        "democratic republic of congo": "Africa", "congo democratic republic": "Africa", "congo republic": "Africa",
        "sao tome and principe": "Africa", "sao tome principe": "Africa", "reunion": "Africa",
        "russian federation": "Asia", "viet nam": "Asia", "lao pdr": "Asia", "state of palestine": "Asia",
        "palestine": "Asia", "east timor": "Asia", "czechia": "Europe", "north macedonia": "Europe",
        "aland islands": "Europe", "curacao": "North America", "saint barthelemy": "North America",
        "saint martin": "North America", "united states of america": "North America", "us virgin islands": "North America",
        "u s virgin islands": "North America", "falkland islands malvinas": "South America",
        "micronesia federated states of": "Oceania", "pitcairn": "Oceania",
    }
    mapping.update({normalize_name(k): v for k, v in aliases.items()})
    return mapping


GID_CONTINENT_MAP = build_gid_continent_map()
COUNTRY_CONTINENT_MAP = build_country_continent_map()


def infer_continent(country_name, gid=None):
    """优先根据 GID_0 判断大洲；如果没有 GID_0，则退回国家名称匹配。"""
    if gid is not None and not (isinstance(gid, float) and np.isnan(gid)):
        gid_key = str(gid).strip().upper()
        if gid_key in GID_CONTINENT_MAP:
            return GID_CONTINENT_MAP[gid_key]
    key = normalize_name(country_name)
    return COUNTRY_CONTINENT_MAP.get(key, "Unknown")


def build_gadm_country_lookup(countries):
    """
    构建矿区 Country 字段到 GADM 标准国家名的匹配表。
    先进行规范化匹配，再补充常见国家名称别名，避免因为 USA/Russian Federation/Viet Nam 等写法
    导致矿区无法归属到 GADM 国家。
    """
    gadm_names = countries[COUNTRY_FIELD].dropna().astype(str).str.strip().unique().tolist()
    gadm_set = set(gadm_names)

    norm_to_gadm = {}
    for name in gadm_names:
        key = normalize_name(name)
        if key and key not in norm_to_gadm:
            norm_to_gadm[key] = name

    aliases = {
        "usa": "United States",
        "u s a": "United States",
        "us": "United States",
        "u s": "United States",
        "united states of america": "United States",
        "russian federation": "Russia",
        "russia federation": "Russia",
        "viet nam": "Vietnam",
        "lao pdr": "Laos",
        "lao people s democratic republic": "Laos",
        "lao people's democratic republic": "Laos",
        "bolivia plurinational state of": "Bolivia",
        "venezuela bolivarian republic of": "Venezuela",
        "iran islamic republic of": "Iran",
        "iran islamic rep": "Iran",
        "tanzania united republic of": "Tanzania",
        "united republic of tanzania": "Tanzania",
        "democratic republic of congo": "Democratic Republic of the Congo",
        "democratic republic of the congo": "Democratic Republic of the Congo",
        "dr congo": "Democratic Republic of the Congo",
        "drc": "Democratic Republic of the Congo",
        "congo kinshasa": "Democratic Republic of the Congo",
        "congo democratic republic": "Democratic Republic of the Congo",
        "dem rep congo": "Democratic Republic of the Congo",
        "republic of congo": "Republic of Congo",
        "congo brazzaville": "Republic of Congo",
        "cote d ivoire": "Côte d'Ivoire",
        "côte d ivoire": "Côte d'Ivoire",
        "ivory coast": "Côte d'Ivoire",
        "turkiye": "Turkey",
        "türkiye": "Turkey",
        "czechia": "Czech Republic",
        "north macedonia": "Macedonia",
        "eswatini": "Swaziland",
        "brunei darussalam": "Brunei",
        "korea republic of": "South Korea",
        "republic of korea": "South Korea",
        "south korea": "South Korea",
        "korea dem people s republic": "North Korea",
        "democratic people s republic of korea": "North Korea",
        "north korea": "North Korea",
        "moldova republic of": "Moldova",
        "syrian arab republic": "Syria",
        "cape verde": "Cape Verde",
        "cabo verde": "Cape Verde",
        "slovak republic": "Slovakia",
        "myanmar burma": "Myanmar",
        "burma": "Myanmar",
        "palestine": "Palestina",
        "state of palestine": "Palestina",
        "timor leste": "Timor-Leste",
        "east timor": "Timor-Leste",
        "micronesia federated states of": "Micronesia",
        "federated states of micronesia": "Micronesia",
        "sao tome and principe": "São Tomé and Príncipe",
        "sao tome principe": "São Tomé and Príncipe",
        "saint barthelemy": "Saint-Barthélemy",
        "saint martin": "Saint-Martin",
        "curacao": "Curaçao",
        "aland": "Åland",
        "aland islands": "Åland",
    }

    for alias, target in aliases.items():
        if target in gadm_set:
            norm_to_gadm[normalize_name(alias)] = target

    return norm_to_gadm


def standardize_mine_country_names(mine, countries, mine_country_field):
    """
    将铜矿属性表中的国家名统一到 GADM 国家名。
    如果存在无法匹配的国家名，保存清单并停止，防止后续统计发生隐性漏算。
    """
    country_lookup = build_gadm_country_lookup(countries)
    country_cont_map = countries.set_index(COUNTRY_FIELD)[CONTINENT_FIELD].to_dict()

    raw = mine[mine_country_field].astype(str).str.strip()
    raw = raw.replace({"nan": "", "None": "", "NULL": "", "NaN": ""})

    standardized = []
    unmatched = []

    gadm_set = set(countries[COUNTRY_FIELD].astype(str))
    for val in raw:
        if val is None or str(val).strip() == "":
            standardized.append("")
            unmatched.append(val)
            continue

        val_str = str(val).strip()
        key = normalize_name(val_str)
        gadm_country = country_lookup.get(key, "")

        # exact fallback
        if gadm_country == "" and val_str in gadm_set:
            gadm_country = val_str

        standardized.append(gadm_country)
        if gadm_country == "":
            unmatched.append(val_str)

    mine[COUNTRY_FIELD] = standardized
    mine[MINE_COUNTRY_RAW_FIELD] = raw

    unmatched_series = pd.Series(unmatched, dtype="object")
    unmatched_series = unmatched_series[unmatched_series.astype(str).str.strip().str.len() > 0]
    if len(unmatched_series) > 0:
        unmatched_df = (
            unmatched_series.value_counts()
            .rename_axis(MINE_COUNTRY_RAW_FIELD)
            .reset_index(name="mine_count")
        )
        unmatched_path = OUT_DIR / "Unmatched_Mine_Countries.csv"
        unmatched_df.to_csv(unmatched_path, index=False, encoding="utf-8-sig")
        raise ValueError(
            f"有 {len(unmatched_df)} 个矿区国家名称无法匹配到 GADM，已保存清单: {unmatched_path}\n"
            f"请在 build_gadm_country_lookup() 的 aliases 中补充这些国家别名后重新运行。\n"
            f"未匹配示例: {unmatched_df.head(20).to_dict(orient='records')}"
        )

    mine[CONTINENT_FIELD] = mine[COUNTRY_FIELD].map(country_cont_map).fillna("Unknown")
    return mine


# ============================================================
# 工具函数
# ============================================================

def read_vector(path):
    """读取矢量文件，优先使用 pyogrio。"""
    try:
        return gpd.read_file(path, engine="pyogrio")
    except Exception:
        return gpd.read_file(path)


def detect_field(columns, candidates, label):
    """从字段候选列表中识别字段。"""
    for f in candidates:
        if f in columns:
            return f
    raise ValueError(
        f"未找到 {label} 字段。候选字段为 {candidates}，实际字段为：{list(columns)}"
    )


def detect_id_field(columns):
    """从 IUCN shapefile 字段列表中选择物种识别字段。"""
    return detect_field(columns, ID_FIELD_CANDIDATES, "物种识别")


def ensure_crs_and_project(gdf, target_crs, name="vector"):
    """确保矢量数据具有 CRS，并投影到目标 CRS。"""
    if gdf.crs is None:
        warnings.warn(
            f"[CRS] {name} 没有 CRS，已默认设置为 EPSG:4326。请确认原始数据是否确实为经纬度坐标。"
        )
        gdf = gdf.set_crs(epsg=4326)

    if str(gdf.crs) != str(target_crs):
        gdf = gdf.to_crs(target_crs)

    return gdf


def drop_zm(geom):
    """删除 Z/M 维度，只保留二维坐标。ArcGIS Pro 对部分带 Z/M 的 GPKG 面要素较敏感。"""
    if geom is None or geom.is_empty:
        return None
    try:
        return transform(lambda x, y, *args: (x, y), geom)
    except Exception:
        return geom


def extract_polygonal(geom):
    """
    从任意几何中提取 Polygon / MultiPolygon。
    make_valid 后可能产生 GeometryCollection，其中包含线或点；这些会导致 ArcGIS Pro 报形状完整性错误。
    """
    if geom is None or geom.is_empty:
        return None

    geom_type = geom.geom_type

    if geom_type == "Polygon":
        return geom if geom.area > 0 else None

    if geom_type == "MultiPolygon":
        polys = [g for g in geom.geoms if (g is not None and not g.is_empty and g.area > 0)]
        if len(polys) == 0:
            return None
        if len(polys) == 1:
            return polys[0]
        return MultiPolygon(polys)

    if geom_type == "GeometryCollection":
        polys = []
        for part in geom.geoms:
            pg = extract_polygonal(part)
            if pg is None or pg.is_empty:
                continue
            if pg.geom_type == "Polygon":
                polys.append(pg)
            elif pg.geom_type == "MultiPolygon":
                polys.extend([p for p in pg.geoms if p.area > 0])
        if len(polys) == 0:
            return None
        try:
            u = unary_union(polys)
            return extract_polygonal(u)
        except Exception:
            return MultiPolygon(polys) if len(polys) > 1 else polys[0]

    # LineString / Point 等非面要素直接丢弃
    return None


def repair_polygonal_geometry(geom):
    """
    面要素专用修复函数：返回 ArcGIS Pro 更容易接受的二维 Polygon/MultiPolygon。
    """
    if geom is None or geom.is_empty:
        return None

    geom = drop_zm(geom)

    # 先尝试 make_valid
    try:
        geom = make_valid(geom)
    except Exception:
        pass

    geom = extract_polygonal(geom)
    if geom is None or geom.is_empty:
        return None

    # buffer(0) 对部分自交、环方向、微小拓扑问题更稳
    try:
        if (not geom.is_valid) or geom.geom_type == "GeometryCollection":
            geom = geom.buffer(0)
    except Exception:
        try:
            geom = geom.buffer(0)
        except Exception:
            return None

    try:
        geom = make_valid(geom)
    except Exception:
        pass

    geom = extract_polygonal(geom)
    if geom is None or geom.is_empty:
        return None

    try:
        if geom.area <= 0:
            return None
    except Exception:
        return None

    return geom


def fix_geometries(gdf, name="vector"):
    """修复面几何，并删除无法修复为空几何的记录。"""
    before = len(gdf)

    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[~gdf.geometry.is_empty].copy()
    gdf["geometry"] = gdf.geometry.apply(repair_polygonal_geometry)
    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[~gdf.geometry.is_empty].copy()

    after = len(gdf)
    if after < before:
        print(f"  [GEOM ] {name}: 删除空、非面或无法修复几何 {before - after} 条")

    return gdf


def prepare_gdf_for_arcgis(gdf, name="layer", output_crs=OUTPUT_GPKG_CRS):
    """
    写出 GPKG 前的 ArcGIS Pro 兼容性处理：
    1) 只保留二维 Polygon/MultiPolygon；
    2) 再次修复无效几何；
    3) bool 字段转 0/1，避免部分 GIS 读取不稳定；
    4) 可选转为 EPSG:4326，只影响输出显示几何，不影响已计算面积字段。
    """
    out = gdf.copy()

    out = out[out.geometry.notna()].copy()
    out["geometry"] = out.geometry.apply(repair_polygonal_geometry)
    out = out[out.geometry.notna()].copy()
    out = out[~out.geometry.is_empty].copy()

    for col in out.columns:
        if col == "geometry":
            continue
        if str(out[col].dtype) == "bool":
            out[col] = out[col].astype("int16")

    # 对 object 字段做基础清理，避免 None / 过长异常对象
    for col in out.select_dtypes(include=["object"]).columns:
        out[col] = out[col].fillna("").astype(str)

    if output_crs is not None and str(out.crs) != str(output_crs):
        out = out.to_crs(output_crs)
        out["geometry"] = out.geometry.apply(repair_polygonal_geometry)
        out = out[out.geometry.notna()].copy()
        out = out[~out.geometry.is_empty].copy()

    print(f"  [ARCGIS] {name}: 写出前有效面要素 {len(out):,} 条，CRS={out.crs}")
    return out



def clean_species_id(series):
    """清理物种 ID，例如 12345.0 转为 12345。"""
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


def safe_union(geometries):
    """对几何列表进行 union，返回严格空间去重后的联合几何。"""
    valid_geoms = [geom for geom in geometries if geom is not None and not geom.is_empty]

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

    union_geom = repair_polygonal_geometry(union_geom)
    if union_geom is None or union_geom.is_empty:
        return None

    return union_geom


def safe_intersection(geom_a, geom_b):
    """安全计算两个几何的相交。"""
    if geom_a is None or geom_b is None:
        return None
    if geom_a.is_empty or geom_b.is_empty:
        return None

    try:
        inter = geom_a.intersection(geom_b)
    except Exception:
        try:
            inter = make_valid(geom_a).intersection(make_valid(geom_b))
        except Exception:
            return None

    if inter is None or inter.is_empty:
        return None

    inter = repair_polygonal_geometry(inter)
    if inter is None or inter.is_empty:
        return None

    return inter


# ============================================================
# 数据读取与预处理
# ============================================================

def read_and_prepare_countries():
    """读取 GADM 国家边界，并生成统一 Country / Continent 字段。"""
    if not GADM_PATH.exists():
        raise FileNotFoundError(f"GADM 国家边界文件不存在: {GADM_PATH}")

    print("=" * 70)
    print("Step 1: 读取 GADM 国家边界，并统一投影")
    print("=" * 70)

    countries = read_vector(GADM_PATH)
    print(f"  输入 GADM: {GADM_PATH}")
    print(f"  原始国家面数量: {len(countries):,}")
    print(f"  原始 CRS: {countries.crs}")

    country_field = detect_field(countries.columns, GADM_COUNTRY_FIELD_CANDIDATES, "GADM 国家名称")

    # 用户提供的 GC_GADM.xlsx 属性表显示该数据具有 NAME_0 与 GID_0 字段；
    # 如果实际 shp 中也存在 GID_0，则优先用 GID_0 映射大洲，避免国家名写法差异导致 Unknown。
    has_gid = GADM_GID_FIELD in countries.columns
    keep_cols = [country_field, "geometry"] + ([GADM_GID_FIELD] if has_gid else [])
    countries = countries[keep_cols].rename(columns={country_field: COUNTRY_FIELD})
    countries[COUNTRY_FIELD] = countries[COUNTRY_FIELD].astype(str).str.strip()
    countries = countries[countries[COUNTRY_FIELD].str.len() > 0].copy()

    if has_gid:
        countries[GADM_GID_FIELD] = countries[GADM_GID_FIELD].astype(str).str.strip().str.upper()
        gid_lookup = (
            countries[[COUNTRY_FIELD, GADM_GID_FIELD]]
            .dropna()
            .drop_duplicates(subset=[COUNTRY_FIELD])
            .set_index(COUNTRY_FIELD)[GADM_GID_FIELD]
        )
    else:
        gid_lookup = None

    countries = ensure_crs_and_project(countries, TARGET_CRS, name="GADM countries")
    countries = fix_geometries(countries, name="GADM countries")

    # 如果 GADM 同一国家有多个要素，合并为国家级面
    countries = countries.dissolve(by=COUNTRY_FIELD, as_index=False)[[COUNTRY_FIELD, "geometry"]]
    countries = fix_geometries(countries, name="Dissolved GADM countries")

    if gid_lookup is not None:
        countries[GADM_GID_FIELD] = countries[COUNTRY_FIELD].map(gid_lookup).fillna("")
    else:
        countries[GADM_GID_FIELD] = ""

    countries[CONTINENT_FIELD] = countries.apply(
        lambda row: infer_continent(row[COUNTRY_FIELD], row[GADM_GID_FIELD]),
        axis=1
    )

    unmatched = countries[countries[CONTINENT_FIELD] == "Unknown"][COUNTRY_FIELD].tolist()
    if len(unmatched) > 0:
        warnings.warn(
            f"有 {len(unmatched)} 个 GADM 国家未匹配到大洲，将标记为 Unknown。\n"
            f"未匹配国家示例: {unmatched[:30]}\n"
            f"如需彻底消除 Unknown，请检查这些国家的 NAME_0/GID_0 是否不在内置映射表中。"
        )
    else:
        print("  所有 GADM 国家均已匹配到大洲，不会在大洲统计中出现 Unknown。")

    print(f"  处理后国家数: {len(countries):,}")
    print("  大洲分布：")
    print(countries[CONTINENT_FIELD].value_counts(dropna=False).to_string())

    return countries



def read_and_prepare_mines(countries):
    """
    读取铜矿影响范围，创建 mine_id，并根据铜矿属性表中的国家字段确定矿区所属国家。

    重要：
    - 本函数不再根据矿区缓冲区与 GADM 国界的空间位置重新分配国家。
    - Country 字段表示矿区所属国家 / 资源拥有国 / 项目所在国，是后续所有国家归属统计的唯一依据。
    - 缓冲区几何保留完整范围，不按国界裁剪。
    """
    if not MINE_SHP.exists():
        raise FileNotFoundError(f"铜矿影响范围文件不存在: {MINE_SHP}")

    print("\n" + "=" * 70)
    print("Step 2: 读取铜矿影响范围，并按属性字段确定矿区所属国家")
    print("=" * 70)

    mine = read_vector(MINE_SHP)
    print(f"  输入铜矿影响范围: {MINE_SHP}")
    print(f"  原始铜矿要素数: {len(mine):,}")
    print(f"  原始 CRS: {mine.crs}")

    # 自动识别矿区 ID；若没有则用原始行号生成
    mine_id_field = None
    for f in MINE_ID_FIELD_CANDIDATES:
        if f in mine.columns:
            mine_id_field = f
            break

    if mine_id_field is None:
        print("  [MINE_ID] 未找到矿区唯一 ID 字段，将使用原始行号生成 mine_id。")
        mine[MINE_ID_FIELD] = np.arange(1, len(mine) + 1).astype(str)
    else:
        print(f"  [MINE_ID] 使用字段 {mine_id_field} 作为矿区 ID。")
        mine[MINE_ID_FIELD] = mine[mine_id_field].astype(str).str.strip()
        empty_id = (
            mine[MINE_ID_FIELD].isna()
            | (mine[MINE_ID_FIELD].str.len() == 0)
            | (mine[MINE_ID_FIELD].str.lower() == "nan")
        )
        if empty_id.any():
            mine.loc[empty_id, MINE_ID_FIELD] = (mine.loc[empty_id].index + 1).astype(str)

    # 识别矿区所属国家字段
    mine_country_field = None
    for f in MINE_COUNTRY_FIELD_CANDIDATES:
        if f in mine.columns:
            mine_country_field = f
            break

    if mine_country_field is None:
        raise ValueError(
            "铜矿数据中未找到矿区所属国家字段。\n"
            f"候选字段为: {MINE_COUNTRY_FIELD_CANDIDATES}\n"
            f"实际字段为: {list(mine.columns)}\n"
            "由于当前统计口径要求按矿区所属国家归属，不能使用空间落点国家自动替代。"
        )

    print(f"  [COUNTRY] 使用铜矿属性字段 {mine_country_field} 作为矿区所属国家。")

    keep_cols = [MINE_ID_FIELD, mine_country_field, "geometry"]
    mine = mine[keep_cols].copy()

    mine = standardize_mine_country_names(mine, countries, mine_country_field)

    mine = mine[[MINE_ID_FIELD, MINE_COUNTRY_RAW_FIELD, COUNTRY_FIELD, CONTINENT_FIELD, "geometry"]].copy()
    mine = ensure_crs_and_project(mine, TARGET_CRS, name="Copper mine buffer")
    mine = fix_geometries(mine, name="Copper mine buffer")

    # 同一 mine_id 理论上只对应一个所属国家。如果存在同一 mine_id 多国家，停止检查，避免归属混乱。
    multi_country = (
        mine[[MINE_ID_FIELD, COUNTRY_FIELD]]
        .drop_duplicates()
        .groupby(MINE_ID_FIELD)[COUNTRY_FIELD]
        .nunique()
    )
    bad_ids = multi_country[multi_country > 1]
    if len(bad_ids) > 0:
        bad = mine[mine[MINE_ID_FIELD].isin(bad_ids.index)][
            [MINE_ID_FIELD, MINE_COUNTRY_RAW_FIELD, COUNTRY_FIELD]
        ].drop_duplicates()
        bad_path = OUT_DIR / "Mine_ID_Multiple_Countries_Check.csv"
        bad.to_csv(bad_path, index=False, encoding="utf-8-sig")
        raise ValueError(
            f"发现 {len(bad_ids)} 个 mine_id 对应多个所属国家，已保存清单: {bad_path}\n"
            "请先修正矿区 ID 或国家字段后再运行。"
        )

    print(f"  处理后铜矿影响范围数量: {len(mine):,}")
    print(f"  涉及矿区所属国家数: {mine[COUNTRY_FIELD].nunique():,}")
    print("  矿区所属国家前10：")
    print(mine[COUNTRY_FIELD].value_counts().head(10).to_string())

    return mine


def build_mine_country_buffers(mine):
    """
    按矿区属性中的所属国家组织铜矿影响范围。

    返回：
    country_mines:
        每行表示一个矿区完整影响范围及其所属国家，不按国界裁剪。
    country_mine_dissolved:
        按矿区所属国家 dissolve 后的完整铜矿影响范围，用于国家归属口径下的空间去重统计。
    """
    print("\n" + "=" * 70)
    print("Step 3: 按矿区所属国家汇总完整铜矿影响范围，不按国界裁剪")
    print("=" * 70)

    country_mines = mine[[COUNTRY_FIELD, CONTINENT_FIELD, MINE_ID_FIELD, MINE_COUNTRY_RAW_FIELD, "geometry"]].copy()
    country_mines = fix_geometries(country_mines, name="Mine-country attributed full buffers")

    # 同一 mine_id 如有多个面片，按所属国家合并；几何仍为完整缓冲区 union。
    country_mines = country_mines.dissolve(
        by=[COUNTRY_FIELD, CONTINENT_FIELD, MINE_ID_FIELD],
        as_index=False,
        aggfunc={MINE_COUNTRY_RAW_FIELD: "first"}
    )
    country_mines["mine_buffer_area_km2"] = country_mines.geometry.area / 1e6

    print(f"  国家归属口径下矿区记录数: {len(country_mines):,}")
    print(f"  涉及矿区所属国家数: {country_mines[COUNTRY_FIELD].nunique():,}")
    print(f"  涉及唯一铜矿数: {country_mines[MINE_ID_FIELD].nunique():,}")

    print("  按矿区所属国家 dissolve 完整铜矿影响范围，用于后续空间去重统计...")
    country_mine_dissolved = country_mines.dissolve(
        by=[COUNTRY_FIELD, CONTINENT_FIELD],
        as_index=False
    )[[COUNTRY_FIELD, CONTINENT_FIELD, "geometry"]]
    country_mine_dissolved = fix_geometries(
        country_mine_dissolved,
        name="Country-attributed dissolved mine buffers"
    )
    country_mine_dissolved["total_mine_buffer_area_km2"] = country_mine_dissolved.geometry.area / 1e6

    print("  注意：total_mine_buffer_area_km2 为该国名下矿区完整影响范围 union 面积，可能包含跨境部分。")

    return country_mines, country_mine_dissolved


# ============================================================
# 核心函数：处理一个类群 × 濒危等级组合
# ============================================================

def process_one_group(taxon, status, country_mine_dissolved, country_mines):
    """
    处理一个 taxon-status 组合。

    返回：
    1. gdf_intersections:
       每一行表示某矿区所属国家名下，一个物种的栖息地与该国名下铜矿完整影响范围的相交结果。
       字段：Country, Continent, species_id, taxon, status, species_overlap_area_km2, geometry

    2. conflict_mine_pairs:
       字段：Country, Continent, mine_id，表示该所属国家名下该铜矿至少与一种目标物种栖息地相交。
    """
    sub = IUCN_ROOT / f"{taxon}-{status}"
    shp = sub / "data_0.shp"
    interim_path = INTERIM_DIR / f"{taxon}_{status}_country_species_intersections.parquet"
    interim_mine_path = INTERIM_DIR / f"{taxon}_{status}_conflict_mine_pairs.csv"

    if USE_CACHE and interim_path.exists() and interim_mine_path.exists():
        print(f"  [CACHE] 已有中间结果: {interim_path.name}")
        gdf = gpd.read_parquet(interim_path)
        pairs = pd.read_csv(interim_mine_path, dtype={MINE_ID_FIELD: str})
        return gdf, pairs

    if not shp.exists():
        warnings.warn(f"  [SKIP] 数据文件不存在: {shp}")
        return None, None

    t0 = time.time()
    print(f"  [READ ] {taxon}-{status}: {shp.name}")

    habitat = read_vector(shp)
    if habitat.empty:
        warnings.warn(f"  [EMPTY] {shp} 为空")
        return None, None

    id_field = detect_id_field(habitat.columns)

    habitat = ensure_crs_and_project(habitat, TARGET_CRS, name=f"{taxon}-{status}")
    habitat = fix_geometries(habitat, name=f"{taxon}-{status}")

    habitat = habitat[[id_field, "geometry"]].rename(columns={id_field: "species_id"})
    habitat["species_id"] = clean_species_id(habitat["species_id"])
    habitat = habitat[
        habitat["species_id"].notna()
        & (habitat["species_id"].astype(str).str.len() > 0)
        & (habitat["species_id"].astype(str).str.lower() != "nan")
    ].copy()

    if habitat.empty:
        warnings.warn(f"  [EMPTY] {taxon}-{status}: 清理后没有有效物种记录")
        return None, None

    # --------------------------------------------------------
    # A. 国家归属 × 物种口径：使用按矿区所属国家 dissolve 后的完整矿区范围，计算相交几何
    # --------------------------------------------------------
    print(
        f"  [SJOIN] {taxon}-{status}: "
        f"{len(habitat):,} 个栖息地面 -> 预筛选国家归属矿区范围与物种相交关系"
    )

    mine_for_join = country_mine_dissolved[[COUNTRY_FIELD, CONTINENT_FIELD, "geometry"]].copy()

    joined = gpd.sjoin(
        mine_for_join,
        habitat[["species_id", "geometry"]],
        how="inner",
        predicate="intersects"
    )

    if joined.empty:
        print(f"  [EMPTY] {taxon}-{status}: 无国家归属-物种空间交集")
        return None, None

    print(
        f"  [INTER] {taxon}-{status}: "
        f"{len(joined):,} 条候选国家归属-物种关系 -> 精确相交"
    )

    results = []
    result_geoms = []

    dissolved_index = country_mine_dissolved.set_index(COUNTRY_FIELD)
    countries = joined[COUNTRY_FIELD].dropna().unique()

    for idx_country, country in enumerate(countries, 1):
        if idx_country % 10 == 0 or idx_country == 1 or idx_country == len(countries):
            print(f"    {taxon}-{status}: 处理矿区所属国家 {idx_country}/{len(countries)} - {country}")

        try:
            mine_geom = dissolved_index.loc[country, "geometry"]
            continent = dissolved_index.loc[country, CONTINENT_FIELD]
        except Exception:
            continue

        country_joined = joined[joined[COUNTRY_FIELD] == country].copy()
        species_list = country_joined["species_id"].dropna().unique()

        for species_id in species_list:
            rows_sp = country_joined[country_joined["species_id"] == species_id]
            habitat_indices = rows_sp["index_right"].dropna().unique()
            if len(habitat_indices) == 0:
                continue

            species_geoms = habitat.loc[habitat_indices, "geometry"]

            # 对同一物种的多个栖息地面片进行 union，避免同一物种自身面片重叠导致面积重复。
            species_union = safe_union(species_geoms)
            if species_union is None:
                continue

            inter_geom = safe_intersection(mine_geom, species_union)
            if inter_geom is None:
                continue

            area_km2 = inter_geom.area / 1e6
            if area_km2 <= 0:
                continue

            results.append({
                COUNTRY_FIELD: country,
                CONTINENT_FIELD: continent,
                "species_id": species_id,
                "taxon": taxon,
                "status": status,
                "species_overlap_area_km2": area_km2
            })
            result_geoms.append(inter_geom)

    if len(results) == 0:
        print(f"  [EMPTY] {taxon}-{status}: 精确相交后无有效国家-物种结果")
        gdf_intersections = None
    else:
        gdf_intersections = gpd.GeoDataFrame(results, geometry=result_geoms, crs=TARGET_CRS)
        gdf_intersections.to_parquet(interim_path, index=False)

    # --------------------------------------------------------
    # B. 国家归属 × 矿区口径：统计与任意目标物种栖息地相交的铜矿区
    # --------------------------------------------------------
    print(f"  [MINE ] {taxon}-{status}: 统计冲突铜矿区 mine-country 对")
    mine_join = gpd.sjoin(
        country_mines[[COUNTRY_FIELD, CONTINENT_FIELD, MINE_ID_FIELD, "geometry"]],
        habitat[["species_id", "geometry"]],
        how="inner",
        predicate="intersects"
    )

    if mine_join.empty:
        conflict_mine_pairs = None
    else:
        conflict_mine_pairs = (
            mine_join[[COUNTRY_FIELD, CONTINENT_FIELD, MINE_ID_FIELD]]
            .drop_duplicates()
            .copy()
        )
        conflict_mine_pairs[MINE_ID_FIELD] = conflict_mine_pairs[MINE_ID_FIELD].astype(str)
        conflict_mine_pairs.to_csv(interim_mine_path, index=False, encoding="utf-8-sig")

    dt = time.time() - t0
    n_rec = 0 if gdf_intersections is None else len(gdf_intersections)
    n_sp = 0 if gdf_intersections is None else gdf_intersections["species_id"].nunique()
    n_cty = 0 if gdf_intersections is None else gdf_intersections[COUNTRY_FIELD].nunique()
    n_mine_pairs = 0 if conflict_mine_pairs is None else len(conflict_mine_pairs)

    print(
        f"  [OK   ] {taxon}-{status}: "
        f"国家-物种记录 {n_rec:,}，unique species = {n_sp:,}，"
        f"countries = {n_cty:,}，冲突 mine-country 对 = {n_mine_pairs:,}，"
        f"用时 {dt:.1f}s"
    )

    del habitat, joined
    gc.collect()

    return gdf_intersections, conflict_mine_pairs


# ============================================================
# 全球物种尺度严格空间去重面积
# ============================================================

def generate_global_unique_species_overlap(intersections):
    """
    对每个 species_id × taxon × status，将其在所有国家中的相交几何全部 union，
    输出全球物种尺度严格空间去重重叠面积。
    """
    print("\n" + "=" * 70)
    print("额外统计：全球物种尺度严格空间去重重叠面积")
    print("=" * 70)

    global_species_rows = []
    global_species_geoms = []

    species_groups = list(intersections.groupby(["species_id", "taxon", "status"]))

    for i, ((species_id, taxon, status), group) in enumerate(species_groups, 1):
        if i % 200 == 0 or i == 1 or i == len(species_groups):
            print(f"  全球物种去重面积 {i}/{len(species_groups)}: {species_id} - {taxon} - {status}")

        union_geom = safe_union(group.geometry)
        unique_area_km2 = 0.0 if union_geom is None else union_geom.area / 1e6
        country_accumulated_area_km2 = group["species_overlap_area_km2"].sum()

        global_species_rows.append({
            "species_id": species_id,
            "taxon": taxon,
            "status": status,
            "global_unique_overlap_area_km2": unique_area_km2,
            "country_accumulated_overlap_area_km2": country_accumulated_area_km2,
            "overlap_duplicate_area_km2": country_accumulated_area_km2 - unique_area_km2,
            "overlap_duplicate_percent_of_unique": (
                (country_accumulated_area_km2 - unique_area_km2) / unique_area_km2 * 100
                if unique_area_km2 > 0 else 0.0
            ),
            "num_countries": group[COUNTRY_FIELD].nunique(),
            "country_occurrences": len(group),
            "countries": ";".join(sorted(group[COUNTRY_FIELD].dropna().astype(str).unique()))
        })
        global_species_geoms.append(union_geom)

    global_species_unique = pd.DataFrame(global_species_rows)

    print(f"  保存全球物种尺度严格去重面积表: {OUT_BY_SPECIES_GLOBAL_UNIQUE}")
    global_species_unique.to_csv(OUT_BY_SPECIES_GLOBAL_UNIQUE, index=False, encoding="utf-8-sig")

    try:
        global_species_gdf = gpd.GeoDataFrame(
            global_species_unique,
            geometry=global_species_geoms,
            crs=TARGET_CRS
        )
        global_species_gdf = prepare_gdf_for_arcgis(
            global_species_gdf,
            name="global_unique_species_overlap",
            output_crs=OUTPUT_GPKG_CRS
        )
        print(f"  保存全球物种尺度严格去重几何: {OUT_BY_SPECIES_GLOBAL_UNIQUE_GEOM}")
        if OUT_BY_SPECIES_GLOBAL_UNIQUE_GEOM.exists():
            OUT_BY_SPECIES_GLOBAL_UNIQUE_GEOM.unlink()
        global_species_gdf.to_file(
            OUT_BY_SPECIES_GLOBAL_UNIQUE_GEOM,
            layer="global_unique_species_overlap",
            driver="GPKG"
        )
    except Exception as e:
        warnings.warn(f"保存全球物种尺度 GPKG 失败，但 CSV 已保存。错误信息: {e}")

    return global_species_unique


# ============================================================
# 国家级与大洲级汇总
# ============================================================

def build_country_summary(countries, country_mines, country_mine_dissolved, intersections, conflict_mine_pairs):
    """生成按矿区所属国家归属的国家级汇总表和 GPKG。"""
    print("\n" + "=" * 70)
    print("Step 5: 生成国家级汇总 GPKG")
    print("=" * 70)

    # 1. 各矿区所属国家涉及物种数量；同一物种可在不同所属国家重复统计
    species_count = (
        intersections.groupby(COUNTRY_FIELD)["species_id"]
        .nunique()
        .rename("num_species")
    )

    # 2. 国家归属口径下空间冲突总面积：对该所属国家名下所有国家-物种相交几何 union 后求面积
    print("  计算国家归属口径下严格空间去重冲突面积...")
    strict_overlap_rows = []
    strict_overlap_geoms = []

    for i, (country, group) in enumerate(intersections.groupby(COUNTRY_FIELD), 1):
        if i % 10 == 0 or i == 1 or i == intersections[COUNTRY_FIELD].nunique():
            print(f"    严格去重面积: 处理矿区所属国家 {i}/{intersections[COUNTRY_FIELD].nunique()} - {country}")

        union_geom = safe_union(group.geometry)
        area_km2 = 0.0 if union_geom is None else union_geom.area / 1e6
        continent = group[CONTINENT_FIELD].iloc[0] if CONTINENT_FIELD in group.columns else infer_continent(country)

        strict_overlap_rows.append({
            COUNTRY_FIELD: country,
            CONTINENT_FIELD: continent,
            "overlap_area_km2": area_km2
        })
        strict_overlap_geoms.append(union_geom)

    overlap_df = pd.DataFrame(strict_overlap_rows)
    overlap_gdf = gpd.GeoDataFrame(overlap_df, geometry=strict_overlap_geoms, crs=TARGET_CRS)

    # 3. 每国冲突铜矿区数量：同一个 mine_id 在同一国家只计一次
    if conflict_mine_pairs is None or len(conflict_mine_pairs) == 0:
        conflict_mine_count = pd.Series(dtype=int, name="conflict_mine_count")
    else:
        conflict_mine_count = (
            conflict_mine_pairs.drop_duplicates([COUNTRY_FIELD, MINE_ID_FIELD])
            .groupby(COUNTRY_FIELD)[MINE_ID_FIELD]
            .nunique()
            .rename("conflict_mine_count")
        )

    # 4. 每个矿区所属国家名下全部铜矿影响范围的矿区数量与 union 面积
    all_mine_count = (
        country_mines.drop_duplicates([COUNTRY_FIELD, MINE_ID_FIELD])
        .groupby(COUNTRY_FIELD)[MINE_ID_FIELD]
        .nunique()
        .rename("total_mine_count")
    )

    total_mine_buffer_area = (
        country_mine_dissolved.set_index(COUNTRY_FIELD)["total_mine_buffer_area_km2"]
        .rename("total_mine_buffer_area_km2")
    )

    # 5. 组合到 GADM 国家面上，保留所有国家
    summary = countries[[COUNTRY_FIELD, CONTINENT_FIELD, "geometry"]].copy()
    summary = summary.merge(species_count.reset_index(), on=COUNTRY_FIELD, how="left")
    summary = summary.merge(conflict_mine_count.reset_index(), on=COUNTRY_FIELD, how="left")
    summary = summary.merge(all_mine_count.reset_index(), on=COUNTRY_FIELD, how="left")
    summary = summary.merge(total_mine_buffer_area.reset_index(), on=COUNTRY_FIELD, how="left")
    summary = summary.merge(overlap_df[[COUNTRY_FIELD, "overlap_area_km2"]], on=COUNTRY_FIELD, how="left")

    fill_int_cols = ["num_species", "conflict_mine_count", "total_mine_count"]
    for col in fill_int_cols:
        summary[col] = summary[col].fillna(0).astype(int)

    fill_float_cols = ["total_mine_buffer_area_km2", "overlap_area_km2"]
    for col in fill_float_cols:
        summary[col] = summary[col].fillna(0.0).astype(float)

    summary["overlap_pct_of_mine_buffer"] = np.where(
        summary["total_mine_buffer_area_km2"] > 0,
        summary["overlap_area_km2"] / summary["total_mine_buffer_area_km2"] * 100,
        0.0
    )

    # 保留所有 GADM 国家：
    # - 没有铜矿区的国家：total_mine_count = 0, total_mine_buffer_area_km2 = 0
    # - 有铜矿但没有与目标物种栖息地冲突的国家：conflict_mine_count = 0, overlap_area_km2 = 0, num_species = 0
    summary["has_mine"] = summary["total_mine_count"] > 0
    summary["has_overlap"] = summary["overlap_area_km2"] > 0

    # 字段排序
    summary = summary[[
        COUNTRY_FIELD,
        CONTINENT_FIELD,
        "conflict_mine_count",
        "total_mine_count",
        "overlap_area_km2",
        "num_species",
        "total_mine_buffer_area_km2",
        "overlap_pct_of_mine_buffer",
        "has_mine",
        "has_overlap",
        "geometry"
    ]]

    summary_sorted = summary.drop(columns="geometry").sort_values(
        ["overlap_area_km2", "num_species"], ascending=[False, False]
    )

    print("\n  国家级汇总前10：")
    print(summary_sorted.head(10).to_string(index=False))

    print(f"\n  保存国家级 CSV: {OUT_COUNTRY_SUMMARY_CSV}")
    summary_sorted.to_csv(OUT_COUNTRY_SUMMARY_CSV, index=False, encoding="utf-8-sig")

    print(f"  保存国家级 GPKG: {OUT_COUNTRY_SUMMARY_GPKG}")
    if OUT_COUNTRY_SUMMARY_GPKG.exists():
        OUT_COUNTRY_SUMMARY_GPKG.unlink()

    summary_out = prepare_gdf_for_arcgis(
        summary,
        name="country_summary",
        output_crs=OUTPUT_GPKG_CRS
    )
    summary_out.to_file(
        OUT_COUNTRY_SUMMARY_GPKG,
        layer="country_summary",
        driver="GPKG"
    )

    # 保存国家归属口径下的冲突 union 几何作为第二个图层，仅包含存在冲突的矿区所属国家。
    # 注意：该图层几何可能跨出对应国家边界，因为几何表示该国名下矿区在全球范围内造成的影响区域。
    if len(overlap_gdf) > 0:
        overlap_stats = summary.drop(columns="geometry")
        overlap_union_layer = overlap_gdf.merge(
            overlap_stats.drop(columns=[CONTINENT_FIELD], errors="ignore"),
            on=COUNTRY_FIELD,
            how="left",
            suffixes=("", "_summary")
        )
        # 清理重复字段
        if "overlap_area_km2_summary" in overlap_union_layer.columns:
            overlap_union_layer = overlap_union_layer.drop(columns=["overlap_area_km2_summary"])

        overlap_union_layer = prepare_gdf_for_arcgis(
            overlap_union_layer,
            name="country_overlap_union",
            output_crs=OUTPUT_GPKG_CRS
        )
        overlap_union_layer.to_file(
            OUT_COUNTRY_SUMMARY_GPKG,
            layer="country_overlap_union",
            driver="GPKG"
        )

    return summary, overlap_gdf


def build_continent_summary(summary, intersections, conflict_mine_pairs):
    """生成大洲级面积与物种汇总。"""
    print("\n" + "=" * 70)
    print("Step 6: 大洲级汇总统计")
    print("=" * 70)

    summary_df = pd.DataFrame(summary.drop(columns="geometry"))
    summary_overlap = summary_df[summary_df["overlap_area_km2"] > 0].copy()
    global_overlap_total = summary_overlap["overlap_area_km2"].sum()

    rows = []
    for continent, continent_data in summary_df.groupby(CONTINENT_FIELD):
        countries_in_continent = continent_data[COUNTRY_FIELD].tolist()
        continent_overlap = continent_data["overlap_area_km2"].sum()

        continent_intersections = intersections[intersections[COUNTRY_FIELD].isin(countries_in_continent)]
        num_species_unique = continent_intersections["species_id"].nunique()
        country_species_occurrences = int(continent_data["num_species"].sum())

        if conflict_mine_pairs is not None and len(conflict_mine_pairs) > 0:
            c_pairs = conflict_mine_pairs[conflict_mine_pairs[COUNTRY_FIELD].isin(countries_in_continent)]
            conflict_mine_unique = c_pairs[MINE_ID_FIELD].nunique()
            conflict_mine_country_occurrences = len(c_pairs.drop_duplicates([COUNTRY_FIELD, MINE_ID_FIELD]))
        else:
            conflict_mine_unique = 0
            conflict_mine_country_occurrences = 0

        rows.append({
            CONTINENT_FIELD: continent,
            "overlap_area_km2": continent_overlap,
            "overlap_pct_of_global": (continent_overlap / global_overlap_total * 100) if global_overlap_total > 0 else 0.0,
            "num_countries_total": len(continent_data),
            "num_countries_with_overlap": int((continent_data["overlap_area_km2"] > 0).sum()),
            "num_species_unique": num_species_unique,
            "country_species_occurrences": country_species_occurrences,
            "conflict_mine_unique": conflict_mine_unique,
            "conflict_mine_country_occurrences": conflict_mine_country_occurrences,
            "total_mine_buffer_area_km2": continent_data["total_mine_buffer_area_km2"].sum()
        })

    continent_summary = pd.DataFrame(rows).sort_values("overlap_area_km2", ascending=False)

    print("\n  大洲级汇总：")
    print(continent_summary.to_string(index=False))

    print(f"\n  保存大洲级汇总表: {OUT_CONTINENT_SUMMARY}")
    continent_summary.to_csv(OUT_CONTINENT_SUMMARY, index=False, encoding="utf-8-sig")

    return continent_summary


# ============================================================
# 主流程
# ============================================================

def main():
    countries = read_and_prepare_countries()
    mine = read_and_prepare_mines(countries)
    country_mines, country_mine_dissolved = build_mine_country_buffers(mine)

    print("\n" + "=" * 70)
    print("Step 4: 遍历 15 组物种栖息地，按矿区所属国家计算国家-物种相交与冲突矿区数量")
    print("=" * 70)

    all_intersections = []
    all_conflict_mine_pairs = []

    for taxon in TAXA:
        for status in STATUSES:
            gdf, pairs = process_one_group(taxon, status, country_mine_dissolved, country_mines)
            if gdf is not None and len(gdf) > 0:
                all_intersections.append(gdf)
            if pairs is not None and len(pairs) > 0:
                all_conflict_mine_pairs.append(pairs)

    if len(all_intersections) == 0:
        raise RuntimeError("没有成功处理任何 taxon-status 组合，请检查 IUCN 数据路径与字段。")

    intersections = gpd.GeoDataFrame(
        pd.concat(all_intersections, ignore_index=True),
        geometry="geometry",
        crs=TARGET_CRS
    )

    if len(all_conflict_mine_pairs) == 0:
        conflict_mine_pairs = pd.DataFrame(columns=[COUNTRY_FIELD, CONTINENT_FIELD, MINE_ID_FIELD])
    else:
        conflict_mine_pairs = pd.concat(all_conflict_mine_pairs, ignore_index=True)
        conflict_mine_pairs[MINE_ID_FIELD] = conflict_mine_pairs[MINE_ID_FIELD].astype(str)
        conflict_mine_pairs = conflict_mine_pairs.drop_duplicates([COUNTRY_FIELD, MINE_ID_FIELD])

    print("\n[MERGE] 合并所有 taxon-status 结果")
    print(f"  国家-物种相交记录数: {len(intersections):,}")
    print(f"  涉及物种数: {intersections['species_id'].nunique():,}")
    print(f"  涉及国家数: {intersections[COUNTRY_FIELD].nunique():,}")
    print(f"  冲突 mine-country 对数量: {len(conflict_mine_pairs):,}")
    print(
        f"  国家×物种表中重叠面积累计值: "
        f"{intersections['species_overlap_area_km2'].sum():,.2f} km²"
    )

    # 输出国家 × 物种层面的非几何表
    by_species = intersections.drop(columns="geometry").copy()
    by_species = by_species.sort_values([COUNTRY_FIELD, "taxon", "status", "species_id"])
    print(f"\n  保存国家×物种相交面积表: {OUT_BY_SPECIES}")
    by_species.to_csv(OUT_BY_SPECIES, index=False, encoding="utf-8-sig")

    # 输出国家 × 矿区冲突表
    print(f"  保存国家×冲突矿区表: {OUT_BY_MINE}")
    conflict_mine_pairs.sort_values([COUNTRY_FIELD, MINE_ID_FIELD]).to_csv(
        OUT_BY_MINE,
        index=False,
        encoding="utf-8-sig"
    )

    # 全球物种尺度严格空间去重面积，供后续优先救助模块使用
    global_species_unique = generate_global_unique_species_overlap(intersections)

    # 国家级 GPKG 与 CSV
    summary, overlap_gdf = build_country_summary(
        countries,
        country_mines,
        country_mine_dissolved,
        intersections,
        conflict_mine_pairs
    )

    # 大洲级统计
    continent_summary = build_continent_summary(summary, intersections, conflict_mine_pairs)

    # 全球总览
    print("\n" + "=" * 70)
    print("Step 7: 全球总览与类群/濒危等级统计")
    print("=" * 70)

    # 当前脚本采用“矿区所属国家归属”口径，国家间几何不再天然互斥。
    # 因此：
    # - global_attributed_overlap_area_km2 是各矿区所属国家 overlap_area_km2 的求和，适合国家责任归属统计；
    # - global_spatial_unique_overlap_area_km2 是所有冲突几何在全球范围 union 后的严格空间去重面积，
    #   适合描述真实物理空间暴露规模。
    global_attributed_overlap_area_km2 = summary["overlap_area_km2"].sum()
    global_spatial_union_geom = safe_union(overlap_gdf.geometry) if len(overlap_gdf) > 0 else None
    global_spatial_unique_overlap_area_km2 = (
        0.0 if global_spatial_union_geom is None else global_spatial_union_geom.area / 1e6
    )

    # 为兼容旧版下游代码，保留 global_true_overlap_area_km2 字段；
    # 但在当前国家归属口径下，它等同于 global_attributed_overlap_area_km2，而不是全球物理空间严格去重面积。
    global_true_overlap_area_km2 = global_attributed_overlap_area_km2

    global_affected_species_count = intersections["species_id"].nunique()
    global_conflict_mine_country_occurrences = len(conflict_mine_pairs)
    global_conflict_mine_unique = conflict_mine_pairs[MINE_ID_FIELD].nunique() if len(conflict_mine_pairs) > 0 else 0

    global_overview = pd.DataFrame([{
        "global_true_overlap_area_km2": global_true_overlap_area_km2,
        "global_attributed_overlap_area_km2": global_attributed_overlap_area_km2,
        "global_spatial_unique_overlap_area_km2": global_spatial_unique_overlap_area_km2,
        "global_affected_species_count": global_affected_species_count,
        "num_mine_countries_with_overlap": int((summary["overlap_area_km2"] > 0).sum()),
        "global_conflict_mine_unique": global_conflict_mine_unique,
        "global_conflict_mine_country_occurrences": global_conflict_mine_country_occurrences,
        "global_total_mine_unique": mine[MINE_ID_FIELD].nunique(),
        "global_total_mine_country_occurrences": len(country_mines)
    }])

    print(f"  全球国家归属口径冲突面积累计值: {global_attributed_overlap_area_km2:,.2f} km²")
    print(f"  全球物理空间严格去重冲突面积: {global_spatial_unique_overlap_area_km2:,.2f} km²")
    print(f"  全球受影响物种总数量: {global_affected_species_count:,}")
    print(f"  全球冲突铜矿唯一数量: {global_conflict_mine_unique:,}")
    print(f"  全球冲突 mine-country 归属对数量: {global_conflict_mine_country_occurrences:,}")

    print(f"  保存全球总览表: {OUT_GLOBAL_OVERVIEW}")
    global_overview.to_csv(OUT_GLOBAL_OVERVIEW, index=False, encoding="utf-8-sig")

    affected_by_taxon = (
        intersections.groupby("taxon")["species_id"]
        .nunique()
        .reset_index(name="num_species")
        .sort_values("num_species", ascending=False)
    )
    print(f"  保存不同类群受影响物种数量表: {OUT_GLOBAL_BY_TAXON}")
    affected_by_taxon.to_csv(OUT_GLOBAL_BY_TAXON, index=False, encoding="utf-8-sig")

    affected_by_status = (
        intersections.groupby("status")["species_id"]
        .nunique()
        .reset_index(name="num_species")
    )
    affected_by_status["status"] = pd.Categorical(
        affected_by_status["status"],
        categories=STATUSES,
        ordered=True
    )
    affected_by_status = affected_by_status.sort_values("status").reset_index(drop=True)
    affected_by_status["status"] = affected_by_status["status"].astype(str)

    print(f"  保存不同濒危等级受影响物种数量表: {OUT_GLOBAL_BY_STATUS}")
    affected_by_status.to_csv(OUT_GLOBAL_BY_STATUS, index=False, encoding="utf-8-sig")

    print("\n" + "=" * 70)
    print("完成")
    print("=" * 70)
    print(f"  输出目录: {OUT_DIR}")
    zero_no_mine = int(((summary["total_mine_count"] == 0) & (summary["overlap_area_km2"] == 0) & (summary["num_species"] == 0)).sum())
    zero_no_overlap = int(((summary["total_mine_count"] > 0) & (summary["overlap_area_km2"] == 0) & (summary["num_species"] == 0)).sum())
    print(f"  无所属铜矿国家数（核心统计字段已填 0）: {zero_no_mine:,}")
    print(f"  有所属铜矿但无空间冲突国家数（核心统计字段已填 0）: {zero_no_overlap:,}")
    print(f"  国家级 GPKG: {OUT_COUNTRY_SUMMARY_GPKG.name}")
    print(f"  国家级 CSV: {OUT_COUNTRY_SUMMARY_CSV.name}")
    print(f"  国家×物种表: {OUT_BY_SPECIES.name}")
    print(f"  国家×冲突矿区表: {OUT_BY_MINE.name}")
    print(f"  大洲级汇总表: {OUT_CONTINENT_SUMMARY.name}")
    print(f"  全球物种去重表: {OUT_BY_SPECIES_GLOBAL_UNIQUE.name}")


if __name__ == "__main__":
    main()
