"""Read tabular results and the minimal geographic inputs used in the figures."""
from types import SimpleNamespace
import numpy as np
import pandas as pd
import geopandas as gpd
from common import INPUTS, table


def geography(layer):
    return gpd.read_file(INPUTS / 'main_geography.gpkg', layer=layer)


def lon_formatter(value, pos=None):
    value = int(round(value))
    return f'{abs(value)}°W' if value < 0 else (f'{value}°E' if value > 0 else '0°')


def lat_formatter(value, pos=None):
    value = int(round(value))
    return f'{abs(value)}°S' if value < 0 else (f'{value}°N' if value > 0 else '0°')


def safe_gdf_plot(frame, **kwargs):
    if not frame.empty:
        return frame.plot(**kwargs)


def richness_inputs():
    def load():
        rwr = pd.read_csv(INPUTS / 'figure_1_rwr_points.csv.gz')
        mines = pd.read_csv(INPUTS / 'figure_1_mine_points.csv')
        mines = gpd.GeoDataFrame(mines, geometry=gpd.points_from_xy(mines.longitude, mines.latitude), crs=4326)
        if len(rwr) != 1404410 or rwr.grid_id.duplicated().any():
            raise ValueError('Figure 1a requires the complete 1,404,410-cell grid')
        return rwr, geography('figure_1_land'), geography('figure_1_antarctica'), mines
    codes = ['p0_50', 'p50_75', 'p75_90', 'p90_99', 'p99_100']
    colors = ['#F0F0F0', '#BFD3E6', '#8C96C6', '#F6CD8E', '#CB1C2A']
    return SimpleNamespace(load_plot_data=load, safe_gdf_plot=safe_gdf_plot,
                           CLASS_FIELD='rwr_percentile_class', BIN_ORDER=codes,
                           BIN_COLOR=dict(zip(codes, colors)),
                           lon_formatter=lon_formatter, lat_formatter=lat_formatter)


def overlap_inputs():
    area_bins = [(0, 0, '0', '#000000'), (0, 5000, '(0, 5,000]', '#f7ede7'),
                 (5000, 10000, '(5,000, 10,000]', '#eeb09a'),
                 (10000, 50000, '(10,000, 50,000]', '#e56b53'),
                 (50000, 100000, '(50,000, 100,000]', '#cb1c2a'),
                 (100000, np.inf, '>100,000', '#66000e')]
    species_bins = [(1, 10, '1–10', '#d9f0d3'), (11, 15, '11–15', '#a6dba0'),
                    (16, 20, '16–20', '#5aae61'), (21, 50, '21–50', '#2f854f'),
                    (51, np.inf, '>50', '#145a32')]

    def load():
        world = geography('figure_1_overlap')
        points = geography('figure_1_overlap_points')
        if 'num_species' not in points.columns:
            points = points.merge(world[['Country', 'num_species']].drop_duplicates('Country'),
                                  on='Country', validate='many_to_one')
        points['species_count_class'] = [next(label for lower, upper, label, _ in species_bins
                                             if lower <= value <= upper) for value in points.num_species]
        exposed = table(4, 'Exposed_Species')
        counts = exposed.groupby('Taxonomic group')['IUCN species ID'].nunique()
        taxon = counts.reindex(['Amphibians', 'Birds', 'Fish', 'Mammals', 'Reptiles']).rename('num_species').rename_axis('taxon').reset_index()
        return world, points, taxon

    return SimpleNamespace(load_plot_data=load, AREA_FIELD='overlap_area_km2', AREA_BINS=area_bins,
                           SPECIES_COUNT_BINS=species_bins,
                           get_species_count_color=lambda label: next(c for _, _, lab, c in species_bins if lab == label),
                           lon_formatter=lon_formatter, lat_formatter=lat_formatter)


def country_area_data():
    return table(3, 'Country_Spatial').rename(columns={
        'Mine host country': 'Country', 'Exposed species count': 'num_species',
        'Unioned exposed range area (km²)': 'overlap_area_km2',
        'Influence area in non-zero RWR cells (km²)': 'area_total_rwr_km2',
        'Rwr 0 50 area (km²)': 'area_le50_km2', 'Rwr 50 75 area (km²)': 'area_50_75_km2',
        'Rwr 75 90 area (km²)': 'area_75_90_km2', 'Rwr 90 99 area (km²)': 'area_90_99_km2',
        'Rwr 99 100 area (km²)': 'area_99_100_km2',
        'Classified influence area in top-10% RWR cells (%)': 'pct_ge90',
    })


def format_area(value, pos=None):
    if abs(value) >= 1000000:
        return f'{value / 1000000:.1f}M'
    if abs(value) >= 1000:
        return f'{value / 1000:.0f}k'
    return f'{value:.0f}'


def set_asinh_area_scale(ax, ymax):
    width = 10000.0
    ax.set_yscale('function', functions=(lambda x: np.arcsinh(np.asarray(x, dtype=float)/width),
                                       lambda x: np.sinh(np.asarray(x, dtype=float))*width))
    candidates = np.array([0, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000], dtype=float)
    ticks = candidates[candidates <= ymax * 1.02]
    if not len(ticks) or ticks[0] != 0:
        ticks = np.insert(ticks, 0, 0.)
    if ticks[-1] < ymax * .8:
        magnitude = 10 ** np.floor(np.log10(ymax))
        ticks = np.append(ticks, np.ceil(ymax/magnitude)*magnitude)
    ax.set_yticks(np.unique(ticks))


def country_inputs():
    return SimpleNamespace(country_area_data=country_area_data, COUNTRY_FIELD='Country',
                           OVERLAP_AREA_FIELD='overlap_area_km2', TOTAL_AREA_FIELD='area_total_rwr_km2',
                           PCT_TOP10_FIELD='pct_ge90', MIN_OVERLAP_AREA_KM2=5000.,
                           AREA_COLS=['area_le50_km2', 'area_50_75_km2', 'area_75_90_km2', 'area_90_99_km2', 'area_99_100_km2'],
                           BIN_COLORS=['#f0f0f0', '#bfd3e6', '#8c96c6', '#f6cd8e', '#cb1c2a'],
                           format_area=format_area, set_asinh_area_scale=set_asinh_area_scale)


def classify_three(value, breaks):
    if pd.isna(value):
        return -1
    return 0 if value <= breaks[0] else 1 if value <= breaks[1] else 2


def map_inputs():
    return SimpleNamespace(MAIN_FIXED_AREA_BREAKS_KM2=[1000., 10000.], PERCENT_BREAKS=[20., 50.],
                           MAIN_PALETTE=np.array([['#e8e8e8', '#b5c0da', '#6c83b5'],
                                                  ['#b8d6be', '#90b2b3', '#567994'],
                                                  ['#73ae80', '#5a9178', '#2a5a5b']]),
                           make_area_breaks=lambda values, breaks: np.array(breaks, dtype=float),
                           classify_three=classify_three)
