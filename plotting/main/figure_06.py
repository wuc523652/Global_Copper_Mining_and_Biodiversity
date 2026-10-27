"""Figure 6: registration map and two complete, data-driven alluvial panels."""
import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from common import INPUTS, table, fixed_canvas, axes_from_frame
from plot_inputs import geography
from ownership_panels import draw_bivariate_map
from sankey import draw_alluvial

def flows():
    data = table(7, 'CrossBorder_Flows').rename(columns={
        'Owner registration iso3': 'investor_iso3', 'Owner registration country': 'investor_country',
        'Mine host iso3': 'resource_iso3', 'Mine host country': 'resource_country',
        'Attributed exposure area (km²)': 'attributed_exposure_area_km2', 'Num mines': 'num_mines'})
    owners = table(6, 'Owner_Countries').set_index('Owner registration iso3')['Owner registration continent']
    hosts = table(6, 'Mine_Host_Countries').set_index('Mine host iso3')['Mine host continent']
    data['investor_continent'] = data.investor_iso3.map(owners)
    data['resource_continent'] = data.resource_iso3.map(hosts)
    data['is_domestic'] = data.investor_iso3.eq(data.resource_iso3).astype(int)
    if data.is_domestic.any() or data[['investor_continent', 'resource_continent']].isna().any().any():
        raise ValueError('Cross-border flows or continent metadata are inconsistent')
    return data



def build():
    data = flows()
    registrations = table(6, 'Owner_Countries').rename(columns={
        'Owner registration iso3': 'investor_iso3', 'Total attributed area (km²)': 'total_attributed_area_km2',
        'Cross border share (%)': 'foreign_project_share_pct'})
    # The map uses the fixed set shown in the manuscript; flow panels retain
    # their own explicitly selected registrations and mine hosts.
    selection = pd.read_csv(INPUTS / 'figure_6a_registration_selection.csv').investor_iso3
    registrations = registrations[registrations.investor_iso3.isin(selection)]
    if len(registrations) != len(selection):
        raise ValueError('An approved map registration is missing from Data 6')
    fig, spec = fixed_canvas(6)
    axes = [axes_from_frame(fig, spec, frame, physical=index < 2)
            for index, frame in enumerate(spec['frames_mm'][:3])]
    area_breaks, share_breaks, mapped = draw_bivariate_map(axes[0], geography('figure_6_land'), registrations)
    owner_model = draw_alluvial(axes[1], data, 'owner')
    host_model = draw_alluvial(axes[2], data, 'host')
    expected = {'MEX': 83538, 'CHL': 72262, 'PER': 58660, 'COD': 28236, 'CHN': 22679,
                'BRA': 18677, 'PHL': 14557, 'ZAF': 12925, 'ECU': 10062}
    if {k: round(v) for k, v in host_model['country_totals_km2']['left'].items()} != expected:
        raise ValueError('Figure 6c differs from the approved source totals')
    if round(host_model['country_totals_km2']['right']['HKG']) != 18868:
        raise ValueError('Unexpected HKG subtotal')
    if not np.isclose(host_model['displayed_area_km2'], 321596.277248, atol=0.01, rtol=0):
        raise ValueError('Unexpected Figure 6c area total')
    width, height = spec['canvas_mm']
    for label, bounds in zip('abc', spec['panel_content_bounds_mm']):
        w, h = fig.get_size_inches()*25.4 if label in 'ab' else (width, height)
        fig.text(bounds[0]/w, bounds[3]/h, label, fontsize=7, fontweight='bold', ha='left', va='top')
    return fig, {'mapped_registrations': mapped, 'owner_flows': owner_model, 'host_flows': host_model}


if __name__ == '__main__':
    from common import run
    run(6, build)
