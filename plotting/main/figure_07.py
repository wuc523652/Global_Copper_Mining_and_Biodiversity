"""Figure 7: equity-weighted portfolios and distinct-mine country-pair counts."""
import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from common import ZENODO, table, fixed_canvas, axes_from_frame


def heatmap_counts():
    owners = table(8, 'Owner_Portfolios').rename(columns={
        'Owner registration iso3': 'investor_iso3',
        'Cross border conservation priority unique mines': 'n'})
    owner_order = owners[owners.n.gt(0)].sort_values(['n', 'investor_iso3'], ascending=[False, True]).head(8).investor_iso3.tolist()
    path = ZENODO / '7Ownership_Attribution/Mine_OwnerCountry_DerivedWeights_30km_Primary_HKGIncluded.csv'
    ledger = pd.read_csv(path)
    selected = ledger[ledger.investor_iso3.ne(ledger.resource_iso3)
        & ledger.investor_iso3.isin(owner_order)
        & ledger.ecological_attention_category.ne('Neither')
        & ledger.investor_country_share_pct.gt(0)].copy()
    # Count distinct mines across the selected portfolios, not sums of owner counts.
    hosts = selected.groupby('resource_iso3').property_id.nunique().rename('n').reset_index()
    host_order = hosts.sort_values(['n', 'resource_iso3'], ascending=[False, True]).head(10).resource_iso3.tolist()
    selected['location_node'] = np.where(selected.resource_iso3.isin(host_order), selected.resource_iso3, 'Other')
    matrix = selected.pivot_table(index='investor_iso3', columns='location_node', values='property_id', aggfunc='nunique', fill_value=0)
    order = host_order + (['Other'] if 'Other' in matrix.columns else [])
    matrix = matrix.reindex(index=owner_order, columns=order, fill_value=0).astype(int)
    expected_owners = ['CAN', 'CHN', 'AUS', 'USA', 'JPN', 'GBR', 'CHE', 'HKG']
    expected_hosts = ['CHL', 'PER', 'COD', 'MEX', 'ZAF', 'PHL', 'ARG', 'AUS', 'BRA', 'ECU']
    if owner_order != expected_owners or host_order != expected_hosts:
        raise ValueError('Figure 7b selections differ from the approved complete-data rankings')
    expected = {c: {'CHL': 2, 'COD': 1, 'AUS': 1}.get(c, 0) for c in order}
    if matrix.loc['HKG'].to_dict() != expected:
        raise ValueError('Incorrect HKG mine counts')
    return matrix


def build():
    matrix = heatmap_counts()
    fig, spec = fixed_canvas(7)
    from ownership_panels import draw_portfolios
    from layout_alignment import _external_legends
    portfolio = table(8, 'Owner_Portfolios').rename(columns={
        'Owner registration iso3': 'investor_iso3', 'Owner registration country': 'investor_country',
        'Num unique mines': 'num_unique_mines',
        'Neither portfolio share (%)': 'neither_portfolio_share_percent',
        'High rwr only portfolio share (%)': 'high_rwr_only_portfolio_share_percent',
        'Priority species only portfolio share (%)': 'priority_screened_species_only_portfolio_share_percent',
        'Dual conservation priority portfolio share (%)': 'dual_ecological_attention_portfolio_share_percent',
    })
    # Panel a preserves the portfolio selection in the manuscript.
    codes = ['CHN', 'CAN', 'AUS', 'PER', 'CHL', 'ZAF', 'USA', 'JPN', 'COD', 'GBR', 'MEX', 'CHE']
    selection = pd.DataFrame({'selection_group': ['Fig5a investor portfolio']*len(codes),
                              'selection_rank': range(1, len(codes)+1), 'country_iso3': codes})
    ax_a = axes_from_frame(fig, spec, spec['frames_mm'][0], physical=True)
    if draw_portfolios(ax_a, portfolio, selection) != 12:
        raise ValueError('Figure 7a requires all twelve portfolio rows')
    ax = axes_from_frame(fig, spec, spec['frames_mm'][1])
    cax = axes_from_frame(fig, spec, spec['frames_mm'][2])
    cmap = LinearSegmentedColormap.from_list('mine_counts', ['#fffafa', '#7B2E46'])
    values = matrix.to_numpy()
    image = ax.imshow(values, cmap=cmap, aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(matrix.shape[1]), labels=matrix.columns)
    ax.set_yticks(np.arange(matrix.shape[0]), labels=matrix.index)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = int(values[row, column])
            ax.text(column, row, str(value), ha='center', va='center', fontsize=5.5,
                    color='white' if value > values.max()*0.58 else '#222222')
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.35)
    ax.tick_params(which='minor', bottom=False, left=False)
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label('Mines meeting either criterion (n)')
    ax.set_xlabel('Mine-host country')
    ax.set_ylabel('Owner-registration countries and regions')
    ax.tick_params(axis='both', length=0)
    for label in ax.get_xticklabels():
        label.set_rotation(55)
        label.set_rotation_mode('anchor')
        label.set_ha('right')
    # Position the original external portfolio legend with the same clearance.
    from matplotlib.text import Text
    for text in fig.findobj(match=Text):
        if text.get_fontsize() < 5.5:
            text.set_fontsize(5.5)
    fig.canvas.draw()
    _external_legends(fig, [[ax_a], [ax, cax]], 7)
    width, height = spec['canvas_mm']
    for label, bounds in zip('ab', spec['panel_content_bounds_mm']):
        w, h = fig.get_size_inches()*25.4 if label == 'a' else (width, height)
        fig.text(bounds[0]/w, bounds[3]/h, label, fontsize=7, fontweight='bold', ha='left', va='top')
    return fig, {'owner_order': matrix.index.tolist(), 'host_order': matrix.columns.tolist(),
                 'matrix': matrix.to_dict(orient='index')}


if __name__ == '__main__':
    from common import run
    run(7, build)
