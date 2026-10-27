"""Figure 3: stack exact species counts in all six panels."""
import sys
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import configure_style, table
from priority_panels import GROUP_ORDER, RANGE_ORDER, draw_composition


def composition_counts():
    species = table(2, 'Species')
    exposed = table(4, 'Exposed_Species')
    priority = table(4, 'Priority_Species')
    id_col = 'IUCN species ID'
    for frame, expected in [(species, 5890), (exposed, 891), (priority, 203)]:
        if len(frame) != expected or frame[id_col].duplicated().any():
            raise ValueError('Species counts or unique identifiers do not match the approved data')
    if not set(priority[id_col]) <= set(exposed[id_col]) <= set(species[id_col]):
        raise ValueError('Priority/exposed species are not nested in the species inventory')
    data = species.assign(exposed=species[id_col].isin(exposed[id_col]),
                          priority=species[id_col].isin(priority[id_col]))
    rows = []
    for group in GROUP_ORDER:
        subset = data if group == 'All' else data[data['Taxonomic group'].eq(group)]
        for size in RANGE_ORDER:
            part = subset[subset['Range size class'].eq(size)]
            rows.append({'Group': group, 'Range': size,
                'n_unaffected_species': int((~part.exposed).sum()),
                'n_affected_non_priority_species': int((part.exposed & ~part.priority).sum()),
                'n_priority_species': int(part.priority.sum())})
    return pd.DataFrame(rows)


def build():
    counts = composition_counts()
    columns = ['n_unaffected_species', 'n_affected_non_priority_species', 'n_priority_species']
    expected = [[1668, 0, 65], [3263, 444, 138], [68, 244, 0]]
    if counts[counts.Group.eq('All')][columns].values.tolist() != expected:
        raise ValueError('Figure 3a counts differ from the approved source counts')
    configure_style()
    fig = plt.figure(figsize=(165 / 25.4, 115 / 25.4), dpi=600)
    slot = fig.add_gridspec(1, 1, left=0.09, right=0.975, bottom=0.10, top=0.77)
    axes = draw_composition(fig, slot[0, 0], counts)
    for ax, group in zip(axes, GROUP_ORDER):
        block = counts[counts.Group.eq(group)].set_index('Range').loc[RANGE_ORDER]
        for index in range(3):
            bars = [p for p in ax.patches if hasattr(p, 'get_x')
                    and np.isclose(p.get_x(), index-0.62/2) and p.get_transform() != ax.transAxes]
            actual_top = max(p.get_y()+p.get_height() for p in bars)
            if not np.isclose(actual_top, block.iloc[index][columns].sum()):
                raise ValueError(f'Incorrect stacked count: {group}/{RANGE_ORDER[index]}')
    return fig, {'counts': counts.to_dict('records')}


if __name__ == '__main__':
    from common import run
    run(3, build)
