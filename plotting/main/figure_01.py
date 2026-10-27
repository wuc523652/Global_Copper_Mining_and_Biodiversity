"""Figure 1: two global maps and country-level overlap-area composition."""

import sys
sys.dont_write_bytecode = True

from common import new_figure, require_expected
import global_patterns as panels
from plot_inputs import richness_inputs, overlap_inputs, country_inputs


def build():
    richness = richness_inputs()
    overlap = overlap_inputs()
    country = country_inputs()
    fig = new_figure(1)
    grid = fig.add_gridspec(3, 1, height_ratios=[1.1, 1.1, 0.75],
                           left=0.085, right=0.89, bottom=0.11, top=0.955, hspace=0.25)
    axes = [fig.add_subplot(grid[index, 0]) for index in range(3)]
    qa = {
        "A": panels.draw_richness_map(axes[0], richness),
        "B": panels.draw_overlap_map(axes[1], overlap),
        "C": panels.draw_country_area(axes[2], country),
    }
    require_expected(qa["A"]["mine_count"], 1278, "Figure 1 mine count")
    require_expected(qa["B"]["species_total"], 891, "Figure 1 overlapping species")
    require_expected(qa["C"]["displayed_country_count"], 51, "Figure 1 displayed countries")
    return fig, {"panels": qa, "country_selection": "overlap_area_km2 >= 5000; Global excluded",
                 "panel_C_percentage_denominator": "area_total_rwr_km2"}


if __name__ == '__main__':
    from common import run
    run(1, build)
