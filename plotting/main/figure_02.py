"""Figure 2: global randomization and country background comparisons."""

import sys
sys.dont_write_bytecode = True

from types import SimpleNamespace

from common import ZENODO, new_figure, require_expected
import global_patterns as panels


def build():
    source_dir = ZENODO / "6Randomization"
    # Read frozen results directly; do not import or rerun the randomization program.
    inputs = SimpleNamespace(
        OUT_ITERATIONS=source_dir / "Randomization_Iterations.csv",
        OUT_SUMMARY=source_dir / "Randomization_Summary.csv",
        OUT_COUNTRY_DESC=source_dir / "Country_Background_Summary.csv",
        COUNTRY_FIELD="Country", COUNTRY_MIN_MINES_FOR_FIG=10,
    )
    fig = new_figure(2)
    grid = fig.add_gridspec(1, 2, width_ratios=[0.90, 1.55],
                           left=0.095, right=0.97, bottom=0.34, top=0.93, wspace=0.45)
    qa_a = panels.draw_randomization(fig.add_subplot(grid[0, 0]), inputs)
    qa_b = panels.draw_country_comparison(fig.add_subplot(grid[0, 1]), inputs)
    require_expected(qa_a["iterations"], 9999, "Figure 2 randomizations")
    require_expected(qa_b["displayed_country_count"], 17, "Figure 2 country comparisons")
    return fig, {"panels": {"A": qa_a, "B": qa_b},
                 "interval": "2.5th–97.5th percentiles of country-stratified randomizations",
                 "country_comparisons": "descriptive; no country-level significance tests"}


if __name__ == '__main__':
    from common import run
    run(2, build)
