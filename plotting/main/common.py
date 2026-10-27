"""Shared input paths, publication typography and figure export."""
from pathlib import Path
import argparse
import sys
sys.dont_write_bytecode = True
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.text import Text
from dataclasses import dataclass


@dataclass(frozen=True)
class FigureSpec:
    number: int
    stem: str
    title: str
    width_mm: float
    height_mm: float
    pending_revision: bool = False


FIGURES = {
    1: FigureSpec(1, "Fig1_GlobalPatterns_RangeOverlap",
                  "Global distribution of copper mines and potential impacts on vertebrate habitats",
                  165.0, 215.0),
    2: FigureSpec(2, "Fig2_RarityWeightedRichness_Randomization",
                  "Copper mine occurrence in areas of high rarity-weighted richness",
                  165.0, 90.0),
    3: FigureSpec(3, "Fig3_ConservationPriority_Composition",
                  "Conservation-priority species across vertebrate groups and habitat area classes",
                  165.0, 115.0),
    4: FigureSpec(4, "Fig4_ConservationPriority_Distribution",
                  "Global and regional distributions of conservation-priority species",
                  183.0, 134.0),
    5: FigureSpec(5, "Fig5_MineSpecies_CountryAssociations",
                  "Country-level associations between copper mines and conservation-priority species",
                  200.0, 110.0),
    6: FigureSpec(6, "Fig6_OwnershipTraceability_CrossBorderLinks",
                  "Ownership-based tracing of cross-border ecological responsibility for habitats affected by copper mining",
                  165.0, 250.0, pending_revision=False),
    7: FigureSpec(7, "Fig7_OwnershipPortfolios_Biodiversity",
                  "Conservation priorities within copper mine ownership portfolios",
                  165.0, 110.0),
}

# Working proof layout only: left, bottom, width, height in millimeters.
# Three full-width rows, A then B then C; thresholds reduce the label count.
FIGURE_6_PANELS_MM = {
    "A": (14.0, 187.0, 145.0, 56.0),
    "B": (19.402947362065696, 102.0, 134.1941052758686, 80.0),
    "C": (19.402947362065696, 10.0, 134.1941052758686, 84.0),
}


CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR.parent))
from data_paths import data_root_from_arguments, external_output_root

ZENODO = data_root_from_arguments()
INPUTS = ZENODO / '9Figure_Inputs/Main' if ZENODO is not None else None

STYLE = {'font.family': 'sans-serif',
 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
 'font.size': 6.2,
 'axes.labelsize': 6.5,
 'axes.titlesize': 7.0,
 'xtick.labelsize': 5.8,
 'ytick.labelsize': 5.8,
 'legend.fontsize': 5.8,
 'axes.linewidth': 0.65,
 'xtick.major.width': 0.55,
 'ytick.major.width': 0.55,
 'xtick.major.size': 2.5,
 'ytick.major.size': 2.5,
 'axes.spines.right': False,
 'axes.spines.top': False,
 'axes.unicode_minus': True,
 'legend.frameon': False,
 'svg.fonttype': 'none',
 'pdf.fonttype': 42,
 'ps.fonttype': 42,
 'savefig.facecolor': 'white',
 'figure.facecolor': 'white'}

LAYOUTS = {'6': {'canvas_pixels': [1763, 2781],
       'canvas_mm': [149.27727777777778, 235.51402951388883],
       'panel_content_bounds_mm': [[2.5,
                                    172.71014062499998,
                                    146.77727777777773,
                                    233.01402951388886],
                                   [2.5,
                                    89.58307031249998,
                                    146.77727777777775,
                                    169.710140625],
                                   [2.5,
                                    2.499999999999999,
                                    146.77727777777775,
                                    86.58307031250001]],
       'frames_mm': [[12.377277777777774,
                      173.87430729166667,
                      146.77727777777773,
                      229.87430729166664],
                     [12.377277777777726,
                      89.710140625,
                      146.77727777777775,
                      169.710140625],
                     [12.377277777777726,
                      2.5830703124999967,
                      146.77727777777775,
                      86.5830703125],
                     [26.377277777777767,
                      185.87430729166667,
                      41.37727777777776,
                      200.87430729166664]]},
 '7': {'canvas_pixels': [1881, 1048],
       'canvas_mm': [159.27899561849085, 88.77260684931504],
       'panel_content_bounds_mm': [[2.5,
                                    2.499999999999997,
                                    85.29652585298062,
                                    86.27260684931505],
                                   [90.29652585298061,
                                    2.538412857900415,
                                    156.77899561849082,
                                    86.27260684931504]],
       'frames_mm': [[23.174041666666668,
                      12.606377777777775,
                      85.29652585298062,
                      76.40637777777779],
                     [103.53577585298062,
                      12.606377777777775,
                      145.29023243722446,
                      76.40637777777779],
                     [146.81783450737973,
                      24.728377777777833,
                      148.79563450737973,
                      64.28437777777773]]}}

def configure_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(STYLE)

def new_figure(number: int):
    import matplotlib.pyplot as plt

    configure_style()
    spec = FIGURES[number]
    return plt.figure(figsize=(spec.width_mm / 25.4, spec.height_mm / 25.4), dpi=600)

def panel_label(ax, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label.lower(), transform=ax.transAxes, ha="left", va="top",
            fontsize=8.0, fontweight="bold", clip_on=False, zorder=100,
            gid="panel-label")

def fix_geographic_aspect(ax, extent, *, global_map=False):
    """Lock map scales after plotting, independent of allocated panel shape.

    Global maps use equal angular scales (Plate Carree); regional maps retain
    the original central-latitude cosine correction. This is a display rule,
    not a claim of globally equal ground distances or equal mapped areas.
    """
    import math

    xmin, xmax, ymin, ymax = map(float, extent)
    if not (xmax > xmin and ymax > ymin):
        raise ValueError("Geographic limits must have positive spans")
    latitude = (ymin + ymax) / 2
    aspect = 1.0 if global_map else 1.0 / max(math.cos(math.radians(latitude)), 0.25)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect(aspect, adjustable="box", anchor="C")
    ax._geographic_scale_contract = {
        "kind": "global_equal_angular_scale" if global_map else "regional_central_latitude_scale",
        "expected_y_over_x": aspect,
        "extent": [xmin, xmax, ymin, ymax],
        "central_latitude": latitude,
    }

def require_columns(frame, columns, label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}")

def read_csv(path: Path, columns=(), numeric=()):
    """Load existing tabular results; missing or nonfinite measures are errors."""
    import numpy as np
    import pandas as pd

    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    frame = pd.read_csv(path, encoding="utf-8-sig")
    require_columns(frame, (*columns, *numeric), path.name)
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not np.isfinite(frame[column].to_numpy(dtype=float)).all():
            raise ValueError(f"Nonfinite values in {path.name}: {column}")
    return frame

def require_expected(actual, expected, label: str) -> None:
    if actual != expected:
        raise ValueError(f"{label}: expected {expected}, found {actual}; review the data version")

def table(number, sheet):
    frame = pd.read_excel(ZENODO / '8Final_Results' / f'Supplementary_Data_{number}.xlsx', sheet_name=sheet)
    # Normalize published headers only for internal field lookup.
    frame.columns = [str(c).replace('country or region', 'country') for c in frame.columns]
    return frame


def fixed_canvas(number):
    spec = LAYOUTS[str(number)]
    configure_style()
    width, height = spec['canvas_pixels']
    fig = plt.figure(figsize=((width+1e-6)/300, (height+1e-6)/300), dpi=300)
    return fig, spec


def axes_from_frame(fig, spec, frame, *, physical=False):
    width, height = fig.get_size_inches()*25.4 if physical else spec['canvas_mm']
    x0, y0, x1, y1 = frame
    return fig.add_axes([x0/width, y0/height, (x1-x0)/width, (y1-y0)/height])


def run(number, builder):
    parser = argparse.ArgumentParser(description=f'Reproduce manuscript Figure {number}.')
    parser.add_argument('--data-root', type=Path, help='Dataset root; overrides COPPER_DATA_ROOT and bundled data/')
    parser.add_argument('--output-dir', type=Path, required=True, help='Directory outside the code and data packages')
    parser.add_argument('--format', choices=['png', 'pdf', 'svg'], default='png')
    args = parser.parse_args()
    if ZENODO is None:
        parser.error('Supply --data-root or set COPPER_DATA_ROOT')
    try:
        args.output_dir = external_output_root(args.output_dir, ZENODO)
    except ValueError as error:
        parser.error(str(error))
    fig, _ = builder()
    for text in fig.findobj(match=Text):
        if text.get_fontsize() < 5.5:
            text.set_fontsize(5.5)
    if number not in (6, 7):
        from layout_alignment import align_figure
        align_figure(fig, number)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / f'Figure_{number}.{args.format}'
    fig.savefig(path, dpi=300, bbox_inches=None)
    if number == 5 and args.format == 'png':
        from PIL import Image
        with Image.open(path) as image:
            rgb = image.convert('RGB')
        rgb.save(path, dpi=(300, 300))
        rgb.close()
    plt.close(fig)
    print(path)
