# Reproduce all manuscript figures

The plotting code regenerates all 7 main figures and 9 supplementary figures
from the complete research dataset in `data/`. No previously rendered manuscript
image is used to construct a data panel. Numerical values, display selections,
colours, dimensions and panel layouts are retained from the supplied figures.

## Run

Install the repository's root `requirements.txt` in the documented environment.
By default, input comes from the code package's `data/`, containing
`8Final_Results` and `9Figure_Inputs`. To use a separately stored dataset,
set `--data-root` to that complete data directory.

From the repository root:

```text
python plotting/render_all.py --output-root /path/to/reproduced_figures
```

On Windows, paths containing spaces must be quoted. `COPPER_DATA_ROOT` may be
set instead of passing `--data-root`; the priority is explicit argument,
environment variable, then the code package's `data/` directory.
An output directory is always required and must be outside both downloaded
packages. Existing figures of the same names in that output directory will
be replaced; input data and supplied manuscript figures are not changed.

```text
python plotting/render_all.py --data-root /path/to/data --output-root /path/to/figures --figures main
python plotting/render_all.py --data-root /path/to/data --output-root /path/to/figures --figures supplementary
python plotting/render_all.py --data-root /path/to/data --output-root /path/to/figures --figures 1 7 S2 S9 --formats png pdf svg
```

The output contains `Main/Figure_1.png`–`Figure_7.png`,
`Supplementary/Supplementary_Figure_1.png`–`Supplementary_Figure_9.png`, and
`figure_render_report.json` with per-figure timings and return codes. Selection
is sequential and each figure runs in a separate Python process so settings
and memory do not leak between plots. A failed figure stops the run and leaves
a report with `complete: false`; inspect the error, then rerun that selection.
Main figures are rebuilt separately for each requested export format.

## Inputs and interpretation

| Data directory | Used for |
|---|---|
| `8Final_Results/Supplementary_Data_1.xlsx`–`9.xlsx` | Published numerical results and figure source tables |
| `9Figure_Inputs/Main` | Main-figure display boundaries, coordinates, region definitions and icons |
| `9Figure_Inputs/Supplementary` | Supplementary-figure display boundaries |
| `3Grid_and_Biodiversity` | Complete biodiversity grid and priority-species grid |
| `6Randomization` | Frozen 9,999-iteration randomization results |
| `7Ownership_Attribution` | Mine–owner ledger used for Figure 7 |

Plotting reads the archived final results and does not rerun the global spatial
analysis. The repository's analysis workflow is a separate step for rebuilding
and checking those results. See the main and supplementary README files for
panel-level inputs and selection rules. The fixed region extents, display
selections, and diagram text are presentation specifications, not new estimates.

Supplementary Figure 1 reads all 1,404,410 grid cells and has the largest memory
and time requirements. No sampling is introduced to reduce its workload.

## Typography and rendering environment

The original fonts are Arial and Microsoft YaHei for the main figures;
Times New Roman for Supplementary Figures 1–6; DejaVu Sans for Supplementary
Figure 7; and Arial for Supplementary Figures 8–9. Install these fonts separately
under their applicable terms. They are not redistributed here. Matplotlib uses
available fallbacks if a requested font is missing, which can change line breaks,
canvas crops and text spacing. Matplotlib, Pillow, FreeType and font versions can
also affect pixel-level agreement even when all scientific values are identical.

All original PNG export resolutions are retained: main figures 300 dpi;
Supplementary Figures 1 and 6–8 400 dpi; Figures 2–5 300 dpi; Figure 9 600 dpi.
PDF and SVG exports preserve editable text. The cartographic and species-range
inputs retain their third-party source conditions; consult the data-package
documentation before redistribution.
