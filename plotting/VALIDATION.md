# Figure reproduction validation — 18 October 2026

All seven main figures and all nine supplementary figures were regenerated
successfully using only this code package and the specified Zenodo data root.
Outputs were written to a separate validation directory. The original figure
and table files were not modified. A single sequential PNG run took about
236 seconds in the validation environment; hardware and storage affect timing.

The 29 newly archived drawing assets were checked with SHA-256 against the
original inputs: 29/29 were byte-identical. All nine archived source workbooks
were also byte-identical to those used by the original plotting scripts.

## Comparison against the supplied PNGs

All 16 reproduced images have exactly the same pixel dimensions as their
supplied counterparts. Twelve images have identical pixel arrays. The following
four differ:

| Figure | Differing pixels | Total pixels | Largest channel difference (0–255) |
|---|---:|---:|---:|
| Main 6 | 31,310 | 4,902,903 | 255 |
| Main 7 | 9,862 | 1,971,288 | 255 |
| Supplementary 4 | 7,960 | 2,804,368 | 2 |
| Supplementary 5 | 25,100 | 3,985,597 | 254 |

For each of these four figures, the original, unmodified plotting script was
also rerun against its original input paths in the same environment. Its PNG
was byte-for-byte identical to the new package's PNG. Thus the four differences
already exist between the supplied PNGs and outputs from the supplied source
code under this environment; they were not introduced by the packaging or input
path changes. Their historical cause has not been established. Scientific data,
screening rules, numerical checks, selected categories, colours and layouts
were not changed to force a pixel match.

Main Figure 3 and Supplementary Figure 9 have identical pixels but different
PNG file bytes, which does not change their visual contents. Other matching
figures are also byte-identical.

## Export and runtime checks

- Main Figure 2 and Supplementary Figure 9 were additionally exported as PDF
  and SVG. PDF text audits found minimum text sizes of 5.5 pt and 6.5 pt,
  respectively, with no text runs below 5 pt. These are representative vector
  export checks, not a full journal typography audit of every figure.
- The shared entry point stops on a failed subprocess and records per-figure
  completion and timings in `figure_render_report.json` in the output directory.
- The retained source-data assertions passed, including species counts,
  randomization iterations, ownership rankings and buffer-distance diagnostics.
- Supplementary Figure 1 emits the original geographic-CRS centroid warning
  and a non-finite-centroid warning during drawing. Its complete rendered image
  is byte-identical to the supplied PNG. The centroids are display positions;
  this drawing function does not calculate exposure areas.
- The source preflight cannot identify dynamically selected export suffixes
  in the shared export function; actual PNG/PDF/SVG generation was checked
  instead of changing the existing export or scientific drawing code.

Validation environment: Windows, Python 3.13.3, NumPy 2.2.5, pandas 2.2.3,
Matplotlib 3.10.1, GeoPandas 1.0.1, pyogrio 0.10.0, Shapely 2.1.0,
pyproj 3.7.1, openpyxl 3.1.5, and Pillow 11.1.0. The original requested fonts
were available. Reproduction on other environments should preserve numerical
content but can differ in font metrics, rasterization or image metadata.
