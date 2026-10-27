# Local data layout validation — 2026-10-27

This update copies the complete local companion dataset into the code package's `data/` directory and makes it the default input for the unified runner, raw-analysis runner, statistical reanalysis and figure entry points. Explicit external data paths remain supported. The original dataset directory is retained. Scientific algorithms, parameters and data bytes are unchanged.

## Completed checks

| Check | Result and scope |
| --- | --- |
| Complete copy | 229 files, 3,050,985,840 bytes, including the original data README and manifest. Both source and destination pass verification of the 228 listed files. The manifests are byte-identical. |
| Manifest SHA-256 | `e0b1d90250b4eb8a6c94f4d465672b4cccff314cf0c46efb6f97cd5bd510ee48` |
| Default unified input | `python reproduce.py --mode check` succeeds without `--data-root`. |
| Raw-pipeline preflight | All 17 stages' raw files, shapefile sidecars and planned prerequisites pass without `--data-root`; no full-pipeline computation is claimed. |
| Real raw-input execution | Stage 10 reruns from `data/10Raw_Inputs`, producing the 1,278-mine master, 1,491-row ownership ledger and 1,278-mine coverage table. All three CSVs match the archived checkpoints after the existing display-name normalization, with absolute and relative numerical tolerance `1e-10`. |
| Real checkpoint reanalysis | `reproduction/recompute.py --steps priority` uses the default `data/`; all 13 comparisons pass, including the 203-species priority selection. |
| Real plotting input | `plotting/render_all.py --figures 2` uses the default `data/` and successfully generates the main Figure 2 PNG. No pixel-comparison claim is made for this check. |
| Path regression tests | All 10 checks in `test_local_data_layout.py` pass: default inputs from another working directory, explicit path/environment overrides, raw-pipeline dispatch and protection against writing into bundled inputs or source directories. Mock dispatch tests do not execute the full raw pipeline. |
| Existing comparator tests | All 3 tests in `reproduction/test_comparison.py` pass. |

Machine-readable evidence is in `local_data_layout_20261027.json`. Outputs and execution logs remain in the local audit workspace, outside the code package and dataset.

## Environment and limits

Real-data checks used the existing Windows/Python 3.13.3 validation environment backed by the `geo` Conda environment. Launching its Python directly initially failed to locate PROJ/GDAL resources. After explicitly setting `PROJ_DATA` and `GDAL_DATA` to that same environment's resource directories, stage 10 succeeded in a separate output directory. No global environment configuration or scientific code was changed to resolve this. Activate the appropriate Conda environment before normal execution.

The entire 17-stage raw pipeline was **not** rerun. The complete statistical suite and all 16 figures were also **not** rerun for this path-only update; their historical validation remains documented in `README.md` and the 2026-10-18 reports. The `raw` route does not automatically export newly calculated results into every formatted supplementary workbook or rewire all plots to those new tables.

The complete `data/` directory remains ignored by ordinary Git and is intended for matching GitHub Release assets. No Git commit, push, Release publication, Zenodo upload or new DOI was created by this update. Original third-party data conditions remain applicable.

## Subsequent removal of optional examples

On 2026-10-27, the complete `demo/` directory (four files) was removed from the distribution and retained in a separate local backup. No scientific entry point depends on those files. After removal, the 228-file data integrity check, all-stage input preflight and all 10 path regression tests passed again. Documentation and the code inventory were refreshed. Scientific algorithms and the complete dataset were unchanged; no remote publication was performed.
