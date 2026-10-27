# Validation of the 2026-10-18 local release

For the later 2026-10-27 copy into `data/` and default-path update, see
[`local_data_layout_20261027.md`](local_data_layout_20261027.md). The historical
execution results below retain their original scope and date.

This directory records completed executions, not a promise that every supported platform or the full raw pipeline has been run. Personal absolute paths in copied JSON reports are replaced by `<validation-workspace>`, `<zenodo-data-root>` and `<code-root>`; numerical results are unchanged.

| Execution | Actual result and scope |
| --- | --- |
| Unified `reproduce.py --mode all` | Exit 0 after verifying 228 deposited files, recomputing all seven groups and rendering all 16 figures. The final PNGs match the first validated run byte for byte; see `unified_run.json`. |
| Seven reanalysis groups | PASS, 50 checks; 24 directly compare declared fields with final Supplementary Data 3-9. See `calculation_comparison.json` and `reproduction/publication_checks.py`. |
| Randomization | All 9,999 seeded iterations and both test summaries recomputed and compared. |
| Priority rules | Primary 203-species membership, alternative rules and diagnostics; all five distance-specific priority sets checked exactly. |
| Ownership and portfolios | Attribution and ecological portfolios independently recomputed. Unknown equity remains unattributed. |
| Complete case | 959-mine subset selected again, country/global exposure unions recalculated from archived geometry, followed by fresh attribution. |
| Distance diagnostics | Five buffer geometries re-unioned; area increments and priority-species efficiency recomputed. Underlying species intersection areas are archived inputs in this shorter path. |
| Figure rendering | 7 main and 9 supplementary PNGs ran successfully; dimensions all match. 12 images pixel-identical to archived PNGs. |
| Four remaining PNGs | Figures 6, 7, S4 and S5 exactly match runs of the unmodified original drawing sources in the same environment. Differences from previously archived PNGs were not introduced by migration; the historical cause is not established. |
| Vector figures | PNG/PDF/SVG entry points exercised for main Figure 2 and supplementary Figure 9. |
| Raw-input stage 10 | Full mine master, ownership ledger and coverage regenerated from the copied original inputs. All three tables match checkpoints after the existing HKG display-name normalization; see `raw_ownership_comparison.json`. |
| Raw-IUCN stage 17 | Actual spatial intersections rerun for all 203 priority species, using the archived base grid and validated priority list. All 1,404,063 grid rows and 7,783 association records match; 4,771 cells are nonzero. This does not rerun stages 01 or 07. |
| Full 17-stage preflight | All required raw files, shapefile components and ordered prerequisites present. |
| Failure checks | Missing/changed/unlisted data are rejected; output directories overlapping code/data are rejected. Stage 07 propagates scientific failures. Missing expected numerical columns and other numerical-comparison regressions are tested in `reproduction/test_comparison.py`. |

The full original-input workflow was **not completed end to end** during this packaging update. A trial of the complete raw route was stopped during the expensive base-grid construction stage. Its partial files are not used as validated results and are outside the release. Therefore this release establishes successful statistical reanalysis, all-figure rendering, raw ownership preparation and priority-grid recalculation; it does not establish a fresh full-global rerun of all IUCN intersections, RWR geometry calculations or all buffer-distance spatial intersections.

The shorter calculation route reads archived spatial checkpoints. Drawing reads the archived final workbooks, not the newly generated numerical CSVs. The reports specify which workbook columns were checked; they do not claim byte-identical regeneration of every Excel workbook, its formatting or its data dictionary.

Environment: Windows, Python 3.13.3 and library versions in `environment_versions.json`. The validation environment inherits installed geospatial libraries; this update did not test a new online installation or other operating systems. The numerical comparator uses declared tolerances (typically absolute `1e-7` and relative `1e-10`); see each JSON check. Species sets and association/grid comparisons are exact as documented. Existing HKG labels are standardized without changing codes or values.

Code review and Python review were completed. Type checks for the five new orchestration/calculation modules and three comparison regression tests passed. A Figure 6 help-path error and a suppressed stage-07 exception were corrected; no scientific formula was changed.
