# Recompute the manuscript's numerical results

Run from the code-package root, after installing its scientific dependencies:

```console
python reproduction/recompute.py --output-root /path/to/new_results
```

The default input is the code package's `data/` directory. Use `--data-root`
to select a different complete dataset. This route uses real research
checkpoints. Use a new or empty output directory outside
both the code package and data. Input files
are never overwritten. The command exits nonzero if a calculation or comparison
fails. `--steps spatial randomization priority ownership portfolios distance complete-case`
selects individual tasks; required preceding tasks are included automatically.

This is the **archived-spatial-checkpoint route**, intended for a practical,
repeatable regeneration of statistical results. It is not a claim that original
IUCN range polygons and all spatial intersections have been processed again.
Use the package's raw-input pipeline for the original spatial operations.

## What is recalculated

| Task | Computation | Inputs |
|---|---|---|
| `spatial` | RWR quantile boundaries, bin classification, country/global area sums and percentages; affected-species counts | Full global grid values, cell-level intersection areas, country/species intersections, archived country/global spatial-union areas |
| `randomization` | All 9,999 country-stratified draws, observed/background statistics, null summaries, empirical probabilities | Ordered background-grid values and mine-grid observations; seed 20260728 |
| `priority` | Exposure fractions, main thresholds, 203-species selection, alternative-rule and separate cap sensitivity, boundary diagnostics | Species global range areas, globally deduplicated exposure areas, country/species links |
| `ownership` | Country-balanced mine attribution, investor/host/bilateral summaries, unknown-equity accounting | Row-level normalized ownership ledger, coverage, and spatial exposure tables |
| `portfolios` | Mine priority categories, mine/species/country association summaries, equity-weighted country and bilateral portfolios | Recomputed priority list, positive-area mine/species links, mine-grid observations, row-level ownership |
| `distance` | Screening and ownership calculations at 10/20/30/40/50 km; priority-set and rank stability; marginal screening yield | Each distance's spatial exposure areas; five buffer GPKGs are freshly unioned for area/yield diagnostics |
| `complete-case` | 959-mine selection, new country/global unions of those mines' exposure geometries, attribution and rank stability | Mine-exposure geometries and normalized ownership records; no missing-share imputation |

The implementation calls the original `analysis/05`–`analysis/16` functions
where applicable. It does not use the submitted Excel workbooks or
`11Reanalysis_Checkpoints/expected/` as scientific inputs. Those files are read
only after recalculation to compare results. The optional 1,000,000-km² cap
check in step 08 is enabled; manuscript main and alternative classifications
are also checked explicitly using their **uncapped** thresholds.

## Outputs and verification

- Analysis-specific CSV tables and QA records under directories numbered 05–16.
- `publication_numeric_tables/`: recalculated columns using the submitted
  supplementary workbooks' column labels. These are numeric verification
  tables, not replacements for the workbooks' formatting or data dictionaries.
- `comparison_report.json`: each checked table, exact identifier-set comparison,
  numeric columns, tolerance, largest difference, and PASS/FAIL. Usual tolerances
  are absolute `1e-7` and relative `1e-10`; randomization iterations use absolute
  `1e-10`. Species sets and table row counts must match exactly. This does not
  apply a general percentage-error allowance.

The report checks recomputed rows against the archived scientific tables and
the submitted Supplementary Data 3–9 workbooks, including all 9,999 iterations,
all 203 priority species, five distance-specific species sets, 78 owner
registrations, 77 host countries, cross-border pairs and portfolio values.
The species and mine inventories (Supplementary Data 1–2) are foundational
inputs here; re-derivation of their geometry/metadata belongs to the raw route.

Plot rendering uses the archived publication tables so the submitted figure
layout is retained. It does not silently replace those plotting inputs with
new CSV outputs. Successful calculation comparisons establish the numerical
agreement for the columns listed in the report; unlisted columns are not
asserted to have been checked.

The provenance manifest in `11Reanalysis_Checkpoints/` distinguishes derived
spatial inputs from comparison-only expected outputs. Range areas, primary
spatial unions, cell intersections and mine/species intersection links remain
archived computational starting points in this route. Recalculating them from
raw licensed geometries is a separate and more expensive workflow.
