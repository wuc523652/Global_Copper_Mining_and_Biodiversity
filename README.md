# Global Copper Mining and Biodiversity

Research code for **Identifying conservation priorities to protect threatened vertebrates from global copper mining**.

The 16 ordered scripts cover spatial exposure, conservation-priority screening, ownership attribution and sensitivity analyses. The study uses 1,278 copper mines and a frozen cohort of 5,890 threatened vertebrate species. The primary 30-km analysis identifies 891 exposed species and 203 priority species.

This script-based repository includes small runnable examples. Full study data are distributed separately. Figure styling, manuscript production and historical alternative scripts are outside its scope.

## Installation and quick start

Use 64-bit Python 3.13.3, the version used for repository validation. No GPU is needed. From the repository root:

```text
python -m venv .venv
```

Activate with `.venv\Scripts\activate` in Windows Command Prompt, `.venv/Scripts/Activate.ps1` in PowerShell, or `source .venv/bin/activate` on Linux/macOS. Then:

```text
python -m pip install -r requirements.txt
python demo/run_demo.py --output ../copper_demo/numerical
python demo/run_spatial_demo.py --output ../copper_demo/spatial
python run_analysis.py --list
```

Both demos write a JSON summary with `"status": "PASS"` and run without network access after installation. Keep outputs in separate folders as shown above.

**Tested environment:** Windows 10 build 19045, 64-bit Python 3.13.3. A fresh virtual environment was installed from `requirements.txt`; dependency consistency, both demos and Excel read/write were verified. Each demo completed in approximately two seconds including interpreter startup. Allow about 5–15 minutes for installation depending on connectivity; this is a planning estimate, not a benchmark. Linux/macOS and full-global-run resource requirements have not been validated. Global polygon operations need substantially more memory, disk space and time than the examples.

## Repository contents

```text
README.md                   Purpose, installation, data and execution instructions
LICENSE                     MIT license for the code
requirements.txt            Pinned Python dependencies
.gitignore                  Excludes environments, inputs, outputs and caches
run_analysis.py             List, check, prepare or execute analysis stages
analysis_manifest.json      Ordered stages, dependencies and paths
analysis/
  _paths.py                 Configurable input and output paths
  01_*.py ... 16_*.py        Scientific analysis scripts
demo/
  run_demo.py               Derived-species and synthetic numerical checks
  run_spatial_demo.py       Small synthetic spatial example
  data/
    exposed_species.csv
    expected_priority_ids.csv
```

## Analysis order

Run stages through `run_analysis.py`. Independent stages need not run serially. `analysis_manifest.json` records exact prerequisites. All output folders are relative to `--work-root`.

| Stage | Analysis | Output folder | Related material |
| --- | --- | --- | --- |
| 01 | Build the global 10-km equal-area land grid | `01_grid` | Spatial framework |
| 02 | Calculate area-based rarity-weighted richness (RWR) | `02_rwr` | Biodiversity and exposure analyses |
| 03 | Calculate threatened-species richness | `03_species_richness` | Supplementary Fig. 1 |
| 04 | Calculate mine-host and globally deduplicated species overlap | `04_mine_species_overlap` | Supplementary Data 3 |
| 05 | Summarize RWR classes within overlap unions | `05_rwr_overlap` | Supplementary Data 3 |
| 06 | Run country-stratified mine-location randomization | `06_randomization` | Supplementary Data 3 |
| 07 | Calculate species exposure and screen priority species | `07_priority_species` | Supplementary Data 2 and 4 |
| 08 | Compare primary and alternative screening rules | `08_priority_robustness` | Supplementary Data 4 |
| 09 | Calculate screening margins and boundary diagnostics | `09_priority_diagnostics` | Supplementary Data 4 |
| 10 | Prepare mine identifiers, ownership ledger and coverage | `10_mine_ownership` | Attribution inputs |
| 11 | Calculate mine, host-country and global exposure unions | `11_mine_exposure` | Exposure and attribution analyses |
| 12 | Attribute balanced exposure by reported ownership | `12_ownership_attribution` | Supplementary Data 6 and 7 |
| 13 | Link mines, priority species and ecological categories | `13_mine_species_links` | Supplementary Data 5 and 8 |
| 14 | Summarize ownership-linked conservation portfolios | `14_conservation_portfolios` | Supplementary Data 8 |
| 15 | Repeat analyses at 10/20/30/40/50-km buffer distances | `15_buffer_sensitivity` | Supplementary Data 9 |
| 16 | Run complete-case ownership sensitivity | `16_ownership_sensitivity` | Supplementary Data 6 and 7 |

The numerical methods and study parameters are retained; output directories are standardized. Some output filenames retain historical figure or table numbers. The mapping above gives the current supplementary numbering.

## Study data and input layout

The companion dataset is **Global Copper Mining and Biodiversity**, being deposited in Zenodo. Its DOI will be linked here once the record is available.

That package contains mine footprints, the species inventory, grid biodiversity values, buffers, exposure tables, randomization inputs/results, derived ownership attribution and Supplementary Data 1–9. These support inspection and reuse of the results. They do not replace all raw inputs needed to reconstruct the complete spatial workflow.

For full reanalysis, obtain study-matched inputs from the original providers. Arrange them under one directory supplied as `--data-root`:

| Relative path | Required content |
| --- | --- |
| `Global_GADM/GC_GADM.shp` | Country boundaries with `NAME_0` and `GID_0`; the study input contains 256 features. |
| `Global_Copper_Mine/Global_Copper_Mine/Global_Copper_Mine-allinfo.shp` | The 1,278 footprints with `Property_I`, `Property_N`, `Country` and geometry. |
| `Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp` | Matching 30-km zones, identifiers and metadata; alternatively generate using `--prepare-buffer`. |
| `IUCNredlist/202502/{taxon}-{status}/data_0.shp` | Fifteen frozen range layers with `ID_NO` and `SCI_NAME`. Taxa: `Amphibians`, `Birds`, `Fish`, `Mammals`, `Reptiles`; statuses: `CR`, `EN`, `VU`. |
| `Global_Copper_Mine/Excel/Global_Copper_Mine-regionrenew.xlsx` | First sheet: `FID`, `Property_I`, `Property_N`, `Country`, `List_of_Ow`. |
| `Global_Copper_Mine/Excel/Copper_mine_all_info.xlsx` | `Sheet1`: `PROP_ID`, `OWNER_LIST`, from the primary ownership snapshot. |
| `Global_Copper_Mine/Excel/Copper_mine_ownerlist-countryrenew.xlsx` | `Sheet1`: `OWNER_NAME`, `OWNER_COUNTRY`. |

Each shapefile needs its `.shp`, `.shx`, `.dbf` and `.prj` files; retain `.cpg` if supplied. Keep identifiers, crosswalks and relevant input ordering consistent. A current provider download may differ from the frozen cohort. The code does not infer historical threat classifications from current online records.

Sources and access:

- Mine footprints: Li et al. (2025), *Machine learning-enhanced monitoring of global copper mining areas*, [dataset version 1](https://doi.org/10.6084/m9.figshare.28680863.v1), with the study selection and stable identifiers.
- Species ranges: [IUCN Red List](https://www.iucnredlist.org/resources/spatial-data-download), February 2025 local snapshot. The cohort contains mining-and-quarrying-threat-listed CR/EN/VU vertebrates.
- Country boundaries: [GADM](https://gadm.org/data.html).
- Ownership: [S&P Global Market Intelligence](https://www.spglobal.com/marketintelligence/), with the study's owner-country crosswalk.

Original IUCN/GADM geometries and commercial ownership records are not bundled; obtain them under provider conditions. The exact historical GADM release, ownership-snapshot date and original acquisition dates have not been established in this repository. Dates are not inferred from file timestamps. A blank primary-snapshot owner list is not replaced with an older snapshot.

## Examples and expected results

### Numerical example

`demo/run_demo.py` calls selected production functions from stages 06, 07, 08 and 12.

| Result | Expected |
| --- | --- |
| Exposed / primary-priority species | 891 / 203 |
| Alternative-priority / shared species | 247 / 202 |
| Synthetic host-country union | 150 km² |
| Synthetic attributed / unattributed area | 120 / 30 km² |

It also checks threshold equality, agreement between two primary-rule implementations, HKG retention, unknown equity, country balancing, seeded sampling and the empirical probability correction.

`exposed_species.csv` contains actual derived results for 891 species: `species_id`, `scientific_name`, `taxon`, `status`, `total_habitat_area_km2`, `total_overlap_area_km2`, `Impact_Percent` and `Is_Copper_Affected`. Areas are km²; exposure percentage is on a 0–100 scale. Full available precision is retained. `expected_priority_ids.csv` contains the 203 expected identifiers. These correspond to Supplementary Data 4 and contain no range geometries. Ownership and randomization examples use fictitious inputs generated in memory.

### Spatial example

`demo/run_spatial_demo.py` generates synthetic polygons and calls production functions from stages 15, 04 and 02; external spatial files are unnecessary.

| Result | Expected |
| --- | --- |
| A 1-km-square footprint with a 1-km buffer | Approximately 8.13654849 km² |
| Sum of two overlapping zones / deduplicated union | 8 / 6 km² |
| Separate overlap areas for two species | 4 and 4 km² |
| RWR values for three cells | 0.5, 5/6 and 1/3 |

One species extends outside the grid, checking use of its full range as the denominator. The example writes two CSVs and a JSON summary; temporary files are removed. These are fictitious geometries, not study observations.

## Running the study workflow

Use separate input and output directories. Replace example paths with your own:

```text
python run_analysis.py --check --data-root ../study_inputs --work-root ../study_results
```

This checks files, shapefile sidecars and intermediate prerequisites, not geometry quality, cohort membership or numerical correctness. Files produced by earlier selected stages are checked again before each dependent stage runs.

If the precomputed primary buffer is unavailable:

```text
python run_analysis.py --prepare-buffer --data-root ../study_inputs --work-root ../study_results
```

Preparation uses the study's per-feature local azimuthal-equidistant buffer function. It writes `prepared_inputs/primary_buffer.gpkg` under the output root, which then takes precedence over the raw buffer shapefile. An existing prepared buffer is not overwritten.

Run all stages or resume selected stages:

```text
python run_analysis.py --all --data-root ../study_inputs --work-root ../study_results
python run_analysis.py --steps 08 09 --data-root ../study_inputs --work-root ../study_results
```

Resuming 08/09 requires `07_priority_species/Species_PrioritySpecies_Full.csv` in the output directory. Missing dependencies are not executed automatically. The optional `--capped-check --steps 08` runs the additional 1,000,000-km² target-area cap check; the default comparison is uncapped.

Use a fresh output folder after changing inputs or parameters: selected stages can replace outputs, and inherited caches do not all check input content. Direct script execution can use `NC_DATA_ROOT` and `NC_WORK_ROOT`, but the runner provides the documented preflight checks.

## Validation and scientific conventions

Checks include a fresh dependency installation, both examples, module imports, syntax, manifest/path consistency and re-execution of 08/09 from the full 5,890-species derived table. Their ten CSV outputs match the source-code outputs exactly. Calculation logic and numeric parameters were preserved during path reorganization.

The complete global grid, IUCN geometry processing, global intersections and all sensitivity workflows have not been rerun end to end for this repository. A passing example is not full-study reproduction.

- RWR sums each species' within-cell area divided by its global range area; it is dimensionless.
- Primary buffers are 30 km; sensitivity distances are 10, 20, 30, 40 and 50 km.
- Full randomization uses 9,999 iterations, seed 20260728, and sampling with replacement within 77 mine-host countries.
- Spatial overlap measures potential exposure, not observed habitat loss.
- Attribution retains HKG and leaves unknown shares unattributed; known shares are not rescaled to 100%. Inherited `NamedDisplay` outputs omit HKG only from selected named displays.
- The 959-mine complete-case analysis is a sensitivity check, not imputation.

## License and citation

Source code and synthetic examples use the [MIT License](LICENSE). This does not relicense included species-derived CSVs or external data, which retain applicable source-data conditions, including [IUCN terms](https://www.iucnredlist.org/terms/terms-of-use).

Please cite the accompanying manuscript, the exact code version or commit used, and the companion Zenodo dataset. Author-level citation metadata and persistent identifiers will be added when finalized; no paper DOI or code-archive DOI is assigned here. Usage questions can be raised through the repository's Issues page.
