# Global Copper Mining and Biodiversity

Code for **Global conservation priorities for vertebrates threatened by copper mining**.

This package reproduces the 7 main figures and 9 supplementary figures and recomputes the scientific analyses using the complete research dataset in `data/`. It includes both a shorter reanalysis from archived intermediate data and a 17-stage workflow starting from the frozen original spatial and ownership inputs. No submission-materials `Figure` or `Table` directory is needed.

## Start here

Use 64-bit Python 3.13.3 and the versions in `requirements.txt`. No GPU is needed. Create an environment, activate it, and install the dependencies:

```text
python -m venv .venv
python -m pip install -r requirements.txt
```

Activation: `.venv\Scripts\Activate.ps1` in PowerShell, `.venv\Scripts\activate` in Windows Command Prompt, or `source .venv/bin/activate` on Linux/macOS. Run the installation command after activation. For Conda environments, activate the environment before running commands so that GDAL/PROJ can find their matching resource files.

In the combined local package, `data/` already contains a complete copy of the original companion dataset, including raw inputs, checkpoints, plotting assets and final tables. Default input paths are resolved relative to the installed code, not the terminal's current directory. To use a separately stored dataset, pass `--data-root /path/to/copper-data`. Quote paths containing spaces.

For GitHub distribution, code belongs in the repository and the complete data archives belong in the matching Release assets. Extract all data assets into the code package's `data/` directory, so that `data/deposit_manifest.json` and `data/10Raw_Inputs/` exist. Ordinary Git and GitHub's automatic source-code ZIP do not include this ignored data directory. This local reorganization does not publish a Release or certify third-party redistribution rights. Keep generated results outside the code package and its data.

```text
python reproduce.py --mode check
python reproduce.py --output-root ../copper-results --mode all
```

The first command verifies every deposited file against `deposit_manifest.json`. The second independently recomputes statistics and redraws all 16 figures. `all` uses archived intermediate spatial data; it does not repeat every global polygon operation. Recomputed numerical outputs and comparison reports are written under `calculations/`; the figures under `figures/` read the archived final Excel tables and plotting assets. They do not silently replace the submitted tables with newly computed values.

| Mode | Work performed |
| --- | --- |
| `check` | Verify deposited file names, sizes and SHA-256 hashes; no scientific calculation |
| `calculations` | Recompute statistics from archived intermediate inputs and compare with archived references; see `reproduction/README.md` for exact coverage |
| `figures` | Render 7 main and 9 supplementary figures from the dataset alone |
| `all` | Run `calculations`, then `figures` |
| `raw` | Run all 17 spatial/statistical analysis stages from `10Raw_Inputs`; substantially more time and memory required |

For full spatial reanalysis:

```text
python reproduce.py --output-root ../copper-raw-results --mode raw
```

The `raw` route uses real inputs in `data/10Raw_Inputs`. It executes the spatial/statistical stages; it does not automatically rebuild the formatted supplementary workbooks or redraw figures from newly computed tables. Figure rendering continues to read archived final tables until an explicit, validated table-export step replaces them. The complete 17-stage route has not been validated end to end; see the validation section.

Use a fresh output directory when changing inputs, software versions or parameters. The original scientific scripts can replace their own generated outputs and some retain intermediate caches. The runner never requires writing into the dataset.

## Figures

```text
python plotting/render_all.py --output-root ../copper-figures --figures all
python plotting/render_all.py --output-root ../copper-figures-vector --figures 2 S9 --formats png pdf svg
```

`--figures` accepts `all`, `main`, `supplementary`, or individual identifiers such as `1 S2`. Arial, Microsoft YaHei, Times New Roman and DejaVu Sans are used by the original figures. Fonts are not redistributed. A different font or rendering library can alter layout and pixels without changing the underlying data. See `plotting/README.md` and `plotting/VALIDATION.md`.

## Files and data contracts

```text
data/                     Complete real dataset, including 10Raw_Inputs; release assets
reproduce.py              One entry point, defaulting to this package's data/
requirements.txt          Pinned numerical, spatial, Excel and plotting libraries
analysis_manifest.json    Ordered scientific steps and explicit prerequisites
run_analysis.py           Run/check selected original-data stages
analysis/01_*.py ... 17_*  Spatial and statistical analyses
reproduction/             Independent reanalysis from deposited intermediates
plotting/main/            Main figures 1-7 and reusable plotting modules
plotting/supplementary/   Supplementary figures 1-9
plotting/render_all.py     Unified plotting entry point
validation/               Reports from the local validation run
tools/data_manifest.py    Dataset integrity checker and maintainer inventory builder
tools/package_raw_inputs.py  Maintainer-only copying/fingerprinting of frozen raw files
```

The `data/` directory contains original-layout raw inputs in `10Raw_Inputs`, reusable checkpoints in `11Reanalysis_Checkpoints`, all original plotting assets in `9Figure_Inputs`, and the nine unmodified final supplementary workbooks in `8Final_Results`. All data files, the data README and its manifest are preserved byte for byte from the original local Zenodo package. Its README retains the historical deposit description; the current local layout and intended GitHub download arrangement are documented here. The legacy grid split folders are retained for compatibility; the workflow uses `3Grid_and_Biodiversity`.

Identifiers must be joined by `mine_id`/`property_id`, `species_id` and `grid_id`, not by row number. Preserve the background `grid_order` and country ordering for exact seeded randomization. Some reference files are explicitly retained only for validation; recomputation does not treat those expected answers as its calculation inputs. See the checkpoint manifest and reanalysis documentation for provenance.

## Full analysis order

| Step | Calculation |
| --- | --- |
| 01 | Global 10-km equal-area land grid |
| 02 | Area-weighted rarity-weighted richness (RWR) |
| 03 | Threatened-species richness |
| 04 | Mine-host and globally deduplicated species overlap |
| 05 | Country/global RWR classes within exposure unions |
| 06 | Country-stratified mine-location randomization, 9,999 iterations |
| 07 | Species exposure and the primary priority rule |
| 08 | Alternative priority rules and boundary cases |
| 09 | Screening margins and boundary diagnostics |
| 10 | Stable mine identifiers, reported ownership ledger and coverage |
| 11 | Mine, host-country and global exposure unions |
| 12 | Ownership attribution with country-overlap balancing, including HKG |
| 13 | Mine-priority-species associations and ecological categories |
| 14 | Equity-weighted ecological portfolios |
| 15 | 10/20/30/40/50-km buffer-distance sensitivity |
| 16 | Complete-case ownership sensitivity, 959 mines |
| 17 | Priority-species grid associations, counts and dominant taxon |

List prerequisites or resume selected stages:

```text
python run_analysis.py --list
python run_analysis.py --check --work-root ../copper-raw-results/raw_analysis
python run_analysis.py --steps 08 09 --work-root ../copper-raw-results/raw_analysis
```

The runner defaults to `data/` and also accepts an explicit complete dataset or `10Raw_Inputs` folder. Dependencies must already exist when resuming; they are not run implicitly. `--capped-check --steps 08` additionally runs the historical 1,000,000-km2 cap sensitivity. `--prepare-buffer` can rebuild the primary 30-km buffer from the footprints; the frozen primary buffer is supplied for matching the submitted analysis.

## Scientific conventions and expected study totals

The frozen cohort contains 1,278 mines and 5,890 threatened vertebrates. The primary 30-km analysis exposes 891 species and selects 203 priority species. The alternative rule selects 247 species, with 202 shared. The base grid has 1,404,410 cells; priority-species mapping retains 1,404,063 valid WGS84 cells, 4,771 nonzero cells and 7,783 unique grid-species pairs. The 347 excluded cells must not be filled with zero.

RWR is the sum of within-cell species range area divided by each species' full global range area. Exposure is potential spatial exposure, not observed habitat loss. Randomization uses seed 20260728, 9,999 iterations and within-country sampling with replacement across 77 mine-host countries. Unknown ownership remains unattributed; reported shares are not rescaled to 100%. Complete-case results are a sensitivity analysis, not imputation. Primary and sensitivity 30-km geometry results are kept distinct.

## Validation and limits

Local validation for the 2026-10-18 update uses Windows and Python 3.13.3 with the pinned scientific-library versions. It uses an isolated environment inheriting installed geospatial libraries; this update does not claim a newly downloaded clean installation or verification on Linux/macOS. See `validation/README.md` for the actual completed runs and any remaining limitations. A hash check or successful small demo is not a full-study scientific validation.

All 16 PNG figures were executed. All dimensions match existing figures; 12 PNGs match pixel for pixel. For figures 6, 7, S4 and S5, re-running the unmodified original plotting sources in the same environment produces exactly the new package output; the differences from archived PNGs therefore predate this migration. Their historical cause is not established. PDF/SVG export was exercised for figures 2 and S9. No scientific values, thresholds or selection rules were changed to improve a comparison.

## Access, licensing and citation

Local layout update: 2026-10-27. The complete dataset was copied into `data/`, entry points were configured to use it by default, and the optional examples were subsequently removed from this distribution. The manuscript title and funding retain the 2026-10-20 metadata. The full statistical and figure validation above remains the 2026-10-18 run; see `validation/local_data_layout_20261027.md` for the narrower checks of this path update. No remote repository, Release or Zenodo record was updated by this local operation.

Code repository: [https://github.com/wuc523652/Global_Copper_Mining_and_Biodiversity](https://github.com/wuc523652/Global_Copper_Mining_and_Biodiversity). The intended data download location is the corresponding [GitHub Release](https://github.com/wuc523652/Global_Copper_Mining_and_Biodiversity/releases); upload/publication of these updated local files remains pending. The manuscript also records the **reserved, unpublished Zenodo DOI [10.5281/zenodo.23277107](https://doi.org/10.5281/zenodo.23277107)**. That identifier does not establish public data availability and is not required to run this local package. No new DOI, release tag or publication status is asserted here.

The manuscript author list and final author-based citations will be added after author confirmation. Do not infer authors or a publication year from this local package date.

The code uses the existing MIT license. That license does not cover third-party data, maps, flags, icons or fonts. The local raw-input copy includes IUCN range geometries, GADM boundaries and S&P ownership records under their original conditions; inclusion for local reproduction is not certification that they may be publicly redistributed. The same consideration applies to derived data. See the dataset README and `10Raw_Inputs/README.md` before choosing public, restricted-access or provider-obtained distribution. Removing restricted files changes the self-contained reproduction claim and requires an updated manifest and access instructions.

Please cite the manuscript, the exact code release and dataset version, together with the original data sources. The February 2025 IUCN snapshot is known; historical GADM release and ownership-snapshot acquisition dates have not been inferred from file timestamps. File hashes identify the exact bytes used.

## Study funding

- National Natural Science Foundation of China (72274187, 72534006).
- Deep Earth Probe and Mineral Resources Exploration-National Science and Technology Major Project (2025ZD1007007).
- Key R&D Program of China (2024YFC3712900).
