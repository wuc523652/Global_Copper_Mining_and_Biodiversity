# Main-figure plotting code

`figure_01.py`–`figure_07.py` reproduce all seven manuscript figures, including every panel. Other Python files contain shared plotting, layout and input functions.

Install the repository's root `requirements.txt`. From the repository root, run:

```text
python plotting/main/figure_01.py --data-root /path/to/data --output-dir /path/to/figures/Main
python plotting/main/figure_06.py --data-root /path/to/data --output-dir /path/to/figures/Main
```

Each script writes its corresponding `Figure_N.png` to the required output directory, replacing that file if it already exists. The output must be outside the code and data packages. Optional editable exports: `--format pdf` or `--format svg`. `COPPER_DATA_ROOT` can replace the `--data-root` argument.

Inputs default to the code package's `data/`; the optional `--data-root` argument takes precedence over `COPPER_DATA_ROOT` and that default. Inputs include published workbooks in `8Final_Results`, analysis outputs in their numbered directories, and map, region and icon inputs in `9Figure_Inputs/Main`. No sibling Figure or Table directory is needed. No previous figure image is used to construct a panel. See `../README.md` for the unified all-figure entry point.

| Figure | Main numeric inputs |
|---|---|
| 1 | Supplementary Data 3–4; complete richness-grid and mine-point plotting coordinates |
| 2 | Zenodo `6Randomization` results |
| 3 | Supplementary Data 2 and 4 species inventories |
| 4 | Zenodo priority-species grid; regional species habitats and lists |
| 5 | Supplementary Data 5, `Country_Summary` |
| 6 | Supplementary Data 6–7; the manuscript map selection |
| 7 | Supplementary Data 8 and the Zenodo mine–owner ledger |

Figure 3 uses continuous stacks of exact source counts. Figure 6 retains its manuscript map and five-registration flow selection, while panel c uses the complete bilateral-flow table. Figure 7 retains the twelve portfolios in panel a; panel b uses the complete-data top-eight registration and top-ten host selections.

Use Arial for the main figure text and Microsoft YaHei for the alluvial labels to reproduce the supplied typography. Available fallback fonts are used if either font is absent. Different font versions or rendering-library versions can slightly alter text spacing.
