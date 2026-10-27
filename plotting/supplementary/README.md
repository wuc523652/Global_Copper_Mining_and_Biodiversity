# Supplementary Figures 1–9

`render_supplementary.py` is the entry point for all nine supplementary figures. The nine `supplementary_N.py` files contain their drawing functions. The code reads the accompanying source data and does not rerun the spatial analyses.

## Run

Use Python 3.10 or later. From the repository root:

```text
python -m pip install -r requirements.txt
python plotting/supplementary/render_supplementary.py --data-root /path/to/data --output /path/to/figures/Supplementary
```

The command regenerates all nine PNGs in the required output directory, replacing files with the same names. The output must be outside the code and data packages. To select figures or export additional formats:

```text
python plotting/supplementary/render_supplementary.py --data-root /path/to/data --figures 2 5 9 --output /path/to/preview
python plotting/supplementary/render_supplementary.py --data-root /path/to/data --figures 6 --output /path/to/preview --formats png pdf svg
```

Inputs default to the code package's `data/`; no original project folder or sibling Table/Figure folder is required. An explicit `--data-root` overrides `COPPER_DATA_ROOT`, which overrides the local default. Figure 1 reads the full 1,404,410-cell grid and therefore uses more memory and time than the other figures. Original canvas sizes, colour scales and display selection rules are retained. Figure 1 and Figures 6–8 use 400 dpi, Figures 2–5 use 300 dpi, and Figure 9 uses 600 dpi.

Times New Roman is used for Figures 1–6, DejaVu Sans for Figure 7, and Arial for Figures 8–9. Install these fonts to reproduce the supplied typography. Minor raster differences can occur across library and font versions.

## Source data

| Figure | Input within the materials package |
|---|---|
| 1 | `3Grid_and_Biodiversity/Global_10km_Grid_Biodiversity.gpkg`; `9Figure_Inputs/Supplementary/world_s1.gpkg` |
| 2 | `8Final_Results/Supplementary_Data_3.xlsx`: `Country_Spatial` |
| 3 | Data 3: `Randomization_Countries`, `Randomization_Iterations`, `Randomization_Global` |
| 4 | Data 4: `Exposed_Species` |
| 5 | Data 4: `Rule_Comparison`; the main and alternative threshold curves are calculated in `supplementary_5.py` |
| 6 | Data 6: `Mine_Host_Countries`, `Owner_Countries`; `9Figure_Inputs/Supplementary/world_s6.gpkg` |
| 7 | Data 6: `Completeness_Scenarios`, `Owner_Rank_Comparison`; Data 7: `Pair_Rank_Comparison` |
| 8 | Data 9: `Core_Results`, `Set_Rank_Stability`, `Distance_Diagnostics` |
| 9 | Methods workflow defined as editable text, boxes and arrows in `supplementary_9.py` |

The two GeoPackages are necessary GADM display boundaries in EPSG:4326, simplified with topology preserved at the existing 0.15° (Figure 1) and 0.12° (Figure 6) tolerances. They are drawing inputs, not previously rendered figures, and are not used to calculate areas. The original data terms apply to their redistribution.

Figure 1 retains the display cap of 30 species. Figure 2 includes all 77 mine-host countries; Figure 3 displays the 17 countries with at least ten mines and all 9,999 randomization iterations. Figures 4–5 use all 891 affected species. Figure 6 ranks all registrations and hosts with positive cross-border attributed area. Figure 7 selects the top 20 using the supplied main ranks. Figure 8 uses all five buffer distances. No registration-specific exclusion is applied.

Excel display headers are mapped explicitly by the entry point. ISO identifiers and the supplied ranks and correlations are preserved. Default execution writes only the selected figure files, without audit folders, log files or local Python bytecode caches.
