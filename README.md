# Gümüşler Pb–Zn–Cu Mineral Prospectivity Mapping

Reproducibility code for the validation-oriented machine-learning mineral prospectivity mapping workflow developed for the Gümüşler Pb–Zn–Cu system, Biga Peninsula, NW Türkiye.

The study addresses mineral prospectivity mapping under severe occurrence-data scarcity (five documented mineral occurrences) by integrating Pb, Zn and Cu ordinary-kriging surfaces with lithology and distances to faults, silicified zones and brecciated zones. Random Forest (RF), Extreme Gradient Boosting (XGB) and radial Support Vector Machine (SVM) are evaluated using repeated leave-one-occurrence-out validation and independently sampled pseudo-background realizations. The final workflow also evaluates operational threshold sensitivity and between-run RF prediction variability.

> **Repository status:** pre-submission reproducibility draft. The analysis code is included; source-data redistribution permissions and a clean-clone execution test must be completed before public release.

## Key methodological features

- Five documented Pb–Zn–Cu occurrences.
- Seven predictors: `Pb_OK`, `Zn_OK`, `Cu_OK`, lithology, distance to faults, distance to silicified zones, and distance to brecciated zones.
- Five leave-one-occurrence-out folds × 30 repeats = 150 outer runs per algorithm.
- Fold-specific pseudo-background selection.
- Held-out occurrence excluded from model fitting and training-side preprocessing/tuning.
- Common independent pseudo-background validation realizations for RF, XGB and SVM.
- Threshold-free and threshold-dependent model comparison.
- Study-specific threshold-sensitivity analysis.
- Final RF ensemble mean and between-run standard deviation on a 25 × 25 m, 9,600-cell grid.
- Reproducible publication figures generated directly from the final RF rasters.

## R requirements

Use a current R installation compatible with the listed packages. RStudio is optional but convenient.

| Package | Purpose in this workflow |
|---|---|
| `sp` | legacy spatial classes used by geostatistical scripts |
| `gstat` | variogram modeling and ordinary kriging |
| `openxlsx` | Excel summary-table output |
| `raster` | raster production in the OK workflow |
| `terra` | modern raster/vector grid handling |
| `sf` | vector spatial data and CRS handling |
| `ranger` | Random Forest |
| `xgboost` | Extreme Gradient Boosting |
| `e1071` | radial SVM |
| `pROC` | ROC-AUC evaluation |
| `ggplot2` | publication-quality map generation |
| `scales` | publication-figure axis and label helpers |

Install/check them with:

```r
source("00_install_packages.R")
```

Before running the full workflow:

```r
source("01_check_repository.R")
```

For the archival release:

```r
writeLines(capture.output(sessionInfo()), "docs/sessionInfo.txt")
```

## Directory structure

```text
.
├── 00_install_packages.R
├── 01_check_repository.R
├── run_pipeline.R
├── docs/
├── Project/
│   ├── Data/
│   ├── Raster/
│   ├── Vector/
│   ├── Models/
│   ├── Tables/
│   ├── Figures/
│   └── QAQC/
└── .gitignore
```

The `Project/` layout is intentionally retained because the final analysis scripts use repository-relative paths such as `Project/Data/...`.

## Required source/prepared inputs

```text
Project/Data/Suppl_material.csv
Project/Raster/reference_grid_25m.tif
Project/Raster/lithology_25m.tif
Project/Raster/dist_fault_25m.tif
Project/Raster/dist_silicified_25m.tif
Project/Raster/dist_brecciated_25m.tif
Project/Vector/occurrences.shp (+ shapefile sidecars)
```

Additional geological vectors can be supplied for the final map overlays.

## Reproducing the analysis

From the repository root:

```r
source("run_pipeline.R")
```

The master pipeline follows the verified order:

1. geochemical variogram modeling and OK cross-validation;
2. Pb–Zn–Cu OK raster production;
3. full-grid MPM table construction;
4. initial/global pseudo-background candidate-pool audit;
5. LOOCV fold/repeat plan creation;
6. RF training and evaluation;
7. XGBoost training and evaluation;
8. SVM training and evaluation;
9. common independent validation pseudo-background generation;
10. RF/XGB/SVM independent-validation prediction;
11. RF/XGB/SVM full-grid prediction;
12. algorithm comparison;
13. RF Median + 2MAD validation;
14. RF threshold-sensitivity summary;
15. final RF-P90 scientific products;
16. final publication figures.

`Prepare_Fold_Data.R` and `Prepare_Model_Matrix.R` are helper modules and are not standalone analysis stages.

For exact script order, purpose, inputs, outputs, and expected checkpoints, see:

**[`docs/PIPELINE_GUIDE.md`](docs/PIPELINE_GUIDE.md)**

## Expected manuscript-level checks

- RF ROC-AUC ≈ 0.915.
- XGB ROC-AUC ≈ 0.858.
- SVM ROC-AUC ≈ 0.710.
- RF P90 held-out occurrence capture ≈ 84.67%.
- RF P90 pseudo-background rejection ≈ 91.09%.
- P90 prospective domain = 10.0% of the 9,600-cell grid.
- Final RF P90 threshold ≈ 0.08310891.
- Distance to silicified zones is the highest mean grouped RF permutation-importance predictor.

These values are verification targets, not hard-coded acceptance criteria.

## Final RF and publication outputs

The final scientific-product script generates:

- `RF_Final_Continuous_Probability.tif`
- `RF_Final_Probability_SD.tif`
- `RF_Final_P90_Prospectivity.tif`

The publication-figure script then generates:

- `Figure_RF_Continuous_Probability`
- `Figure_RF_P90_Prospectivity`
- `Figure_RF_Ensemble_Uncertainty`

The file name `Figure_RF_Ensemble_Uncertainty` is retained for manuscript consistency; the mapped quantity is the cell-wise standard deviation of RF predictions across the 150 accepted runs and is interpreted as between-run prediction variability.

## Data and pseudo-background interpretation

Pseudo-background cells are locations without documented mineralization, not confirmed barren sites. Classification metrics therefore quantify discrimination between documented occurrences and sampled pseudo-background under the study design.

## Code provenance and cleanup

This repository contains the manuscript-facing analysis pipeline only. CrewAI orchestration files, agent prompts, knowledge-base files, exploratory/obsolete script versions, temporary logs and manuscript-writing artifacts are intentionally excluded.

The original `Summarize_RF_Threshold_Sensitivity_v2.R` was retained as the canonical `Summarize_RF_Threshold_Sensitivity.R`. Only the canonical publication-figure script, `Create_Publication_RF_Figures.R`, is retained.

## Citation

Citation metadata will be added after the manuscript bibliographic details are final.

## License

No software license has been selected in this pre-submission draft.
