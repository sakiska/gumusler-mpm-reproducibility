# Gümüşler MPM Reproducibility Pipeline

This guide gives the exact verified execution order for the Gümüşler Pb–Zn–Cu mineral prospectivity mapping workflow.

## 0. Repository check

```r
source("00_install_packages.R")
source("01_check_repository.R")
```

Expected: required packages, code files, and core inputs are available.

## 1. Variogram modeling and OK cross-validation

```r
source("Project/Models/Run_Variogram_OK_Pb_Zn_Cu.R")
```

Purpose: compare raw and log10 representations, fit candidate variogram models, and select final Pb, Zn, and Cu models using ordinary-kriging LOOCV.

Expected:
- Pb: Raw + Exponential, range ≈ 70.21 m
- Zn: Raw + Exponential, range ≈ 120.79 m
- Cu: Raw + Gaussian, range ≈ 109.77 m
- RMSE ≈ 165.37, 34.70, and 125.50 ppm

## 2. Produce OK rasters

```r
source("Project/Models/Run_OK_Rasters_Pb_Zn_Cu.R")
```

Expected:
- 9,600 cells
- 100% coverage
- Pb ≈ 3.24–1491.54 ppm
- Zn ≈ 16.94–317.86 ppm
- Cu ≈ 9.98–485.88 ppm

## 3. Build the full-grid MPM table

```r
source("Project/Models/Create_MPM_FullGrid.R")
```

Expected:
- 9,600 cells
- 7 predictors
- 5 positive cells
- 9,595 unlabeled cells
- no predictor NA values

## 4. Initial/global pseudo-background candidate-pool audit

```r
source("Project/Models/Create_PseudoBackground.R")
```

Expected:
- 64 cells inside/on the initial 50 m exclusion boundary
- 9,536 candidate cells
- minimum candidate distance >50 m
- QA/QC PASS

> The global pool is an audit product. During model fitting, pseudo-background eligibility is reconstructed separately within each fold using only the four training occurrences; the held-out occurrence is excluded before fold-specific background generation.

## 5. Create the LOOCV fold/repeat plan

```r
source("Project/Models/Create_LOOCV_Folds.R")
```

Expected:
- 5 folds
- 4 training occurrences/fold
- 1 held-out occurrence/fold
- 30 repeats/fold
- 100 pseudo-background cells/repeat
- 150 runs/algorithm
- 450 total planned model runs
- QA/QC PASS

### Helper scripts

`Prepare_Fold_Data.R` defines `prepare_fold_data()`.

`Prepare_Model_Matrix.R` defines `prepare_model_matrix()`.

Do not run either as a standalone analysis step.

## 6. Random Forest

```r
source("Project/Models/Train_RF_LOOCV.R")
source("Project/Models/Evaluate_RF.R")
```

Expected: 150 successful runs, 0 failures, mean held-out probability ≈ 0.1145, top grouped predictor = Dist_Silicified.

## 7. XGBoost

```r
source("Project/Models/Train_XGBoost_LOOCV.R")
source("Project/Models/Evaluate_XGBoost.R")
```

Expected: 150 successful runs, 0 failures, mean held-out probability ≈ 0.1830, top grouped predictor = Dist_Silicified.

## 8. Support Vector Machine

```r
source("Project/Models/Train_SVM_LOOCV.R")
source("Project/Models/Evaluate_SVM.R")
```

Expected: 150 successful runs, most frequently selected Cost = 1 and Gamma = 0.01.

## 9. Common independent validation pseudo-background

```r
source("Project/Models/Create_Independent_Validation_Background.R")
```

Expected:
- 150 runs
- 100 validation pseudo-background cells/run
- 15,000 total validation pseudo-background rows
- no training-validation overlap
- all validation cells >50 m from training occurrences
- QA/QC PASS

## 10. Independent-validation prediction

```r
source("Project/Models/Predict_RF_Independent_Validation.R")
source("Project/Models/Predict_XGB_Independent_Validation.R")
source("Project/Models/Predict_SVM_Independent_Validation.R")
```

Expected for each algorithm:
- 150 held-out occurrence rows
- 15,000 pseudo-background rows
- 15,150 total prediction rows
- no training-validation CellID overlap
- QA/QC PASS

## 11. Full-grid prediction

```r
source("Project/Models/Predict_RF_FullGrid.R")
source("Project/Models/Predict_XGBoost_FullGrid.R")
source("Project/Models/Predict_SVM_FullGrid.R")
```

These three steps must complete before model comparison.

## 12. Compare RF, XGB, and SVM

```r
source("Project/Models/Compare_RF_XGB_SVM.R")
```

Expected:

| Metric | RF | XGB | SVM |
|---|---:|---:|---:|
| ROC-AUC | 0.915 | 0.858 | 0.710 |
| Average precision | 0.060 | 0.038 | 0.024 |
| Brier score | 0.0126 | 0.0381 | 0.0122 |
| P90 recall | 0.847 | 0.407 | 0.367 |
| P90 pseudo-background rejection | 0.911 | 0.911 | 0.910 |
| P90 balanced accuracy | 0.879 | 0.659 | 0.638 |

Expected RF ensemble thresholds:
- P90 ≈ 0.08310891
- P95 ≈ 0.15440079
- P99 ≈ 0.29136142
- Mean + 2SD ≈ 0.14599507
- scaled Median + 2MAD ≈ 0.03997547

## 13. RF Median + 2MAD validation

```r
source("Project/Models/Evaluate_RF_MedianMAD_Threshold.R")
```

Expected: recall = 1.000, pseudo-background rejection ≈ 0.7457, QA/QC PASS.

## 14. RF threshold sensitivity

```r
source("Project/Models/Summarize_RF_Threshold_Sensitivity.R")
```

| Threshold | Prospective grid | Recall | Pseudo-background rejection |
|---|---:|---:|---:|
| P90 | 10.00% | 0.8467 | 0.9109 |
| P95 | 5.00% | 0.1267 | 0.9585 |
| P99 | 1.00% | 0.0000 | 0.9923 |
| Mean + 2SD | 5.33% | 0.1800 | 0.9562 |
| scaled Median + 2MAD | 21.10% | 1.0000 | 0.7457 |

## 15. Final RF-P90 scientific products

```r
source("Project/Models/Create_Final_RF_P90_Prospectivity_Map.R")
```

Expected:
- P90 ≈ 0.08310891
- 9,600 cells
- 960 prospective cells
- 10.0% prospective area
- 150 RF predictions/cell
- continuous, SD, and P90 rasters written
- QA/QC PASS

## 16. Publication figures

```r
source("Project/Models/Create_Publication_RF_Figures.R")
```

Expected:
- EPSG:32635 / UTM Zone 35N
- 25 m resolution
- continuous map generated
- binary P90 map generated
- prediction-variability map generated
- publication outputs written
- QA/QC PASS

The file name `Figure_RF_Ensemble_Uncertainty` is retained for manuscript consistency. The mapped quantity is the cell-wise standard deviation across the 150 RF predictions and should be interpreted as between-run prediction variability.

## 17. Save software environment

```r
writeLines(capture.output(sessionInfo()), "docs/sessionInfo.txt")
```
