# Manuscript-to-code map

| Manuscript component | Main script(s) |
|---|---|
| Geochemical variograms and ordinary-kriging cross-validation | `Project/Models/Run_Variogram_OK_Pb_Zn_Cu.R` |
| Pb/Zn/Cu ordinary-kriging surfaces | `Project/Models/Run_OK_Rasters_Pb_Zn_Cu.R` |
| 25 m predictor grid and occurrence-cell assignment | `Project/Models/Create_MPM_FullGrid.R` |
| Global pseudo-background candidate-pool audit | `Project/Models/Create_PseudoBackground.R` |
| Repeated leave-one-occurrence-out fold/repeat plan | `Project/Models/Create_LOOCV_Folds.R` |
| Fold-specific training and pseudo-background construction | `Project/Models/Prepare_Fold_Data.R` |
| Training-only encoding, imputation, and optional scaling | `Project/Models/Prepare_Model_Matrix.R` |
| RF repeated LOOCV | `Project/Models/Train_RF_LOOCV.R` |
| RF evaluation | `Project/Models/Evaluate_RF.R` |
| RF independent-validation predictions | `Project/Models/Predict_RF_Independent_Validation.R` |
| RF full-grid repeated predictions | `Project/Models/Predict_RF_FullGrid.R` |
| XGBoost repeated LOOCV | `Project/Models/Train_XGBoost_LOOCV.R` |
| XGBoost evaluation | `Project/Models/Evaluate_XGBoost.R` |
| XGBoost independent-validation predictions | `Project/Models/Predict_XGB_Independent_Validation.R` |
| XGBoost full-grid repeated predictions | `Project/Models/Predict_XGBoost_FullGrid.R` |
| SVM repeated LOOCV | `Project/Models/Train_SVM_LOOCV.R` |
| SVM evaluation | `Project/Models/Evaluate_SVM.R` |
| SVM independent-validation predictions | `Project/Models/Predict_SVM_Independent_Validation.R` |
| SVM full-grid repeated predictions | `Project/Models/Predict_SVM_FullGrid.R` |
| Independent pseudo-background validation set | `Project/Models/Create_Independent_Validation_Background.R` |
| RF/XGB/SVM comparison and standard threshold assessment | `Project/Models/Compare_RF_XGB_SVM.R` |
| Scaled Median + 2MAD RF evaluation | `Project/Models/Evaluate_RF_MedianMAD_Threshold.R` |
| Final RF threshold-sensitivity summary | `Project/Models/Summarize_RF_Threshold_Sensitivity.R` |
| Final RF/P90 raster products | `Project/Models/Create_Final_RF_P90_Prospectivity_Map.R` |
| Publication-ready RF probability, P90, and prediction-variability figures | `Project/Models/Create_Publication_RF_Figures.R` |

## Leakage-control implementation

The central leakage-control logic is implemented in `Prepare_Fold_Data.R` and `Prepare_Model_Matrix.R`, with additional run-level checks in the algorithm-specific training and validation scripts.

The study contains five documented mineral occurrences. For each LOOCV fold, one occurrence is held out for testing and the remaining four are used as positive training labels. Thirty pseudo-background repeats are evaluated for every fold, producing 5 × 30 = 150 fold-repeat analyses per algorithm.

The held-out occurrence is excluded from model fitting. Fold-specific training pseudo-background selection is constructed using the training occurrences, and the workflow applies explicit CellID and distance checks to prevent inappropriate overlap. Predictor preprocessing is learned from the training data and then applied to the held-out occurrence. For SVM, numeric scaling is therefore based on training data rather than the test occurrence.

Algorithm-specific model fitting and tuning are performed within the corresponding RF, XGBoost, and SVM training scripts.

## Independent validation

`Create_Independent_Validation_Background.R` constructs a validation pseudo-background set that is separate from the pseudo-background cells used for model training. Each of the 150 fold-repeat runs receives 100 independent validation pseudo-background cells, producing 15,000 pseudo-background validation rows.

The algorithm-specific `Predict_*_Independent_Validation.R` scripts reconstruct predictions for both the held-out occurrences and independent pseudo-background cells and perform QA/QC checks for overlap, probability validity, reproducibility, and run status.

This separation is important because pseudo-background samples are not treated as confirmed mineral absences. Performance statistics therefore quantify discrimination between held-out known occurrences and independently sampled pseudo-background/unlabeled locations under the implemented sampling design.

## Full-grid prediction and final model products

The `Predict_*_FullGrid.R` scripts generate repeated predictions across the complete 9,600-cell reference grid. `Compare_RF_XGB_SVM.R` compares the three algorithms and evaluates the standard threshold alternatives.

For the selected RF model, `Evaluate_RF_MedianMAD_Threshold.R` and `Summarize_RF_Threshold_Sensitivity.R` provide the final threshold-sensitivity assessment. The retained P90 threshold captures 127 of 150 held-out occurrence predictions (84.67% recall), rejects 91.09% of independent pseudo-background predictions, and classifies 10.0% of the reference grid as prospective.

`Create_Final_RF_P90_Prospectivity_Map.R` creates the final continuous RF ensemble, prediction-variability/uncertainty, and binary P90 raster products. `Create_Publication_RF_Figures.R` converts the final products into the publication-ready map figures using the EPSG:32635 projected coordinate system.

## Scope

This document maps the manuscript-facing scientific workflow to the release code. It should be read together with `PIPELINE_GUIDE.md`, `CODE_AUDIT_NOTES.md`, and the QA/QC outputs. Final manuscript numbers should be checked directly against the archived release outputs before publication.
