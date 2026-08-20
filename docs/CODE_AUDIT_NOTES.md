# Code audit notes for repository preparation

This repository contains the manuscript-facing reproducibility workflow for the Gümüşler mineral prospectivity mapping (MPM) study. The public repository is intentionally focused on the scientific analysis and excludes the separate CrewAI-assisted manuscript-development environment.

## Included

The repository includes the R scripts and supporting files required for the manuscript-facing workflow:

- Pb, Zn, and Cu variogram analysis and ordinary-kriging cross-validation.
- Ordinary-kriging predictor-surface generation.
- Construction of the common 25 m MPM reference grid and occurrence-cell assignments.
- Pseudo-background candidate-pool auditing.
- Repeated leave-one-occurrence-out (LOOCV) fold construction.
- Fold-specific training and pseudo-background preparation.
- Training-only predictor encoding and, where required, scaling.
- Random Forest (RF), XGBoost (XGB), and support vector machine (SVM) fitting.
- Independent pseudo-background validation.
- Full-grid repeated prediction and ensemble summaries.
- RF/XGB/SVM model comparison.
- RF threshold-sensitivity analysis.
- Final RF-P90 prospectivity products.
- Publication-ready RF prospectivity and prediction-variability figures.
- QA/QC records, relevant result summaries, final rasters, and publication figures retained for transparent inspection.

## Validation design represented by the release code

The analysis uses five documented mineral occurrences. Each LOOCV fold holds out one occurrence and uses the remaining four occurrences as positive training labels. Thirty pseudo-background repeats are performed for each of the five folds, producing 150 fold-repeat analyses per algorithm.

Pseudo-background cells are treated as unlabeled/background samples rather than confirmed mineral absences. Fold-specific preparation prevents the held-out occurrence from entering model fitting. The training pseudo-background selection is constructed with the training occurrences separated from the held-out occurrence, and the workflow contains explicit CellID- and distance-based QA/QC checks.

Independent validation pseudo-background samples are generated separately from the training pseudo-background samples. The validation workflow contains 100 independent pseudo-background cells per fold-repeat run, yielding 15,000 pseudo-background validation rows, plus 150 held-out-occurrence predictions for each algorithm.

## Verified execution status

The repository pre-flight check was run under R 4.6.1 and returned `PASS`. The check confirmed the required repository files, analysis scripts, core input files, occurrence-shapefile sidecars, and required R packages.

The RF, XGBoost, and SVM repeated-LOOCV workflows each completed 150 runs with no failed runs in the reported execution. Their independent-validation prediction workflows also completed all 150 runs and passed the implemented QA/QC checks.

For the selected RF model, the final threshold-sensitivity workflow identified P90 as the retained final threshold method. The reported P90 results are:

- held-out occurrence recall/capture: 84.67% (127/150);
- pseudo-background specificity/rejection: 91.09%;
- prospective grid area: 10.0%;
- final RF P90 probability threshold: 0.08310891.

The final RF-P90 raster-production script reported 9,600 grid cells, 960 prospective cells, and successful creation of the continuous, uncertainty/variability, and binary P90 raster products. The publication-figure script also passed its QA/QC checks.

## R packages verified by the repository pre-flight check

The checked environment reported the following packages:

- `sp` 2.2.3
- `gstat` 2.1.6
- `openxlsx` 4.2.8.1
- `raster` 3.6.32
- `terra` 1.9.34
- `sf` 1.1.2
- `ranger` 0.18.0
- `xgboost` 3.2.1.1
- `e1071` 1.7.17
- `pROC` 1.19.0.1
- `ggplot2` 4.0.3
- `scales` 1.4.0

The complete R environment record is stored in `docs/sessionInfo.txt`.

## Excluded from the public scientific workflow

The repository does not require the CrewAI project setup, agent definitions, prompts, handoff records, manuscript-synthesis outputs, or unrelated manuscript-development material. Temporary, cache, session, editor, operating-system, and backup files are excluded through `.gitignore`.

## Remaining checks before archival release

Before creating the archival release, the authors should:

1. perform a fresh-clone execution in a clean R environment;
2. verify every numerical value quoted in the final manuscript directly against the release outputs, including the final ROC-AUC/model-comparison statistics;
3. compare regenerated publication products with the exact submitted manuscript figures;
4. confirm redistribution rights and attribution requirements for source data;
5. perform a final privacy/credential/path audit; and
6. add final manuscript metadata and the selected repository license.

These notes document the verified state of the repository preparation; they do not substitute for the final fresh-clone reproducibility test.
