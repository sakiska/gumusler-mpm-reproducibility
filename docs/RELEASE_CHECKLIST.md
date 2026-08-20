# Public-release checklist

## Scientific reproducibility

- [x] Run the repository pre-flight check successfully with `01_check_repository.R`.
- [x] Save `sessionInfo()` to `docs/sessionInfo.txt`.
- [x] Confirm the repeated LOOCV design contains 5 folds and 30 pseudo-background repeats, giving 150 fold-repeat analyses per algorithm.
- [x] Confirm 150 successful RF LOOCV runs.
- [x] Confirm 150 successful XGBoost LOOCV runs.
- [x] Confirm 150 successful SVM LOOCV runs.
- [x] Confirm independent pseudo-background validation contains 100 validation cells per fold-repeat run (15,000 pseudo-background validation rows in total).
- [x] Confirm the independent validation QA/QC checks pass, including no training-validation background overlap, no training-occurrence overlap, no held-out-occurrence overlap, and the required occurrence-distance exclusions.
- [x] Confirm independent validation predictions are reproduced for RF, XGBoost, and SVM.
- [x] Confirm final RF threshold-sensitivity analysis reproduces the P90 held-out occurrence capture of 84.67%, pseudo-background specificity/rejection of 91.09%, and a prospective grid area of 10.0%.
- [x] Confirm the final RF-P90 prospectivity raster and publication figures are generated successfully.
- [ ] Confirm the final manuscript ROC-AUC value and all other manuscript-reported model-comparison statistics directly against the final release outputs before submission/archiving.
- [ ] Compare the final regenerated raster/figure products with the exact figures included in the submitted manuscript.

## Repository hygiene

- [x] Remove `.RData` and `.Rhistory`.
- [x] Configure `.gitignore` to exclude session, cache, editor, operating-system, temporary, and backup files without broadly excluding scientific inputs or outputs.
- [x] Confirm required repository files, analysis scripts, packages, and core inputs are detected by `01_check_repository.R`.
- [x] Retain QA/QC records, summary tables, final rasters, and publication figures needed for transparent inspection of the analysis.
- [ ] Perform a final check for private absolute paths, credentials, tokens, temporary files, CrewAI prompts/outputs, and unpublished material not intended for public release.

## Data, metadata, and publication

- [ ] Confirm redistribution rights and appropriate attribution for the 1983 MTA-derived geochemical and geological source data.
- [ ] Add the exact manuscript title, authors, journal, and DOI when these details are final/available.
- [ ] Add a repository license selected by the authors; no license should be assumed automatically.
- [ ] Run the repository from a fresh clone in a clean R environment as the final external-style reproducibility test.
- [ ] Create a versioned GitHub release (for example, `v1.0.0-manuscript`).
- [ ] Optionally archive the exact GitHub release with Zenodo and cite the resulting DOI in the manuscript.
