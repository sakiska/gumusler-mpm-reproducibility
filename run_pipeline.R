# Gümüşler MPM reproducibility pipeline
# Run from the repository root. The order mirrors the verified manuscript workflow.

run_script <- function(path) {
  if (!file.exists(path)) stop("Required pipeline script not found: ", path)
  cat("\n============================================================\n")
  cat("Running:", path, "\n")
  cat("============================================================\n")
  source(path, echo = FALSE, chdir = FALSE)
}

# A. Geochemical spatial modeling
run_script("Project/Models/Run_Variogram_OK_Pb_Zn_Cu.R")
run_script("Project/Models/Run_OK_Rasters_Pb_Zn_Cu.R")

# B. Common MPM grid and occurrence assignment
run_script("Project/Models/Create_MPM_FullGrid.R")

# Global candidate-pool audit only. Actual ML pseudo-background eligibility
# is reconstructed separately within each LOOCV fold from training occurrences.
run_script("Project/Models/Create_PseudoBackground.R")

# C. Repeated leave-one-occurrence-out design
run_script("Project/Models/Create_LOOCV_Folds.R")

# Prepare_Fold_Data.R and Prepare_Model_Matrix.R are helper modules sourced
# downstream; they are not standalone pipeline stages.

# D. Random Forest
run_script("Project/Models/Train_RF_LOOCV.R")
run_script("Project/Models/Evaluate_RF.R")

# E. XGBoost
run_script("Project/Models/Train_XGBoost_LOOCV.R")
run_script("Project/Models/Evaluate_XGBoost.R")

# F. Support Vector Machine
run_script("Project/Models/Train_SVM_LOOCV.R")
run_script("Project/Models/Evaluate_SVM.R")

# G. Common independent pseudo-background validation set
run_script("Project/Models/Create_Independent_Validation_Background.R")
run_script("Project/Models/Predict_RF_Independent_Validation.R")
run_script("Project/Models/Predict_XGB_Independent_Validation.R")
run_script("Project/Models/Predict_SVM_Independent_Validation.R")

# H. Full-grid predictions for all accepted runs
run_script("Project/Models/Predict_RF_FullGrid.R")
run_script("Project/Models/Predict_XGBoost_FullGrid.R")
run_script("Project/Models/Predict_SVM_FullGrid.R")

# I. Algorithm comparison and RF threshold sensitivity
run_script("Project/Models/Compare_RF_XGB_SVM.R")
run_script("Project/Models/Evaluate_RF_MedianMAD_Threshold.R")
run_script("Project/Models/Summarize_RF_Threshold_Sensitivity.R")

# J. Final RF P90 scientific products
run_script("Project/Models/Create_Final_RF_P90_Prospectivity_Map.R")

# K. Manuscript/publication figures
run_script("Project/Models/Create_Publication_RF_Figures.R")

cat("\n============================================================\n")
cat("GUMUSLER MPM REPRODUCIBILITY PIPELINE COMPLETED\n")
cat("============================================================\n")
