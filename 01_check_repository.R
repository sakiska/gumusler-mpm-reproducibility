# Pre-flight check for the Gümüşler MPM reproducibility repository.
# Run from the repository root:
#   source("01_check_repository.R")

required_packages <- c(
  "sp", "gstat", "openxlsx", "raster", "terra", "sf",
  "ranger", "xgboost", "e1071", "pROC", "ggplot2", "scales"
)

required_root_files <- c(
  "README.md",
  "run_pipeline.R",
  "00_install_packages.R",
  "docs/PIPELINE_GUIDE.md"
)

required_code <- c(
  "Project/Models/Run_Variogram_OK_Pb_Zn_Cu.R",
  "Project/Models/Run_OK_Rasters_Pb_Zn_Cu.R",
  "Project/Models/Create_MPM_FullGrid.R",
  "Project/Models/Create_PseudoBackground.R",
  "Project/Models/Create_LOOCV_Folds.R",
  "Project/Models/Prepare_Fold_Data.R",
  "Project/Models/Prepare_Model_Matrix.R",

  "Project/Models/Train_RF_LOOCV.R",
  "Project/Models/Evaluate_RF.R",
  "Project/Models/Predict_RF_Independent_Validation.R",
  "Project/Models/Predict_RF_FullGrid.R",

  "Project/Models/Train_XGBoost_LOOCV.R",
  "Project/Models/Evaluate_XGBoost.R",
  "Project/Models/Predict_XGB_Independent_Validation.R",
  "Project/Models/Predict_XGBoost_FullGrid.R",

  "Project/Models/Train_SVM_LOOCV.R",
  "Project/Models/Evaluate_SVM.R",
  "Project/Models/Predict_SVM_Independent_Validation.R",
  "Project/Models/Predict_SVM_FullGrid.R",

  "Project/Models/Create_Independent_Validation_Background.R",
  "Project/Models/Compare_RF_XGB_SVM.R",
  "Project/Models/Evaluate_RF_MedianMAD_Threshold.R",
  "Project/Models/Summarize_RF_Threshold_Sensitivity.R",
  "Project/Models/Create_Final_RF_P90_Prospectivity_Map.R",
  "Project/Models/Create_Publication_RF_Figures.R"
)

required_inputs <- c(
  "Project/Data/Suppl_material.csv",
  "Project/Raster/reference_grid_25m.tif",
  "Project/Raster/lithology_25m.tif",
  "Project/Raster/dist_fault_25m.tif",
  "Project/Raster/dist_silicified_25m.tif",
  "Project/Raster/dist_brecciated_25m.tif",
  "Project/Vector/occurrences.shp"
)

required_occurrence_sidecars <- c(
  "Project/Vector/occurrences.dbf",
  "Project/Vector/occurrences.shx",
  "Project/Vector/occurrences.prj"
)

optional_map_overlays <- c(
  "Project/Vector/faults.shp",
  "Project/Vector/silicifiedZones.shp",
  "Project/Vector/brecciatedZones.shp"
)

status <- list()

cat("============================================================\n")
cat("GUMUSLER MPM REPOSITORY PRE-FLIGHT CHECK\n")
cat("============================================================\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "\n\n")

cat("Package availability:\n")
pkg_ok <- logical(length(required_packages))
for (i in seq_along(required_packages)) {
  pkg <- required_packages[i]
  pkg_ok[i] <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (pkg_ok[i]) as.character(utils::packageVersion(pkg)) else "MISSING"
  cat(sprintf("  %-12s %s\n", pkg, ver))
}

cat("\nRepository files:\n")
root_ok <- file.exists(required_root_files)
for (i in seq_along(required_root_files)) {
  cat(sprintf("  %-75s %s\n",
              required_root_files[i],
              if (root_ok[i]) "OK" else "MISSING"))
}

cat("\nAnalysis code files:\n")
code_ok <- file.exists(required_code)
for (i in seq_along(required_code)) {
  cat(sprintf("  %-75s %s\n",
              required_code[i],
              if (code_ok[i]) "OK" else "MISSING"))
}

cat("\nCore input files:\n")
input_ok <- file.exists(required_inputs)
for (i in seq_along(required_inputs)) {
  cat(sprintf("  %-75s %s\n",
              required_inputs[i],
              if (input_ok[i]) "OK" else "MISSING"))
}

cat("\nOccurrence shapefile sidecars:\n")
sidecar_ok <- file.exists(required_occurrence_sidecars)
for (i in seq_along(required_occurrence_sidecars)) {
  cat(sprintf("  %-75s %s\n",
              required_occurrence_sidecars[i],
              if (sidecar_ok[i]) "OK" else "MISSING"))
}

cat("\nOptional final-map overlays:\n")
overlay_ok <- file.exists(optional_map_overlays)
for (i in seq_along(optional_map_overlays)) {
  cat(sprintf("  %-75s %s\n",
              optional_map_overlays[i],
              if (overlay_ok[i]) "AVAILABLE" else "NOT PROVIDED"))
}

all_required_ok <- all(pkg_ok) &&
  all(root_ok) &&
  all(code_ok) &&
  all(input_ok) &&
  all(sidecar_ok)

cat("\n============================================================\n")
if (all_required_ok) {
  cat("PRE-FLIGHT STATUS: PASS\n")
  cat("All required packages, repository files, analysis scripts, and core inputs are available.\n")
} else {
  cat("PRE-FLIGHT STATUS: FAIL\n")
  cat("One or more required items are missing. Resolve the items marked MISSING before running run_pipeline.R.\n")
}
cat("============================================================\n")

cat("\nNotes:\n")
cat("- Optional map overlays do not prevent the scientific pipeline from running.\n")
cat("- ESRI Shapefile inputs require matching .dbf, .shx, and .prj sidecars.\n")
cat("- The publication-figure script additionally uses ggplot2 and scales.\n")
cat("- Before archiving the final release, save sessionInfo() to docs/sessionInfo.txt.\n")

invisible(all_required_ok)
