# ============================================================
# Compare_RF_XGB_SVM.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Compare RF, XGB, and SVM using identical independent validation observations and threshold-free and P90 metrics.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("pROC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------
project_dir <- "Project"
models_dir <- file.path(project_dir, "Models")
qaqc_dir <- file.path(project_dir, "QAQC")
tables_dir <- file.path(project_dir, "Tables")
figures_dir <- file.path(project_dir, "Figures")
comparison_output_dir <- file.path(models_dir, "Model_Comparison")
dir.create(comparison_output_dir, recursive = TRUE, showWarnings = FALSE)
# ============================================================
# ============================================================
rf_fullgrid_summary_file <- file.path(
  models_dir, "RF_Results", "FullGrid", "RF_FullGrid_Ensemble_Summary.csv"
)
rf_threshold_summary_file <- file.path(
  models_dir, "RF_Results", "FullGrid", "RF_Threshold_Capture_Summary.csv"
)
rf_heldout_capture_file <- file.path(
  models_dir, "RF_Results", "FullGrid", "RF_HeldOut_Threshold_Capture.csv"
)
xgb_fullgrid_summary_file <- file.path(
  models_dir, "XGB_Results", "FullGrid", "XGB_FullGrid_Ensemble_Summary.csv"
)
xgb_threshold_summary_file <- file.path(
  models_dir, "XGB_Results", "FullGrid", "XGB_Threshold_Capture_Summary.csv"
)
xgb_heldout_capture_file <- file.path(
  models_dir, "XGB_Results", "FullGrid", "XGB_HeldOut_Threshold_Capture.csv"
)
svm_fullgrid_summary_file <- file.path(
  models_dir, "SVM_Results", "FullGrid", "SVM_FullGrid_Ensemble_Summary.csv"
)
svm_threshold_summary_file <- file.path(
  models_dir, "SVM_Results", "FullGrid", "SVM_Threshold_Capture_Summary.csv"
)
svm_heldout_capture_file <- file.path(
  models_dir, "SVM_Results", "FullGrid", "SVM_HeldOut_Threshold_Capture.csv"
)
input_files <- data.frame(
  Algorithm = c("RF", "RF", "RF", "XGB", "XGB", "XGB", "SVM", "SVM", "SVM"),
  File_Type = rep(c("FullGrid Ensemble", "Threshold Summary", "HeldOut Capture"), 3),
  File = c(
    rf_fullgrid_summary_file,
    rf_threshold_summary_file,
    rf_heldout_capture_file,
    xgb_fullgrid_summary_file,
    xgb_threshold_summary_file,
    xgb_heldout_capture_file,
    svm_fullgrid_summary_file,
    svm_threshold_summary_file,
    svm_heldout_capture_file
  ),
  stringsAsFactors = FALSE
)
input_files$Exists <- file.exists(input_files$File)
print(input_files)
if (!all(input_files$Exists)) {
  print(input_files[!input_files$Exists, , drop = FALSE])
  stop(
    "One or more required model-comparison input files are missing.",
    call. = FALSE
  )
}
message("")
message("============================================")
message("MODEL COMPARISON INPUT AUDIT")
message("============================================")
message("All required RF, XGB and SVM files found.")
message("Stage 1 PASS")
# ============================================================
# ============================================================
rf_fullgrid_summary <- read.csv(rf_fullgrid_summary_file, stringsAsFactors = FALSE)
rf_threshold_summary <- read.csv(rf_threshold_summary_file, stringsAsFactors = FALSE)
rf_heldout_capture <- read.csv(rf_heldout_capture_file, stringsAsFactors = FALSE)
xgb_fullgrid_summary <- read.csv(xgb_fullgrid_summary_file, stringsAsFactors = FALSE)
xgb_threshold_summary <- read.csv(xgb_threshold_summary_file, stringsAsFactors = FALSE)
xgb_heldout_capture <- read.csv(xgb_heldout_capture_file, stringsAsFactors = FALSE)
svm_fullgrid_summary <- read.csv(svm_fullgrid_summary_file, stringsAsFactors = FALSE)
svm_threshold_summary <- read.csv(svm_threshold_summary_file, stringsAsFactors = FALSE)
svm_heldout_capture <- read.csv(svm_heldout_capture_file, stringsAsFactors = FALSE)
standardize_fullgrid <- function(df, algorithm) {
  old_names <- c(
    paste0(algorithm, "_Mean_Probability"),
    paste0(algorithm, "_SD_Probability"),
    paste0(algorithm, "_Min_Probability"),
    paste0(algorithm, "_Max_Probability"),
    paste0(algorithm, "_Prediction_n")
  )
  new_names <- c(
    "Mean_Probability",
    "SD_Probability",
    "Min_Probability",
    "Max_Probability",
    "Prediction_n"
  )
  for (i in seq_along(old_names)) {
    names(df)[names(df) == old_names[i]] <- new_names[i]
  }
  df$Algorithm <- algorithm
  df
}
rf_fullgrid_std <- standardize_fullgrid(rf_fullgrid_summary, "RF")
xgb_fullgrid_std <- standardize_fullgrid(xgb_fullgrid_summary, "XGB")
svm_fullgrid_std <- standardize_fullgrid(svm_fullgrid_summary, "SVM")
structure_qaqc <- data.frame(
  Check = c(
    "Standardized FullGrid structures identical",
    "Same CellID",
    "Same X coordinates",
    "Same Y coordinates",
    "Threshold Summary structures identical",
    "HeldOut Capture structures identical"
  ),
  Result = c(
    identical(names(rf_fullgrid_std), names(xgb_fullgrid_std)) &&
      identical(names(rf_fullgrid_std), names(svm_fullgrid_std)),
    identical(rf_fullgrid_std$CellID, xgb_fullgrid_std$CellID) &&
      identical(rf_fullgrid_std$CellID, svm_fullgrid_std$CellID),
    isTRUE(all.equal(rf_fullgrid_std$X, xgb_fullgrid_std$X)) &&
      isTRUE(all.equal(rf_fullgrid_std$X, svm_fullgrid_std$X)),
    isTRUE(all.equal(rf_fullgrid_std$Y, xgb_fullgrid_std$Y)) &&
      isTRUE(all.equal(rf_fullgrid_std$Y, svm_fullgrid_std$Y)),
    identical(names(rf_threshold_summary), names(xgb_threshold_summary)) &&
      identical(names(rf_threshold_summary), names(svm_threshold_summary)),
    identical(names(rf_heldout_capture), names(xgb_heldout_capture)) &&
      identical(names(rf_heldout_capture), names(svm_heldout_capture))
  ),
  stringsAsFactors = FALSE
)
print(structure_qaqc)
if (!all(structure_qaqc$Result)) {
  stop("Stage 2 structure QA/QC failed.", call. = FALSE)
}
message("")
message("============================================")
message("MODEL COMPARISON STRUCTURE AUDIT")
message("============================================")
message("FullGrid structures successfully standardized.")
message("RF, XGB and SVM use the same 9600-cell grid.")
message("Threshold Summary structures are identical.")
message("HeldOut Capture structures are identical.")
message("Stage 2 PASS")
# ============================================================
# ============================================================
summarize_model_grid <- function(df, algorithm_name) {
  data.frame(
    Algorithm = algorithm_name,
    Mean_of_MeanProbability = mean(df$Mean_Probability, na.rm = TRUE),
    SD_of_MeanProbability = sd(df$Mean_Probability, na.rm = TRUE),
    Median_MeanProbability = median(df$Mean_Probability, na.rm = TRUE),
    Min_MeanProbability = min(df$Mean_Probability, na.rm = TRUE),
    Max_MeanProbability = max(df$Mean_Probability, na.rm = TRUE),
    Mean_Uncertainty = mean(df$SD_Probability, na.rm = TRUE),
    Median_Uncertainty = median(df$SD_Probability, na.rm = TRUE),
    Max_Uncertainty = max(df$SD_Probability, na.rm = TRUE),
    Mean_Prediction_Range = mean(
      df$Max_Probability - df$Min_Probability,
      na.rm = TRUE
    ),
    Median_Prediction_Range = median(
      df$Max_Probability - df$Min_Probability,
      na.rm = TRUE
    ),
    Min_Prediction_n = min(df$Prediction_n, na.rm = TRUE),
    Max_Prediction_n = max(df$Prediction_n, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
model_grid_summary <- rbind(
  summarize_model_grid(rf_fullgrid_std, "RF"),
  summarize_model_grid(xgb_fullgrid_std, "XGB"),
  summarize_model_grid(svm_fullgrid_std, "SVM")
)
print(model_grid_summary)
prediction_count_audit <- model_grid_summary[
  ,
  c("Algorithm", "Min_Prediction_n", "Max_Prediction_n"),
  drop = FALSE
]
write.csv(
  model_grid_summary,
  file.path(comparison_output_dir, "RF_XGB_SVM_FullGrid_Summary.csv"),
  row.names = FALSE
)
write.csv(
  prediction_count_audit,
  file.path(comparison_output_dir, "RF_XGB_SVM_Prediction_Count_Audit.csv"),
  row.names = FALSE
)
message("")
message("============================================")
message("FULLGRID DISTRIBUTION AND UNCERTAINTY AUDIT")
message("============================================")
message("Model-level grid summaries created.")
message("Stage 3 PASS")
# ============================================================
# ============================================================
rf_threshold_summary$Algorithm <- "RF"
xgb_threshold_summary$Algorithm <- "XGB"
svm_threshold_summary$Algorithm <- "SVM"
threshold_comparison <- rbind(
  rf_threshold_summary,
  xgb_threshold_summary,
  svm_threshold_summary
)
threshold_comparison$Threshold_CV <- ifelse(
  threshold_comparison$Mean_Threshold > 0,
  threshold_comparison$SD_Threshold / threshold_comparison$Mean_Threshold,
  NA_real_
)
threshold_comparison$Capture_Percent <- 100 * threshold_comparison$Capture_Rate
threshold_comparison$Algorithm <- factor(
  threshold_comparison$Algorithm,
  levels = c("RF", "XGB", "SVM")
)
threshold_comparison$Threshold_Method <- factor(
  threshold_comparison$Threshold_Method,
  levels = c("P90", "P95", "P99", "Mean_Plus_2SD")
)
threshold_comparison <- threshold_comparison[
  order(
    threshold_comparison$Threshold_Method,
    threshold_comparison$Algorithm
  ),
  ,
  drop = FALSE
]
threshold_comparison$Algorithm <- as.character(threshold_comparison$Algorithm)
threshold_comparison$Threshold_Method <- as.character(threshold_comparison$Threshold_Method)
p90_comparison <- threshold_comparison[
  threshold_comparison$Threshold_Method == "P90",
  ,
  drop = FALSE
]
p90_comparison <- p90_comparison[
  order(-p90_comparison$Capture_Rate),
  ,
  drop = FALSE
]
cat("\n============================================\n")
cat("THRESHOLD CAPTURE COMPARISON\n")
cat("============================================\n")
print(threshold_comparison)
cat("\n============================================\n")
cat("P90 MODEL COMPARISON\n")
cat("============================================\n")
print(p90_comparison)
threshold_qaqc <- data.frame(
  Check = c(
    "Three algorithms present",
    "Four threshold methods present",
    "All models have 150 runs",
    "Capture rates within 0-1",
    "No missing threshold values"
  ),
  Result = c(
    length(unique(threshold_comparison$Algorithm)) == 3,
    length(unique(threshold_comparison$Threshold_Method)) == 4,
    all(threshold_comparison$Run_n == 150),
    all(
      threshold_comparison$Capture_Rate >= 0 &
        threshold_comparison$Capture_Rate <= 1
    ),
    all(
      is.finite(threshold_comparison$Mean_Threshold) &
        is.finite(threshold_comparison$SD_Threshold)
    )
  ),
  stringsAsFactors = FALSE
)
print(threshold_qaqc)
if (!all(threshold_qaqc$Result)) {
  stop("Stage 4 threshold QA/QC failed.", call. = FALSE)
}
write.csv(
  threshold_comparison,
  file.path(comparison_output_dir, "RF_XGB_SVM_Threshold_Comparison.csv"),
  row.names = FALSE
)
write.csv(
  p90_comparison,
  file.path(comparison_output_dir, "RF_XGB_SVM_P90_Comparison.csv"),
  row.names = FALSE
)
message("")
message("============================================")
message("THRESHOLD CAPTURE AND STABILITY AUDIT")
message("============================================")
message("Threshold comparison completed successfully.")
message("Stage 4 PASS")
# ============================================================
# ============================================================
rf_validation_file <- file.path(
  models_dir,
  "RF_Results",
  "Independent_Validation",
  "RF_Independent_Validation_Predictions.csv"
)
xgb_validation_file <- file.path(
  models_dir,
  "XGB_Results",
  "Independent_Validation",
  "XGB_Independent_Validation_Predictions.csv"
)
svm_validation_file <- file.path(
  models_dir,
  "SVM_Results",
  "Independent_Validation",
  "SVM_Independent_Validation_Predictions.csv"
)
validation_files <- c(
  RF = rf_validation_file,
  XGB = xgb_validation_file,
  SVM = svm_validation_file
)
if (any(!file.exists(validation_files))) {
  stop(
    paste(
      "Missing independent-validation file(s):",
      paste(
        validation_files[!file.exists(validation_files)],
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}
rf_val <- read.csv(rf_validation_file, stringsAsFactors = FALSE)
xgb_val <- read.csv(xgb_validation_file, stringsAsFactors = FALSE)
svm_val <- read.csv(svm_validation_file, stringsAsFactors = FALSE)
required_validation_columns <- c(
  "Algorithm",
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "CellID",
  "Validation_Role",
  "Observed_Class",
  "Probability"
)
check_required_columns <- function(df, algorithm_name) {
  missing_columns <- setdiff(required_validation_columns, names(df))
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        algorithm_name,
        " validation table missing columns: ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
check_required_columns(rf_val, "RF")
check_required_columns(xgb_val, "XGB")
check_required_columns(svm_val, "SVM")
key_columns <- c(
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "CellID",
  "Validation_Role",
  "Observed_Class"
)
rf_keys <- rf_val[, key_columns, drop = FALSE]
xgb_keys <- xgb_val[, key_columns, drop = FALSE]
svm_keys <- svm_val[, key_columns, drop = FALSE]
validation_alignment_qaqc <- data.frame(
  Check = c(
    "RF and XGB validation samples identical",
    "RF and SVM validation samples identical",
    "RF has 15150 rows",
    "XGB has 15150 rows",
    "SVM has 15150 rows",
    "Each model has 150 positives",
    "Each model has 15000 pseudo-background negatives"
  ),
  Result = c(
    identical(rf_keys, xgb_keys),
    identical(rf_keys, svm_keys),
    nrow(rf_val) == 15150L,
    nrow(xgb_val) == 15150L,
    nrow(svm_val) == 15150L,
    all(
      c(
        sum(rf_val$Observed_Class == 1),
        sum(xgb_val$Observed_Class == 1),
        sum(svm_val$Observed_Class == 1)
      ) == 150L
    ),
    all(
      c(
        sum(rf_val$Observed_Class == 0),
        sum(xgb_val$Observed_Class == 0),
        sum(svm_val$Observed_Class == 0)
      ) == 15000L
    )
  ),
  stringsAsFactors = FALSE
)
print(validation_alignment_qaqc)
if (!all(validation_alignment_qaqc$Result)) {
  stop(
    "Independent-validation sample alignment QA/QC failed.",
    call. = FALSE
  )
}
calculate_logloss <- function(observed, probability) {
  epsilon <- 1e-15
  probability <- pmin(pmax(probability, epsilon), 1 - epsilon)
  -mean(
    observed * log(probability) +
      (1 - observed) * log(1 - probability)
  )
}
calculate_brier <- function(observed, probability) {
  mean((probability - observed)^2)
}
calculate_average_precision <- function(observed, probability) {
  ordering <- order(probability, decreasing = TRUE)
  y <- observed[ordering]
  cumulative_tp <- cumsum(y == 1)
  cumulative_fp <- cumsum(y == 0)
  precision <- cumulative_tp / (cumulative_tp + cumulative_fp)
  positive_positions <- which(y == 1)
  if (length(positive_positions) == 0) {
    return(NA_real_)
  }
  mean(precision[positive_positions])
}
evaluate_validation_model <- function(df, algorithm_name) {
  observed <- as.integer(df$Observed_Class)
  probability <- as.numeric(df$Probability)
  roc_object <- pROC::roc(
    response = observed,
    predictor = probability,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  positive_mean_probability <- mean(probability[observed == 1])
  background_mean_probability <- mean(probability[observed == 0])
  data.frame(
    Algorithm = algorithm_name,
    ROC_AUC = as.numeric(pROC::auc(roc_object)),
    Average_Precision = calculate_average_precision(observed, probability),
    Brier_Score = calculate_brier(observed, probability),
    Log_Loss = calculate_logloss(observed, probability),
    Mean_Positive_Probability = positive_mean_probability,
    Mean_Background_Probability = background_mean_probability,
    Mean_Probability_Separation =
      positive_mean_probability - background_mean_probability,
    stringsAsFactors = FALSE
  )
}
validation_metrics <- rbind(
  evaluate_validation_model(rf_val, "RF"),
  evaluate_validation_model(xgb_val, "XGB"),
  evaluate_validation_model(svm_val, "SVM")
)
cat("\n============================================\n")
cat("INDEPENDENT VALIDATION - THRESHOLD FREE\n")
cat("============================================\n")
print(validation_metrics)
write.csv(
  validation_metrics,
  file.path(
    comparison_output_dir,
    "RF_XGB_SVM_Independent_Validation_ThresholdFree.csv"
  ),
  row.names = FALSE
)
write.csv(
  validation_alignment_qaqc,
  file.path(
    comparison_output_dir,
    "RF_XGB_SVM_Independent_Validation_Alignment_QAQC.csv"
  ),
  row.names = FALSE
)
message("")
message("============================================")
message("INDEPENDENT VALIDATION THRESHOLD-FREE AUDIT")
message("============================================")
message("RF, XGB and SVM evaluated on identical validation samples.")
message("Stage 5A PASS")
# ============================================================
# ============================================================
rf_heldout_capture$Algorithm <- "RF"
xgb_heldout_capture$Algorithm <- "XGB"
svm_heldout_capture$Algorithm <- "SVM"
threshold_cols <- c(
  "Algorithm",
  "Run_ID",
  "Threshold_P90",
  "Threshold_P95",
  "Threshold_P99",
  "Threshold_Mean_Plus_2SD"
)
run_thresholds <- rbind(
  rf_heldout_capture[, threshold_cols, drop = FALSE],
  xgb_heldout_capture[, threshold_cols, drop = FALSE],
  svm_heldout_capture[, threshold_cols, drop = FALSE]
)
if (nrow(run_thresholds) != 450L) {
  stop(
    paste0(
      "Expected 450 algorithm-run threshold rows; found ",
      nrow(run_thresholds),
      "."
    ),
    call. = FALSE
  )
}
if (anyDuplicated(run_thresholds[, c("Algorithm", "Run_ID")]) > 0) {
  stop(
    "Duplicate Algorithm + Run_ID combinations in threshold table.",
    call. = FALSE
  )
}
validation_all <- rbind(rf_val, xgb_val, svm_val)
validation_thresholded <- merge(
  validation_all,
  run_thresholds,
  by = c("Algorithm", "Run_ID"),
  all.x = TRUE,
  sort = FALSE
)
if (nrow(validation_thresholded) != nrow(validation_all)) {
  stop("Threshold merge changed validation row count.", call. = FALSE)
}
threshold_columns <- c(
  "Threshold_P90",
  "Threshold_P95",
  "Threshold_P99",
  "Threshold_Mean_Plus_2SD"
)
if (anyNA(validation_thresholded[, threshold_columns, drop = FALSE])) {
  stop("Missing run-specific thresholds after merge.", call. = FALSE)
}
calculate_confusion_metrics <- function(
    observed,
    probability,
    threshold
) {
  predicted <- ifelse(probability >= threshold, 1L, 0L)
  TP <- sum(observed == 1L & predicted == 1L)
  FN <- sum(observed == 1L & predicted == 0L)
  TN <- sum(observed == 0L & predicted == 0L)
  FP <- sum(observed == 0L & predicted == 1L)
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  f1 <- if (
    is.finite(precision) &&
      is.finite(recall) &&
      (precision + recall) > 0
  ) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  balanced_accuracy <- if (
    is.finite(recall) &&
      is.finite(specificity)
  ) {
    (recall + specificity) / 2
  } else {
    NA_real_
  }
  background_false_positive_rate <- if ((FP + TN) > 0) {
    FP / (FP + TN)
  } else {
    NA_real_
  }
  data.frame(
    TP = TP,
    FN = FN,
    TN = TN,
    FP = FP,
    Recall = recall,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Balanced_Accuracy = balanced_accuracy,
    Predicted_Positive_Rate = mean(predicted == 1L),
    Background_False_Positive_Rate = background_false_positive_rate,
    stringsAsFactors = FALSE
  )
}
algorithm_names <- c("RF", "XGB", "SVM")
threshold_map <- c(
  P90 = "Threshold_P90",
  P95 = "Threshold_P95",
  P99 = "Threshold_P99",
  Mean_Plus_2SD = "Threshold_Mean_Plus_2SD"
)
classification_results <- list()
result_counter <- 1L
for (algorithm_name in algorithm_names) {
  model_data <- validation_thresholded[
    validation_thresholded$Algorithm == algorithm_name,
    ,
    drop = FALSE
  ]
  for (threshold_name in names(threshold_map)) {
    threshold_column <- unname(threshold_map[[threshold_name]])
    metric_result <- calculate_confusion_metrics(
      observed = as.integer(model_data$Observed_Class),
      probability = as.numeric(model_data$Probability),
      threshold = as.numeric(model_data[[threshold_column]])
    )
    metric_result$Algorithm <- algorithm_name
    metric_result$Threshold_Method <- threshold_name
    classification_results[[result_counter]] <- metric_result
    result_counter <- result_counter + 1L
  }
}
classification_metrics <- do.call(rbind, classification_results)
rownames(classification_metrics) <- NULL
classification_metrics <- classification_metrics[
  ,
  c(
    "Algorithm",
    "Threshold_Method",
    "TP",
    "FN",
    "TN",
    "FP",
    "Recall",
    "Specificity",
    "Precision",
    "F1",
    "Balanced_Accuracy",
    "Predicted_Positive_Rate",
    "Background_False_Positive_Rate"
  ),
  drop = FALSE
]
classification_metrics$Algorithm <- factor(
  classification_metrics$Algorithm,
  levels = c("RF", "XGB", "SVM")
)
classification_metrics$Threshold_Method <- factor(
  classification_metrics$Threshold_Method,
  levels = c("P90", "P95", "P99", "Mean_Plus_2SD")
)
classification_metrics <- classification_metrics[
  order(
    classification_metrics$Threshold_Method,
    classification_metrics$Algorithm
  ),
  ,
  drop = FALSE
]
classification_metrics$Algorithm <- as.character(
  classification_metrics$Algorithm
)
classification_metrics$Threshold_Method <- as.character(
  classification_metrics$Threshold_Method
)
cat("\n============================================\n")
cat("INDEPENDENT VALIDATION - CLASSIFICATION\n")
cat("============================================\n")
print(classification_metrics)
p90_classification <- classification_metrics[
  classification_metrics$Threshold_Method == "P90",
  ,
  drop = FALSE
]
p90_classification <- p90_classification[
  order(-p90_classification$Balanced_Accuracy),
  ,
  drop = FALSE
]
cat("\n============================================\n")
cat("P90 INDEPENDENT VALIDATION COMPARISON\n")
cat("============================================\n")
print(p90_classification)
classification_qaqc <- data.frame(
  Check = c(
    "12 model-threshold combinations",
    "Each combination contains 150 positives",
    "Each combination contains 15000 pseudo-background negatives",
    "All Recall values within 0-1",
    "All Specificity values within 0-1",
    "All Precision values within 0-1 where defined",
    "All F1 values within 0-1 where defined",
    "All Balanced Accuracy values within 0-1"
  ),
  Result = c(
    nrow(classification_metrics) == 12L,
    all(classification_metrics$TP + classification_metrics$FN == 150L),
    all(classification_metrics$TN + classification_metrics$FP == 15000L),
    all(
      classification_metrics$Recall >= 0 &
        classification_metrics$Recall <= 1
    ),
    all(
      classification_metrics$Specificity >= 0 &
        classification_metrics$Specificity <= 1
    ),
    all(
      is.na(classification_metrics$Precision) |
        (
          classification_metrics$Precision >= 0 &
            classification_metrics$Precision <= 1
        )
    ),
    all(
      is.na(classification_metrics$F1) |
        (
          classification_metrics$F1 >= 0 &
            classification_metrics$F1 <= 1
        )
    ),
    all(
      classification_metrics$Balanced_Accuracy >= 0 &
        classification_metrics$Balanced_Accuracy <= 1
    )
  ),
  stringsAsFactors = FALSE
)
print(classification_qaqc)
if (!all(classification_qaqc$Result)) {
  stop(
    "Independent-validation classification QA/QC failed.",
    call. = FALSE
  )
}
write.csv(
  classification_metrics,
  file.path(
    comparison_output_dir,
    "RF_XGB_SVM_Independent_Validation_Classification.csv"
  ),
  row.names = FALSE
)
write.csv(
  p90_classification,
  file.path(
    comparison_output_dir,
    "RF_XGB_SVM_P90_Independent_Validation.csv"
  ),
  row.names = FALSE
)
write.csv(
  classification_qaqc,
  file.path(
    comparison_output_dir,
    "RF_XGB_SVM_Independent_Validation_Classification_QAQC.csv"
  ),
  row.names = FALSE
)
message("")
message("============================================")
message("INDEPENDENT VALIDATION CLASSIFICATION AUDIT")
message("============================================")
message("Run-specific thresholds applied successfully.")
message("Stage 5B PASS")
# ============================================================
# ============================================================
message("")
message("============================================")
message("RF - XGB - SVM COMPARISON COMPLETE")
message("============================================")
message("Stage 1  : PASS")
message("Stage 2  : PASS")
message("Stage 3  : PASS")
message("Stage 4  : PASS")
message("Stage 5A : PASS")
message("Stage 5B : PASS")
message("")
message("Outputs saved to: ", comparison_output_dir)
# ============================================================
# Stage 6A:
# ============================================================
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_probability <- as.numeric(
  rf_fullgrid_summary$RF_Mean_Probability
)
if (
  length(rf_probability) != 9600L
) {
  stop(
    "RF full-grid probability count is not 9600."
  )
}
if (
  anyNA(rf_probability)
) {
  stop(
    "RF full-grid probabilities contain missing values."
  )
}
if (
  any(
    rf_probability < 0 |
      rf_probability > 1
  )
) {
  stop(
    "RF full-grid probabilities outside 0-1."
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_mean <- mean(
  rf_probability
)
rf_median <- median(
  rf_probability
)
rf_sd <- sd(
  rf_probability
)
rf_mad_raw <- median(
  abs(
    rf_probability -
      rf_median
  )
)
rf_mad_scaled <- mad(
  rf_probability,
  center = rf_median,
  constant = 1.4826
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_n <- length(
  rf_probability
)
rf_centered <- rf_probability -
  rf_mean
rf_moment2 <- mean(
  rf_centered^2
)
rf_moment3 <- mean(
  rf_centered^3
)
rf_moment4 <- mean(
  rf_centered^4
)
rf_skewness <- rf_moment3 /
  (
    rf_moment2^(3 / 2)
  )
rf_excess_kurtosis <- rf_moment4 /
  (
    rf_moment2^2
  ) -
  3
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_quantiles <- quantile(
  rf_probability,
  probs = c(
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99
  ),
  names = FALSE,
  type = 7
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_threshold_mean_2sd <-
  rf_mean +
  2 * rf_sd
rf_threshold_median_2mad_raw <-
  rf_median +
  2 * rf_mad_raw
rf_threshold_median_2mad_scaled <-
  rf_median +
  2 * rf_mad_scaled
rf_threshold_p90 <- quantile(
  rf_probability,
  probs = 0.90,
  names = FALSE,
  type = 7
)
rf_threshold_p95 <- quantile(
  rf_probability,
  probs = 0.95,
  names = FALSE,
  type = 7
)
rf_threshold_p99 <- quantile(
  rf_probability,
  probs = 0.99,
  names = FALSE,
  type = 7
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_distribution_summary <- data.frame(
  Statistic = c(
    "N",
    "Mean",
    "Median",
    "SD",
    "MAD_Raw",
    "MAD_Scaled_1.4826",
    "Skewness",
    "Excess_Kurtosis",
    "Minimum",
    "P01",
    "P05",
    "P10",
    "P25",
    "P50",
    "P75",
    "P90",
    "P95",
    "P99",
    "Maximum"
  ),
  Value = c(
    rf_n,
    rf_mean,
    rf_median,
    rf_sd,
    rf_mad_raw,
    rf_mad_scaled,
    rf_skewness,
    rf_excess_kurtosis,
    min(rf_probability),
    rf_quantiles[1],
    rf_quantiles[2],
    rf_quantiles[3],
    rf_quantiles[4],
    rf_quantiles[5],
    rf_quantiles[6],
    rf_quantiles[7],
    rf_quantiles[8],
    rf_quantiles[9],
    max(rf_probability)
  ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_distribution_thresholds <- data.frame(
  Threshold_Method = c(
    "P90",
    "P95",
    "P99",
    "Mean_Plus_2SD",
    "Median_Plus_2MAD_Raw",
    "Median_Plus_2MAD_Scaled"
  ),
  Threshold = c(
    rf_threshold_p90,
    rf_threshold_p95,
    rf_threshold_p99,
    rf_threshold_mean_2sd,
    rf_threshold_median_2mad_raw,
    rf_threshold_median_2mad_scaled
  ),
  stringsAsFactors = FALSE
)
rf_distribution_thresholds$Grid_Cell_n <-
  vapply(
    rf_distribution_thresholds$Threshold,
    function(current_threshold) {
      sum(
        rf_probability >=
          current_threshold
      )
    },
    integer(1)
  )
rf_distribution_thresholds$Grid_Percent <-
  100 *
  rf_distribution_thresholds$Grid_Cell_n /
  length(
    rf_probability
  )
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("RF FULL-GRID DISTRIBUTION AUDIT\n")
cat("============================================\n")
print(
  rf_distribution_summary,
  row.names = FALSE
)
cat("\n")
cat("============================================\n")
cat("RF CANDIDATE THRESHOLDS\n")
cat("============================================\n")
print(
  rf_distribution_thresholds,
  row.names = FALSE
)
# ------------------------------------------------------------
# 62. SAVE RESULTS
# ------------------------------------------------------------
write.csv(
  rf_distribution_summary,
  file.path(
    comparison_output_dir,
    "RF_FullGrid_Probability_Distribution.csv"
  ),
  row.names = FALSE
)
write.csv(
  rf_distribution_thresholds,
  file.path(
    comparison_output_dir,
    "RF_FullGrid_Candidate_Thresholds.csv"
  ),
  row.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
png(
  filename = file.path(
    comparison_output_dir,
    "RF_FullGrid_Probability_Distribution.png"
  ),
  width = 1800,
  height = 1200,
  res = 180
)
hist(
  rf_probability,
  breaks = 50,
  main = "RF Full-Grid Prospectivity Probability Distribution",
  xlab = "Mean RF probability",
  ylab = "Grid-cell frequency"
)
abline(
  v = rf_threshold_p90,
  lty = 2,
  lwd = 2
)
abline(
  v = rf_threshold_mean_2sd,
  lty = 3,
  lwd = 2
)
abline(
  v = rf_threshold_median_2mad_scaled,
  lty = 4,
  lwd = 2
)
legend(
  "topright",
  legend = c(
    "P90",
    "Mean + 2SD",
    "Median + 2MAD (scaled)"
  ),
  lty = c(
    2,
    3,
    4
  ),
  lwd = 2,
  bty = "n"
)
dev.off()
# ------------------------------------------------------------
# 64. QA/QC
# ------------------------------------------------------------
rf_distribution_qaqc <- data.frame(
  Check = c(
    "9600 RF grid probabilities",
    "No missing probabilities",
    "All probabilities within 0-1",
    "Mean is finite",
    "Median is finite",
    "SD is finite",
    "Raw MAD is finite",
    "Scaled MAD is finite",
    "Skewness is finite",
    "Excess kurtosis is finite",
    "All candidate thresholds are finite"
  ),
  Result = c(
    length(rf_probability) == 9600L,
    !anyNA(rf_probability),
    all(
      rf_probability >= 0 &
        rf_probability <= 1
    ),
    is.finite(rf_mean),
    is.finite(rf_median),
    is.finite(rf_sd),
    is.finite(rf_mad_raw),
    is.finite(rf_mad_scaled),
    is.finite(rf_skewness),
    is.finite(rf_excess_kurtosis),
    all(
      is.finite(
        rf_distribution_thresholds$Threshold
      )
    )
  ),
  stringsAsFactors = FALSE
)
cat("\n")
print(
  rf_distribution_qaqc,
  row.names = FALSE
)
if (
  !all(
    rf_distribution_qaqc$Result
  )
) {
  stop(
    "RF probability distribution QA/QC failed."
  )
}
write.csv(
  rf_distribution_qaqc,
  file.path(
    comparison_output_dir,
    "RF_FullGrid_Probability_Distribution_QAQC.csv"
  ),
  row.names = FALSE
)
message("")
message("============================================")
message("RF FULL-GRID DISTRIBUTION AUDIT COMPLETE")
message("============================================")
message("Stage 6A PASS")
