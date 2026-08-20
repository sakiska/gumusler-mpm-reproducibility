# ============================================================
# Evaluate_RF_MedianMAD_Threshold.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Evaluate the scaled Median + 2MAD operational threshold for the selected Random Forest model.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("ranger")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
# 2. SHARED MODULES
# ------------------------------------------------------------
source("Project/Models/Prepare_Fold_Data.R")
source("Project/Models/Prepare_Model_Matrix.R")
# ------------------------------------------------------------
# ------------------------------------------------------------
predictor_columns <- c(
  "Pb_OK",
  "Zn_OK",
  "Cu_OK",
  "Lithology",
  "Dist_Fault",
  "Dist_Silicified",
  "Dist_Brecciated"
)
categorical_predictors <- c("Lithology")
num_trees <- 1000L
mad_constant <- 1.4826
# ------------------------------------------------------------
# ------------------------------------------------------------
full_grid_file <- "Project/Data/MPM_FullGrid_25m.csv"
run_results_file <- paste0(
  "Project/Models/RF_Results/",
  "RF_LOOCV_Run_Results.csv"
)
validation_prediction_file <- paste0(
  "Project/Models/RF_Results/",
  "Independent_Validation/",
  "RF_Independent_Validation_Predictions.csv"
)
output_directory <- paste0(
  "Project/Models/RF_Results/",
  "MedianMAD_Validation"
)
dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  "Project/QAQC",
  recursive = TRUE,
  showWarnings = FALSE
)
run_threshold_file <- file.path(
  output_directory,
  "RF_MedianMAD_Run_Thresholds.csv"
)
run_metrics_file <- file.path(
  output_directory,
  "RF_MedianMAD_Run_Validation_Metrics.csv"
)
overall_metrics_file <- file.path(
  output_directory,
  "RF_MedianMAD_Overall_Validation_Metrics.csv"
)
error_log_file <- file.path(
  output_directory,
  "RF_MedianMAD_Error_Log.csv"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "RF_MedianMAD_Validation_QAQC.txt"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "RF_MedianMAD_Validation_Methodology.txt"
)
# ------------------------------------------------------------
# 5. HELPERS
# ------------------------------------------------------------
assert_required_columns <- function(
    data,
    required_columns,
    data_name
) {
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        data_name,
        " is missing required columns: ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
build_aligned_full_grid_matrix <- function(
    full_grid,
    predictor_columns,
    categorical_predictors,
    training_column_names
) {
  predictor_data <- full_grid[
    ,
    predictor_columns,
    drop = FALSE
  ]
  numeric_predictors <- setdiff(
    predictor_columns,
    categorical_predictors
  )
  for (column_name in numeric_predictors) {
    predictor_data[[column_name]] <- as.numeric(
      predictor_data[[column_name]]
    )
  }
  for (column_name in categorical_predictors) {
    predictor_data[[column_name]] <- factor(
      predictor_data[[column_name]]
    )
  }
  full_matrix <- stats::model.matrix(
    ~ . - 1,
    data = predictor_data,
    na.action = stats::na.pass
  )
  full_matrix <- as.data.frame(
    full_matrix,
    check.names = FALSE
  )
  missing_columns <- setdiff(
    training_column_names,
    names(full_matrix)
  )
  if (length(missing_columns) > 0) {
    for (column_name in missing_columns) {
      full_matrix[[column_name]] <- 0
    }
  }
  extra_columns <- setdiff(
    names(full_matrix),
    training_column_names
  )
  if (length(extra_columns) > 0) {
    full_matrix <- full_matrix[
      ,
      setdiff(names(full_matrix), extra_columns),
      drop = FALSE
    ]
  }
  full_matrix <- full_matrix[
    ,
    training_column_names,
    drop = FALSE
  ]
  full_matrix
}
calculate_median_mad_threshold <- function(
    probability,
    mad_constant = 1.4826
) {
  grid_median <- stats::median(
    probability
  )
  mad_raw <- stats::median(
    abs(
      probability - grid_median
    )
  )
  mad_scaled <- mad_constant * mad_raw
  threshold <- min(
    1,
    grid_median + 2 * mad_scaled
  )
  data.frame(
    Grid_Median = grid_median,
    Grid_MAD_Raw = mad_raw,
    Grid_MAD_Scaled = mad_scaled,
    MAD_Constant = mad_constant,
    Threshold_Median_Plus_2MAD = threshold,
    stringsAsFactors = FALSE
  )
}
calculate_confusion_metrics <- function(
    observed,
    probability,
    threshold
) {
  predicted <- ifelse(
    probability >= threshold,
    1L,
    0L
  )
  TP <- sum(observed == 1L & predicted == 1L)
  FN <- sum(observed == 1L & predicted == 0L)
  TN <- sum(observed == 0L & predicted == 0L)
  FP <- sum(observed == 0L & predicted == 1L)
  recall <- if ((TP + FN) > 0) {
    TP / (TP + FN)
  } else {
    NA_real_
  }
  specificity <- if ((TN + FP) > 0) {
    TN / (TN + FP)
  } else {
    NA_real_
  }
  precision <- if ((TP + FP) > 0) {
    TP / (TP + FP)
  } else {
    NA_real_
  }
  f1 <- if (
    is.finite(precision) &&
      is.finite(recall) &&
      (precision + recall) > 0
  ) {
    2 * precision * recall /
      (precision + recall)
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
    Background_False_Positive_Rate =
      if ((FP + TN) > 0) FP / (FP + TN) else NA_real_,
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  full_grid_file,
  run_results_file,
  validation_prediction_file
)
missing_files <- input_files[
  !file.exists(input_files)
]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input file(s):\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}
full_grid <- read.csv(
  full_grid_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
run_results <- read.csv(
  run_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
validation_predictions <- read.csv(
  validation_prediction_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
assert_required_columns(
  full_grid,
  c("CellID", "X", "Y", predictor_columns),
  "Full-grid table"
)
assert_required_columns(
  run_results,
  c(
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "Held_Out_CellID",
    "Best_mtry",
    "Best_min_node_size",
    "Occurrence_Class_Weight",
    "RF_Seed",
    "QAQC_Status"
  ),
  "RF run-results table"
)
assert_required_columns(
  validation_predictions,
  c(
    "Algorithm",
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "CellID",
    "Validation_Role",
    "Observed_Class",
    "Probability"
  ),
  "RF independent-validation predictions"
)
if (nrow(run_results) != 150L) {
  stop(
    paste0(
      "Expected 150 RF runs; found ",
      nrow(run_results),
      "."
    ),
    call. = FALSE
  )
}
if (any(run_results$QAQC_Status != "PASS")) {
  stop(
    "RF run-results contain failed runs.",
    call. = FALSE
  )
}
if (anyDuplicated(run_results$Run_ID) > 0) {
  stop(
    "RF run-results contain duplicated Run_ID values.",
    call. = FALSE
  )
}
if (anyDuplicated(full_grid$CellID) > 0) {
  stop(
    "Full-grid contains duplicated CellID values.",
    call. = FALSE
  )
}
if (
  nrow(validation_predictions) != 15150L ||
    sum(validation_predictions$Observed_Class == 1L) != 150L ||
    sum(validation_predictions$Observed_Class == 0L) != 15000L
) {
  stop(
    "Independent-validation prediction structure is invalid.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 8. STORAGE
# ------------------------------------------------------------
threshold_rows <- vector(
  "list",
  nrow(run_results)
)
run_metric_rows <- vector(
  "list",
  nrow(run_results)
)
error_rows <- list()
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("RF MEDIAN + 2MAD VALIDATION STARTED\n")
cat("============================================\n")
cat("Runs:", nrow(run_results), "\n")
cat("Grid cells:", nrow(full_grid), "\n\n")
for (run_i in seq_len(nrow(run_results))) {
  current_run <- run_results[
    run_i,
    ,
    drop = FALSE
  ]
  run_id <- current_run$Run_ID[1]
  cat(
    "[",
    run_i,
    "/",
    nrow(run_results),
    "] ",
    run_id,
    "\n",
    sep = ""
  )
  current_output <- tryCatch({
    # --------------------------------------------------------
    # --------------------------------------------------------
    fold_data <- prepare_fold_data(
      run_id = run_id
    )
    model_matrix <- prepare_model_matrix(
      fold_data = fold_data,
      predictor_columns = predictor_columns,
      categorical_predictors = categorical_predictors,
      scale_numeric = FALSE
    )
    x_train <- as.data.frame(
      model_matrix$x_train,
      check.names = FALSE
    )
    training_rf_data <- x_train
    training_rf_data$Class <-
      model_matrix$y_train_factor
    # --------------------------------------------------------
    # --------------------------------------------------------
    full_grid_matrix <- build_aligned_full_grid_matrix(
      full_grid = full_grid,
      predictor_columns = predictor_columns,
      categorical_predictors = categorical_predictors,
      training_column_names = names(x_train)
    )
    if (
      anyNA(full_grid_matrix) ||
        any(
          !is.finite(
            as.matrix(full_grid_matrix)
          )
        )
    ) {
      stop(
        "Full-grid model matrix contains invalid values.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # 9.3 Reconstruct accepted RF model
    # --------------------------------------------------------
    rf_class_weights <- c(
      PseudoBackground = 1,
      Occurrence = as.numeric(
        current_run$Occurrence_Class_Weight[1]
      )
    )
    final_model <- ranger::ranger(
      formula = Class ~ .,
      data = training_rf_data,
      num.trees = num_trees,
      mtry = as.integer(
        current_run$Best_mtry[1]
      ),
      min.node.size = as.integer(
        current_run$Best_min_node_size[1]
      ),
      probability = TRUE,
      class.weights = rf_class_weights,
      importance = "none",
      seed = as.integer(
        current_run$RF_Seed[1]
      ),
      num.threads = 1,
      write.forest = TRUE,
      verbose = FALSE
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    prediction <- predict(
      final_model,
      data = full_grid_matrix
    )$predictions
    if (!"Occurrence" %in% colnames(prediction)) {
      stop(
        "Occurrence probability column was not returned.",
        call. = FALSE
      )
    }
    full_probability <- as.numeric(
      prediction[, "Occurrence"]
    )
    if (
      length(full_probability) != nrow(full_grid) ||
        any(!is.finite(full_probability)) ||
        any(full_probability < 0) ||
        any(full_probability > 1)
    ) {
      stop(
        "Invalid full-grid probability vector.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    threshold_row <- calculate_median_mad_threshold(
      probability = full_probability,
      mad_constant = mad_constant
    )
    threshold_row$Algorithm <- "RF"
    threshold_row$Run_ID <- run_id
    threshold_row$Repeat <- current_run$Repeat[1]
    threshold_row$Fold_ID <- current_run$Fold_ID[1]
    threshold_row$Fold_Number <- current_run$Fold_Number[1]
    threshold_row <- threshold_row[
      ,
      c(
        "Algorithm",
        "Run_ID",
        "Repeat",
        "Fold_ID",
        "Fold_Number",
        "Grid_Median",
        "Grid_MAD_Raw",
        "Grid_MAD_Scaled",
        "MAD_Constant",
        "Threshold_Median_Plus_2MAD"
      ),
      drop = FALSE
    ]
    # --------------------------------------------------------
    # --------------------------------------------------------
    current_validation <- validation_predictions[
      validation_predictions$Run_ID == run_id,
      ,
      drop = FALSE
    ]
    if (
      nrow(current_validation) != 101L ||
        sum(current_validation$Observed_Class == 1L) != 1L ||
        sum(current_validation$Observed_Class == 0L) != 100L
    ) {
      stop(
        "Run-specific independent-validation structure is invalid.",
        call. = FALSE
      )
    }
    metric_row <- calculate_confusion_metrics(
      observed = as.integer(
        current_validation$Observed_Class
      ),
      probability = as.numeric(
        current_validation$Probability
      ),
      threshold = threshold_row$Threshold_Median_Plus_2MAD[1]
    )
    metric_row$Algorithm <- "RF"
    metric_row$Run_ID <- run_id
    metric_row$Repeat <- current_run$Repeat[1]
    metric_row$Fold_ID <- current_run$Fold_ID[1]
    metric_row$Fold_Number <- current_run$Fold_Number[1]
    metric_row$Threshold_Method <- "Median_Plus_2MAD_Scaled"
    metric_row$Threshold <-
      threshold_row$Threshold_Median_Plus_2MAD[1]
    metric_row <- metric_row[
      ,
      c(
        "Algorithm",
        "Run_ID",
        "Repeat",
        "Fold_ID",
        "Fold_Number",
        "Threshold_Method",
        "Threshold",
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
    list(
      threshold = threshold_row,
      metric = metric_row
    )
  }, error = function(e) {
    list(
      error = data.frame(
        Algorithm = "RF",
        Run_ID = run_id,
        Repeat = current_run$Repeat[1],
        Fold_ID = current_run$Fold_ID[1],
        Error_Message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  })
  if ("error" %in% names(current_output)) {
    error_rows[[length(error_rows) + 1L]] <-
      current_output$error
    cat(
      "   FAIL: ",
      current_output$error$Error_Message,
      "\n",
      sep = ""
    )
  } else {
    threshold_rows[[run_i]] <-
      current_output$threshold
    run_metric_rows[[run_i]] <-
      current_output$metric
    cat("   PASS\n")
  }
}
# ------------------------------------------------------------
# 10. ERROR LOG
# ------------------------------------------------------------
error_df <- if (length(error_rows) > 0) {
  do.call(rbind, error_rows)
} else {
  data.frame(
    Algorithm = character(0),
    Run_ID = character(0),
    Repeat = integer(0),
    Fold_ID = character(0),
    Error_Message = character(0),
    stringsAsFactors = FALSE
  )
}
write.csv(
  error_df,
  error_log_file,
  row.names = FALSE,
  na = ""
)
if (nrow(error_df) > 0) {
  stop(
    paste0(
      nrow(error_df),
      " Median+2MAD RF run(s) failed. See: ",
      error_log_file
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
run_thresholds <- do.call(
  rbind,
  threshold_rows
)
run_metrics <- do.call(
  rbind,
  run_metric_rows
)
rownames(run_thresholds) <- NULL
rownames(run_metrics) <- NULL
# ------------------------------------------------------------
#
# ------------------------------------------------------------
TP_total <- sum(run_metrics$TP)
FN_total <- sum(run_metrics$FN)
TN_total <- sum(run_metrics$TN)
FP_total <- sum(run_metrics$FP)
recall_total <- TP_total / (TP_total + FN_total)
specificity_total <- TN_total / (TN_total + FP_total)
precision_total <- if ((TP_total + FP_total) > 0) {
  TP_total / (TP_total + FP_total)
} else {
  NA_real_
}
f1_total <- if (
  is.finite(precision_total) &&
    is.finite(recall_total) &&
    (precision_total + recall_total) > 0
) {
  2 * precision_total * recall_total /
    (precision_total + recall_total)
} else {
  NA_real_
}
balanced_accuracy_total <- (
  recall_total + specificity_total
) / 2
overall_metrics <- data.frame(
  Algorithm = "RF",
  Threshold_Method = "Median_Plus_2MAD_Scaled",
  TP = TP_total,
  FN = FN_total,
  TN = TN_total,
  FP = FP_total,
  Recall = recall_total,
  Specificity = specificity_total,
  Precision = precision_total,
  F1 = f1_total,
  Balanced_Accuracy = balanced_accuracy_total,
  Predicted_Positive_Rate =
    (TP_total + FP_total) /
    (TP_total + FN_total + TN_total + FP_total),
  Background_False_Positive_Rate =
    FP_total / (FP_total + TN_total),
  Mean_Threshold = mean(
    run_thresholds$Threshold_Median_Plus_2MAD
  ),
  SD_Threshold = stats::sd(
    run_thresholds$Threshold_Median_Plus_2MAD
  ),
  Threshold_CV =
    stats::sd(
      run_thresholds$Threshold_Median_Plus_2MAD
    ) /
    mean(
      run_thresholds$Threshold_Median_Plus_2MAD
    ),
  Minimum_Threshold = min(
    run_thresholds$Threshold_Median_Plus_2MAD
  ),
  Maximum_Threshold = max(
    run_thresholds$Threshold_Median_Plus_2MAD
  ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# 13. QA/QC
# ------------------------------------------------------------
qaqc <- data.frame(
  Check = c(
    "150 successful RF runs",
    "150 run-specific Median+2MAD thresholds",
    "150 run-level validation metric rows",
    "Total positives equal 150",
    "Total pseudo-background negatives equal 15000",
    "All thresholds finite and within 0-1",
    "All MAD values non-negative",
    "Recall within 0-1",
    "Specificity within 0-1",
    "Balanced Accuracy within 0-1",
    "No run errors"
  ),
  Result = c(
    nrow(run_thresholds) == 150L,
    length(run_thresholds$Threshold_Median_Plus_2MAD) == 150L,
    nrow(run_metrics) == 150L,
    (TP_total + FN_total) == 150L,
    (TN_total + FP_total) == 15000L,
    all(
      is.finite(run_thresholds$Threshold_Median_Plus_2MAD) &
        run_thresholds$Threshold_Median_Plus_2MAD >= 0 &
        run_thresholds$Threshold_Median_Plus_2MAD <= 1
    ),
    all(
      run_thresholds$Grid_MAD_Raw >= 0 &
        run_thresholds$Grid_MAD_Scaled >= 0
    ),
    is.finite(recall_total) &&
      recall_total >= 0 &&
      recall_total <= 1,
    is.finite(specificity_total) &&
      specificity_total >= 0 &&
      specificity_total <= 1,
    is.finite(balanced_accuracy_total) &&
      balanced_accuracy_total >= 0 &&
      balanced_accuracy_total <= 1,
    nrow(error_df) == 0L
  ),
  stringsAsFactors = FALSE
)
print(qaqc)
if (!all(qaqc$Result)) {
  stop(
    "RF Median+2MAD validation QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------
write.csv(
  run_thresholds,
  run_threshold_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  run_metrics,
  run_metrics_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  overall_metrics,
  overall_metrics_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 15. QA/QC + METHODOLOGY NOTES
# ------------------------------------------------------------
qaqc_lines <- c(
  "RF MEDIAN + 2MAD VALIDATION QA/QC",
  "============================================",
  "",
  paste(qaqc$Check, qaqc$Result, sep = ": "),
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  qaqc_file
)
methodology_lines <- c(
  "RF MEDIAN + 2MAD VALIDATION METHODOLOGY",
  "============================================",
  "",
  paste(
    "For each of the 150 accepted RF outer-LOOCV runs,",
    "the model was reconstructed using the original",
    "fold-specific training data, selected hyperparameters,",
    "class weights and RF seed."
  ),
  "",
  paste(
    "Each reconstructed model predicted the complete 25 m grid."
  ),
  "",
  paste(
    "A run-specific robust threshold was calculated as",
    "Median + 2 * MAD_scaled, where MAD_scaled equals",
    "1.4826 times the median absolute deviation from the",
    "run-specific grid median."
  ),
  "",
  paste(
    "The held-out occurrence and independent validation",
    "pseudo-background cells were not used to calculate",
    "the threshold."
  ),
  "",
  paste(
    "Each run-specific threshold was then applied to the",
    "existing leakage-free independent validation predictions",
    "for the same run."
  ),
  "",
  paste(
    "Pseudo-background cells are not confirmed absences;",
    "classification metrics quantify discrimination against",
    "pseudo-background."
  )
)
writeLines(
  methodology_lines,
  methodology_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("RF MEDIAN + 2MAD VALIDATION COMPLETE\n")
cat("============================================\n")
cat("Successful runs           :", nrow(run_thresholds), "\n")
cat("MAD scaling constant      :", mad_constant, "\n")
cat(
  "Mean run threshold        :",
  overall_metrics$Mean_Threshold,
  "\n"
)
cat(
  "Threshold CV              :",
  overall_metrics$Threshold_CV,
  "\n"
)
cat("\n")
cat("OVERALL INDEPENDENT VALIDATION METRICS\n")
print(overall_metrics)
cat("\n")
print(qaqc)
cat("\n")
cat("RF MEDIAN + 2MAD VALIDATION: PASS\n")
