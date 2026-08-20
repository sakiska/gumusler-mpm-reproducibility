# ============================================================
# Predict_RF_FullGrid.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Apply accepted Random Forest models to the complete modeling grid.
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
      paste(missing_packages, collapse = ", "),
      ". Install with install.packages()."
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 2. Shared modules
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
num_trees <- 1000
rf_base_seed <- 5000
# ------------------------------------------------------------
# ------------------------------------------------------------
full_grid_file <- "Project/Data/MPM_FullGrid_25m.csv"
run_results_file <- paste0(
  "Project/Models/RF_Results/",
  "RF_LOOCV_Run_Results.csv"
)
output_directory <- paste0(
  "Project/Models/RF_Results/",
  "FullGrid"
)
dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
full_grid_summary_file <- paste0(
  output_directory,
  "/RF_FullGrid_Ensemble_Summary.csv"
)
run_threshold_file <- paste0(
  output_directory,
  "/RF_FullGrid_Run_Thresholds.csv"
)
held_out_capture_file <- paste0(
  output_directory,
  "/RF_HeldOut_Threshold_Capture.csv"
)
threshold_summary_file <- paste0(
  output_directory,
  "/RF_Threshold_Capture_Summary.csv"
)
error_log_file <- paste0(
  output_directory,
  "/RF_FullGrid_Error_Log.csv"
)
qaqc_file <- "Project/QAQC/RF_FullGrid_QAQC.txt"
methodology_file <- "Project/QAQC/RF_FullGrid_Methodology.txt"
# ------------------------------------------------------------
# ------------------------------------------------------------
assert_required_columns <- function(data, required_columns, data_name) {
  missing_columns <- setdiff(required_columns, names(data))
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
calculate_thresholds <- function(probability) {
  mean_probability <- mean(probability)
  sd_probability <- stats::sd(probability)
  mean_plus_2sd <- min(
    1,
    mean_probability + 2 * sd_probability
  )
  data.frame(
    Grid_Mean = mean_probability,
    Grid_SD = sd_probability,
    Threshold_P90 = as.numeric(
      stats::quantile(
        probability,
        probs = 0.90,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_P95 = as.numeric(
      stats::quantile(
        probability,
        probs = 0.95,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_P99 = as.numeric(
      stats::quantile(
        probability,
        probs = 0.99,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_Mean_Plus_2SD = mean_plus_2sd,
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  full_grid_file,
  run_results_file
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
    "Held_Out_Occurrence_Probability",
    "Best_mtry",
    "Best_min_node_size",
    "Occurrence_Class_Weight",
    "RF_Seed",
    "QAQC_Status"
  ),
  "RF run-results table"
)
if (anyDuplicated(full_grid$CellID) > 0) {
  stop(
    "Full-grid table contains duplicated CellID values.",
    call. = FALSE
  )
}
if (anyDuplicated(run_results$Run_ID) > 0) {
  stop(
    "RF run-results table contains duplicated Run_ID values.",
    call. = FALSE
  )
}
if (any(run_results$QAQC_Status != "PASS")) {
  stop(
    "RF run-results table contains failed runs.",
    call. = FALSE
  )
}
if (nrow(run_results) != 150) {
  stop(
    paste0(
      "Expected 150 accepted RF runs; found ",
      nrow(run_results),
      "."
    ),
    call. = FALSE
  )
}
if (anyNA(full_grid[, c("CellID", "X", "Y", predictor_columns)])) {
  stop(
    "Full-grid table contains missing values in required fields.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
grid_n <- nrow(full_grid)
run_n <- nrow(run_results)
probability_sum <- numeric(grid_n)
probability_sum_squares <- numeric(grid_n)
probability_min <- rep(Inf, grid_n)
probability_max <- rep(-Inf, grid_n)
prediction_count <- integer(grid_n)
all_run_thresholds <- vector("list", run_n)
all_held_out_capture <- vector("list", run_n)
all_errors <- list()
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\nRF full-grid prediction started.\n")
cat("Runs:", run_n, "\n")
cat("Grid cells:", grid_n, "\n\n")
for (run_i in seq_len(run_n)) {
  current_run <- run_results[run_i, , drop = FALSE]
  current_run_id <- current_run$Run_ID[1]
  cat(
    "[",
    run_i,
    "/",
    run_n,
    "] ",
    current_run_id,
    "\n",
    sep = ""
  )
  current_output <- tryCatch(
    {
      fold_data <- prepare_fold_data(
        run_id = current_run_id
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
      training_rf_data$Class <- model_matrix$y_train_factor
      full_grid_matrix <- build_aligned_full_grid_matrix(
        full_grid = full_grid,
        predictor_columns = predictor_columns,
        categorical_predictors = categorical_predictors,
        training_column_names = names(x_train)
      )
      if (anyNA(full_grid_matrix)) {
        stop(
          "Full-grid model matrix contains missing values.",
          call. = FALSE
        )
      }
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
      probability <- as.numeric(
        prediction[, "Occurrence"]
      )
      if (
        length(probability) != grid_n ||
        any(!is.finite(probability)) ||
        any(probability < 0) ||
        any(probability > 1)
      ) {
        stop(
          "Invalid full-grid probability vector.",
          call. = FALSE
        )
      }
      threshold_row <- calculate_thresholds(
        probability
      )
      threshold_row$Algorithm <- "RF"
      threshold_row$Run_ID <- current_run_id
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
          "Grid_Mean",
          "Grid_SD",
          "Threshold_P90",
          "Threshold_P95",
          "Threshold_P99",
          "Threshold_Mean_Plus_2SD"
        ),
        drop = FALSE
      ]
      held_out_probability <- as.numeric(
        current_run$
          Held_Out_Occurrence_Probability[1]
      )
      capture_row <- data.frame(
        Algorithm = "RF",
        Run_ID = current_run_id,
        Repeat = current_run$Repeat[1],
        Fold_ID = current_run$Fold_ID[1],
        Fold_Number = current_run$Fold_Number[1],
        Held_Out_CellID =
          current_run$Held_Out_CellID[1],
        Held_Out_Probability =
          held_out_probability,
        Threshold_P90 =
          threshold_row$Threshold_P90,
        Captured_P90 =
          held_out_probability >=
          threshold_row$Threshold_P90,
        Threshold_P95 =
          threshold_row$Threshold_P95,
        Captured_P95 =
          held_out_probability >=
          threshold_row$Threshold_P95,
        Threshold_P99 =
          threshold_row$Threshold_P99,
        Captured_P99 =
          held_out_probability >=
          threshold_row$Threshold_P99,
        Threshold_Mean_Plus_2SD =
          threshold_row$
            Threshold_Mean_Plus_2SD,
        Captured_Mean_Plus_2SD =
          held_out_probability >=
          threshold_row$
            Threshold_Mean_Plus_2SD,
        stringsAsFactors = FALSE
      )
      list(
        probability = probability,
        threshold = threshold_row,
        capture = capture_row
      )
    },
    error = function(e) {
      list(
        error = data.frame(
          Algorithm = "RF",
          Run_ID = current_run_id,
          Repeat = current_run$Repeat[1],
          Fold_ID = current_run$Fold_ID[1],
          Error_Message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      )
    }
  )
  if ("error" %in% names(current_output)) {
    all_errors[[length(all_errors) + 1]] <-
      current_output$error
    cat(
      "   FAIL: ",
      current_output$error$Error_Message,
      "\n",
      sep = ""
    )
  } else {
    probability <- current_output$probability
    probability_sum <- probability_sum + probability
    probability_sum_squares <-
      probability_sum_squares + probability^2
    probability_min <- pmin(
      probability_min,
      probability
    )
    probability_max <- pmax(
      probability_max,
      probability
    )
    prediction_count <- prediction_count + 1L
    all_run_thresholds[[run_i]] <-
      current_output$threshold
    all_held_out_capture[[run_i]] <-
      current_output$capture
    cat("   PASS\n")
  }
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (length(all_errors) > 0) {
  error_results <- do.call(
    rbind,
    all_errors
  )
  write.csv(
    error_results,
    error_log_file,
    row.names = FALSE,
    na = ""
  )
  stop(
    paste0(
      length(all_errors),
      " RF full-grid run(s) failed. See: ",
      error_log_file
    ),
    call. = FALSE
  )
} else {
  error_results <- data.frame(
    Algorithm = character(0),
    Run_ID = character(0),
    Repeat = integer(0),
    Fold_ID = character(0),
    Error_Message = character(0),
    stringsAsFactors = FALSE
  )
  write.csv(
    error_results,
    error_log_file,
    row.names = FALSE,
    na = ""
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (any(prediction_count != run_n)) {
  stop(
    "At least one grid cell does not have 150 predictions.",
    call. = FALSE
  )
}
mean_probability <- probability_sum / prediction_count
variance_probability <- (
  probability_sum_squares -
    (probability_sum^2 / prediction_count)
) / (prediction_count - 1)
variance_probability <- pmax(
  variance_probability,
  0
)
sd_probability <- sqrt(
  variance_probability
)
full_grid_summary <- data.frame(
  CellID = full_grid$CellID,
  X = full_grid$X,
  Y = full_grid$Y,
  RF_Mean_Probability = mean_probability,
  RF_SD_Probability = sd_probability,
  RF_Min_Probability = probability_min,
  RF_Max_Probability = probability_max,
  RF_Prediction_n = prediction_count,
  stringsAsFactors = FALSE
)
write.csv(
  full_grid_summary,
  full_grid_summary_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 11. Run thresholds and held-out capture
# ------------------------------------------------------------
run_thresholds <- do.call(
  rbind,
  all_run_thresholds
)
held_out_capture <- do.call(
  rbind,
  all_held_out_capture
)
write.csv(
  run_thresholds,
  run_threshold_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  held_out_capture,
  held_out_capture_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 12. Threshold-capture summary
# ------------------------------------------------------------
threshold_names <- c(
  "P90",
  "P95",
  "P99",
  "Mean_Plus_2SD"
)
capture_columns <- c(
  P90 = "Captured_P90",
  P95 = "Captured_P95",
  P99 = "Captured_P99",
  Mean_Plus_2SD =
    "Captured_Mean_Plus_2SD"
)
threshold_columns <- c(
  P90 = "Threshold_P90",
  P95 = "Threshold_P95",
  P99 = "Threshold_P99",
  Mean_Plus_2SD =
    "Threshold_Mean_Plus_2SD"
)
threshold_summary_list <- vector(
  "list",
  length(threshold_names)
)
for (threshold_i in seq_along(threshold_names)) {
  threshold_name <- threshold_names[threshold_i]
  capture_column <- capture_columns[[threshold_name]]
  threshold_column <- threshold_columns[[threshold_name]]
  capture_values <- held_out_capture[[capture_column]]
  threshold_values <- held_out_capture[[threshold_column]]
  threshold_summary_list[[threshold_i]] <- data.frame(
    Algorithm = "RF",
    Threshold_Method = threshold_name,
    Run_n = nrow(held_out_capture),
    Captured_Run_n = sum(capture_values),
    Capture_Rate = mean(capture_values),
    Mean_Threshold = mean(threshold_values),
    SD_Threshold = stats::sd(threshold_values),
    Minimum_Threshold = min(threshold_values),
    Maximum_Threshold = max(threshold_values),
    stringsAsFactors = FALSE
  )
}
threshold_summary <- do.call(
  rbind,
  threshold_summary_list
)
write.csv(
  threshold_summary,
  threshold_summary_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 13. QA/QC
# ------------------------------------------------------------
successful_run_n <- nrow(run_thresholds)
invalid_summary_probability_n <- sum(
  !is.finite(
    full_grid_summary$RF_Mean_Probability
  ) |
    full_grid_summary$RF_Mean_Probability < 0 |
    full_grid_summary$RF_Mean_Probability > 1
)
invalid_summary_sd_n <- sum(
  !is.finite(
    full_grid_summary$RF_SD_Probability
  ) |
    full_grid_summary$RF_SD_Probability < 0
)
capture_record_n <- nrow(
  held_out_capture
)
all_checks_passed <- all(
  successful_run_n == 150,
  capture_record_n == 150,
  nrow(full_grid_summary) == nrow(full_grid),
  all(full_grid_summary$RF_Prediction_n == 150),
  invalid_summary_probability_n == 0,
  invalid_summary_sd_n == 0,
  nrow(error_results) == 0
)
overall_qaqc_status <- ifelse(
  all_checks_passed,
  "PASS",
  "FAIL"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "RANDOM FOREST FULL-GRID PREDICTION METHODOLOGY",
  "",
  paste0("Accepted RF runs: ", successful_run_n),
  paste0("Grid cells: ", nrow(full_grid_summary)),
  paste0("Trees per model: ", num_trees),
  "",
  paste(
    "Each of the 150 accepted LOOCV RF models was",
    "reconstructed using the same fold-specific training data,",
    "selected hyperparameters, class weights and RF seed."
  ),
  paste(
    "Each reconstructed model predicted every valid 25 m",
    "full-grid cell."
  ),
  paste(
    "Cell-level ensemble mean, standard deviation, minimum",
    "and maximum occurrence probabilities were calculated."
  ),
  paste(
    "P90, P95, P99 and Mean + 2 SD thresholds were calculated",
    "separately from each run's full-grid probability distribution."
  ),
  paste(
    "Threshold capture was evaluated using only the genuinely",
    "held-out occurrence probability from the corresponding run."
  ),
  paste(
    "Ensemble map probabilities at known occurrence cells were",
    "not used as validation values because most ensemble models",
    "were trained with those occurrences."
  )
)
writeLines(
  methodology_lines,
  methodology_file,
  useBytes = TRUE
)
qaqc_lines <- c(
  "RANDOM FOREST FULL-GRID QAQC",
  "",
  paste0("Expected runs: 150"),
  paste0("Successful runs: ", successful_run_n),
  paste0("Failed runs: ", nrow(error_results)),
  paste0("Grid cells: ", nrow(full_grid_summary)),
  paste0(
    "Cells with 150 predictions: ",
    sum(full_grid_summary$RF_Prediction_n == 150)
  ),
  paste0(
    "Invalid mean probabilities: ",
    invalid_summary_probability_n
  ),
  paste0(
    "Invalid probability SD values: ",
    invalid_summary_sd_n
  ),
  paste0(
    "Held-out capture records: ",
    capture_record_n
  ),
  paste0(
    "Overall QAQC status: ",
    overall_qaqc_status
  )
)
writeLines(
  qaqc_lines,
  qaqc_file,
  useBytes = TRUE
)
# ------------------------------------------------------------
# 15. Console summary
# ------------------------------------------------------------
cat("\n========================================\n")
cat("RF full-grid prediction completed.\n")
cat("Successful runs:", successful_run_n, "\n")
cat("Grid cells:", nrow(full_grid_summary), "\n")
cat("Overall QAQC:", overall_qaqc_status, "\n")
cat("Full-grid summary:", full_grid_summary_file, "\n")
cat("Run thresholds:", run_threshold_file, "\n")
cat("Held-out capture:", held_out_capture_file, "\n")
cat("Threshold summary:", threshold_summary_file, "\n")
cat("========================================\n")
