# ============================================================
# Predict_XGBoost_FullGrid.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Apply accepted XGBoost models to the complete modeling grid.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("xgboost")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop(paste0("Missing package(s): ", paste(missing_packages, collapse = ", "), ". Install with install.packages()."), call. = FALSE)
source("Project/Models/Prepare_Fold_Data.R")
source("Project/Models/Prepare_Model_Matrix.R")
predictor_columns <- c("Pb_OK", "Zn_OK", "Cu_OK", "Lithology", "Dist_Fault", "Dist_Silicified", "Dist_Brecciated")
categorical_predictors <- c("Lithology")
full_grid_file <- "Project/Data/MPM_FullGrid_25m.csv"
run_results_file <- "Project/Models/XGB_Results/XGB_LOOCV_Run_Results.csv"
output_directory <- "Project/Models/XGB_Results/FullGrid"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create("Project/QAQC", recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(output_directory, "XGB_FullGrid_Ensemble_Summary.csv")
threshold_file <- file.path(output_directory, "XGB_FullGrid_Run_Thresholds.csv")
capture_file <- file.path(output_directory, "XGB_HeldOut_Threshold_Capture.csv")
capture_summary_file <- file.path(output_directory, "XGB_Threshold_Capture_Summary.csv")
error_file <- file.path(output_directory, "XGB_FullGrid_Error_Log.csv")
if (!file.exists(full_grid_file)) stop(paste0("Missing file: ", full_grid_file), call. = FALSE)
if (!file.exists(run_results_file)) stop(paste0("Missing file: ", run_results_file), call. = FALSE)
full_grid <- read.csv(full_grid_file, stringsAsFactors = FALSE, check.names = FALSE)
run_results <- read.csv(run_results_file, stringsAsFactors = FALSE, check.names = FALSE)
required_grid_columns <- c("CellID", "X", "Y", predictor_columns)
missing_grid_columns <- setdiff(required_grid_columns, names(full_grid))
if (length(missing_grid_columns) > 0) stop(paste0("Full grid missing columns: ", paste(missing_grid_columns, collapse = ", ")), call. = FALSE)
required_run_columns <- c("Run_ID", "Repeat", "Fold_ID", "Fold_Number", "Held_Out_CellID", "Held_Out_Occurrence_Probability", "eta", "max_depth", "min_child_weight", "subsample", "colsample_bytree", "gamma", "scale_pos_weight", "Best_nrounds", "XGB_Seed", "QAQC_Status")
missing_run_columns <- setdiff(required_run_columns, names(run_results))
if (length(missing_run_columns) > 0) stop(paste0("Run results missing columns: ", paste(missing_run_columns, collapse = ", ")), call. = FALSE)
if (nrow(run_results) != 150L || any(run_results$QAQC_Status != "PASS")) stop("Expected 150 accepted XGBoost runs.", call. = FALSE)
if (anyDuplicated(full_grid$CellID) > 0) stop("Full grid contains duplicated CellID values.", call. = FALSE)
if (anyNA(full_grid[, required_grid_columns])) stop("Full grid contains missing required values.", call. = FALSE)
build_full_grid_matrix <- function(full_grid, training_column_names) {
  predictor_data <- full_grid[, predictor_columns, drop = FALSE]
  numeric_predictors <- setdiff(predictor_columns, categorical_predictors)
  for (column_name in numeric_predictors) predictor_data[[column_name]] <- as.numeric(predictor_data[[column_name]])
  for (column_name in categorical_predictors) predictor_data[[column_name]] <- factor(predictor_data[[column_name]])
  matrix_data <- as.data.frame(stats::model.matrix(~ . - 1, data = predictor_data, na.action = stats::na.pass), check.names = FALSE)
  missing_columns <- setdiff(training_column_names, names(matrix_data))
  if (length(missing_columns) > 0) for (column_name in missing_columns) matrix_data[[column_name]] <- 0
  matrix_data <- matrix_data[, training_column_names, drop = FALSE]
  as.matrix(matrix_data)
}
calculate_thresholds <- function(probability) {
  mean_probability <- mean(probability)
  sd_probability <- sd(probability)
  data.frame(
    Grid_Mean = mean_probability,
    Grid_SD = sd_probability,
    Threshold_P90 = as.numeric(quantile(probability, 0.90, names = FALSE, type = 7)),
    Threshold_P95 = as.numeric(quantile(probability, 0.95, names = FALSE, type = 7)),
    Threshold_P99 = as.numeric(quantile(probability, 0.99, names = FALSE, type = 7)),
    Threshold_Mean_Plus_2SD = min(1, mean_probability + 2 * sd_probability),
    stringsAsFactors = FALSE
  )
}
grid_n <- nrow(full_grid)
run_n <- nrow(run_results)
probability_sum <- numeric(grid_n)
probability_sum_squares <- numeric(grid_n)
probability_min <- rep(Inf, grid_n)
probability_max <- rep(-Inf, grid_n)
prediction_count <- integer(grid_n)
threshold_rows <- vector("list", run_n)
capture_rows <- vector("list", run_n)
error_rows <- list()
cat("\nXGBoost full-grid prediction started.\n")
cat("Runs:", run_n, "\n")
cat("Grid cells:", grid_n, "\n\n")
for (run_i in seq_len(run_n)) {
  current_run <- run_results[run_i, , drop = FALSE]
  run_id <- current_run$Run_ID[1]
  cat("[", run_i, "/", run_n, "] ", run_id, "\n", sep = "")
  current_result <- tryCatch({
    fold_data <- prepare_fold_data(run_id = run_id)
    model_matrix <- prepare_model_matrix(
      fold_data = fold_data,
      predictor_columns = predictor_columns,
      categorical_predictors = categorical_predictors,
      scale_numeric = FALSE
    )
    x_train <- as.matrix(model_matrix$x_train)
    y_train <- as.numeric(model_matrix$y_train_numeric)
    full_matrix <- build_full_grid_matrix(full_grid, colnames(x_train))
    if (anyNA(full_matrix)) stop("Full-grid model matrix contains missing values.", call. = FALSE)
    dtrain <- xgboost::xgb.DMatrix(x_train, label = y_train)
    dfull <- xgboost::xgb.DMatrix(full_matrix)
    params <- list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      eta = current_run$eta[1],
      max_depth = as.integer(current_run$max_depth[1]),
      min_child_weight = current_run$min_child_weight[1],
      subsample = current_run$subsample[1],
      colsample_bytree = current_run$colsample_bytree[1],
      gamma = current_run$gamma[1],
      scale_pos_weight = current_run$scale_pos_weight[1],
      tree_method = "hist",
      seed = as.integer(current_run$XGB_Seed[1]),
      nthread = 1
    )
    model <- xgboost::xgb.train(params = params, data = dtrain, nrounds = as.integer(current_run$Best_nrounds[1]), verbose = 0)
    probability <- as.numeric(predict(model, dfull))
    if (length(probability) != grid_n || any(!is.finite(probability)) || any(probability < 0) || any(probability > 1)) stop("Invalid full-grid probability vector.", call. = FALSE)
    threshold_row <- calculate_thresholds(probability)
    threshold_row$Algorithm <- "XGBoost"
    threshold_row$Run_ID <- run_id
    threshold_row$Repeat <- current_run$Repeat[1]
    threshold_row$Fold_ID <- current_run$Fold_ID[1]
    threshold_row$Fold_Number <- current_run$Fold_Number[1]
    threshold_row <- threshold_row[, c("Algorithm", "Run_ID", "Repeat", "Fold_ID", "Fold_Number", "Grid_Mean", "Grid_SD", "Threshold_P90", "Threshold_P95", "Threshold_P99", "Threshold_Mean_Plus_2SD")]
    held_out_probability <- current_run$Held_Out_Occurrence_Probability[1]
    capture_row <- data.frame(
      Algorithm = "XGBoost",
      Run_ID = run_id,
      Repeat = current_run$Repeat[1],
      Fold_ID = current_run$Fold_ID[1],
      Fold_Number = current_run$Fold_Number[1],
      Held_Out_CellID = current_run$Held_Out_CellID[1],
      Held_Out_Probability = held_out_probability,
      Threshold_P90 = threshold_row$Threshold_P90,
      Captured_P90 = held_out_probability >= threshold_row$Threshold_P90,
      Threshold_P95 = threshold_row$Threshold_P95,
      Captured_P95 = held_out_probability >= threshold_row$Threshold_P95,
      Threshold_P99 = threshold_row$Threshold_P99,
      Captured_P99 = held_out_probability >= threshold_row$Threshold_P99,
      Threshold_Mean_Plus_2SD = threshold_row$Threshold_Mean_Plus_2SD,
      Captured_Mean_Plus_2SD = held_out_probability >= threshold_row$Threshold_Mean_Plus_2SD,
      stringsAsFactors = FALSE
    )
    list(probability = probability, threshold = threshold_row, capture = capture_row)
  }, error = function(e) {
    list(error = data.frame(Algorithm = "XGBoost", Run_ID = run_id, Repeat = current_run$Repeat[1], Fold_ID = current_run$Fold_ID[1], Error_Message = conditionMessage(e), stringsAsFactors = FALSE))
  })
  if ("error" %in% names(current_result)) {
    error_rows[[length(error_rows) + 1L]] <- current_result$error
    cat("   FAIL: ", current_result$error$Error_Message, "\n", sep = "")
  } else {
    probability <- current_result$probability
    probability_sum <- probability_sum + probability
    probability_sum_squares <- probability_sum_squares + probability^2
    probability_min <- pmin(probability_min, probability)
    probability_max <- pmax(probability_max, probability)
    prediction_count <- prediction_count + 1L
    threshold_rows[[run_i]] <- current_result$threshold
    capture_rows[[run_i]] <- current_result$capture
    cat("   PASS\n")
  }
}
error_df <- if (length(error_rows) > 0) do.call(rbind, error_rows) else data.frame(Algorithm = character(0), Run_ID = character(0), Repeat = integer(0), Fold_ID = character(0), Error_Message = character(0), stringsAsFactors = FALSE)
write.csv(error_df, error_file, row.names = FALSE)
if (nrow(error_df) > 0) stop(paste0(nrow(error_df), " full-grid run(s) failed. See ", error_file), call. = FALSE)
if (any(prediction_count != run_n)) stop("At least one grid cell does not have 150 predictions.", call. = FALSE)
mean_probability <- probability_sum / prediction_count
variance_probability <- (probability_sum_squares - probability_sum^2 / prediction_count) / (prediction_count - 1)
variance_probability <- pmax(variance_probability, 0)
sd_probability <- sqrt(variance_probability)
full_grid_summary <- data.frame(
  CellID = full_grid$CellID,
  X = full_grid$X,
  Y = full_grid$Y,
  XGB_Mean_Probability = mean_probability,
  XGB_SD_Probability = sd_probability,
  XGB_Min_Probability = probability_min,
  XGB_Max_Probability = probability_max,
  XGB_Prediction_n = prediction_count,
  stringsAsFactors = FALSE
)
run_thresholds <- do.call(rbind, threshold_rows)
held_out_capture <- do.call(rbind, capture_rows)
write.csv(full_grid_summary, summary_file, row.names = FALSE)
write.csv(run_thresholds, threshold_file, row.names = FALSE)
write.csv(held_out_capture, capture_file, row.names = FALSE)
threshold_names <- c("P90", "P95", "P99", "Mean_Plus_2SD")
capture_columns <- c(P90 = "Captured_P90", P95 = "Captured_P95", P99 = "Captured_P99", Mean_Plus_2SD = "Captured_Mean_Plus_2SD")
threshold_columns <- c(P90 = "Threshold_P90", P95 = "Threshold_P95", P99 = "Threshold_P99", Mean_Plus_2SD = "Threshold_Mean_Plus_2SD")
summary_rows <- vector("list", length(threshold_names))
for (i in seq_along(threshold_names)) {
  threshold_name <- threshold_names[i]
  capture_values <- held_out_capture[[capture_columns[[threshold_name]]]]
  threshold_values <- held_out_capture[[threshold_columns[[threshold_name]]]]
  summary_rows[[i]] <- data.frame(
    Algorithm = "XGBoost",
    Threshold_Method = threshold_name,
    Run_n = nrow(held_out_capture),
    Captured_Run_n = sum(capture_values),
    Capture_Rate = mean(capture_values),
    Mean_Threshold = mean(threshold_values),
    SD_Threshold = sd(threshold_values),
    Minimum_Threshold = min(threshold_values),
    Maximum_Threshold = max(threshold_values),
    stringsAsFactors = FALSE
  )
}
threshold_summary <- do.call(rbind, summary_rows)
write.csv(threshold_summary, capture_summary_file, row.names = FALSE)
overall_qaqc <- if (nrow(full_grid_summary) == grid_n && all(full_grid_summary$XGB_Prediction_n == 150L) && nrow(run_thresholds) == 150L && nrow(held_out_capture) == 150L && nrow(error_df) == 0L) "PASS" else "FAIL"
writeLines(c(
  "XGBOOST FULL-GRID QAQC",
  "",
  paste0("Expected runs: 150"),
  paste0("Successful runs: ", nrow(run_thresholds)),
  paste0("Failed runs: ", nrow(error_df)),
  paste0("Grid cells: ", nrow(full_grid_summary)),
  paste0("Cells with 150 predictions: ", sum(full_grid_summary$XGB_Prediction_n == 150L)),
  paste0("Held-out capture records: ", nrow(held_out_capture)),
  paste0("Overall QAQC status: ", overall_qaqc)
), "Project/QAQC/XGB_FullGrid_QAQC.txt")
writeLines(c(
  "XGBOOST FULL-GRID METHODOLOGY",
  "",
  "Each accepted outer-LOOCV XGBoost model was reconstructed using its fold-specific training data, selected hyperparameters and stored random seed.",
  "Each model predicted all valid 25 m grid cells.",
  "Cell-level ensemble mean, standard deviation, minimum and maximum probabilities were calculated.",
  "P90, P95, P99 and Mean + 2 SD thresholds were calculated separately from each run's full-grid prediction distribution.",
  "Threshold capture was evaluated only with the corresponding genuinely held-out occurrence probability."
), "Project/QAQC/XGB_FullGrid_Methodology.txt")
cat("\n========================================\n")
cat("XGBoost full-grid prediction completed.\n")
cat("Successful runs:", nrow(run_thresholds), "\n")
cat("Grid cells:", nrow(full_grid_summary), "\n")
cat("Overall QAQC:", overall_qaqc, "\n")
cat("Full-grid summary:", summary_file, "\n")
cat("Run thresholds:", threshold_file, "\n")
cat("Held-out capture:", capture_file, "\n")
cat("Threshold summary:", capture_summary_file, "\n")
cat("========================================\n")
if (overall_qaqc != "PASS") stop("XGBoost full-grid QAQC failed.", call. = FALSE)
