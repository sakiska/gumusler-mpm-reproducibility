# ============================================================
# Train_XGBoost_LOOCV.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Train and tune XGBoost models within the repeated outer validation framework.
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
output_directory <- "Project/Models/XGB_Results"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create("Project/QAQC", recursive = TRUE, showWarnings = FALSE)
run_results_file <- file.path(output_directory, "XGB_LOOCV_Run_Results.csv")
tuning_results_file <- file.path(output_directory, "XGB_LOOCV_Tuning_Results.csv")
importance_results_file <- file.path(output_directory, "XGB_LOOCV_Variable_Importance.csv")
error_log_file <- file.path(output_directory, "XGB_LOOCV_Error_Log.csv")
qaqc_file <- "Project/QAQC/XGB_LOOCV_QAQC.txt"
methodology_file <- "Project/QAQC/XGB_LOOCV_Methodology.txt"
repeat_values <- seq_len(30)
fold_values <- seq_len(5)
expected_run_n <- length(repeat_values) * length(fold_values)
base_seed <- 7000L
inner_nfold <- 4L
early_stopping_rounds <- 25L
max_nrounds <- 500L
parameter_grid <- expand.grid(
  eta = c(0.03, 0.10),
  max_depth = c(2L, 4L),
  min_child_weight = c(1, 3),
  subsample = 0.80,
  colsample_bytree = 0.80,
  gamma = 0,
  stringsAsFactors = FALSE
)
extract_positive_probability <- function(prediction_object) {
  as.numeric(prediction_object)
}
safe_numeric <- function(x) {
  as.numeric(as.character(x))
}
run_results <- vector("list", expected_run_n)
tuning_results <- list()
importance_results <- list()
error_results <- list()
run_counter <- 0L
cat("\nXGBoost LOOCV training started.\n")
cat("Expected runs:", expected_run_n, "\n\n")
for (repeat_i in repeat_values) {
  for (fold_i in fold_values) {
    run_counter <- run_counter + 1L
    run_id <- sprintf("Repeat_%02d_Fold_%02d", repeat_i, fold_i)
    cat("[", run_counter, "/", expected_run_n, "] ", run_id, "\n", sep = "")
    result <- tryCatch({
      fold_data <- prepare_fold_data(run_id = run_id)
      model_matrix <- prepare_model_matrix(
        fold_data = fold_data,
        predictor_columns = predictor_columns,
        categorical_predictors = categorical_predictors,
        scale_numeric = FALSE
      )
      x_train <- as.matrix(model_matrix$x_train)
      x_test <- as.matrix(model_matrix$x_test)
      y_train <- safe_numeric(model_matrix$y_train_numeric)
      if (length(y_train) != nrow(x_train)) stop("Training response length does not match x_train.", call. = FALSE)
      if (nrow(x_test) != 1L) stop("Each LOOCV run must contain exactly one held-out occurrence.", call. = FALSE)
      if (anyNA(x_train) || anyNA(x_test) || anyNA(y_train)) stop("Model matrix contains missing values.", call. = FALSE)
      if (!all(y_train %in% c(0, 1))) stop("Training response must be coded as 0/1.", call. = FALSE)
      positive_n <- sum(y_train == 1)
      negative_n <- sum(y_train == 0)
      if (positive_n < 2L || negative_n < 2L) stop("Insufficient class counts for XGBoost training.", call. = FALSE)
      scale_pos_weight <- negative_n / positive_n
      dtrain <- xgboost::xgb.DMatrix(data = x_train, label = y_train)
      dtest <- xgboost::xgb.DMatrix(data = x_test)
      best_score <- Inf
      best_row <- NULL
      best_nrounds <- NA_integer_
      for (grid_i in seq_len(nrow(parameter_grid))) {
        grid_row <- parameter_grid[grid_i, , drop = FALSE]
        current_seed <- base_seed + repeat_i * 100L + fold_i * 10L + grid_i
        params <- list(
          objective = "binary:logistic",
          eval_metric = "logloss",
          eta = grid_row$eta,
          max_depth = as.integer(grid_row$max_depth),
          min_child_weight = grid_row$min_child_weight,
          subsample = grid_row$subsample,
          colsample_bytree = grid_row$colsample_bytree,
          gamma = grid_row$gamma,
          scale_pos_weight = scale_pos_weight,
          tree_method = "hist",
          seed = current_seed,
          nthread = 1
        )
        cv_fit <- xgboost::xgb.cv(
          params = params,
          data = dtrain,
          nrounds = max_nrounds,
          nfold = min(inner_nfold, positive_n),
          stratified = TRUE,
          early_stopping_rounds = early_stopping_rounds,
          maximize = FALSE,
          verbose = 0,
          prediction = FALSE
        )
	evaluation_log <- cv_fit$evaluation_log
	best_iteration <- which.min(
   		evaluation_log$test_logloss_mean
	)
	best_test_logloss <-
    		evaluation_log$test_logloss_mean[
        		best_iteration
   		 ]
        tuning_results[[length(tuning_results) + 1L]] <- data.frame(
          Algorithm = "XGBoost",
          Run_ID = run_id,
          Repeat = repeat_i,
          Fold_ID = sprintf("Fold_%02d", fold_i),
          Grid_ID = grid_i,
          eta = grid_row$eta,
          max_depth = grid_row$max_depth,
          min_child_weight = grid_row$min_child_weight,
          subsample = grid_row$subsample,
          colsample_bytree = grid_row$colsample_bytree,
          gamma = grid_row$gamma,
          scale_pos_weight = scale_pos_weight,
          Best_Iteration = best_iteration,
          Inner_CV_Logloss = best_test_logloss,
          stringsAsFactors = FALSE
        )
        if (is.finite(best_test_logloss) && best_test_logloss < best_score) {
          best_score <- best_test_logloss
          best_row <- grid_row
          best_nrounds <- best_iteration
        }
      }
      if (is.null(best_row) || !is.finite(best_score) || is.na(best_nrounds)) stop("No valid XGBoost tuning result was obtained.", call. = FALSE)
      final_seed <- base_seed + repeat_i * 100L + fold_i
      final_params <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = best_row$eta,
        max_depth = as.integer(best_row$max_depth),
        min_child_weight = best_row$min_child_weight,
        subsample = best_row$subsample,
        colsample_bytree = best_row$colsample_bytree,
        gamma = best_row$gamma,
        scale_pos_weight = scale_pos_weight,
        tree_method = "hist",
        seed = final_seed,
        nthread = 1
      )
      final_model <- xgboost::xgb.train(
        params = final_params,
        data = dtrain,
        nrounds = best_nrounds,
        verbose = 0
      )
      held_out_probability <- extract_positive_probability(predict(final_model, dtest))
      if (length(held_out_probability) != 1L || !is.finite(held_out_probability) || held_out_probability < 0 || held_out_probability > 1) stop("Invalid held-out occurrence probability.", call. = FALSE)
      importance <- xgboost::xgb.importance(feature_names = colnames(x_train), model = final_model)
      if (nrow(importance) > 0) {
        importance$Algorithm <- "XGBoost"
        importance$Run_ID <- run_id
        importance$Repeat <- repeat_i
        importance$Fold_ID <- sprintf("Fold_%02d", fold_i)
        importance_results[[length(importance_results) + 1L]] <- importance[, c("Algorithm", "Run_ID", "Repeat", "Fold_ID", "Feature", "Gain", "Cover", "Frequency")]
      }
      held_out_cell_id <- if ("Held_Out_CellID" %in% names(fold_data$run_metadata)) fold_data$run_metadata$Held_Out_CellID[1] else fold_data$test_data$CellID[1]
      run_results[[run_counter]] <- data.frame(
        Algorithm = "XGBoost",
        Run_ID = run_id,
        Repeat = repeat_i,
        Fold_ID = sprintf("Fold_%02d", fold_i),
        Fold_Number = fold_i,
        Held_Out_CellID = held_out_cell_id,
        Held_Out_Occurrence_Probability = held_out_probability,
        eta = best_row$eta,
        max_depth = best_row$max_depth,
        min_child_weight = best_row$min_child_weight,
        subsample = best_row$subsample,
        colsample_bytree = best_row$colsample_bytree,
        gamma = best_row$gamma,
        scale_pos_weight = scale_pos_weight,
        Best_nrounds = best_nrounds,
        Inner_CV_Logloss = best_score,
        XGB_Seed = final_seed,
        QAQC_Status = "PASS",
        stringsAsFactors = FALSE
      )
      cat("   PASS\n")
      TRUE
    }, error = function(e) {
      error_results[[length(error_results) + 1L]] <<- data.frame(
        Algorithm = "XGBoost",
        Run_ID = run_id,
        Repeat = repeat_i,
        Fold_ID = sprintf("Fold_%02d", fold_i),
        Error_Message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      cat("   FAIL: ", conditionMessage(e), "\n", sep = "")
      FALSE
    })
  }
}
successful_results <- run_results[!vapply(run_results, is.null, logical(1))]
run_results_df <- if (length(successful_results) > 0) do.call(rbind, successful_results) else data.frame()
tuning_results_df <- if (length(tuning_results) > 0) do.call(rbind, tuning_results) else data.frame()
importance_results_df <- if (length(importance_results) > 0) do.call(rbind, importance_results) else data.frame()
error_results_df <- if (length(error_results) > 0) do.call(rbind, error_results) else data.frame(Algorithm = character(0), Run_ID = character(0), Repeat = integer(0), Fold_ID = character(0), Error_Message = character(0), stringsAsFactors = FALSE)
write.csv(run_results_df, run_results_file, row.names = FALSE, na = "")
write.csv(tuning_results_df, tuning_results_file, row.names = FALSE, na = "")
write.csv(importance_results_df, importance_results_file, row.names = FALSE, na = "")
write.csv(error_results_df, error_log_file, row.names = FALSE, na = "")
overall_qaqc <- if (nrow(run_results_df) == expected_run_n && nrow(error_results_df) == 0 && all(run_results_df$QAQC_Status == "PASS")) "PASS" else "FAIL"
writeLines(c(
  "XGBOOST LOOCV QAQC",
  "",
  paste0("Expected runs: ", expected_run_n),
  paste0("Successful runs: ", nrow(run_results_df)),
  paste0("Failed runs: ", nrow(error_results_df)),
  paste0("Overall QAQC status: ", overall_qaqc)
), qaqc_file)
writeLines(c(
  "XGBOOST LOOCV METHODOLOGY",
  "",
  "The same repeated LOOCV folds, fold-specific pseudo-background samples and predictor matrix used for RF were retained.",
  "Hyperparameters were selected independently within each outer run using stratified inner cross-validation on that run's training data only.",
  "The final model was refitted on the complete outer-training set and evaluated only at the genuinely held-out occurrence.",
  "Pseudo-background locations were not treated as verified absences in final performance reporting."
), methodology_file)
cat("\n========================================\n")
cat("XGBoost LOOCV training completed.\n")
cat("Successful runs:", nrow(run_results_df), "\n")
cat("Failed runs:", nrow(error_results_df), "\n")
cat("Overall QAQC:", overall_qaqc, "\n")
cat("Run results:", run_results_file, "\n")
cat("========================================\n")
if (overall_qaqc != "PASS") stop("XGBoost LOOCV QAQC failed. Review the error log.", call. = FALSE)
