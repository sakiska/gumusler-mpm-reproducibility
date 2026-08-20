# ============================================================
# Train_RF_LOOCV.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Train and tune Random Forest models within the repeated outer validation framework.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c(
  "ranger"
)
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
      paste(
        missing_packages,
        collapse = ", "
      ),
      ". Install with: install.packages(c(",
      paste0(
        '"',
        missing_packages,
        '"',
        collapse = ", "
      ),
      "))"
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 2. Source shared modules
# ------------------------------------------------------------
source(
  "Project/Models/Prepare_Fold_Data.R"
)
source(
  "Project/Models/Prepare_Model_Matrix.R"
)
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
categorical_predictors <- c(
  "Lithology"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
num_trees <- 1000
min_node_size_candidates <- c(
  1,
  3,
  5
)
rf_base_seed <- 5000
importance_method <- "permutation"
save_models <- FALSE
# ------------------------------------------------------------
# ------------------------------------------------------------
run_plan_file <- paste0(
  "Project/Tables/",
  "LOOCV_Fold_Repeat_Plan.csv"
)
results_directory <- paste0(
  "Project/Models/",
  "RF_Results"
)
model_directory <- paste0(
  results_directory,
  "/Saved_Models"
)
dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
if (save_models) {
  dir.create(
    model_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}
run_results_file <- paste0(
  results_directory,
  "/RF_LOOCV_Run_Results.csv"
)
tuning_results_file <- paste0(
  results_directory,
  "/RF_LOOCV_Tuning_Results.csv"
)
importance_results_file <- paste0(
  results_directory,
  "/RF_LOOCV_Variable_Importance.csv"
)
matrix_audit_file <- paste0(
  results_directory,
  "/RF_LOOCV_Model_Matrix_Audit.csv"
)
error_log_file <- paste0(
  results_directory,
  "/RF_LOOCV_Error_Log.csv"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "RF_LOOCV_Methodology.txt"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "RF_LOOCV_QAQC.txt"
)
# ------------------------------------------------------------
# 6. Read run plan
# ------------------------------------------------------------
if (!file.exists(run_plan_file)) {
  stop(
    paste0(
      "Run plan file not found: ",
      run_plan_file
    ),
    call. = FALSE
  )
}
run_plan <- read.csv(
  run_plan_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_run_columns <- c(
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "Background_Seed"
)
missing_run_columns <- setdiff(
  required_run_columns,
  names(run_plan)
)
if (length(missing_run_columns) > 0) {
  stop(
    paste0(
      "Missing columns in run plan: ",
      paste(
        missing_run_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}
if (anyDuplicated(run_plan$Run_ID) > 0) {
  stop(
    "Run plan contains duplicated Run_ID values.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
all_run_results <- list()
all_tuning_results <- list()
all_importance_results <- list()
all_matrix_audits <- list()
all_error_results <- list()
# ------------------------------------------------------------
# ------------------------------------------------------------
total_runs <- nrow(run_plan)
cat(
  "\nRandom Forest LOOCV started.\n"
)
cat(
  "Total runs:",
  total_runs,
  "\n\n"
)
for (run_i in seq_len(total_runs)) {
  current_run_id <- run_plan$Run_ID[run_i]
  cat(
    "[",
    run_i,
    "/",
    total_runs,
    "] ",
    current_run_id,
    "\n",
    sep = ""
  )
  current_result <- tryCatch(
    {
      # ------------------------------------------------------
      # ------------------------------------------------------
      fold_data <- prepare_fold_data(
        run_id = current_run_id
      )
      # ------------------------------------------------------
      #
      # ------------------------------------------------------
      model_matrix <- prepare_model_matrix(
        fold_data = fold_data,
        predictor_columns = predictor_columns,
        categorical_predictors =
          categorical_predictors,
        scale_numeric = FALSE
      )
      # ------------------------------------------------------
      # ------------------------------------------------------
      x_train <- as.data.frame(
        model_matrix$x_train,
        check.names = FALSE
      )
      x_test <- as.data.frame(
        model_matrix$x_test,
        check.names = FALSE
      )
      training_rf_data <- x_train
      training_rf_data$Class <-
        model_matrix$y_train_factor
      # ------------------------------------------------------
      # ------------------------------------------------------
      predictor_n <- ncol(x_train)
      mtry_candidates <- unique(
        pmax(
          1,
          pmin(
            predictor_n,
            c(
              1,
              floor(
                sqrt(predictor_n)
              ),
              ceiling(
                predictor_n / 3
              ),
              ceiling(
                predictor_n / 2
              ),
              predictor_n
            )
          )
        )
      )
      mtry_candidates <- sort(
        mtry_candidates
      )
      # ------------------------------------------------------
      #
      #   100 pseudo-background
      #   4 occurrence
      #
      # ------------------------------------------------------
      class_counts <- table(
        model_matrix$y_train_factor
      )
      occurrence_weight <- as.numeric(
        class_counts["PseudoBackground"] /
          class_counts["Occurrence"]
      )
      rf_class_weights <- c(
        PseudoBackground = 1,
        Occurrence = occurrence_weight
      )
      # ------------------------------------------------------
      # ------------------------------------------------------
      tuning_grid <- expand.grid(
        mtry = mtry_candidates,
        min.node.size =
          min_node_size_candidates,
        stringsAsFactors = FALSE
      )
      tuning_grid$OOB_Prediction_Error <- NA_real_
      tuning_grid$RF_Seed <- NA_integer_
      for (
        tuning_i in
        seq_len(nrow(tuning_grid))
      ) {
        current_mtry <-
          tuning_grid$mtry[tuning_i]
        current_min_node_size <-
          tuning_grid$min.node.size[tuning_i]
        current_rf_seed <-
          rf_base_seed +
          run_i * 100 +
          tuning_i
        tuning_model <- ranger::ranger(
          formula = Class ~ .,
          data = training_rf_data,
          num.trees = num_trees,
          mtry = current_mtry,
          min.node.size =
            current_min_node_size,
          probability = TRUE,
          class.weights =
            rf_class_weights,
          importance = "none",
          seed = current_rf_seed,
          num.threads = 1,
          write.forest = FALSE,
          verbose = FALSE
        )
        tuning_grid$OOB_Prediction_Error[
          tuning_i
        ] <- tuning_model$prediction.error
        tuning_grid$RF_Seed[
          tuning_i
        ] <- current_rf_seed
      }
      # ------------------------------------------------------
      #
      # ------------------------------------------------------
      tuning_grid <- tuning_grid[
        order(
          tuning_grid$OOB_Prediction_Error,
          tuning_grid$mtry,
          tuning_grid$min.node.size
        ),
        ,
        drop = FALSE
      ]
      rownames(tuning_grid) <- NULL
      best_mtry <-
        tuning_grid$mtry[1]
      best_min_node_size <-
        tuning_grid$min.node.size[1]
      best_oob_error <-
        tuning_grid$OOB_Prediction_Error[1]
      # ------------------------------------------------------
      # ------------------------------------------------------
      final_rf_seed <-
        rf_base_seed +
        run_i
      final_rf_model <- ranger::ranger(
        formula = Class ~ .,
        data = training_rf_data,
        num.trees = num_trees,
        mtry = best_mtry,
        min.node.size =
          best_min_node_size,
        probability = TRUE,
        class.weights =
          rf_class_weights,
        importance =
          importance_method,
        seed = final_rf_seed,
        num.threads = 1,
        write.forest = TRUE,
        verbose = FALSE
      )
      # ------------------------------------------------------
      # ------------------------------------------------------
      held_out_prediction <- predict(
        final_rf_model,
        data = x_test
      )$predictions
      if (
        !"Occurrence" %in%
          colnames(held_out_prediction)
      ) {
        stop(
          paste0(
            "Occurrence probability column not found for ",
            current_run_id,
            "."
          ),
          call. = FALSE
        )
      }
      held_out_probability <-
        as.numeric(
          held_out_prediction[
            1,
            "Occurrence"
          ]
        )
      # ------------------------------------------------------
      # ------------------------------------------------------
      model_file <- NA_character_
      if (save_models) {
        model_file <- paste0(
          model_directory,
          "/RF_",
          current_run_id,
          ".rds"
        )
        saveRDS(
          final_rf_model,
          model_file
        )
      }
      # ------------------------------------------------------
      # 8.11. Run-level result
      # ------------------------------------------------------
      run_result <- data.frame(
        Algorithm = "RF",
        Run_ID = current_run_id,
        Repeat =
          run_plan$Repeat[run_i],
        Fold_ID =
          run_plan$Fold_ID[run_i],
        Fold_Number =
          run_plan$Fold_Number[run_i],
        Background_Seed =
          run_plan$Background_Seed[run_i],
        Held_Out_CellID =
          fold_data$test_data$CellID[1],
        Held_Out_X =
          fold_data$test_data$X[1],
        Held_Out_Y =
          fold_data$test_data$Y[1],
        Held_Out_Observed_Class =
          fold_data$test_data$Class[1],
        Held_Out_Occurrence_Probability =
          held_out_probability,
        Predictor_n =
          predictor_n,
        Training_Row_n =
          nrow(x_train),
        Training_Occurrence_n =
          sum(
            model_matrix$y_train_numeric == 1
          ),
        Training_Background_n =
          sum(
            model_matrix$y_train_numeric == 0
          ),
        Best_mtry =
          best_mtry,
        Best_min_node_size =
          best_min_node_size,
        Num_Trees =
          num_trees,
        Occurrence_Class_Weight =
          occurrence_weight,
        Best_Tuning_OOB_Error =
          best_oob_error,
        Final_Model_OOB_Error =
          final_rf_model$prediction.error,
        RF_Seed =
          final_rf_seed,
        Model_File =
          model_file,
        QAQC_Status = "PASS",
        stringsAsFactors = FALSE
      )
      # ------------------------------------------------------
      # ------------------------------------------------------
      tuning_result <- tuning_grid
      tuning_result$Algorithm <- "RF"
      tuning_result$Run_ID <-
        current_run_id
      tuning_result$Repeat <-
        run_plan$Repeat[run_i]
      tuning_result$Fold_ID <-
        run_plan$Fold_ID[run_i]
      tuning_result$Selected <-
        seq_len(
          nrow(tuning_result)
        ) == 1
      tuning_result <- tuning_result[
        ,
        c(
          "Algorithm",
          "Run_ID",
          "Repeat",
          "Fold_ID",
          "mtry",
          "min.node.size",
          "OOB_Prediction_Error",
          "RF_Seed",
          "Selected"
        )
      ]
      # ------------------------------------------------------
      # ------------------------------------------------------
      importance_values <-
        final_rf_model$variable.importance
      importance_result <- data.frame(
        Algorithm = "RF",
        Run_ID = current_run_id,
        Repeat =
          run_plan$Repeat[run_i],
        Fold_ID =
          run_plan$Fold_ID[run_i],
        Predictor =
          names(importance_values),
        Importance =
          as.numeric(importance_values),
        stringsAsFactors = FALSE
      )
      # ------------------------------------------------------
      # ------------------------------------------------------
      matrix_audit <-
        model_matrix$audit
      matrix_audit$Algorithm <- "RF"
      matrix_audit <- matrix_audit[
        ,
        c(
          "Algorithm",
          setdiff(
            names(matrix_audit),
            "Algorithm"
          )
        ),
        drop = FALSE
      ]
      # ------------------------------------------------------
      # 8.15. Return current run
      # ------------------------------------------------------
      list(
        run_result =
          run_result,
        tuning_result =
          tuning_result,
        importance_result =
          importance_result,
        matrix_audit =
          matrix_audit
      )
    },
    error = function(e) {
      error_result <- data.frame(
        Algorithm = "RF",
        Run_ID = current_run_id,
        Repeat =
          run_plan$Repeat[run_i],
        Fold_ID =
          run_plan$Fold_ID[run_i],
        Fold_Number =
          run_plan$Fold_Number[run_i],
        Error_Message =
          conditionMessage(e),
        QAQC_Status = "FAIL",
        stringsAsFactors = FALSE
      )
      list(
        error_result =
          error_result
      )
    }
  )
  # ----------------------------------------------------------
  # 8.16. Store current run output
  # ----------------------------------------------------------
  if (
    "error_result" %in%
      names(current_result)
  ) {
    all_error_results[[length(all_error_results) + 1]] <- current_result$error_result
    cat(
      "   FAIL: ",
      current_result$error_result$Error_Message,
      "\n",
      sep = ""
    )
  } else {
    all_run_results[[length(all_run_results) + 1]] <- current_result$run_result
    all_tuning_results[[length(all_tuning_results) + 1]] <- current_result$tuning_result
    all_importance_results[[length(all_importance_results) + 1]] <- current_result$importance_result
    all_matrix_audits[[length(all_matrix_audits) + 1]] <- current_result$matrix_audit
    cat(
      "   PASS | Held-out probability: ",
      round(
        current_result$run_result$
          Held_Out_Occurrence_Probability,
        4
      ),
      "\n",
      sep = ""
    )
  }
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (length(all_run_results) == 0) {
  stop(
    paste(
      "No RF run completed successfully.",
      "Check the error log."
    ),
    call. = FALSE
  )
}
run_results <- do.call(
  rbind,
  all_run_results
)
tuning_results <- do.call(
  rbind,
  all_tuning_results
)
importance_results <- do.call(
  rbind,
  all_importance_results
)
matrix_audits <- do.call(
  rbind,
  all_matrix_audits
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (length(all_error_results) > 0) {
  error_results <- do.call(
    rbind,
    all_error_results
  )
} else {
  error_results <- data.frame(
    Algorithm = character(0),
    Run_ID = character(0),
    Repeat = integer(0),
    Fold_ID = character(0),
    Fold_Number = integer(0),
    Error_Message = character(0),
    QAQC_Status = character(0),
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
successful_run_n <- nrow(
  run_results
)
failed_run_n <- nrow(
  error_results
)
expected_run_n <- nrow(
  run_plan
)
duplicate_successful_run_n <- sum(
  duplicated(
    run_results$Run_ID
  )
)
missing_run_ids <- setdiff(
  run_plan$Run_ID,
  run_results$Run_ID
)
invalid_probability_n <- sum(
  is.na(
    run_results$
      Held_Out_Occurrence_Probability
  ) |
    run_results$
      Held_Out_Occurrence_Probability < 0 |
    run_results$
      Held_Out_Occurrence_Probability > 1
)
matrix_qaqc_fail_n <- sum(
  matrix_audits$QAQC_Status != "PASS"
)
all_runs_passed <- all(
  successful_run_n == expected_run_n,
  failed_run_n == 0,
  duplicate_successful_run_n == 0,
  length(missing_run_ids) == 0,
  invalid_probability_n == 0,
  matrix_qaqc_fail_n == 0
)
overall_qaqc_status <- ifelse(
  all_runs_passed,
  "PASS",
  "FAIL"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  run_results,
  run_results_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  tuning_results,
  tuning_results_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  importance_results,
  importance_results_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  matrix_audits,
  matrix_audit_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  error_results,
  error_log_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "RANDOM FOREST LOOCV METHODOLOGY",
  "",
  paste0(
    "Expected runs: ",
    expected_run_n
  ),
  paste0(
    "Number of trees: ",
    num_trees
  ),
  paste0(
    "Minimum node size candidates: ",
    paste(
      min_node_size_candidates,
      collapse = ", "
    )
  ),
  paste0(
    "Predictors: ",
    paste(
      predictor_columns,
      collapse = ", "
    )
  ),
  paste0(
    "Categorical predictors: ",
    paste(
      categorical_predictors,
      collapse = ", "
    )
  ),
  "",
  paste(
    "Each run used four training occurrences,",
    "one held-out occurrence and 100 fold-specific",
    "pseudo-background cells."
  ),
  paste(
    "Pseudo-background cells were generated using only",
    "the training occurrence buffer."
  ),
  paste(
    "Known occurrence cells were excluded from",
    "pseudo-background selection."
  ),
  paste(
    "All model-matrix transformations were learned",
    "from the training data."
  ),
  paste(
    "Random Forest hyperparameters were selected using",
    "training-set out-of-bag prediction error."
  ),
  paste(
    "The held-out occurrence was not used in",
    "background generation, preprocessing, tuning or",
    "model fitting."
  )
)
writeLines(
  methodology_lines,
  methodology_file,
  useBytes = TRUE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc_lines <- c(
  "RANDOM FOREST LOOCV QAQC",
  "",
  paste0(
    "Expected runs: ",
    expected_run_n
  ),
  paste0(
    "Successful runs: ",
    successful_run_n
  ),
  paste0(
    "Failed runs: ",
    failed_run_n
  ),
  paste0(
    "Duplicated successful Run_ID values: ",
    duplicate_successful_run_n
  ),
  paste0(
    "Missing Run_ID values: ",
    length(missing_run_ids)
  ),
  paste0(
    "Invalid held-out probabilities: ",
    invalid_probability_n
  ),
  paste0(
    "Failed model-matrix audits: ",
    matrix_qaqc_fail_n
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
cat(
  "\n========================================\n"
)
cat(
  "Random Forest LOOCV completed.\n"
)
cat(
  "Expected runs:",
  expected_run_n,
  "\n"
)
cat(
  "Successful runs:",
  successful_run_n,
  "\n"
)
cat(
  "Failed runs:",
  failed_run_n,
  "\n"
)
cat(
  "Overall QAQC:",
  overall_qaqc_status,
  "\n"
)
cat(
  "Run results:",
  run_results_file,
  "\n"
)
cat(
  "Tuning results:",
  tuning_results_file,
  "\n"
)
cat(
  "Variable importance:",
  importance_results_file,
  "\n"
)
cat(
  "Error log:",
  error_log_file,
  "\n"
)
cat(
  "========================================\n"
)
