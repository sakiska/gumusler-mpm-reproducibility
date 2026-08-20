# ============================================================
# Predict_XGB_Independent_Validation.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Generate XGBoost predictions for the independent held-out occurrence and pseudo-background validation sets.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c(
  "xgboost"
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
      )
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 2. SHARED MODULES
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
run_results_file <- paste0(
  "Project/Models/XGB_Results/",
  "XGB_LOOCV_Run_Results.csv"
)
validation_file <- paste0(
  "Project/Tables/",
  "Independent_Validation_PseudoBackground.csv"
)
# ------------------------------------------------------------
# 5. OUTPUT PATHS
# ------------------------------------------------------------
output_directory <- paste0(
  "Project/Models/XGB_Results/",
  "Independent_Validation"
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
prediction_file <- file.path(
  output_directory,
  "XGB_Independent_Validation_Predictions.csv"
)
run_summary_file <- file.path(
  output_directory,
  "XGB_Independent_Validation_Run_Summary.csv"
)
error_file <- file.path(
  output_directory,
  "XGB_Independent_Validation_Error_Log.csv"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "XGB_Independent_Validation_QAQC.txt"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "XGB_Independent_Validation_Methodology.txt"
)
# ------------------------------------------------------------
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
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
}
safe_numeric <- function(x) {
  as.numeric(
    as.character(x)
  )
}
# ------------------------------------------------------------
#
# ------------------------------------------------------------
build_validation_matrix <- function(
    validation_data,
    model_matrix
) {
  predictor_data <- validation_data[
    ,
    predictor_columns,
    drop = FALSE
  ]
  numeric_predictors <-
    model_matrix$numeric_predictors
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  for (
    current_predictor in numeric_predictors
  ) {
    values <- as.numeric(
      predictor_data[[current_predictor]]
    )
    values[
      is.infinite(values)
    ] <- NA_real_
    training_median <-
      model_matrix$
        preprocessing$
        numeric_imputation_values[
          current_predictor
        ]
    if (
      length(training_median) != 1 ||
        !is.finite(training_median)
    ) {
      stop(
        paste0(
          "Invalid training-derived imputation value for ",
          current_predictor,
          "."
        ),
        call. = FALSE
      )
    }
    values[
      is.na(values)
    ] <- training_median
    predictor_data[[current_predictor]] <-
      values
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  for (
    current_predictor in
      categorical_predictors
  ) {
    values <- as.character(
      predictor_data[[current_predictor]]
    )
    values <- trimws(values)
    values[
      is.na(values) |
        values == ""
    ] <- "Missing"
    allowed_levels <-
  model_matrix$
    preprocessing$
    factor_level_map[[current_predictor]]
    unknown_level <-
      model_matrix$
        preprocessing$
        unknown_factor_level
    unseen_values <- !values %in%
      allowed_levels
    if (
      any(unseen_values)
    ) {
      values[
        unseen_values
      ] <- unknown_level
    }
    predictor_data[[current_predictor]] <-
      factor(
        values,
        levels = allowed_levels
      )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  matrix_data <- stats::model.matrix(
    ~ . - 1,
    data = predictor_data,
    na.action = stats::na.pass
  )
  matrix_data <- as.data.frame(
    matrix_data,
    check.names = FALSE
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  final_columns <-
    model_matrix$final_matrix_columns
  missing_columns <- setdiff(
    final_columns,
    names(matrix_data)
  )
  if (
    length(missing_columns) > 0
  ) {
    for (
      column_name in missing_columns
    ) {
      matrix_data[[column_name]] <- 0
    }
  }
  extra_columns <- setdiff(
    names(matrix_data),
    final_columns
  )
  if (
    length(extra_columns) > 0
  ) {
    matrix_data <- matrix_data[
      ,
      setdiff(
        names(matrix_data),
        extra_columns
      ),
      drop = FALSE
    ]
  }
  matrix_data <- matrix_data[
    ,
    final_columns,
    drop = FALSE
  ]
  matrix_data <- as.matrix(
    matrix_data
  )
  if (
    anyNA(matrix_data) ||
      any(
        !is.finite(
          matrix_data
        )
      )
  ) {
    stop(
      "Validation model matrix contains invalid values.",
      call. = FALSE
    )
  }
  matrix_data
}
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  run_results_file,
  validation_file
)
missing_files <- input_files[
  !file.exists(
    input_files
  )
]
if (
  length(missing_files) > 0
) {
  stop(
    paste0(
      "Missing input file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}
run_results <- read.csv(
  run_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
validation_data <- read.csv(
  validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
assert_required_columns(
  run_results,
  c(
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "Held_Out_CellID",
    "Held_Out_Occurrence_Probability",
    "eta",
    "max_depth",
    "min_child_weight",
    "subsample",
    "colsample_bytree",
    "gamma",
    "scale_pos_weight",
    "Best_nrounds",
    "XGB_Seed",
    "QAQC_Status"
  ),
  "XGBoost run-results table"
)
assert_required_columns(
  validation_data,
  c(
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "CellID",
    predictor_columns
  ),
  "Independent-validation table"
)
if (
  nrow(run_results) != 150L
) {
  stop(
    paste0(
      "Expected 150 XGBoost runs; found ",
      nrow(run_results),
      "."
    ),
    call. = FALSE
  )
}
if (
  any(
    run_results$QAQC_Status != "PASS"
  )
) {
  stop(
    "XGBoost run-results contain failed runs.",
    call. = FALSE
  )
}
if (
  anyDuplicated(
    run_results$Run_ID
  ) > 0
) {
  stop(
    "XGBoost run-results contain duplicated Run_ID values.",
    call. = FALSE
  )
}
validation_count_by_run <- table(
  validation_data$Run_ID
)
if (
  length(validation_count_by_run) != 150L ||
    any(
      validation_count_by_run != 100L
    )
) {
  stop(
    paste(
      "Expected exactly 100 independent validation",
      "pseudo-background cells for each of the 150 runs."
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
prediction_rows <- vector(
  "list",
  nrow(run_results)
)
summary_rows <- vector(
  "list",
  nrow(run_results)
)
error_rows <- list()
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("XGBOOST INDEPENDENT VALIDATION STARTED\n")
cat("============================================\n")
cat(
  "Runs:",
  nrow(run_results),
  "\n\n"
)
for (
  run_i in seq_len(
    nrow(run_results)
  )
) {
  current_run <- run_results[
    run_i,
    ,
    drop = FALSE
  ]
  run_id <-
    current_run$Run_ID[1]
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
  current_result <- tryCatch({
    # --------------------------------------------------------
    # --------------------------------------------------------
    fold_data <- prepare_fold_data(
      run_id = run_id
    )
    model_matrix <- prepare_model_matrix(
      fold_data = fold_data,
      predictor_columns =
        predictor_columns,
      categorical_predictors =
        categorical_predictors,
      scale_numeric = FALSE
    )
    x_train <- as.matrix(
      model_matrix$x_train
    )
    x_test <- as.matrix(
      model_matrix$x_test
    )
    y_train <- safe_numeric(
      model_matrix$y_train_numeric
    )
    if (
      nrow(x_test) != 1L
    ) {
      stop(
        "Outer test set must contain exactly one held-out occurrence.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    current_validation <- validation_data[
      validation_data$Run_ID ==
        run_id,
      ,
      drop = FALSE
    ]
    if (
      nrow(current_validation) !=
        100L
    ) {
      stop(
        paste0(
          "Expected 100 validation pseudo-background cells; found ",
          nrow(current_validation),
          "."
        ),
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    training_overlap_n <- sum(
      current_validation$CellID %in%
        model_matrix$training_cellids
    )
    heldout_overlap_n <- sum(
      current_validation$CellID ==
        current_run$Held_Out_CellID[1]
    )
    if (
      training_overlap_n != 0L
    ) {
      stop(
        "Independent validation cells overlap outer-training CellIDs.",
        call. = FALSE
      )
    }
    if (
      heldout_overlap_n != 0L
    ) {
      stop(
        "Independent validation background contains held-out occurrence CellID.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    validation_matrix <-
      build_validation_matrix(
        validation_data =
          current_validation,
        model_matrix =
          model_matrix
      )
    if (
      !identical(
        colnames(
          validation_matrix
        ),
        colnames(
          x_train
        )
      )
    ) {
      stop(
        "Validation matrix columns do not match outer-training matrix.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    dtrain <- xgboost::xgb.DMatrix(
      data = x_train,
      label = y_train
    )
    dtest <- xgboost::xgb.DMatrix(
      data = x_test
    )
    dvalidation <- xgboost::xgb.DMatrix(
      data = validation_matrix
    )
    # --------------------------------------------------------
    # 10.6 Reconstruct accepted XGBoost model
    # --------------------------------------------------------
    params <- list(
      objective =
        "binary:logistic",
      eval_metric =
        "logloss",
      eta =
        safe_numeric(
          current_run$eta[1]
        ),
      max_depth =
        as.integer(
          current_run$max_depth[1]
        ),
      min_child_weight =
        safe_numeric(
          current_run$
            min_child_weight[1]
        ),
      subsample =
        safe_numeric(
          current_run$subsample[1]
        ),
      colsample_bytree =
        safe_numeric(
          current_run$
            colsample_bytree[1]
        ),
      gamma =
        safe_numeric(
          current_run$gamma[1]
        ),
      scale_pos_weight =
        safe_numeric(
          current_run$
            scale_pos_weight[1]
        ),
      tree_method =
        "hist",
      seed =
        as.integer(
          current_run$XGB_Seed[1]
        ),
      nthread = 1
    )
    final_model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds =
        as.integer(
          current_run$
            Best_nrounds[1]
        ),
      verbose = 0
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    held_out_probability <-
      as.numeric(
        predict(
          final_model,
          dtest
        )
      )
    if (
      length(
        held_out_probability
      ) != 1L ||
        !is.finite(
          held_out_probability
        ) ||
        held_out_probability < 0 ||
        held_out_probability > 1
    ) {
      stop(
        "Invalid reconstructed held-out occurrence probability.",
        call. = FALSE
      )
    }
    stored_held_out_probability <-
      safe_numeric(
        current_run$
          Held_Out_Occurrence_Probability[1]
      )
    reproducibility_difference <- abs(
      held_out_probability -
        stored_held_out_probability
    )
    reproducibility_pass <-
      reproducibility_difference <=
        1e-8
    if (
      !reproducibility_pass
    ) {
      stop(
        paste0(
          "Held-out probability reproducibility failed. ",
          "Stored = ",
          stored_held_out_probability,
          "; reconstructed = ",
          held_out_probability,
          "; difference = ",
          reproducibility_difference
        ),
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    validation_probability <-
      as.numeric(
        predict(
          final_model,
          dvalidation
        )
      )
    if (
      length(
        validation_probability
      ) != 100L ||
        any(
          !is.finite(
            validation_probability
          )
        ) ||
        any(
          validation_probability < 0 |
            validation_probability > 1
        )
    ) {
      stop(
        "Invalid independent validation probability vector.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # 10.9 Create held-out occurrence row
    # --------------------------------------------------------
    positive_row <- data.frame(
      Algorithm =
        "XGB",
      Run_ID =
        run_id,
      Repeat =
        current_run$Repeat[1],
      Fold_ID =
        current_run$Fold_ID[1],
      Fold_Number =
        current_run$Fold_Number[1],
      CellID =
        current_run$Held_Out_CellID[1],
      Validation_Role =
        "Held_Out_Occurrence",
      Observed_Class =
        1L,
      Probability =
        held_out_probability,
      stringsAsFactors = FALSE
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    background_rows <- data.frame(
      Algorithm =
        "XGB",
      Run_ID =
        run_id,
      Repeat =
        current_run$Repeat[1],
      Fold_ID =
        current_run$Fold_ID[1],
      Fold_Number =
        current_run$Fold_Number[1],
      CellID =
        current_validation$CellID,
      Validation_Role =
        "Independent_PseudoBackground",
      Observed_Class =
        0L,
      Probability =
        validation_probability,
      stringsAsFactors = FALSE
    )
    combined_predictions <- rbind(
      positive_row,
      background_rows
    )
    # --------------------------------------------------------
    # 10.11 Run summary
    # --------------------------------------------------------
    run_summary <- data.frame(
      Algorithm =
        "XGB",
      Run_ID =
        run_id,
      Repeat =
        current_run$Repeat[1],
      Fold_ID =
        current_run$Fold_ID[1],
      Held_Out_CellID =
        current_run$Held_Out_CellID[1],
      Stored_HeldOut_Probability =
        stored_held_out_probability,
      Reconstructed_HeldOut_Probability =
        held_out_probability,
      Reproducibility_Difference =
        reproducibility_difference,
      Reproducibility_PASS =
        reproducibility_pass,
      Validation_Background_n =
        length(
          validation_probability
        ),
      Mean_Validation_Background_Probability =
        mean(
          validation_probability
        ),
      Median_Validation_Background_Probability =
        median(
          validation_probability
        ),
      Maximum_Validation_Background_Probability =
        max(
          validation_probability
        ),
      Training_Validation_Overlap_n =
        training_overlap_n,
      Validation_HeldOut_Overlap_n =
        heldout_overlap_n,
      QAQC_Status =
        "PASS",
      stringsAsFactors = FALSE
    )
    list(
      predictions =
        combined_predictions,
      summary =
        run_summary
    )
  }, error = function(e) {
    list(
      error = data.frame(
        Algorithm =
          "XGB",
        Run_ID =
          run_id,
        Repeat =
          current_run$Repeat[1],
        Fold_ID =
          current_run$Fold_ID[1],
        Error_Message =
          conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  })
  # ----------------------------------------------------------
  # 10.12 Store run result
  # ----------------------------------------------------------
  if (
    "error" %in%
      names(current_result)
  ) {
    error_rows[[length(error_rows) + 1L]] <- current_result$error <- current_result$error
    cat(
      "   FAIL: ",
      current_result$error$
        Error_Message,
      "\n",
      sep = ""
    )
  } else {
    prediction_rows[[run_i]] <-
      current_result$predictions
    summary_rows[[run_i]] <-
      current_result$summary
    cat(
      "   PASS\n"
    )
  }
}
# ------------------------------------------------------------
# 11. ERROR LOG
# ------------------------------------------------------------
error_df <- if (
  length(error_rows) > 0
) {
  do.call(
    rbind,
    error_rows
  )
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
  error_file,
  row.names = FALSE,
  na = ""
)
if (
  nrow(error_df) > 0
) {
  stop(
    paste0(
      nrow(error_df),
      " XGBoost independent-validation run(s) failed. ",
      "See: ",
      error_file
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
prediction_df <- do.call(
  rbind,
  prediction_rows
)
rownames(
  prediction_df
) <- NULL
run_summary_df <- do.call(
  rbind,
  summary_rows
)
rownames(
  run_summary_df
) <- NULL
# ------------------------------------------------------------
# 13. GLOBAL QA/QC
# ------------------------------------------------------------
global_qaqc <- data.frame(
  Check = c(
    "150 successful XGBoost runs",
    "15150 validation prediction rows",
    "150 held-out occurrence predictions",
    "15000 independent pseudo-background predictions",
    "No training-validation CellID overlap",
    "No validation-held-out occurrence overlap",
    "All reconstructed held-out probabilities reproduce originals",
    "All probabilities within 0-1",
    "All run QAQC statuses PASS"
  ),
  Result = c(
    nrow(run_summary_df) ==
      150L,
    nrow(prediction_df) ==
      15150L,
    sum(
      prediction_df$
        Observed_Class == 1L
    ) == 150L,
    sum(
      prediction_df$
        Observed_Class == 0L
    ) == 15000L,
    all(
      run_summary_df$
        Training_Validation_Overlap_n ==
        0L
    ),
    all(
      run_summary_df$
        Validation_HeldOut_Overlap_n ==
        0L
    ),
    all(
      run_summary_df$
        Reproducibility_PASS
    ),
    all(
      is.finite(
        prediction_df$Probability
      ) &
        prediction_df$Probability >= 0 &
        prediction_df$Probability <= 1
    ),
    all(
      run_summary_df$
        QAQC_Status == "PASS"
    )
  ),
  stringsAsFactors = FALSE
)
print(
  global_qaqc
)
if (
  !all(
    global_qaqc$Result
  )
) {
  stop(
    "XGBoost independent-validation global QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 14. SAVE OUTPUT TABLES
# ------------------------------------------------------------
write.csv(
  prediction_df,
  prediction_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  run_summary_df,
  run_summary_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc_lines <- c(
  "XGBOOST INDEPENDENT VALIDATION QA/QC",
  "============================================",
  "",
  paste(
    global_qaqc$Check,
    global_qaqc$Result,
    sep = ": "
  ),
  "",
  "Leakage-control conditions:",
  "",
  "- Fold-specific outer-training data were reconstructed using the original repeated LOOCV run definition.",
  "- Independent validation pseudo-background cells were not included in model fitting.",
  "- Validation data were not used for preprocessing, hyperparameter selection, early stopping or model configuration.",
  "- Only preprocessing rules learned from the outer-training data were applied to validation cells.",
  "- Training and validation CellIDs were explicitly checked for overlap.",
  "- The held-out occurrence remained outside model fitting.",
  "- Previously selected hyperparameters, Best_nrounds and XGB_Seed were reused.",
  "- Reconstructed held-out probabilities were compared against original stored LOOCV probabilities.",
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  qaqc_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "XGBOOST INDEPENDENT VALIDATION METHODOLOGY",
  "============================================",
  "",
  paste(
    "Each of the 150 accepted XGBoost outer-LOOCV models",
    "was reconstructed using its original fold-specific",
    "training dataset."
  ),
  "",
  paste(
    "The previously selected hyperparameters, number",
    "of boosting rounds and XGBoost random seed were",
    "reused. Hyperparameter tuning was not repeated."
  ),
  "",
  paste(
    "All preprocessing applied to independent",
    "validation cells was derived exclusively from",
    "the corresponding outer-training dataset."
  ),
  "",
  paste(
    "For each run, the genuinely held-out occurrence",
    "and 100 independent pseudo-background cells were",
    "predicted."
  ),
  "",
  paste(
    "The validation pseudo-background cells were not",
    "used during fitting, preprocessing, tuning, early",
    "stopping or threshold estimation."
  ),
  "",
  paste(
    "Pseudo-background cells are not verified mineral",
    "absences. Therefore subsequent confusion-matrix",
    "metrics quantify occurrence-versus-pseudo-background",
    "discrimination."
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
cat("XGBOOST INDEPENDENT VALIDATION COMPLETE\n")
cat("============================================\n")
cat(
  "Successful runs              :",
  nrow(run_summary_df),
  "\n"
)
cat(
  "Held-out occurrence rows     :",
  sum(
    prediction_df$
      Observed_Class == 1L
  ),
  "\n"
)
cat(
  "Pseudo-background rows       :",
  sum(
    prediction_df$
      Observed_Class == 0L
  ),
  "\n"
)
cat(
  "Total prediction rows        :",
  nrow(prediction_df),
  "\n"
)
cat(
  "Max reproducibility error    :",
  max(
    run_summary_df$
      Reproducibility_Difference
  ),
  "\n"
)
cat("\n")
print(
  global_qaqc
)
cat("\n")
cat(
  "XGBOOST INDEPENDENT VALIDATION: PASS\n"
)
