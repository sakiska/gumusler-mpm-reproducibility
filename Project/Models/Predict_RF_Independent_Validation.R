# ============================================================
# Predict_RF_Independent_Validation.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Generate Random Forest predictions for the independent held-out occurrence and pseudo-background validation sets.
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
  "Project/Models/RF_Results/",
  "RF_LOOCV_Run_Results.csv"
)
validation_file <- paste0(
  "Project/Tables/",
  "Independent_Validation_PseudoBackground.csv"
)
# ------------------------------------------------------------
# 5. OUTPUTS
# ------------------------------------------------------------
output_directory <- paste0(
  "Project/Models/RF_Results/",
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
  "RF_Independent_Validation_Predictions.csv"
)
run_summary_file <- file.path(
  output_directory,
  "RF_Independent_Validation_Run_Summary.csv"
)
error_file <- file.path(
  output_directory,
  "RF_Independent_Validation_Error_Log.csv"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "RF_Independent_Validation_QAQC.txt"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "RF_Independent_Validation_Methodology.txt"
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
build_aligned_prediction_matrix <- function(
    new_data,
    training_column_names
) {
  predictor_data <- new_data[
    ,
    predictor_columns,
    drop = FALSE
  ]
  numeric_predictors <- setdiff(
    predictor_columns,
    categorical_predictors
  )
  for (
    column_name in numeric_predictors
  ) {
    predictor_data[[column_name]] <-
      as.numeric(
        predictor_data[[column_name]]
      )
  }
  for (
    column_name in categorical_predictors
  ) {
    predictor_data[[column_name]] <-
      factor(
        predictor_data[[column_name]]
      )
  }
  matrix_data <- stats::model.matrix(
    ~ . - 1,
    data = predictor_data,
    na.action = stats::na.pass
  )
  matrix_data <- as.data.frame(
    matrix_data,
    check.names = FALSE
  )
  missing_columns <- setdiff(
    training_column_names,
    names(matrix_data)
  )
  if (length(missing_columns) > 0) {
    for (
      column_name in missing_columns
    ) {
      matrix_data[[column_name]] <- 0
    }
  }
  extra_columns <- setdiff(
    names(matrix_data),
    training_column_names
  )
  if (length(extra_columns) > 0) {
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
    training_column_names,
    drop = FALSE
  ]
  matrix_data
}
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  run_results_file,
  validation_file
)
missing_files <- input_files[
  !file.exists(input_files)
]
if (length(missing_files) > 0) {
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
    "Best_mtry",
    "Best_min_node_size",
    "Num_Trees",
    "Occurrence_Class_Weight",
    "RF_Seed",
    "QAQC_Status"
  ),
  "RF run-results table"
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
      "Expected 150 RF runs; found ",
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
    "RF run-results table contains failed runs.",
    call. = FALSE
  )
}
if (
  anyDuplicated(
    run_results$Run_ID
  ) > 0
) {
  stop(
    "RF run-results contain duplicated Run_ID values.",
    call. = FALSE
  )
}
validation_count_by_run <- table(
  validation_data$Run_ID
)
if (
  length(validation_count_by_run) != 150L ||
    any(validation_count_by_run != 100L)
) {
  stop(
    paste(
      "Expected exactly 100 independent validation",
      "pseudo-background cells for every one of the",
      "150 runs."
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
cat("RF INDEPENDENT VALIDATION STARTED\n")
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
  current_result <- tryCatch({
    # --------------------------------------------------------
    # --------------------------------------------------------
    fold_data <- prepare_fold_data(
      run_id = run_id
    )
    model_matrix <- prepare_model_matrix(
      fold_data = fold_data,
      predictor_columns = predictor_columns,
      categorical_predictors =
        categorical_predictors,
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
    current_validation <- validation_data[
      validation_data$Run_ID == run_id,
      ,
      drop = FALSE
    ]
    if (
      nrow(current_validation) != 100L
    ) {
      stop(
        paste0(
          "Expected 100 validation backgrounds; found ",
          nrow(current_validation),
          "."
        ),
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    if (
      any(
        current_validation$CellID %in%
          model_matrix$training_cellids
      )
    ) {
      stop(
        paste(
          "Independent validation cells overlap",
          "outer-training cells."
        ),
        call. = FALSE
      )
    }
    if (
      any(
        current_validation$CellID ==
          current_run$Held_Out_CellID[1]
      )
    ) {
      stop(
        paste(
          "Independent validation background contains",
          "the held-out occurrence CellID."
        ),
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    validation_matrix <-
      build_aligned_prediction_matrix(
        new_data = current_validation,
        training_column_names =
          names(x_train)
      )
    if (
      anyNA(validation_matrix) ||
        any(
          !is.finite(
            as.matrix(
              validation_matrix
            )
          )
        )
    ) {
      stop(
        "Validation matrix contains invalid values.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # 10.5 Reconstruct accepted RF model
    # --------------------------------------------------------
    rf_class_weights <- c(
      PseudoBackground = 1,
      Occurrence =
        as.numeric(
          current_run$
            Occurrence_Class_Weight[1]
        )
    )
    final_model <- ranger::ranger(
      formula = Class ~ .,
      data = training_rf_data,
      num.trees =
        as.integer(
          current_run$Num_Trees[1]
        ),
      mtry =
        as.integer(
          current_run$Best_mtry[1]
        ),
      min.node.size =
        as.integer(
          current_run$
            Best_min_node_size[1]
        ),
      probability = TRUE,
      class.weights =
        rf_class_weights,
      importance = "none",
      seed =
        as.integer(
          current_run$RF_Seed[1]
        ),
      num.threads = 1,
      write.forest = TRUE,
      verbose = FALSE
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    held_out_matrix <- as.data.frame(
      model_matrix$x_test,
      check.names = FALSE
    )
    held_out_prediction <- predict(
      final_model,
      data = held_out_matrix
    )$predictions
    if (
      !"Occurrence" %in%
        colnames(
          held_out_prediction
        )
    ) {
      stop(
        paste(
          "Occurrence probability column missing",
          "from held-out prediction."
        ),
        call. = FALSE
      )
    }
    held_out_probability <-
      as.numeric(
        held_out_prediction[
          ,
          "Occurrence"
        ]
      )
    stored_held_out_probability <-
      as.numeric(
        current_run$
          Held_Out_Occurrence_Probability[1]
      )
    reproducibility_difference <- abs(
      held_out_probability -
        stored_held_out_probability
    )
    reproducibility_pass <-
      reproducibility_difference <= 1e-10
    if (!reproducibility_pass) {
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
    validation_prediction <- predict(
      final_model,
      data = validation_matrix
    )$predictions
    if (
      !"Occurrence" %in%
        colnames(
          validation_prediction
        )
    ) {
      stop(
        "Occurrence probability column missing.",
        call. = FALSE
      )
    }
    validation_probability <-
      as.numeric(
        validation_prediction[
          ,
          "Occurrence"
        ]
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
        "Invalid validation probability vector.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    positive_row <- data.frame(
      Algorithm = "RF",
      Run_ID = run_id,
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
      Observed_Class = 1L,
      Probability =
        held_out_probability,
      stringsAsFactors = FALSE
    )
    # --------------------------------------------------------
    # 10.9 Create pseudo-background rows
    # --------------------------------------------------------
    background_rows <- data.frame(
      Algorithm = "RF",
      Run_ID = run_id,
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
      Observed_Class = 0L,
      Probability =
        validation_probability,
      stringsAsFactors = FALSE
    )
    combined_predictions <- rbind(
      positive_row,
      background_rows
    )
    # --------------------------------------------------------
    # 10.10 Run summary
    # --------------------------------------------------------
    run_summary <- data.frame(
      Algorithm = "RF",
      Run_ID = run_id,
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
      QAQC_Status = "PASS",
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
        Algorithm = "RF",
        Run_ID = run_id,
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
  if (
    "error" %in%
      names(current_result)
  ) {
    error_rows[[length(error_rows) + 1L]] <- current_result$error
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
      " RF independent-validation run(s) failed. ",
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
    "150 successful RF runs",
    "15150 validation prediction rows",
    "150 held-out occurrence predictions",
    "15000 independent pseudo-background predictions",
    "No training-validation CellID overlap detected",
    "All reconstructed held-out probabilities reproduce originals",
    "All probabilities within 0-1",
    "All run QAQC statuses PASS"
  ),
  Result = c(
    nrow(run_summary_df) == 150L,
    nrow(prediction_df) == 15150L,
    sum(
      prediction_df$
        Observed_Class == 1L
    ) == 150L,
    sum(
      prediction_df$
        Observed_Class == 0L
    ) == 15000L,
    TRUE,
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
    "RF independent-validation global QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 14. SAVE OUTPUTS
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
# 15. QAQC REPORT
# ------------------------------------------------------------
qaqc_lines <- c(
  "RF INDEPENDENT VALIDATION QA/QC",
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
  "- Fold-specific training data were reconstructed using the original LOOCV run definition.",
  "- Validation pseudo-background cells were never included in model fitting.",
  "- Validation pseudo-background cells were never used in preprocessing or tuning.",
  "- Validation CellIDs were checked against outer-training CellIDs.",
  "- The held-out occurrence remained outside model fitting.",
  "- The reconstructed held-out probability was required to reproduce the original stored LOOCV probability.",
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  qaqc_file
)
# ------------------------------------------------------------
# 16. METHODOLOGY REPORT
# ------------------------------------------------------------
methodology_lines <- c(
  "RF INDEPENDENT VALIDATION METHODOLOGY",
  "============================================",
  "",
  paste(
    "Each of the 150 accepted RF outer-LOOCV models",
    "was reconstructed from its original fold-specific",
    "training dataset."
  ),
  "",
  paste(
    "Previously selected RF hyperparameters and the",
    "stored RF random seed were reused. No tuning was",
    "repeated and the validation data had no influence",
    "on model configuration."
  ),
  "",
  paste(
    "For each run, one genuinely held-out occurrence",
    "and 100 independent pseudo-background cells were",
    "predicted."
  ),
  "",
  paste(
    "The held-out probability was recalculated only",
    "to verify exact model reconstruction and was",
    "compared against the probability stored during",
    "the original LOOCV workflow."
  ),
  "",
  paste(
    "Pseudo-background cells are not confirmed",
    "mineralization absences. Subsequent confusion",
    "matrix metrics therefore quantify discrimination",
    "between known occurrences and pseudo-background."
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
cat("RF INDEPENDENT VALIDATION COMPLETE\n")
cat("============================================\n")
cat(
  "Successful runs              :",
  nrow(run_summary_df),
  "\n"
)
cat(
  "Held-out occurrence rows     :",
  sum(
    prediction_df$Observed_Class == 1
  ),
  "\n"
)
cat(
  "Pseudo-background rows       :",
  sum(
    prediction_df$Observed_Class == 0
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
  "RF INDEPENDENT VALIDATION: PASS\n"
)
