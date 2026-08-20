# ============================================================
# Predict_SVM_Independent_Validation.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Generate SVM predictions for the independent held-out occurrence and pseudo-background validation sets.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("e1071")
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
# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------
run_results_file <- paste0(
  "Project/Models/SVM_Results/",
  "SVM_LOOCV_Run_Results.csv"
)
validation_file <- paste0(
  "Project/Tables/",
  "Independent_Validation_PseudoBackground.csv"
)
output_directory <- paste0(
  "Project/Models/SVM_Results/",
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
  "SVM_Independent_Validation_Predictions.csv"
)
run_summary_file <- file.path(
  output_directory,
  "SVM_Independent_Validation_Run_Summary.csv"
)
error_file <- file.path(
  output_directory,
  "SVM_Independent_Validation_Error_Log.csv"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "SVM_Independent_Validation_QAQC.txt"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "SVM_Independent_Validation_Methodology.txt"
)
# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------
safe_numeric <- function(x) {
  as.numeric(as.character(x))
}
extract_occurrence_probability <- function(
    prediction_object
) {
  probability_matrix <- attr(
    prediction_object,
    "probabilities"
  )
  if (is.null(probability_matrix)) {
    stop(
      "SVM prediction did not return class probabilities.",
      call. = FALSE
    )
  }
  probability_matrix <- as.matrix(
    probability_matrix
  )
  if (
    !"Occurrence" %in%
      colnames(probability_matrix)
  ) {
    stop(
      "Occurrence probability column was not found.",
      call. = FALSE
    )
  }
  probability <- as.numeric(
    probability_matrix[, "Occurrence"]
  )
  if (
    length(probability) == 0L ||
      anyNA(probability) ||
      any(!is.finite(probability)) ||
      any(probability < 0) ||
      any(probability > 1)
  ) {
    stop(
      "Invalid SVM occurrence probabilities.",
      call. = FALSE
    )
  }
  probability
}
fit_probability_svm <- function(
    x,
    y_factor,
    cost,
    gamma,
    class_weights,
    seed
) {
  set.seed(seed)
  e1071::svm(
    x = x,
    y = y_factor,
    type = "C-classification",
    kernel = "radial",
    cost = cost,
    gamma = gamma,
    class.weights = class_weights,
    probability = TRUE,
    scale = FALSE,
    fitted = FALSE,
    cachesize = 200,
    tolerance = 0.001
  )
}
build_validation_matrix <- function(
    validation_data,
    model_matrix
) {
  preprocessing <- model_matrix$preprocessing
  predictor_data <- validation_data[
    ,
    predictor_columns,
    drop = FALSE
  ]
  numeric_predictors <-
    model_matrix$numeric_predictors
  for (column_name in numeric_predictors) {
    values <- safe_numeric(
      predictor_data[[column_name]]
    )
    values[is.infinite(values)] <- NA_real_
    imputation_value <-
      preprocessing$numeric_imputation_values[[column_name]]
    center_value <-
      preprocessing$numeric_center_values[[column_name]]
    scale_value <-
      preprocessing$numeric_scale_values[[column_name]]
    if (
      length(imputation_value) != 1L ||
        !is.finite(imputation_value)
    ) {
      stop(
        paste0(
          "Invalid imputation value for ",
          column_name
        ),
        call. = FALSE
      )
    }
    if (
      length(center_value) != 1L ||
        !is.finite(center_value)
    ) {
      stop(
        paste0(
          "Invalid center value for ",
          column_name
        ),
        call. = FALSE
      )
    }
    if (
      length(scale_value) != 1L ||
        !is.finite(scale_value) ||
        scale_value <= 0
    ) {
      stop(
        paste0(
          "Invalid scale value for ",
          column_name
        ),
        call. = FALSE
      )
    }
    values[is.na(values)] <- imputation_value
    values <- (
      values - center_value
    ) / scale_value
    predictor_data[[column_name]] <- values
  }
  for (
    column_name in
      model_matrix$categorical_predictors
  ) {
    values <- as.character(
      predictor_data[[column_name]]
    )
    values <- trimws(values)
    values[
      is.na(values) |
        values == ""
    ] <- "Missing"
    training_levels <-
      preprocessing$factor_level_map[[column_name]]
    unknown_level <-
      preprocessing$unknown_factor_level
    if (
      length(training_levels) == 0L ||
        !unknown_level %in% training_levels
    ) {
      stop(
        paste0(
          "Invalid factor-level map for ",
          column_name
        ),
        call. = FALSE
      )
    }
    unseen_values <- !values %in%
      training_levels
    values[unseen_values] <-
      unknown_level
    predictor_data[[column_name]] <-
      factor(
        values,
        levels = training_levels
      )
  }
  encoded_data <- stats::model.matrix(
    ~ . - 1,
    data = predictor_data,
    na.action = stats::na.pass
  )
  encoded_data <- as.data.frame(
    encoded_data,
    check.names = FALSE
  )
  training_columns <-
    model_matrix$final_matrix_columns
  missing_columns <- setdiff(
    training_columns,
    names(encoded_data)
  )
  if (length(missing_columns) > 0) {
    for (column_name in missing_columns) {
      encoded_data[[column_name]] <- 0
    }
  }
  extra_columns <- setdiff(
    names(encoded_data),
    training_columns
  )
  if (length(extra_columns) > 0) {
    encoded_data <- encoded_data[
      ,
      setdiff(
        names(encoded_data),
        extra_columns
      ),
      drop = FALSE
    ]
  }
  encoded_data <- encoded_data[
    ,
    training_columns,
    drop = FALSE
  ]
  validation_matrix <- as.matrix(
    encoded_data
  )
  storage.mode(validation_matrix) <- "double"
  if (
    nrow(validation_matrix) !=
      nrow(validation_data)
  ) {
    stop(
      "Validation matrix row count is invalid.",
      call. = FALSE
    )
  }
  if (
    anyNA(validation_matrix) ||
      any(!is.finite(validation_matrix))
  ) {
    stop(
      "Validation matrix contains invalid values.",
      call. = FALSE
    )
  }
  if (
    !identical(
      colnames(validation_matrix),
      training_columns
    )
  ) {
    stop(
      "Validation and training matrix columns differ.",
      call. = FALSE
    )
  }
  validation_matrix
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!file.exists(run_results_file)) {
  stop(
    paste0(
      "Missing file: ",
      run_results_file
    ),
    call. = FALSE
  )
}
if (!file.exists(validation_file)) {
  stop(
    paste0(
      "Missing file: ",
      validation_file
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
required_run_columns <- c(
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "Held_Out_CellID",
  "Held_Out_Occurrence_Probability",
  "Cost",
  "Gamma",
  "Class_Weight_PseudoBackground",
  "Class_Weight_Occurrence",
  "Scaling_Applied",
  "SVM_Seed",
  "QAQC_Status"
)
missing_run_columns <- setdiff(
  required_run_columns,
  names(run_results)
)
if (length(missing_run_columns) > 0) {
  stop(
    paste0(
      "Run results missing columns: ",
      paste(
        missing_run_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}
required_validation_columns <- c(
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "CellID",
  predictor_columns
)
missing_validation_columns <- setdiff(
  required_validation_columns,
  names(validation_data)
)
if (length(missing_validation_columns) > 0) {
  stop(
    paste0(
      "Validation table missing columns: ",
      paste(
        missing_validation_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}
if (
  nrow(run_results) != 150L ||
    any(run_results$QAQC_Status != "PASS")
) {
  stop(
    "Expected 150 accepted SVM runs.",
    call. = FALSE
  )
}
if (
  any(!run_results$Scaling_Applied)
) {
  stop(
    "At least one SVM run was trained without scaling.",
    call. = FALSE
  )
}
if (
  anyDuplicated(
    run_results$Run_ID
  ) > 0
) {
  stop(
    "Duplicated Run_ID values in SVM results.",
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
    "Each run must contain exactly 100 validation pseudo-background cells.",
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
cat("SVM INDEPENDENT VALIDATION STARTED\n")
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
      predictor_columns =
        predictor_columns,
      categorical_predictors =
        categorical_predictors,
      scale_numeric = TRUE
    )
    x_train <- as.matrix(
      model_matrix$x_train
    )
    x_test <- as.matrix(
      model_matrix$x_test
    )
    y_train_factor <-
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
        "Expected 100 validation pseudo-background cells.",
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
        "Validation CellIDs overlap outer-training CellIDs.",
        call. = FALSE
      )
    }
    if (
      heldout_overlap_n != 0L
    ) {
      stop(
        "Validation background contains held-out occurrence CellID.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    validation_matrix <- build_validation_matrix(
      validation_data =
        current_validation,
      model_matrix =
        model_matrix
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    class_weights <- c(
      PseudoBackground =
        safe_numeric(
          current_run$
            Class_Weight_PseudoBackground[1]
        ),
      Occurrence =
        safe_numeric(
          current_run$
            Class_Weight_Occurrence[1]
        )
    )
    if (
      any(!is.finite(class_weights)) ||
        any(class_weights <= 0)
    ) {
      stop(
        "Invalid stored SVM class weights.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # Reconstruct accepted SVM model
    # --------------------------------------------------------
    model <- fit_probability_svm(
      x = x_train,
      y_factor = y_train_factor,
      cost =
        safe_numeric(
          current_run$Cost[1]
        ),
      gamma =
        safe_numeric(
          current_run$Gamma[1]
        ),
      class_weights =
        class_weights,
      seed =
        as.integer(
          current_run$SVM_Seed[1]
        )
    )
    # --------------------------------------------------------
    # --------------------------------------------------------
    held_out_prediction <- predict(
      model,
      x_test,
      probability = TRUE
    )
    held_out_probability <-
      extract_occurrence_probability(
        held_out_prediction
      )
    if (
      length(held_out_probability) != 1L
    ) {
      stop(
        "Invalid held-out probability length.",
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
      model,
      validation_matrix,
      probability = TRUE
    )
    validation_probability <-
      extract_occurrence_probability(
        validation_prediction
      )
    if (
      length(validation_probability) != 100L
    ) {
      stop(
        "Invalid validation probability vector length.",
        call. = FALSE
      )
    }
    # --------------------------------------------------------
    # --------------------------------------------------------
    positive_row <- data.frame(
      Algorithm = "SVM",
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
    # Background rows
    # --------------------------------------------------------
    background_rows <- data.frame(
      Algorithm = "SVM",
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
    # Summary row
    # --------------------------------------------------------
    run_summary <- data.frame(
      Algorithm = "SVM",
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
        length(validation_probability),
      Mean_Validation_Background_Probability =
        mean(validation_probability),
      Median_Validation_Background_Probability =
        median(validation_probability),
      Maximum_Validation_Background_Probability =
        max(validation_probability),
      Training_Validation_Overlap_n =
        training_overlap_n,
      Validation_HeldOut_Overlap_n =
        heldout_overlap_n,
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
        Algorithm = "SVM",
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
    error_rows[[length(error_rows) + 1L]] <-
      current_result$error
    cat(
      "   FAIL: ",
      current_result$error$Error_Message,
      "\n",
      sep = ""
    )
  } else {
    prediction_rows[[run_i]] <-
      current_result$predictions
    summary_rows[[run_i]] <-
      current_result$summary
    cat("   PASS\n")
  }
}
# ------------------------------------------------------------
# 8. Error log
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
      " SVM independent-validation run(s) failed. ",
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
rownames(prediction_df) <- NULL
run_summary_df <- do.call(
  rbind,
  summary_rows
)
rownames(run_summary_df) <- NULL
# ------------------------------------------------------------
# 10. Global QA/QC
# ------------------------------------------------------------
global_qaqc <- data.frame(
  Check = c(
    "150 successful SVM runs",
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
    nrow(run_summary_df) == 150L,
    nrow(prediction_df) == 15150L,
    sum(
      prediction_df$Observed_Class == 1L
    ) == 150L,
    sum(
      prediction_df$Observed_Class == 0L
    ) == 15000L,
    all(
      run_summary_df$
        Training_Validation_Overlap_n == 0L
    ),
    all(
      run_summary_df$
        Validation_HeldOut_Overlap_n == 0L
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
print(global_qaqc)
if (
  !all(global_qaqc$Result)
) {
  stop(
    "SVM independent-validation global QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 11. Save outputs
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
# 12. QA/QC report
# ------------------------------------------------------------
qaqc_lines <- c(
  "SVM INDEPENDENT VALIDATION QA/QC",
  "============================================",
  "",
  paste(
    global_qaqc$Check,
    global_qaqc$Result,
    sep = ": "
  ),
  "",
  "Leakage-control conditions:",
  "- Numeric imputation, centering and scaling were learned exclusively from the corresponding outer-training set.",
  "- Independent validation pseudo-background cells were never used in model fitting or hyperparameter tuning.",
  "- SVM internal scaling remained disabled because scaling was performed explicitly using outer-training parameters.",
  "- Training and validation CellIDs were explicitly checked for overlap.",
  "- The held-out occurrence remained outside model fitting.",
  "- Stored Cost, Gamma, class weights and SVM seed were reused.",
  "- Reconstructed held-out probabilities were required to reproduce the original LOOCV probabilities.",
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  qaqc_file
)
# ------------------------------------------------------------
# 13. Methodology
# ------------------------------------------------------------
methodology_lines <- c(
  "SVM INDEPENDENT VALIDATION METHODOLOGY",
  "============================================",
  "",
  paste(
    "Each of the 150 accepted radial-SVM outer-LOOCV",
    "models was reconstructed using its original",
    "fold-specific training dataset."
  ),
  "",
  paste(
    "Numeric predictors in the independent validation",
    "set were imputed, centered and scaled using",
    "parameters learned exclusively from the",
    "corresponding outer-training set."
  ),
  "",
  paste(
    "Stored Cost, Gamma, class weights and SVM random",
    "seed were reused. Hyperparameter tuning was not",
    "repeated."
  ),
  "",
  paste(
    "For each run, one genuinely held-out occurrence",
    "and 100 independent pseudo-background cells were",
    "predicted."
  ),
  "",
  paste(
    "Pseudo-background cells are not verified mineral",
    "absences. Subsequent confusion-matrix metrics",
    "therefore quantify occurrence-versus-pseudo-",
    "background discrimination."
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
cat("SVM INDEPENDENT VALIDATION COMPLETE\n")
cat("============================================\n")
cat(
  "Successful runs              :",
  nrow(run_summary_df),
  "\n"
)
cat(
  "Held-out occurrence rows     :",
  sum(
    prediction_df$Observed_Class == 1L
  ),
  "\n"
)
cat(
  "Pseudo-background rows       :",
  sum(
    prediction_df$Observed_Class == 0L
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
print(global_qaqc)
cat("\n")
cat("SVM INDEPENDENT VALIDATION: PASS\n")
