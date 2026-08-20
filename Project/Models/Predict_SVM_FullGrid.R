# ============================================================
# Predict_SVM_FullGrid.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Apply accepted SVM models to the complete modeling grid.
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
      paste(missing_packages, collapse = ", "),
      ". Install with install.packages()."
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
full_grid_file <- "Project/Data/MPM_FullGrid_25m.csv"
run_results_file <- paste0(
  "Project/Models/SVM_Results/",
  "SVM_LOOCV_Run_Results.csv"
)
output_directory <- "Project/Models/SVM_Results/FullGrid"
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
summary_file <- file.path(
  output_directory,
  "SVM_FullGrid_Ensemble_Summary.csv"
)
threshold_file <- file.path(
  output_directory,
  "SVM_FullGrid_Run_Thresholds.csv"
)
capture_file <- file.path(
  output_directory,
  "SVM_HeldOut_Threshold_Capture.csv"
)
capture_summary_file <- file.path(
  output_directory,
  "SVM_Threshold_Capture_Summary.csv"
)
error_file <- file.path(
  output_directory,
  "SVM_FullGrid_Error_Log.csv"
)
qaqc_file <- "Project/QAQC/SVM_FullGrid_QAQC.txt"
methodology_file <- paste0(
  "Project/QAQC/",
  "SVM_FullGrid_Methodology.txt"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!file.exists(full_grid_file)) {
  stop(
    paste0("Missing file: ", full_grid_file),
    call. = FALSE
  )
}
if (!file.exists(run_results_file)) {
  stop(
    paste0("Missing file: ", run_results_file),
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
required_grid_columns <- c(
  "CellID",
  "X",
  "Y",
  predictor_columns
)
missing_grid_columns <- setdiff(
  required_grid_columns,
  names(full_grid)
)
if (length(missing_grid_columns) > 0) {
  stop(
    paste0(
      "Full grid missing columns: ",
      paste(missing_grid_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}
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
      paste(missing_run_columns, collapse = ", ")
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
if (any(!run_results$Scaling_Applied)) {
  stop(
    "At least one SVM run was trained without scaling.",
    call. = FALSE
  )
}
if (anyDuplicated(full_grid$CellID) > 0) {
  stop(
    "Full grid contains duplicated CellID values.",
    call. = FALSE
  )
}
if (anyNA(full_grid[, required_grid_columns])) {
  stop(
    "Full grid contains missing required values.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
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
      paste0(
        "Occurrence probability column was not found. ",
        "Available columns: ",
        paste(
          colnames(probability_matrix),
          collapse = ", "
        )
      ),
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
build_full_grid_matrix <- function(
    full_grid,
    model_matrix
) {
  preprocessing <- model_matrix$preprocessing
  predictor_data <- full_grid[
    ,
    predictor_columns,
    drop = FALSE
  ]
  numeric_predictors <- model_matrix$numeric_predictors
  for (column_name in numeric_predictors) {
    values <- safe_numeric(
      predictor_data[[column_name]]
    )
    values[is.infinite(values)] <- NA_real_
    imputation_value <- preprocessing$
      numeric_imputation_values[[column_name]]
    center_value <- preprocessing$
      numeric_center_values[[column_name]]
    scale_value <- preprocessing$
      numeric_scale_values[[column_name]]
    if (
      length(imputation_value) != 1L ||
        !is.finite(imputation_value)
    ) {
      stop(
        paste0(
          "Invalid imputation value for ",
          column_name,
          "."
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
          column_name,
          "."
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
          column_name,
          "."
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
    values <- clean_factor_values(
      predictor_data[[column_name]]
    )
    training_levels <- preprocessing$
      factor_level_map[[column_name]]
    unknown_level <- preprocessing$
      unknown_factor_level
    if (
      length(training_levels) == 0L ||
        !unknown_level %in% training_levels
    ) {
      stop(
        paste0(
          "Invalid factor-level map for ",
          column_name,
          "."
        ),
        call. = FALSE
      )
    }
    unseen_values <- !values %in% training_levels
    values[unseen_values] <- unknown_level
    predictor_data[[column_name]] <- factor(
      values,
      levels = training_levels
    )
  }
  encoded_data <- stats::model.matrix(
    object = ~ . - 1,
    data = predictor_data,
    na.action = stats::na.pass
  )
  encoded_data <- as.data.frame(
    encoded_data,
    check.names = FALSE
  )
  training_columns <- model_matrix$final_matrix_columns
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
      setdiff(names(encoded_data), extra_columns),
      drop = FALSE
    ]
  }
  encoded_data <- encoded_data[
    ,
    training_columns,
    drop = FALSE
  ]
  full_matrix <- as.matrix(
    encoded_data
  )
  storage.mode(full_matrix) <- "double"
  if (
    nrow(full_matrix) != nrow(full_grid) ||
      ncol(full_matrix) != length(training_columns)
  ) {
    stop(
      "Full-grid matrix dimensions are invalid.",
      call. = FALSE
    )
  }
  if (
    anyNA(full_matrix) ||
      any(!is.finite(full_matrix))
  ) {
    stop(
      "Full-grid model matrix contains invalid values.",
      call. = FALSE
    )
  }
  if (!identical(
    colnames(full_matrix),
    training_columns
  )) {
    stop(
      "Full-grid and training matrix columns are not identical.",
      call. = FALSE
    )
  }
  full_matrix
}
calculate_thresholds <- function(
    probability
) {
  mean_probability <- mean(probability)
  sd_probability <- stats::sd(probability)
  data.frame(
    Grid_Mean = mean_probability,
    Grid_SD = sd_probability,
    Threshold_P90 = as.numeric(
      stats::quantile(
        probability,
        0.90,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_P95 = as.numeric(
      stats::quantile(
        probability,
        0.95,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_P99 = as.numeric(
      stats::quantile(
        probability,
        0.99,
        names = FALSE,
        type = 7
      )
    ),
    Threshold_Mean_Plus_2SD = min(
      1,
      mean_probability + 2 * sd_probability
    ),
    stringsAsFactors = FALSE
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
threshold_rows <- vector("list", run_n)
capture_rows <- vector("list", run_n)
error_rows <- list()
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\nSVM full-grid prediction started.\n")
cat("Runs:", run_n, "\n")
cat("Grid cells:", grid_n, "\n\n")
for (run_i in seq_len(run_n)) {
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
    run_n,
    "] ",
    run_id,
    "\n",
    sep = ""
  )
  current_result <- tryCatch({
    fold_data <- prepare_fold_data(
      run_id = run_id
    )
    model_matrix <- prepare_model_matrix(
      fold_data = fold_data,
      predictor_columns = predictor_columns,
      categorical_predictors = categorical_predictors,
      scale_numeric = TRUE
    )
    x_train <- as.matrix(
      model_matrix$x_train
    )
    y_train_factor <- model_matrix$y_train_factor
    if (
      ncol(x_train) !=
        length(model_matrix$final_matrix_columns)
    ) {
      stop(
        "Training matrix column count is inconsistent.",
        call. = FALSE
      )
    }
    full_matrix <- build_full_grid_matrix(
      full_grid = full_grid,
      model_matrix = model_matrix
    )
    class_weights <- c(
      PseudoBackground = safe_numeric(
        current_run$
          Class_Weight_PseudoBackground[1]
      ),
      Occurrence = safe_numeric(
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
    model <- fit_probability_svm(
      x = x_train,
      y_factor = y_train_factor,
      cost = safe_numeric(
        current_run$Cost[1]
      ),
      gamma = safe_numeric(
        current_run$Gamma[1]
      ),
      class_weights = class_weights,
      seed = as.integer(
        current_run$SVM_Seed[1]
      )
    )
    full_prediction <- predict(
      model,
      full_matrix,
      probability = TRUE
    )
    probability <- extract_occurrence_probability(
      full_prediction
    )
    if (length(probability) != grid_n) {
      stop(
        "Full-grid probability vector has invalid length.",
        call. = FALSE
      )
    }
    threshold_row <- calculate_thresholds(
      probability
    )
    threshold_row$Algorithm <- "SVM"
    threshold_row$Run_ID <- run_id
    threshold_row$Repeat <- current_run$Repeat[1]
    threshold_row$Fold_ID <- current_run$Fold_ID[1]
    threshold_row$Fold_Number <-
      current_run$Fold_Number[1]
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
    held_out_probability <- safe_numeric(
      current_run$
        Held_Out_Occurrence_Probability[1]
    )
    capture_row <- data.frame(
      Algorithm = "SVM",
      Run_ID = run_id,
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
  }, error = function(e) {
    list(
      error = data.frame(
        Algorithm = "SVM",
        Run_ID = run_id,
        Repeat = current_run$Repeat[1],
        Fold_ID = current_run$Fold_ID[1],
        Error_Message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  })
  if ("error" %in% names(current_result)) {
    error_rows[[length(error_rows) + 1L]] <-
      current_result$error
    cat(
      "   FAIL: ",
      current_result$error$Error_Message,
      "\n",
      sep = ""
    )
  } else {
    probability <- current_result$probability
    probability_sum <- (
      probability_sum + probability
    )
    probability_sum_squares <- (
      probability_sum_squares +
        probability^2
    )
    probability_min <- pmin(
      probability_min,
      probability
    )
    probability_max <- pmax(
      probability_max,
      probability
    )
    prediction_count <- (
      prediction_count + 1L
    )
    threshold_rows[[run_i]] <-
      current_result$threshold
    capture_rows[[run_i]] <-
      current_result$capture
    cat("   PASS\n")
  }
}
# ------------------------------------------------------------
# 6. Error log and successful-run QA/QC
# ------------------------------------------------------------
error_df <- if (length(error_rows) > 0) {
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
if (nrow(error_df) > 0) {
  stop(
    paste0(
      nrow(error_df),
      " full-grid run(s) failed. See ",
      error_file
    ),
    call. = FALSE
  )
}
if (any(prediction_count != run_n)) {
  stop(
    "At least one grid cell does not have 150 predictions.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 7. Ensemble summary
# ------------------------------------------------------------
mean_probability <- (
  probability_sum / prediction_count
)
variance_probability <- (
  probability_sum_squares -
    probability_sum^2 / prediction_count
) / (
  prediction_count - 1
)
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
  SVM_Mean_Probability =
    mean_probability,
  SVM_SD_Probability =
    sd_probability,
  SVM_Min_Probability =
    probability_min,
  SVM_Max_Probability =
    probability_max,
  SVM_Prediction_n =
    prediction_count,
  stringsAsFactors = FALSE
)
run_thresholds <- do.call(
  rbind,
  threshold_rows
)
held_out_capture <- do.call(
  rbind,
  capture_rows
)
write.csv(
  full_grid_summary,
  summary_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  run_thresholds,
  threshold_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  held_out_capture,
  capture_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 8. Threshold-capture summary
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
summary_rows <- vector(
  "list",
  length(threshold_names)
)
for (i in seq_along(threshold_names)) {
  threshold_name <- threshold_names[i]
  capture_values <- held_out_capture[[capture_columns[[threshold_name]]]]
  threshold_values <- held_out_capture[[threshold_columns[[threshold_name]]]]
  summary_rows[[i]] <- data.frame(
    Algorithm = "SVM",
    Threshold_Method =
      threshold_name,
    Run_n = nrow(held_out_capture),
    Captured_Run_n =
      sum(capture_values),
    Capture_Rate =
      mean(capture_values),
    Mean_Threshold =
      mean(threshold_values),
    SD_Threshold =
      stats::sd(threshold_values),
    Minimum_Threshold =
      min(threshold_values),
    Maximum_Threshold =
      max(threshold_values),
    stringsAsFactors = FALSE
  )
}
threshold_summary <- do.call(
  rbind,
  summary_rows
)
write.csv(
  threshold_summary,
  capture_summary_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
invalid_mean_n <- sum(
  !is.finite(
    full_grid_summary$SVM_Mean_Probability
  ) |
    full_grid_summary$SVM_Mean_Probability < 0 |
    full_grid_summary$SVM_Mean_Probability > 1
)
invalid_sd_n <- sum(
  !is.finite(
    full_grid_summary$SVM_SD_Probability
  ) |
    full_grid_summary$SVM_SD_Probability < 0
)
overall_qaqc <- if (
  nrow(full_grid_summary) == grid_n &&
    all(
      full_grid_summary$SVM_Prediction_n ==
        150L
    ) &&
    nrow(run_thresholds) == 150L &&
    nrow(held_out_capture) == 150L &&
    nrow(error_df) == 0L &&
    invalid_mean_n == 0L &&
    invalid_sd_n == 0L
) {
  "PASS"
} else {
  "FAIL"
}
writeLines(
  c(
    "SVM FULL-GRID QAQC",
    "",
    "Expected runs: 150",
    paste0(
      "Successful runs: ",
      nrow(run_thresholds)
    ),
    paste0(
      "Failed runs: ",
      nrow(error_df)
    ),
    paste0(
      "Grid cells: ",
      nrow(full_grid_summary)
    ),
    paste0(
      "Cells with 150 predictions: ",
      sum(
        full_grid_summary$
          SVM_Prediction_n == 150L
      )
    ),
    paste0(
      "Held-out capture records: ",
      nrow(held_out_capture)
    ),
    paste0(
      "Invalid ensemble means: ",
      invalid_mean_n
    ),
    paste0(
      "Invalid ensemble SD values: ",
      invalid_sd_n
    ),
    paste0(
      "Overall QAQC status: ",
      overall_qaqc
    )
  ),
  qaqc_file
)
writeLines(
  c(
    "SVM FULL-GRID METHODOLOGY",
    "",
    paste(
      "Each accepted outer-LOOCV radial SVM model",
      "was reconstructed using its fold-specific",
      "training data, selected Cost and Gamma, stored",
      "class weights and random seed."
    ),
    paste(
      "Numeric full-grid predictors were imputed and",
      "standardized using only the corresponding",
      "outer-training medians, means and standard",
      "deviations."
    ),
    paste(
      "Categorical levels and one-hot encoded model",
      "columns were reproduced from the corresponding",
      "outer-training preprocessing object."
    ),
    paste(
      "Each reconstructed model predicted all valid",
      "25 m grid cells."
    ),
    paste(
      "Cell-level ensemble mean, standard deviation,",
      "minimum and maximum occurrence probabilities",
      "were calculated from 150 predictions."
    ),
    paste(
      "P90, P95, P99 and Mean + 2 SD thresholds were",
      "calculated separately from each run's full-grid",
      "probability distribution."
    ),
    paste(
      "Threshold capture was evaluated only with the",
      "corresponding genuinely held-out occurrence",
      "probability."
    )
  ),
  methodology_file
)
# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------
cat("\n========================================\n")
cat("SVM full-grid prediction completed.\n")
cat(
  "Successful runs:",
  nrow(run_thresholds),
  "\n"
)
cat(
  "Grid cells:",
  nrow(full_grid_summary),
  "\n"
)
cat(
  "Overall QAQC:",
  overall_qaqc,
  "\n"
)
cat(
  "Full-grid summary:",
  summary_file,
  "\n"
)
cat(
  "Run thresholds:",
  threshold_file,
  "\n"
)
cat(
  "Held-out capture:",
  capture_file,
  "\n"
)
cat(
  "Threshold summary:",
  capture_summary_file,
  "\n"
)
cat("========================================\n")
if (overall_qaqc != "PASS") {
  stop(
    "SVM full-grid QAQC failed.",
    call. = FALSE
  )
}
