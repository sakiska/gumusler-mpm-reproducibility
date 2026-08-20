# ============================================================
# Train_SVM_LOOCV.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Train and tune radial SVM models within the repeated outer validation framework.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("e1071")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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
output_directory <- "Project/Models/SVM_Results"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create("Project/QAQC", recursive = TRUE, showWarnings = FALSE)
run_results_file <- file.path(
  output_directory,
  "SVM_LOOCV_Run_Results.csv"
)
tuning_results_file <- file.path(
  output_directory,
  "SVM_LOOCV_Tuning_Results.csv"
)
error_log_file <- file.path(
  output_directory,
  "SVM_LOOCV_Error_Log.csv"
)
qaqc_file <- "Project/QAQC/SVM_LOOCV_QAQC.txt"
methodology_file <- "Project/QAQC/SVM_LOOCV_Methodology.txt"
repeat_values <- seq_len(30)
fold_values <- seq_len(5)
expected_run_n <- length(repeat_values) * length(fold_values)
base_seed <- 9000L
inner_nfold <- 4L
parameter_grid <- expand.grid(
  cost = c(0.1, 1, 10, 100),
  gamma = c(0.01, 0.1, 1),
  stringsAsFactors = FALSE
)
safe_numeric <- function(x) {
  as.numeric(as.character(x))
}
binary_logloss <- function(observed, probability) {
  epsilon <- 1e-15
  probability <- pmin(
    pmax(as.numeric(probability), epsilon),
    1 - epsilon
  )
  observed <- safe_numeric(observed)
  -mean(
    observed * log(probability) +
      (1 - observed) * log(1 - probability)
  )
}
extract_occurrence_probability <- function(prediction_object) {
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
  probability_matrix <- as.matrix(probability_matrix)
  if (!"Occurrence" %in% colnames(probability_matrix)) {
    stop(
      paste0(
        "Occurrence probability column was not found. Available columns: ",
        paste(colnames(probability_matrix), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  as.numeric(
    probability_matrix[, "Occurrence"]
  )
}
make_stratified_folds <- function(y, k, seed) {
  y <- safe_numeric(y)
  positive_indices <- which(y == 1)
  negative_indices <- which(y == 0)
  k <- min(
    as.integer(k),
    length(positive_indices),
    length(negative_indices)
  )
  if (k < 2L) {
    stop(
      "Insufficient class counts for stratified inner CV.",
      call. = FALSE
    )
  }
  set.seed(seed)
  positive_indices <- sample(positive_indices)
  negative_indices <- sample(negative_indices)
  fold_assignment <- integer(length(y))
  fold_assignment[positive_indices] <- rep(
    seq_len(k),
    length.out = length(positive_indices)
  )
  fold_assignment[negative_indices] <- rep(
    seq_len(k),
    length.out = length(negative_indices)
  )
  lapply(
    seq_len(k),
    function(fold_i) {
      which(fold_assignment == fold_i)
    }
  )
}
make_class_weights <- function(y_factor) {
  class_counts <- table(y_factor)
  if (
    !"PseudoBackground" %in% names(class_counts) ||
      !"Occurrence" %in% names(class_counts)
  ) {
    stop(
      "Both classes are required for class weighting.",
      call. = FALSE
    )
  }
  c(
    PseudoBackground = 1,
    Occurrence =
      as.numeric(class_counts["PseudoBackground"]) /
      as.numeric(class_counts["Occurrence"])
  )
}
fit_svm_model <- function(
    x,
    y,
    cost,
    gamma,
    class_weights,
    seed
) {
  set.seed(seed)
  e1071::svm(
    x = x,
    y = y,
    type = "C-classification",
    kernel = "radial",
    cost = cost,
    gamma = gamma,
    class.weights = class_weights,
    probability = TRUE,
    scale = FALSE,
    fitted = FALSE
  )
}
run_results <- vector("list", expected_run_n)
tuning_results <- list()
error_results <- list()
run_counter <- 0L
cat("\nSVM LOOCV training started.\n")
cat("Expected runs:", expected_run_n, "\n\n")
for (repeat_i in repeat_values) {
  for (fold_i in fold_values) {
    run_counter <- run_counter + 1L
    run_id <- sprintf(
      "Repeat_%02d_Fold_%02d",
      repeat_i,
      fold_i
    )
    cat(
      "[",
      run_counter,
      "/",
      expected_run_n,
      "] ",
      run_id,
      "\n",
      sep = ""
    )
    result <- tryCatch({
      fold_data <- prepare_fold_data(
        run_id = run_id
      )
      model_matrix <- prepare_model_matrix(
        fold_data = fold_data,
        predictor_columns = predictor_columns,
        categorical_predictors = categorical_predictors,
        scale_numeric = TRUE
      )
      x_train <- as.matrix(model_matrix$x_train)
      x_test <- as.matrix(model_matrix$x_test)
      y_train_numeric <- safe_numeric(
        model_matrix$y_train_numeric
      )
      y_train_factor <- model_matrix$y_train_factor
      if (length(y_train_numeric) != nrow(x_train)) {
        stop(
          "Training response length does not match x_train.",
          call. = FALSE
        )
      }
      if (nrow(x_test) != 1L) {
        stop(
          "Each LOOCV run must contain exactly one held-out occurrence.",
          call. = FALSE
        )
      }
      if (
        anyNA(x_train) ||
          anyNA(x_test) ||
          anyNA(y_train_numeric)
      ) {
        stop(
          "Model matrix contains missing values.",
          call. = FALSE
        )
      }
      if (
        any(!is.finite(x_train)) ||
          any(!is.finite(x_test))
      ) {
        stop(
          "Model matrix contains non-finite values.",
          call. = FALSE
        )
      }
      if (!all(y_train_numeric %in% c(0, 1))) {
        stop(
          "Training response must be coded as 0/1.",
          call. = FALSE
        )
      }
      positive_n <- sum(y_train_numeric == 1)
      negative_n <- sum(y_train_numeric == 0)
      if (positive_n < 2L || negative_n < 2L) {
        stop(
          "Insufficient class counts for SVM training.",
          call. = FALSE
        )
      }
      inner_seed <- (
        base_seed +
          repeat_i * 100L +
          fold_i * 10L
      )
      inner_folds <- make_stratified_folds(
        y = y_train_numeric,
        k = inner_nfold,
        seed = inner_seed
      )
      best_score <- Inf
      best_row <- NULL
      for (grid_i in seq_len(nrow(parameter_grid))) {
        grid_row <- parameter_grid[
          grid_i,
          ,
          drop = FALSE
        ]
        cv_observed <- numeric(0)
        cv_probability <- numeric(0)
        for (inner_fold_i in seq_along(inner_folds)) {
          validation_indices <- inner_folds[[inner_fold_i]]
          training_indices <- setdiff(
            seq_len(nrow(x_train)),
            validation_indices
          )
          inner_x_train <- x_train[
            training_indices,
            ,
            drop = FALSE
          ]
          inner_x_validation <- x_train[
            validation_indices,
            ,
            drop = FALSE
          ]
          inner_y_train <- droplevels(
            y_train_factor[training_indices]
          )
          inner_y_validation <- y_train_numeric[
            validation_indices
          ]
          if (length(unique(inner_y_train)) != 2L) {
            stop(
              "An inner-training fold does not contain both classes.",
              call. = FALSE
            )
          }
          inner_class_weights <- make_class_weights(
            inner_y_train
          )
          current_seed <- (
            base_seed +
              repeat_i * 10000L +
              fold_i * 1000L +
              grid_i * 10L +
              inner_fold_i
          )
          inner_model <- fit_svm_model(
            x = inner_x_train,
            y = inner_y_train,
            cost = grid_row$cost,
            gamma = grid_row$gamma,
            class_weights = inner_class_weights,
            seed = current_seed
          )
          inner_prediction <- predict(
            inner_model,
            inner_x_validation,
            probability = TRUE
          )
          inner_probability <- extract_occurrence_probability(
            inner_prediction
          )
          cv_observed <- c(
            cv_observed,
            inner_y_validation
          )
          cv_probability <- c(
            cv_probability,
            inner_probability
          )
        }
        current_logloss <- binary_logloss(
          observed = cv_observed,
          probability = cv_probability
        )
        tuning_results[[length(tuning_results) + 1L]] <- data.frame(
          Algorithm = "SVM",
          Run_ID = run_id,
          Repeat = repeat_i,
          Fold_ID = sprintf("Fold_%02d", fold_i),
          Grid_ID = grid_i,
          Kernel = "radial",
          Cost = grid_row$cost,
          Gamma = grid_row$gamma,
          Inner_nfold = length(inner_folds),
          Inner_CV_Logloss = current_logloss,
          stringsAsFactors = FALSE
        )
        if (
          is.finite(current_logloss) &&
            current_logloss < best_score
        ) {
          best_score <- current_logloss
          best_row <- grid_row
        }
      }
      if (
        is.null(best_row) ||
          !is.finite(best_score)
      ) {
        stop(
          "No valid SVM tuning result was obtained.",
          call. = FALSE
        )
      }
      final_class_weights <- make_class_weights(
        y_train_factor
      )
      final_seed <- (
        base_seed +
          repeat_i * 100L +
          fold_i
      )
      final_model <- fit_svm_model(
        x = x_train,
        y = y_train_factor,
        cost = best_row$cost,
        gamma = best_row$gamma,
        class_weights = final_class_weights,
        seed = final_seed
      )
      held_out_prediction <- predict(
        final_model,
        x_test,
        probability = TRUE
      )
      held_out_probability <- extract_occurrence_probability(
        held_out_prediction
      )
      if (
        length(held_out_probability) != 1L ||
          !is.finite(held_out_probability) ||
          held_out_probability < 0 ||
          held_out_probability > 1
      ) {
        stop(
          "Invalid held-out occurrence probability.",
          call. = FALSE
        )
      }
      held_out_cell_id <- if (
        "Held_Out_CellID" %in%
          names(fold_data$run_metadata)
      ) {
        fold_data$run_metadata$Held_Out_CellID[1]
      } else {
        fold_data$test_data$CellID[1]
      }
      run_results[[run_counter]] <- data.frame(
        Algorithm = "SVM",
        Run_ID = run_id,
        Repeat = repeat_i,
        Fold_ID = sprintf("Fold_%02d", fold_i),
        Fold_Number = fold_i,
        Held_Out_CellID = held_out_cell_id,
        Held_Out_Occurrence_Probability =
          held_out_probability,
        Kernel = "radial",
        Cost = best_row$cost,
        Gamma = best_row$gamma,
        Class_Weight_PseudoBackground =
          unname(
            final_class_weights["PseudoBackground"]
          ),
        Class_Weight_Occurrence =
          unname(
            final_class_weights["Occurrence"]
          ),
        Inner_nfold = length(inner_folds),
        Inner_CV_Logloss = best_score,
        Scaling_Applied =
          model_matrix$preprocessing$scaling_applied,
        SVM_Seed = final_seed,
        QAQC_Status = "PASS",
        stringsAsFactors = FALSE
      )
      cat("   PASS\n")
      TRUE
    }, error = function(e) {
      error_results[[length(error_results) + 1L]] <<- data.frame(
        Algorithm = "SVM",
        Run_ID = run_id,
        Repeat = repeat_i,
        Fold_ID = sprintf("Fold_%02d", fold_i),
        Error_Message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      cat(
        "   FAIL: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )
      FALSE
    })
  }
}
successful_results <- run_results[
  !vapply(run_results, is.null, logical(1))
]
run_results_df <- if (
  length(successful_results) > 0
) {
  do.call(rbind, successful_results)
} else {
  data.frame()
}
tuning_results_df <- if (
  length(tuning_results) > 0
) {
  do.call(rbind, tuning_results)
} else {
  data.frame()
}
error_results_df <- if (
  length(error_results) > 0
) {
  do.call(rbind, error_results)
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
  run_results_df,
  run_results_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  tuning_results_df,
  tuning_results_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  error_results_df,
  error_log_file,
  row.names = FALSE,
  na = ""
)
overall_qaqc <- if (
  nrow(run_results_df) == expected_run_n &&
    nrow(error_results_df) == 0 &&
    all(run_results_df$QAQC_Status == "PASS") &&
    all(run_results_df$Scaling_Applied)
) {
  "PASS"
} else {
  "FAIL"
}
writeLines(
  c(
    "SVM LOOCV QAQC",
    "",
    paste0("Expected runs: ", expected_run_n),
    paste0("Successful runs: ", nrow(run_results_df)),
    paste0("Failed runs: ", nrow(error_results_df)),
    paste0("Overall QAQC status: ", overall_qaqc)
  ),
  qaqc_file
)
writeLines(
  c(
    "SVM LOOCV METHODOLOGY",
    "",
    paste(
      "The same repeated LOOCV folds, fold-specific",
      "pseudo-background samples and predictor matrix",
      "used for RF and XGBoost were retained."
    ),
    paste(
      "Numeric predictors were standardized using",
      "parameters learned exclusively from each",
      "outer-training set."
    ),
    paste(
      "Cost and gamma were selected independently",
      "within each outer run using stratified",
      "inner cross-validation and binary logloss."
    ),
    paste(
      "The final radial SVM was refitted on the",
      "complete outer-training set and evaluated",
      "only at the genuinely held-out occurrence."
    ),
    paste(
      "Pseudo-background locations were not treated",
      "as verified absences in final performance reporting."
    )
  ),
  methodology_file
)
cat("\n========================================\n")
cat("SVM LOOCV training completed.\n")
cat("Successful runs:", nrow(run_results_df), "\n")
cat("Failed runs:", nrow(error_results_df), "\n")
cat("Overall QAQC:", overall_qaqc, "\n")
cat("Run results:", run_results_file, "\n")
cat("========================================\n")
if (overall_qaqc != "PASS") {
  stop(
    "SVM LOOCV QAQC failed. Review the error log.",
    call. = FALSE
  )
}
