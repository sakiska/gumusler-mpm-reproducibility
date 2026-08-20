# ============================================================
# Prepare_Model_Matrix.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Build fold-specific model matrices while learning transformations from training data only.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

assert_columns_exist <- function(
    data,
    columns,
    data_name
) {
  missing_columns <- setdiff(
    columns,
    names(data)
  )
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        data_name,
        " missing columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
}
# ------------------------------------------------------------
# ------------------------------------------------------------
replace_infinite_with_na <- function(x) {
  if (!is.numeric(x)) {
    return(x)
  }
  x[is.infinite(x)] <- NA_real_
  return(x)
}
# ------------------------------------------------------------
# ------------------------------------------------------------
clean_factor_values <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[
    is.na(x) |
      x == ""
  ] <- "Missing"
  return(x)
}
# ------------------------------------------------------------
# ------------------------------------------------------------
prepare_model_matrix <- function(
    fold_data,
    predictor_columns,
    categorical_predictors = character(0),
    scale_numeric = FALSE,
    remove_zero_variance = TRUE,
    unknown_factor_level = "Unknown"
) {
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (!is.list(fold_data)) {
    stop(
      "fold_data bir liste olmali.",
      call. = FALSE
    )
  }
  required_list_elements <- c(
    "training_data",
    "test_data",
    "audit"
  )
  missing_elements <- setdiff(
    required_list_elements,
    names(fold_data)
  )
  if (length(missing_elements) > 0) {
    stop(
      paste0(
        "Missing elements in the fold_data list: ",
        paste(
          missing_elements,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  training_data <- fold_data$training_data
  test_data <- fold_data$test_data
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (length(predictor_columns) == 0) {
    stop(
      "En az bir predictor belirtilmelidir.",
      call. = FALSE
    )
  }
  if (anyDuplicated(predictor_columns) > 0) {
    stop(
      "Duplicate column names were found in the predictor list.",
      call. = FALSE
    )
  }
  categorical_predictors <- unique(
    categorical_predictors
  )
  invalid_categorical <- setdiff(
    categorical_predictors,
    predictor_columns
  )
  if (length(invalid_categorical) > 0) {
    stop(
      paste0(
        "Kategorik olarak belirtilen ancak predictor ",
        " contains columns that were not found: ",
        paste(
          invalid_categorical,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  assert_columns_exist(
    training_data,
    c(
      predictor_columns,
      "Class",
      "CellID"
    ),
    "Training verisi"
  )
  assert_columns_exist(
    test_data,
    c(
      predictor_columns,
      "Class",
      "CellID"
    ),
    "Test verisi"
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  forbidden_predictors <- c(
    "Class",
    "Class_Label",
    "CellID",
    "Run_ID",
    "Fold_ID",
    "Fold_Number",
    "Repeat",
    "Background_Seed",
    "Buffer_m",
    "Data_Role",
    "Source_Type",
    "Occurrence_Row",
    "Inside_Reference_Mask",
    "Background_Eligible",
    "Min_Dist_to_Occurrence_m",
    "Min_Dist_to_Training_Occurrence_m"
  )
  forbidden_used <- intersect(
    predictor_columns,
    forbidden_predictors
  )
  if (length(forbidden_used) > 0) {
    stop(
      paste0(
        "Columns that must not enter the model were found in the predictor ",
        "olarak secilmis: ",
        paste(
          forbidden_used,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # 4.4. Class kontrolu
  # ----------------------------------------------------------
  if (!all(
    training_data$Class %in% c(
      0,
      1
    )
  )) {
    stop(
      "Training Class values must contain only 0 and 1.",
      call. = FALSE
    )
  }
  if (!all(
    test_data$Class %in% c(
      0,
      1
    )
  )) {
    stop(
      "Test Class values must contain only 0 and 1.",
      call. = FALSE
    )
  }
  if (length(unique(training_data$Class)) != 2) {
    stop(
      paste(
        "Training verisinde hem occurrence hem de",
        "pseudo-background sinifi bulunmali."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  train_predictors <- training_data[
    ,
    predictor_columns,
    drop = FALSE
  ]
  test_predictors <- test_data[
    ,
    predictor_columns,
    drop = FALSE
  ]
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  numeric_predictors <- setdiff(
    predictor_columns,
    categorical_predictors
  )
  non_numeric_unmarked <- numeric_predictors[
    !vapply(
      train_predictors[
        ,
        numeric_predictors,
        drop = FALSE
      ],
      is.numeric,
      logical(1)
    )
  ]
  if (length(non_numeric_unmarked) > 0) {
    stop(
      paste0(
        "The following predictors are not numeric and ",
        "categorical_predictors listesine eklenmemis: ",
        paste(
          non_numeric_unmarked,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  numeric_imputation_values <- numeric(0)
  numeric_center_values <- numeric(0)
  numeric_scale_values <- numeric(0)
  if (length(numeric_predictors) > 0) {
    for (current_predictor in numeric_predictors) {
      train_values <- replace_infinite_with_na(
        train_predictors[[current_predictor]]
      )
      test_values <- replace_infinite_with_na(
        test_predictors[[current_predictor]]
      )
      if (all(is.na(train_values))) {
        stop(
          paste0(
            "Training verisinde predictor tamamen NA: ",
            current_predictor
          ),
          call. = FALSE
        )
      }
      training_median <- median(
        train_values,
        na.rm = TRUE
      )
      numeric_imputation_values[
        current_predictor
      ] <- training_median
      train_values[
        is.na(train_values)
      ] <- training_median
      test_values[
        is.na(test_values)
      ] <- training_median
      if (scale_numeric) {
        training_center <- mean(
          train_values
        )
        training_scale <- sd(
          train_values
        )
        if (
          is.na(training_scale) ||
            training_scale == 0
        ) {
          training_scale <- 1
        }
        numeric_center_values[
          current_predictor
        ] <- training_center
        numeric_scale_values[
          current_predictor
        ] <- training_scale
        train_values <- (
          train_values -
            training_center
        ) / training_scale
        test_values <- (
          test_values -
            training_center
        ) / training_scale
      }
      train_predictors[[current_predictor]] <- train_values
      test_predictors[[current_predictor]] <- test_values
    }
  }
  # ----------------------------------------------------------
  #
  # ----------------------------------------------------------
  factor_level_map <- list()
  if (length(categorical_predictors) > 0) {
    for (
      current_predictor in
      categorical_predictors
    ) {
      train_values <- clean_factor_values(train_predictors[[current_predictor]])
      test_values <- clean_factor_values(test_predictors[[current_predictor]])
      training_levels <- sort(
        unique(train_values)
      )
      training_levels <- unique(
        c(
          training_levels,
          unknown_factor_level
        )
      )
      unseen_test_levels <- setdiff(
        unique(test_values),
        training_levels
      )
      if (length(unseen_test_levels) > 0) {
        test_values[
          test_values %in%
            unseen_test_levels
        ] <- unknown_factor_level
      }
      train_predictors[[current_predictor]] <- factor(train_values, levels = training_levels)
      test_predictors[[current_predictor]] <- factor(test_values, levels = training_levels)
      factor_level_map[[current_predictor]] <- training_levels
    }
  }
  # ----------------------------------------------------------
  #
  # ----------------------------------------------------------
  combined_predictors <- rbind(
    train_predictors,
    test_predictors
  )
  combined_matrix <- model.matrix(
    object = ~ . - 1,
    data = combined_predictors
  )
  training_row_n <- nrow(
    train_predictors
  )
  train_matrix <- combined_matrix[
    seq_len(training_row_n),
    ,
    drop = FALSE
  ]
  test_matrix <- combined_matrix[
    training_row_n +
      seq_len(
        nrow(test_predictors)
      ),
    ,
    drop = FALSE
  ]
  # ----------------------------------------------------------
  #
  # ----------------------------------------------------------
  removed_zero_variance_columns <- character(0)
  if (
    remove_zero_variance &&
      ncol(train_matrix) > 0
  ) {
    zero_variance_flags <- vapply(
      seq_len(
        ncol(train_matrix)
      ),
      function(column_i) {
        length(
          unique(
            train_matrix[
              ,
              column_i
            ]
          )
        ) <= 1
      },
      logical(1)
    )
    removed_zero_variance_columns <-
      colnames(train_matrix)[
        zero_variance_flags
      ]
    if (any(zero_variance_flags)) {
      train_matrix <- train_matrix[
        ,
        !zero_variance_flags,
        drop = FALSE
      ]
      test_matrix <- test_matrix[
        ,
        !zero_variance_flags,
        drop = FALSE
      ]
    }
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (ncol(train_matrix) == 0) {
    stop(
      paste(
        "Model matrisinde kullanilabilir",
        "predictor kalmadi."
      ),
      call. = FALSE
    )
  }
  if (!identical(
    colnames(train_matrix),
    colnames(test_matrix)
  )) {
    stop(
      paste(
        "Training ve test model matrislerinin",
        "columns are not identical."
      ),
      call. = FALSE
    )
  }
  train_na_n <- sum(
    is.na(train_matrix)
  )
  test_na_n <- sum(
    is.na(test_matrix)
  )
  train_nonfinite_n <- sum(
    !is.finite(train_matrix)
  )
  test_nonfinite_n <- sum(
    !is.finite(test_matrix)
  )
  if (
    train_na_n > 0 ||
      test_na_n > 0
  ) {
    stop(
      paste0(
        "NA values were found in the model matrix. Training: ",
        train_na_n,
        "; test: ",
        test_na_n,
        "."
      ),
      call. = FALSE
    )
  }
  if (
    train_nonfinite_n > 0 ||
      test_nonfinite_n > 0
  ) {
    stop(
      paste0(
        "Non-finite values were found in the model matrix. ",
        "Training: ",
        train_nonfinite_n,
        "; test: ",
        test_nonfinite_n,
        "."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  y_train_numeric <- as.integer(
    training_data$Class
  )
  y_train_factor <- factor(
    ifelse(
      y_train_numeric == 1,
      "Occurrence",
      "PseudoBackground"
    ),
    levels = c(
      "PseudoBackground",
      "Occurrence"
    )
  )
  y_test_numeric <- as.integer(
    test_data$Class
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  run_id_value <- if (
    "Run_ID" %in% names(training_data)
  ) {
    unique(training_data$Run_ID)
  } else {
    NA_character_
  }
  if (length(run_id_value) != 1) {
    run_id_value <- NA_character_
  }
  audit <- data.frame(
    Run_ID = run_id_value,
    Training_Row_n = nrow(train_matrix),
    Test_Row_n = nrow(test_matrix),
    Original_Predictor_n =
      length(predictor_columns),
    Numeric_Predictor_n =
      length(numeric_predictors),
    Categorical_Predictor_n =
      length(categorical_predictors),
    Encoded_Column_n =
      ncol(combined_matrix),
    Removed_Zero_Variance_n =
      length(
        removed_zero_variance_columns
      ),
    Final_Model_Column_n =
      ncol(train_matrix),
    Training_NA_n =
      train_na_n,
    Test_NA_n =
      test_na_n,
    Training_Nonfinite_n =
      train_nonfinite_n,
    Test_Nonfinite_n =
      test_nonfinite_n,
    Scaling_Applied =
      scale_numeric,
    Training_Test_Columns_Identical =
      identical(
        colnames(train_matrix),
        colnames(test_matrix)
      ),
    QAQC_Status = "PASS",
    stringsAsFactors = FALSE
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  result <- list(
    x_train = train_matrix,
    x_test = test_matrix,
    y_train_numeric =
      y_train_numeric,
    y_train_factor =
      y_train_factor,
    y_test_numeric =
      y_test_numeric,
    training_cellids =
      training_data$CellID,
    test_cellids =
      test_data$CellID,
    predictor_columns =
      predictor_columns,
    numeric_predictors =
      numeric_predictors,
    categorical_predictors =
      categorical_predictors,
    final_matrix_columns =
      colnames(train_matrix),
    removed_zero_variance_columns =
      removed_zero_variance_columns,
    preprocessing = list(
      numeric_imputation_values =
        numeric_imputation_values,
      numeric_center_values =
        numeric_center_values,
      numeric_scale_values =
        numeric_scale_values,
      factor_level_map =
        factor_level_map,
      scaling_applied =
        scale_numeric,
      unknown_factor_level =
        unknown_factor_level
    ),
    audit = audit
  )
  return(result)
}
# ============================================================
# MODUL SONU
# ============================================================
