# ============================================================
# Evaluate_RF.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Summarize Random Forest held-out validation performance and permutation importance.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

results_directory <- paste0(
  "Project/Models/",
  "RF_Results"
)
run_results_file <- paste0(
  results_directory,
  "/RF_LOOCV_Run_Results.csv"
)
importance_results_file <- paste0(
  results_directory,
  "/RF_LOOCV_Variable_Importance.csv"
)
evaluation_directory <- paste0(
  results_directory,
  "/Evaluation"
)
dir.create(
  evaluation_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
overall_summary_file <- paste0(
  evaluation_directory,
  "/RF_Evaluation_Overall.csv"
)
by_occurrence_file <- paste0(
  evaluation_directory,
  "/RF_Evaluation_By_Occurrence.csv"
)
by_fold_file <- paste0(
  evaluation_directory,
  "/RF_Evaluation_By_Fold.csv"
)
by_repeat_file <- paste0(
  evaluation_directory,
  "/RF_Evaluation_By_Repeat.csv"
)
importance_encoded_file <- paste0(
  evaluation_directory,
  "/RF_Importance_Encoded_Summary.csv"
)
importance_grouped_file <- paste0(
  evaluation_directory,
  "/RF_Importance_Grouped_Summary.csv"
)
qaqc_file <- paste0(
  "Project/QAQC/",
  "RF_Evaluation_QAQC.txt"
)
methodology_file <- paste0(
  "Project/QAQC/",
  "RF_Evaluation_Methodology.txt"
)
# ------------------------------------------------------------
# 2. Plot paths
# ------------------------------------------------------------
histogram_file <- paste0(
  evaluation_directory,
  "/RF_HeldOut_Probability_Histogram.png"
)
fold_boxplot_file <- paste0(
  evaluation_directory,
  "/RF_HeldOut_Probability_By_Fold.png"
)
occurrence_ci_file <- paste0(
  evaluation_directory,
  "/RF_HeldOut_Probability_By_Occurrence_CI95.png"
)
importance_plot_file <- paste0(
  evaluation_directory,
  "/RF_Grouped_Variable_Importance.png"
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
safe_sd <- function(x) {
  x <- x[
    is.finite(x)
  ]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x)
}
safe_se <- function(x) {
  x <- x[
    is.finite(x)
  ]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}
summarise_probability <- function(data) {
  probability <- data$
    Held_Out_Occurrence_Probability
  probability <- probability[
    is.finite(probability)
  ]
  n_value <- length(probability)
  mean_value <- mean(probability)
  sd_value <- safe_sd(probability)
  se_value <- safe_se(probability)
  if (
    n_value > 1 &&
      is.finite(se_value)
  ) {
    t_critical <- stats::qt(
      0.975,
      df = n_value - 1
    )
    ci_lower <- mean_value -
      t_critical * se_value
    ci_upper <- mean_value +
      t_critical * se_value
  } else {
    ci_lower <- NA_real_
    ci_upper <- NA_real_
  }
  data.frame(
    N = n_value,
    Mean_Probability = mean_value,
    Median_Probability =
      stats::median(probability),
    SD_Probability = sd_value,
    SE_Probability = se_value,
    Minimum_Probability = min(probability),
    Q25_Probability =
      as.numeric(
        stats::quantile(
          probability,
          probs = 0.25,
          names = FALSE
        )
      ),
    Q75_Probability =
      as.numeric(
        stats::quantile(
          probability,
          probs = 0.75,
          names = FALSE
        )
      ),
    Maximum_Probability = max(probability),
    CI95_Lower = max(
      0,
      ci_lower,
      na.rm = TRUE
    ),
    CI95_Upper = min(
      1,
      ci_upper,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}
group_predictor_name <- function(x) {
  x <- as.character(x)
  x[
    grepl(
      "^Lithology",
      x
    )
  ] <- "Lithology"
  x
}
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  run_results_file,
  importance_results_file
)
missing_files <- input_files[
  !file.exists(input_files)
]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following input files were not found:\n",
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
importance_results <- read.csv(
  importance_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
assert_required_columns(
  run_results,
  c(
    "Algorithm",
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "Held_Out_CellID",
    "Held_Out_X",
    "Held_Out_Y",
    "Held_Out_Occurrence_Probability",
    "Best_Tuning_OOB_Error",
    "Final_Model_OOB_Error",
    "QAQC_Status"
  ),
  "RF run-results table"
)
assert_required_columns(
  importance_results,
  c(
    "Algorithm",
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Predictor",
    "Importance"
  ),
  "RF variable-importance table"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
expected_run_n <- 150L
run_n <- nrow(
  run_results
)
unique_run_n <- length(
  unique(
    run_results$Run_ID
  )
)
duplicate_run_n <- sum(
  duplicated(
    run_results$Run_ID
  )
)
failed_run_n <- sum(
  run_results$QAQC_Status != "PASS"
)
invalid_probability_n <- sum(
  !is.finite(
    run_results$
      Held_Out_Occurrence_Probability
  ) |
    run_results$
      Held_Out_Occurrence_Probability < 0 |
    run_results$
      Held_Out_Occurrence_Probability > 1
)
invalid_oob_n <- sum(
  !is.finite(
    run_results$
      Final_Model_OOB_Error
  ) |
    run_results$
      Final_Model_OOB_Error < 0 |
    run_results$
      Final_Model_OOB_Error > 1
)
missing_importance_n <- sum(
  !is.finite(
    importance_results$Importance
  )
)
run_id_without_importance <- setdiff(
  run_results$Run_ID,
  unique(
    importance_results$Run_ID
  )
)
importance_unknown_run_id <- setdiff(
  unique(
    importance_results$Run_ID
  ),
  run_results$Run_ID
)
if (
  invalid_probability_n > 0 ||
    failed_run_n > 0
) {
  stop(
    paste(
      "RF evaluation stopped because run-results QA/QC",
      "contains failed runs or invalid probabilities."
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
overall_probability_summary <-
  summarise_probability(
    run_results
  )
overall_summary <- data.frame(
  Algorithm = "RF",
  Expected_Run_n = expected_run_n,
  Successful_Run_n = run_n,
  Unique_Run_n = unique_run_n,
  Unique_Repeat_n =
    length(
      unique(
        run_results$Repeat
      )
    ),
  Unique_Fold_n =
    length(
      unique(
        run_results$Fold_ID
      )
    ),
  Unique_Held_Out_Occurrence_n =
    length(
      unique(
        run_results$Held_Out_CellID
      )
    ),
  overall_probability_summary,
  Mean_Best_Tuning_OOB_Error =
    mean(
      run_results$
        Best_Tuning_OOB_Error
    ),
  SD_Best_Tuning_OOB_Error =
    safe_sd(
      run_results$
        Best_Tuning_OOB_Error
    ),
  Mean_Final_Model_OOB_Error =
    mean(
      run_results$
        Final_Model_OOB_Error
    ),
  SD_Final_Model_OOB_Error =
    safe_sd(
      run_results$
        Final_Model_OOB_Error
    ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrence_keys <- unique(
  run_results[
    ,
    c(
      "Held_Out_CellID",
      "Held_Out_X",
      "Held_Out_Y"
    ),
    drop = FALSE
  ]
)
occurrence_summary_list <- vector(
  "list",
  nrow(occurrence_keys)
)
for (
  occurrence_i in
  seq_len(nrow(occurrence_keys))
) {
  current_cell_id <-
    occurrence_keys$Held_Out_CellID[
      occurrence_i
    ]
  current_data <- run_results[
    run_results$Held_Out_CellID ==
      current_cell_id,
    ,
    drop = FALSE
  ]
  current_summary <-
    summarise_probability(
      current_data
    )
  occurrence_summary_list[[occurrence_i]] <- data.frame(
    Algorithm = "RF",
    Held_Out_CellID =
      current_cell_id,
    Held_Out_X =
      occurrence_keys$Held_Out_X[
        occurrence_i
      ],
    Held_Out_Y =
      occurrence_keys$Held_Out_Y[
        occurrence_i
      ],
    Fold_ID = paste(
      sort(
        unique(
          current_data$Fold_ID
        )
      ),
      collapse = ";"
    ),
    current_summary,
    Mean_Final_Model_OOB_Error =
      mean(
        current_data$
          Final_Model_OOB_Error
      ),
    SD_Final_Model_OOB_Error =
      safe_sd(
        current_data$
          Final_Model_OOB_Error
      ),
    stringsAsFactors = FALSE
  )
}
by_occurrence <- do.call(
  rbind,
  occurrence_summary_list
)
by_occurrence <- by_occurrence[
  order(
    by_occurrence$Held_Out_CellID
  ),
  ,
  drop = FALSE
]
rownames(by_occurrence) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
fold_ids <- sort(
  unique(
    run_results$Fold_ID
  )
)
fold_summary_list <- vector(
  "list",
  length(fold_ids)
)
for (
  fold_i in
  seq_along(fold_ids)
) {
  current_fold_id <-
    fold_ids[fold_i]
  current_data <- run_results[
    run_results$Fold_ID ==
      current_fold_id,
    ,
    drop = FALSE
  ]
  current_summary <-
    summarise_probability(
      current_data
    )
  fold_summary_list[[fold_i]] <- data.frame(
    Algorithm = "RF",
    Fold_ID =
      current_fold_id,
    Fold_Number =
      unique(
        current_data$Fold_Number
      )[1],
    Held_Out_CellID =
      unique(
        current_data$Held_Out_CellID
      )[1],
    current_summary,
    Mean_Final_Model_OOB_Error =
      mean(
        current_data$
          Final_Model_OOB_Error
      ),
    SD_Final_Model_OOB_Error =
      safe_sd(
        current_data$
          Final_Model_OOB_Error
      ),
    stringsAsFactors = FALSE
  )
}
by_fold <- do.call(
  rbind,
  fold_summary_list
)
by_fold <- by_fold[
  order(
    by_fold$Fold_Number
  ),
  ,
  drop = FALSE
]
rownames(by_fold) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
repeat_values <- sort(
  unique(
    run_results$Repeat
  )
)
repeat_summary_list <- vector(
  "list",
  length(repeat_values)
)
for (
  repeat_i in
  seq_along(repeat_values)
) {
  current_repeat <-
    repeat_values[repeat_i]
  current_data <- run_results[
    run_results$Repeat ==
      current_repeat,
    ,
    drop = FALSE
  ]
  current_summary <-
    summarise_probability(
      current_data
    )
  repeat_summary_list[[repeat_i]] <- data.frame(
    Algorithm = "RF",
    Repeat =
      current_repeat,
    current_summary,
    Mean_Final_Model_OOB_Error =
      mean(
        current_data$
          Final_Model_OOB_Error
      ),
    SD_Final_Model_OOB_Error =
      safe_sd(
        current_data$
          Final_Model_OOB_Error
      ),
    stringsAsFactors = FALSE
  )
}
by_repeat <- do.call(
  rbind,
  repeat_summary_list
)
rownames(by_repeat) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
encoded_predictors <- sort(
  unique(
    importance_results$Predictor
  )
)
encoded_importance_list <- vector(
  "list",
  length(encoded_predictors)
)
for (
  predictor_i in
  seq_along(encoded_predictors)
) {
  current_predictor <-
    encoded_predictors[predictor_i]
  current_values <- importance_results$Importance[
    importance_results$Predictor ==
      current_predictor
  ]
  encoded_importance_list[[predictor_i]] <- data.frame(
    Algorithm = "RF",
    Predictor =
      current_predictor,
    Run_n =
      length(current_values),
    Mean_Importance =
      mean(current_values),
    Median_Importance =
      stats::median(current_values),
    SD_Importance =
      safe_sd(current_values),
    Minimum_Importance =
      min(current_values),
    Maximum_Importance =
      max(current_values),
    stringsAsFactors = FALSE
  )
}
importance_encoded <- do.call(
  rbind,
  encoded_importance_list
)
importance_encoded <- importance_encoded[
  order(
    -importance_encoded$Mean_Importance
  ),
  ,
  drop = FALSE
]
importance_encoded$Rank <-
  seq_len(
    nrow(importance_encoded)
  )
importance_encoded <- importance_encoded[
  ,
  c(
    "Algorithm",
    "Rank",
    "Predictor",
    "Run_n",
    "Mean_Importance",
    "Median_Importance",
    "SD_Importance",
    "Minimum_Importance",
    "Maximum_Importance"
  ),
  drop = FALSE
]
rownames(importance_encoded) <- NULL
# ------------------------------------------------------------
#
# ------------------------------------------------------------
importance_results$Original_Predictor <-
  group_predictor_name(
    importance_results$Predictor
  )
grouped_by_run <- stats::aggregate(
  Importance ~
    Algorithm +
    Run_ID +
    Repeat +
    Fold_ID +
    Original_Predictor,
  data = importance_results,
  FUN = sum
)
grouped_predictors <- sort(
  unique(
    grouped_by_run$Original_Predictor
  )
)
grouped_importance_list <- vector(
  "list",
  length(grouped_predictors)
)
for (
  predictor_i in
  seq_along(grouped_predictors)
) {
  current_predictor <-
    grouped_predictors[predictor_i]
  current_values <- grouped_by_run$Importance[
    grouped_by_run$Original_Predictor ==
      current_predictor
  ]
  grouped_importance_list[[predictor_i]] <- data.frame(
    Algorithm = "RF",
    Predictor =
      current_predictor,
    Run_n =
      length(current_values),
    Mean_Importance =
      mean(current_values),
    Median_Importance =
      stats::median(current_values),
    SD_Importance =
      safe_sd(current_values),
    Minimum_Importance =
      min(current_values),
    Maximum_Importance =
      max(current_values),
    stringsAsFactors = FALSE
  )
}
importance_grouped <- do.call(
  rbind,
  grouped_importance_list
)
importance_grouped <- importance_grouped[
  order(
    -importance_grouped$Mean_Importance
  ),
  ,
  drop = FALSE
]
importance_grouped$Rank <-
  seq_len(
    nrow(importance_grouped)
  )
importance_grouped <- importance_grouped[
  ,
  c(
    "Algorithm",
    "Rank",
    "Predictor",
    "Run_n",
    "Mean_Importance",
    "Median_Importance",
    "SD_Importance",
    "Minimum_Importance",
    "Maximum_Importance"
  ),
  drop = FALSE
]
rownames(importance_grouped) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  overall_summary,
  overall_summary_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  by_occurrence,
  by_occurrence_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  by_fold,
  by_fold_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  by_repeat,
  by_repeat_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  importance_encoded,
  importance_encoded_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  importance_grouped,
  importance_grouped_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
grDevices::png(
  filename = histogram_file,
  width = 1800,
  height = 1300,
  res = 200
)
graphics::hist(
  run_results$
    Held_Out_Occurrence_Probability,
  breaks = 20,
  main =
    "RF Held-Out Occurrence Probabilities",
  xlab =
    "Held-out occurrence probability",
  ylab =
    "Number of LOOCV runs",
  xlim = c(
    0,
    1
  )
)
graphics::abline(
  v = mean(
    run_results$
      Held_Out_Occurrence_Probability
  ),
  lwd = 2,
  lty = 2
)
graphics::abline(
  v = stats::median(
    run_results$
      Held_Out_Occurrence_Probability
  ),
  lwd = 2,
  lty = 3
)
graphics::legend(
  "topright",
  legend = c(
    "Mean",
    "Median"
  ),
  lty = c(
    2,
    3
  ),
  lwd = 2,
  bty = "n"
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
grDevices::png(
  filename = fold_boxplot_file,
  width = 1900,
  height = 1300,
  res = 200
)
graphics::boxplot(
  Held_Out_Occurrence_Probability ~
    Fold_ID,
  data = run_results,
  main =
    "RF Held-Out Occurrence Probability by Fold",
  xlab =
    "LOOCV fold",
  ylab =
    "Held-out occurrence probability",
  ylim = c(
    0,
    1
  ),
  las = 2
)
graphics::stripchart(
  Held_Out_Occurrence_Probability ~
    Fold_ID,
  data = run_results,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  cex = 0.55,
  add = TRUE
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
plot_order <- order(
  by_occurrence$Mean_Probability,
  decreasing = FALSE
)
plot_data <- by_occurrence[
  plot_order,
  ,
  drop = FALSE
]
grDevices::png(
  filename = occurrence_ci_file,
  width = 1900,
  height = 1300,
  res = 200
)
graphics::plot(
  x = plot_data$Mean_Probability,
  y = seq_len(
    nrow(plot_data)
  ),
  xlim = c(
    0,
    1
  ),
  ylim = c(
    0.5,
    nrow(plot_data) + 0.5
  ),
  pch = 16,
  yaxt = "n",
  xlab =
    "Mean held-out occurrence probability",
  ylab =
    "Held-out occurrence CellID",
  main =
    "RF Held-Out Occurrence Probability with 95% CI"
)
graphics::axis(
  side = 2,
  at = seq_len(
    nrow(plot_data)
  ),
  labels =
    plot_data$Held_Out_CellID,
  las = 2
)
graphics::segments(
  x0 =
    plot_data$CI95_Lower,
  y0 =
    seq_len(
      nrow(plot_data)
    ),
  x1 =
    plot_data$CI95_Upper,
  y1 =
    seq_len(
      nrow(plot_data)
    )
)
graphics::points(
  x =
    plot_data$Mean_Probability,
  y =
    seq_len(
      nrow(plot_data)
    ),
  pch = 16
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
predictor_order <- importance_grouped$Predictor[
  order(
    importance_grouped$Mean_Importance,
    decreasing = TRUE
  )
]
grouped_by_run$Original_Predictor <- factor(
  grouped_by_run$Original_Predictor,
  levels = rev(predictor_order)
)
grDevices::png(
  filename = importance_plot_file,
  width = 1900,
  height = 1400,
  res = 200
)
graphics::par(
  mar = c(5, 10, 4, 2) + 0.1
)
graphics::boxplot(
  Importance ~ Original_Predictor,
  data = grouped_by_run,
  horizontal = TRUE,
  las = 1,
  outline = TRUE,
  main = "RF Grouped Permutation Importance",
  xlab = "Permutation importance",
  ylab = ""
)
graphics::abline(
  v = 0,
  lty = 2,
  lwd = 1
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrence_repeat_count_ok <- all(
  by_occurrence$N == 30
)
fold_repeat_count_ok <- all(
  by_fold$N == 30
)
repeat_fold_count_ok <- all(
  by_repeat$N == 5
)
grouped_predictor_count_ok <- (
  nrow(importance_grouped) == 7
)
all_checks_passed <- all(
  run_n == expected_run_n,
  unique_run_n == expected_run_n,
  duplicate_run_n == 0,
  failed_run_n == 0,
  invalid_probability_n == 0,
  invalid_oob_n == 0,
  missing_importance_n == 0,
  length(run_id_without_importance) == 0,
  length(importance_unknown_run_id) == 0,
  occurrence_repeat_count_ok,
  fold_repeat_count_ok,
  repeat_fold_count_ok,
  grouped_predictor_count_ok
)
overall_qaqc_status <- ifelse(
  all_checks_passed,
  "PASS",
  "FAIL"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "RANDOM FOREST LOOCV EVALUATION METHODOLOGY",
  "",
  paste0(
    "Evaluation runs: ",
    run_n
  ),
  paste0(
    "Repeats: ",
    length(
      unique(
        run_results$Repeat
      )
    )
  ),
  paste0(
    "LOOCV folds: ",
    length(
      unique(
        run_results$Fold_ID
      )
    )
  ),
  "",
  paste(
    "Evaluation is based on probabilities assigned to",
    "held-out known occurrences."
  ),
  paste(
    "Each known occurrence was held out once per repeat,",
    "producing 30 independent held-out predictions for",
    "each occurrence."
  ),
  paste(
    "Overall, repeat-level, fold-level and occurrence-level",
    "probability summaries were calculated."
  ),
  paste(
    "Occurrence-level 95% confidence intervals use the",
    "Student t distribution across the 30 repeats."
  ),
  paste(
    "Permutation importance was summarized both for encoded",
    "model-matrix columns and for the seven original predictors."
  ),
  paste(
    "Lithology dummy-variable importances were summed within",
    "each run before original-predictor importance was summarized."
  ),
  "",
  paste(
    "ROC-AUC, specificity, precision, recall and F1 were not",
    "calculated because independently verified absence locations",
    "are not available."
  ),
  paste(
    "OOB error is retained as a training diagnostic and must",
    "not be interpreted as independent spatial validation."
  ),
  paste(
    "Full-grid threshold and occurrence-capture evaluation will",
    "be performed after full-grid RF prospectivity prediction."
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
  "RANDOM FOREST LOOCV EVALUATION QAQC",
  "",
  paste0(
    "Expected runs: ",
    expected_run_n
  ),
  paste0(
    "Observed runs: ",
    run_n
  ),
  paste0(
    "Unique Run_ID values: ",
    unique_run_n
  ),
  paste0(
    "Duplicated Run_ID values: ",
    duplicate_run_n
  ),
  paste0(
    "Failed run records: ",
    failed_run_n
  ),
  paste0(
    "Invalid held-out probabilities: ",
    invalid_probability_n
  ),
  paste0(
    "Invalid final OOB errors: ",
    invalid_oob_n
  ),
  paste0(
    "Missing importance values: ",
    missing_importance_n
  ),
  paste0(
    "Run_ID values without importance records: ",
    length(
      run_id_without_importance
    )
  ),
  paste0(
    "Importance records with unknown Run_ID values: ",
    length(
      importance_unknown_run_id
    )
  ),
  paste0(
    "Each occurrence has 30 predictions: ",
    occurrence_repeat_count_ok
  ),
  paste0(
    "Each fold has 30 predictions: ",
    fold_repeat_count_ok
  ),
  paste0(
    "Each repeat has 5 held-out predictions: ",
    repeat_fold_count_ok
  ),
  paste0(
    "Grouped predictor count equals 7: ",
    grouped_predictor_count_ok
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
# 21. Console summary
# ------------------------------------------------------------
cat(
  "\n========================================\n"
)
cat(
  "Random Forest evaluation completed.\n"
)
cat(
  "Runs evaluated:",
  run_n,
  "\n"
)
cat(
  "Mean held-out probability:",
  round(
    overall_summary$Mean_Probability,
    4
  ),
  "\n"
)
cat(
  "Median held-out probability:",
  round(
    overall_summary$Median_Probability,
    4
  ),
  "\n"
)
cat(
  "Top grouped predictor:",
  importance_grouped$Predictor[1],
  "\n"
)
cat(
  "Overall QAQC:",
  overall_qaqc_status,
  "\n"
)
cat(
  "Evaluation directory:",
  evaluation_directory,
  "\n"
)
cat(
  "Overall summary:",
  overall_summary_file,
  "\n"
)
cat(
  "Occurrence summary:",
  by_occurrence_file,
  "\n"
)
cat(
  "Grouped importance:",
  importance_grouped_file,
  "\n"
)
cat(
  "========================================\n"
)
