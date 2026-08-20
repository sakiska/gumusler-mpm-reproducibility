# ============================================================
# Evaluate_SVM.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Summarize SVM held-out validation performance.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

input_directory <- "Project/Models/SVM_Results"
evaluation_directory <- file.path(input_directory, "Evaluation")
dir.create(
  evaluation_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  "Project/QAQC",
  recursive = TRUE,
  showWarnings = FALSE
)
run_results_file <- file.path(
  input_directory,
  "SVM_LOOCV_Run_Results.csv"
)
output_overall_file <- file.path(
  evaluation_directory,
  "SVM_Evaluation_Overall.csv"
)
output_occurrence_file <- file.path(
  evaluation_directory,
  "SVM_Evaluation_By_Occurrence.csv"
)
output_fold_file <- file.path(
  evaluation_directory,
  "SVM_Evaluation_By_Fold.csv"
)
output_repeat_file <- file.path(
  evaluation_directory,
  "SVM_Evaluation_By_Repeat.csv"
)
output_parameter_file <- file.path(
  evaluation_directory,
  "SVM_Hyperparameter_Selection_Summary.csv"
)
histogram_file <- file.path(
  evaluation_directory,
  "SVM_HeldOut_Probability_Histogram.png"
)
fold_boxplot_file <- file.path(
  evaluation_directory,
  "SVM_HeldOut_Probability_By_Fold.png"
)
qaqc_file <- "Project/QAQC/SVM_Evaluation_QAQC.txt"
methodology_file <- paste0(
  "Project/QAQC/",
  "SVM_Evaluation_Methodology.txt"
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
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}
safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}
summarise_probability <- function(data) {
  probability <- data$Held_Out_Occurrence_Probability
  probability <- probability[is.finite(probability)]
  n_value <- length(probability)
  if (n_value == 0L) {
    stop(
      "No finite held-out probabilities were available.",
      call. = FALSE
    )
  }
  mean_value <- mean(probability)
  sd_value <- safe_sd(probability)
  se_value <- safe_se(probability)
  if (
    n_value > 1L &&
      is.finite(se_value)
  ) {
    t_critical <- stats::qt(
      0.975,
      df = n_value - 1L
    )
    ci_lower <- mean_value - t_critical * se_value
    ci_upper <- mean_value + t_critical * se_value
  } else {
    ci_lower <- NA_real_
    ci_upper <- NA_real_
  }
  data.frame(
    Run_n = n_value,
    Mean_Probability = mean_value,
    Median_Probability = stats::median(probability),
    SD_Probability = sd_value,
    SE_Probability = se_value,
    Minimum_Probability = min(probability),
    Q25_Probability = as.numeric(
      stats::quantile(
        probability,
        probs = 0.25,
        names = FALSE
      )
    ),
    Q75_Probability = as.numeric(
      stats::quantile(
        probability,
        probs = 0.75,
        names = FALSE
      )
    ),
    Maximum_Probability = max(probability),
    CI95_Lower = if (is.finite(ci_lower)) max(0, ci_lower) else NA_real_,
    CI95_Upper = if (is.finite(ci_upper)) min(1, ci_upper) else NA_real_,
    stringsAsFactors = FALSE
  )
}
summarise_by_group <- function(
    data,
    group_column,
    output_group_name
) {
  split_data <- split(
    data,
    data[[group_column]],
    drop = TRUE
  )
  summary_list <- lapply(
    names(split_data),
    function(group_value) {
      current_summary <- summarise_probability(
        split_data[[group_value]]
      )
      data.frame(
        Group_Value = group_value,
        current_summary,
        stringsAsFactors = FALSE
      )
    }
  )
  result <- do.call(
    rbind,
    summary_list
  )
  names(result)[1] <- output_group_name
  rownames(result) <- NULL
  result
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
run_results <- read.csv(
  run_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c(
  "Algorithm",
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "Held_Out_CellID",
  "Held_Out_Occurrence_Probability",
  "Cost",
  "Gamma",
  "Inner_nfold",
  "Inner_CV_Logloss",
  "Scaling_Applied",
  "QAQC_Status"
)
assert_required_columns(
  data = run_results,
  required_columns = required_columns,
  data_name = "SVM run results"
)
if (nrow(run_results) != 150L) {
  stop(
    paste0(
      "Expected 150 SVM runs; found ",
      nrow(run_results),
      "."
    ),
    call. = FALSE
  )
}
if (any(run_results$QAQC_Status != "PASS")) {
  stop(
    "SVM run results contain failed QAQC records.",
    call. = FALSE
  )
}
if (any(!run_results$Scaling_Applied)) {
  stop(
    "At least one SVM run was completed without numeric scaling.",
    call. = FALSE
  )
}
probability <- run_results$Held_Out_Occurrence_Probability
invalid_probability <- (
  !is.finite(probability) |
    probability < 0 |
    probability > 1
)
if (any(invalid_probability)) {
  stop(
    "Invalid held-out occurrence probabilities were detected.",
    call. = FALSE
  )
}
if (
  any(!is.finite(run_results$Cost)) ||
    any(run_results$Cost <= 0)
) {
  stop(
    "Invalid selected Cost values were detected.",
    call. = FALSE
  )
}
if (
  any(!is.finite(run_results$Gamma)) ||
    any(run_results$Gamma <= 0)
) {
  stop(
    "Invalid selected Gamma values were detected.",
    call. = FALSE
  )
}
if (
  any(!is.finite(run_results$Inner_CV_Logloss)) ||
    any(run_results$Inner_CV_Logloss < 0)
) {
  stop(
    "Invalid inner-CV logloss values were detected.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
overall_probability <- summarise_probability(
  run_results
)
overall <- data.frame(
  Algorithm = "SVM",
  overall_probability,
  Mean_Inner_CV_Logloss = mean(
    run_results$Inner_CV_Logloss
  ),
  SD_Inner_CV_Logloss = safe_sd(
    run_results$Inner_CV_Logloss
  ),
  Median_Inner_CV_Logloss = stats::median(
    run_results$Inner_CV_Logloss
  ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
by_occurrence <- summarise_by_group(
  data = run_results,
  group_column = "Held_Out_CellID",
  output_group_name = "Held_Out_CellID"
)
by_fold <- summarise_by_group(
  data = run_results,
  group_column = "Fold_ID",
  output_group_name = "Fold_ID"
)
by_repeat <- summarise_by_group(
  data = run_results,
  group_column = "Repeat",
  output_group_name = "Repeat"
)
fold_order <- match(
  by_fold$Fold_ID,
  sprintf("Fold_%02d", seq_len(5))
)
by_fold <- by_fold[
  order(fold_order),
  ,
  drop = FALSE
]
repeat_numeric <- suppressWarnings(
  as.integer(by_repeat$Repeat)
)
by_repeat <- by_repeat[
  order(repeat_numeric),
  ,
  drop = FALSE
]
rownames(by_fold) <- NULL
rownames(by_repeat) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
parameter_counts <- stats::aggregate(
  Run_ID ~ Cost + Gamma,
  data = run_results,
  FUN = length
)
names(parameter_counts)[
  names(parameter_counts) == "Run_ID"
] <- "Selected_Run_n"
parameter_counts$Selected_Percent <- (
  100 * parameter_counts$Selected_Run_n /
    nrow(run_results)
)
parameter_logloss <- stats::aggregate(
  Inner_CV_Logloss ~ Cost + Gamma,
  data = run_results,
  FUN = mean
)
names(parameter_logloss)[
  names(parameter_logloss) == "Inner_CV_Logloss"
] <- "Mean_Inner_CV_Logloss"
parameter_probability <- stats::aggregate(
  Held_Out_Occurrence_Probability ~ Cost + Gamma,
  data = run_results,
  FUN = mean
)
names(parameter_probability)[
  names(parameter_probability) ==
    "Held_Out_Occurrence_Probability"
] <- "Mean_HeldOut_Probability"
parameter_summary <- merge(
  parameter_counts,
  parameter_logloss,
  by = c("Cost", "Gamma"),
  all = TRUE
)
parameter_summary <- merge(
  parameter_summary,
  parameter_probability,
  by = c("Cost", "Gamma"),
  all = TRUE
)
parameter_summary <- parameter_summary[
  order(
    -parameter_summary$Selected_Run_n,
    parameter_summary$Cost,
    parameter_summary$Gamma
  ),
  ,
  drop = FALSE
]
parameter_summary$Selection_Rank <- seq_len(
  nrow(parameter_summary)
)
parameter_summary <- parameter_summary[
  ,
  c(
    "Selection_Rank",
    "Cost",
    "Gamma",
    "Selected_Run_n",
    "Selected_Percent",
    "Mean_Inner_CV_Logloss",
    "Mean_HeldOut_Probability"
  ),
  drop = FALSE
]
rownames(parameter_summary) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  overall,
  output_overall_file,
  row.names = FALSE
)
write.csv(
  by_occurrence,
  output_occurrence_file,
  row.names = FALSE
)
write.csv(
  by_fold,
  output_fold_file,
  row.names = FALSE
)
write.csv(
  by_repeat,
  output_repeat_file,
  row.names = FALSE
)
write.csv(
  parameter_summary,
  output_parameter_file,
  row.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
grDevices::png(
  histogram_file,
  width = 1600,
  height = 1100,
  res = 180
)
graphics::hist(
  probability,
  breaks = 15,
  main = "SVM held-out occurrence probabilities",
  xlab = "Held-out occurrence probability"
)
graphics::abline(
  v = mean(probability),
  lty = 2,
  lwd = 2
)
grDevices::dev.off()
grDevices::png(
  fold_boxplot_file,
  width = 1600,
  height = 1100,
  res = 180
)
graphics::boxplot(
  Held_Out_Occurrence_Probability ~ Fold_ID,
  data = run_results,
  main = "SVM held-out probabilities by fold",
  xlab = "Fold",
  ylab = "Held-out occurrence probability"
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
overall_qaqc <- if (
  nrow(run_results) == 150L &&
    sum(invalid_probability) == 0L &&
    all(run_results$QAQC_Status == "PASS") &&
    all(run_results$Scaling_Applied) &&
    nrow(parameter_summary) > 0L
) {
  "PASS"
} else {
  "FAIL"
}
writeLines(
  c(
    "SVM EVALUATION QAQC",
    "",
    paste0(
      "Runs evaluated: ",
      nrow(run_results)
    ),
    paste0(
      "Invalid probabilities: ",
      sum(invalid_probability)
    ),
    paste0(
      "Runs without scaling: ",
      sum(!run_results$Scaling_Applied)
    ),
    paste0(
      "Selected Cost-Gamma combinations: ",
      nrow(parameter_summary)
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
    "SVM EVALUATION METHODOLOGY",
    "",
    paste(
      "Final evaluation summaries use only genuinely",
      "held-out occurrence probabilities from the",
      "leakage-free repeated LOOCV procedure."
    ),
    paste(
      "Numeric predictors were standardized using",
      "parameters learned exclusively from each",
      "outer-training set."
    ),
    paste(
      "Cost and gamma were selected independently",
      "within each outer run using stratified inner",
      "cross-validation and binary logloss."
    ),
    paste(
      "ROC, precision, recall, specificity and F1",
      "were not calculated at this stage because",
      "pseudo-background cells are not verified absences."
    ),
    paste(
      "Radial SVM does not provide a native variable",
      "importance measure directly comparable with",
      "tree-based RF and XGBoost importance."
    )
  ),
  methodology_file
)
# ------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------
top_parameter <- parameter_summary[1, , drop = FALSE]
cat("\n========================================\n")
cat("SVM evaluation completed.\n")
cat("Runs evaluated:", nrow(run_results), "\n")
cat(
  "Mean held-out probability:",
  round(overall$Mean_Probability, 4),
  "\n"
)
cat(
  "Median held-out probability:",
  round(overall$Median_Probability, 4),
  "\n"
)
cat(
  "Most selected Cost:",
  top_parameter$Cost,
  "\n"
)
cat(
  "Most selected Gamma:",
  top_parameter$Gamma,
  "\n"
)
cat(
  "Selection count:",
  top_parameter$Selected_Run_n,
  "of",
  nrow(run_results),
  "\n"
)
cat("Overall QAQC:", overall_qaqc, "\n")
cat(
  "Evaluation directory:",
  evaluation_directory,
  "\n"
)
cat("========================================\n")
if (overall_qaqc != "PASS") {
  stop(
    "SVM evaluation QAQC failed.",
    call. = FALSE
  )
}
