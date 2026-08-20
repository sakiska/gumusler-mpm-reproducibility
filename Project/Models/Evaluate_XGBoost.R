# ============================================================
# Evaluate_XGBoost.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Summarize XGBoost held-out validation performance.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

input_directory <- "Project/Models/XGB_Results"
evaluation_directory <- file.path(input_directory, "Evaluation")
dir.create(evaluation_directory, recursive = TRUE, showWarnings = FALSE)
dir.create("Project/QAQC", recursive = TRUE, showWarnings = FALSE)
run_results_file <- file.path(input_directory, "XGB_LOOCV_Run_Results.csv")
importance_file <- file.path(input_directory, "XGB_LOOCV_Variable_Importance.csv")
if (!file.exists(run_results_file)) stop(paste0("Missing file: ", run_results_file), call. = FALSE)
if (!file.exists(importance_file)) stop(paste0("Missing file: ", importance_file), call. = FALSE)
run_results <- read.csv(run_results_file, stringsAsFactors = FALSE, check.names = FALSE)
importance <- read.csv(importance_file, stringsAsFactors = FALSE, check.names = FALSE)
required_run_columns <- c("Run_ID", "Repeat", "Fold_ID", "Fold_Number", "Held_Out_CellID", "Held_Out_Occurrence_Probability", "Inner_CV_Logloss", "QAQC_Status")
missing_run_columns <- setdiff(required_run_columns, names(run_results))
if (length(missing_run_columns) > 0) stop(paste0("Run results missing columns: ", paste(missing_run_columns, collapse = ", ")), call. = FALSE)
if (nrow(run_results) != 150L) stop(paste0("Expected 150 runs; found ", nrow(run_results), "."), call. = FALSE)
if (any(run_results$QAQC_Status != "PASS")) stop("Run results contain failed QAQC records.", call. = FALSE)
probability <- run_results$Held_Out_Occurrence_Probability
if (any(!is.finite(probability)) || any(probability < 0) || any(probability > 1)) stop("Invalid held-out probabilities.", call. = FALSE)
overall <- data.frame(
  Algorithm = "XGBoost",
  Run_n = nrow(run_results),
  Mean_HeldOut_Probability = mean(probability),
  SD_HeldOut_Probability = sd(probability),
  Median_HeldOut_Probability = median(probability),
  Minimum_HeldOut_Probability = min(probability),
  Maximum_HeldOut_Probability = max(probability),
  Mean_Inner_CV_Logloss = mean(run_results$Inner_CV_Logloss),
  SD_Inner_CV_Logloss = sd(run_results$Inner_CV_Logloss),
  stringsAsFactors = FALSE
)
summarize_group <- function(data, group_column) {
  groups <- unique(data[[group_column]])
  output <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    group_value <- groups[i]
    subset_data <- data[data[[group_column]] == group_value, , drop = FALSE]
    values <- subset_data$Held_Out_Occurrence_Probability
    output[[i]] <- data.frame(
      Group = as.character(group_value),
      Run_n = nrow(subset_data),
      Mean_Probability = mean(values),
      SD_Probability = sd(values),
      Median_Probability = median(values),
      Minimum_Probability = min(values),
      Maximum_Probability = max(values),
      CI95_Lower = mean(values) - qt(0.975, df = max(1, length(values) - 1)) * sd(values) / sqrt(length(values)),
      CI95_Upper = mean(values) + qt(0.975, df = max(1, length(values) - 1)) * sd(values) / sqrt(length(values)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}
by_occurrence <- summarize_group(run_results, "Held_Out_CellID")
names(by_occurrence)[1] <- "Held_Out_CellID"
by_fold <- summarize_group(run_results, "Fold_ID")
names(by_fold)[1] <- "Fold_ID"
by_repeat <- summarize_group(run_results, "Repeat")
names(by_repeat)[1] <- "Repeat"
group_original_predictor <- function(feature_name) {
  if (grepl("^Lithology", feature_name)) return("Lithology")
  feature_name
}
importance$Original_Predictor <- vapply(importance$Feature, group_original_predictor, character(1))
encoded_features <- unique(importance$Feature)
encoded_summary <- vector("list", length(encoded_features))
for (i in seq_along(encoded_features)) {
  feature_name <- encoded_features[i]
  subset_data <- importance[importance$Feature == feature_name, , drop = FALSE]
  encoded_summary[[i]] <- data.frame(
    Feature = feature_name,
    Run_n = nrow(subset_data),
    Mean_Gain = mean(subset_data$Gain),
    SD_Gain = sd(subset_data$Gain),
    Mean_Cover = mean(subset_data$Cover),
    Mean_Frequency = mean(subset_data$Frequency),
    stringsAsFactors = FALSE
  )
}
encoded_summary <- do.call(rbind, encoded_summary)
encoded_summary <- encoded_summary[order(encoded_summary$Mean_Gain, decreasing = TRUE), , drop = FALSE]
grouped_run <- aggregate(Gain ~ Run_ID + Original_Predictor, data = importance, FUN = sum)
predictors <- unique(grouped_run$Original_Predictor)
grouped_summary <- vector("list", length(predictors))
for (i in seq_along(predictors)) {
  predictor_name <- predictors[i]
  subset_data <- grouped_run[grouped_run$Original_Predictor == predictor_name, , drop = FALSE]
  grouped_summary[[i]] <- data.frame(
    Predictor = predictor_name,
    Run_n = nrow(subset_data),
    Mean_Gain = mean(subset_data$Gain),
    SD_Gain = sd(subset_data$Gain),
    Median_Gain = median(subset_data$Gain),
    Minimum_Gain = min(subset_data$Gain),
    Maximum_Gain = max(subset_data$Gain),
    stringsAsFactors = FALSE
  )
}
grouped_summary <- do.call(rbind, grouped_summary)
grouped_summary <- grouped_summary[order(grouped_summary$Mean_Gain, decreasing = TRUE), , drop = FALSE]
write.csv(overall, file.path(evaluation_directory, "XGB_Evaluation_Overall.csv"), row.names = FALSE)
write.csv(by_occurrence, file.path(evaluation_directory, "XGB_Evaluation_By_Occurrence.csv"), row.names = FALSE)
write.csv(by_fold, file.path(evaluation_directory, "XGB_Evaluation_By_Fold.csv"), row.names = FALSE)
write.csv(by_repeat, file.path(evaluation_directory, "XGB_Evaluation_By_Repeat.csv"), row.names = FALSE)
write.csv(encoded_summary, file.path(evaluation_directory, "XGB_Importance_Encoded_Summary.csv"), row.names = FALSE)
write.csv(grouped_summary, file.path(evaluation_directory, "XGB_Importance_Grouped_Summary.csv"), row.names = FALSE)
png(file.path(evaluation_directory, "XGB_HeldOut_Probability_Histogram.png"), width = 1600, height = 1100, res = 180)
hist(probability, breaks = 15, main = "XGBoost held-out occurrence probabilities", xlab = "Held-out occurrence probability")
abline(v = mean(probability), lty = 2, lwd = 2)
dev.off()
png(file.path(evaluation_directory, "XGB_Grouped_Variable_Importance.png"), width = 1800, height = 1200, res = 180)
barplot(rev(grouped_summary$Mean_Gain), names.arg = rev(grouped_summary$Predictor), horiz = TRUE, las = 1, xlab = "Mean Gain", main = "XGBoost grouped variable importance")
dev.off()
overall_qaqc <- if (nrow(run_results) == 150L && all(is.finite(probability)) && nrow(grouped_summary) > 0) "PASS" else "FAIL"
writeLines(c(
  "XGBOOST EVALUATION QAQC",
  "",
  paste0("Runs evaluated: ", nrow(run_results)),
  paste0("Invalid probabilities: ", sum(!is.finite(probability) | probability < 0 | probability > 1)),
  paste0("Grouped predictors: ", nrow(grouped_summary)),
  paste0("Overall QAQC status: ", overall_qaqc)
), "Project/QAQC/XGB_Evaluation_QAQC.txt")
writeLines(c(
  "XGBOOST EVALUATION METHODOLOGY",
  "",
  "Final performance summaries use only genuinely held-out occurrence probabilities.",
  "ROC, PR, recall, precision, specificity and F1 were not calculated because pseudo-background cells are not verified absences.",
  "Lithology dummy-variable gains were summed within each run before grouped predictor summaries were calculated."
), "Project/QAQC/XGB_Evaluation_Methodology.txt")
cat("\n========================================\n")
cat("XGBoost evaluation completed.\n")
cat("Runs evaluated:", nrow(run_results), "\n")
cat("Mean held-out probability:", round(overall$Mean_HeldOut_Probability, 4), "\n")
cat("Median held-out probability:", round(overall$Median_HeldOut_Probability, 4), "\n")
cat("Top grouped predictor:", grouped_summary$Predictor[1], "\n")
cat("Overall QAQC:", overall_qaqc, "\n")
cat("Evaluation directory:", evaluation_directory, "\n")
cat("========================================\n")
if (overall_qaqc != "PASS") stop("XGBoost evaluation QAQC failed.", call. = FALSE)
