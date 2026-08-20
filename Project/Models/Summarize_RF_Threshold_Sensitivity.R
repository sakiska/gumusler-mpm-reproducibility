# ============================================================
# Summarize_RF_Threshold_Sensitivity.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Compare candidate Random Forest thresholds using occurrence capture, pseudo-background rejection, and spatial selectivity.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

comparison_dir <- "Project/Models/Model_Comparison"
classification_file <- file.path(
  comparison_dir,
  "RF_XGB_SVM_Independent_Validation_Classification.csv"
)
threshold_file <- file.path(
  comparison_dir,
  "RF_XGB_SVM_Threshold_Comparison.csv"
)
ensemble_threshold_file <- file.path(
  comparison_dir,
  "RF_FullGrid_Candidate_Thresholds.csv"
)
mad_file <- paste0(
  "Project/Models/RF_Results/",
  "MedianMAD_Validation/",
  "RF_MedianMAD_Overall_Validation_Metrics.csv"
)
output_file <- file.path(
  comparison_dir,
  "RF_Final_Threshold_Sensitivity.csv"
)
qaqc_file <- file.path(
  comparison_dir,
  "RF_Final_Threshold_Sensitivity_QAQC.csv"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
input_files <- c(
  classification_file,
  threshold_file,
  ensemble_threshold_file,
  mad_file
)
missing_files <- input_files[
  !file.exists(input_files)
]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input file(s):\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
classification <- read.csv(
  classification_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
thresholds <- read.csv(
  threshold_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
ensemble_thresholds <- read.csv(
  ensemble_threshold_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
mad <- read.csv(
  mad_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_class <- classification[
  classification$Algorithm == "RF",
  ,
  drop = FALSE
]
rf_class <- rf_class[
  rf_class$Threshold_Method %in%
    c("P90", "P95", "P99", "Mean_Plus_2SD"),
  ,
  drop = FALSE
]
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_thresholds <- thresholds[
  thresholds$Algorithm == "RF",
  ,
  drop = FALSE
]
rf_thresholds <- rf_thresholds[
  ,
  c(
    "Threshold_Method",
    "Mean_Threshold",
    "SD_Threshold",
    "Threshold_CV",
    "Minimum_Threshold",
    "Maximum_Threshold"
  ),
  drop = FALSE
]
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_standard <- merge(
  rf_class,
  rf_thresholds,
  by = "Threshold_Method",
  all.x = TRUE,
  sort = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
mad_row <- data.frame(
  Threshold_Method = "Median_Plus_2MAD_Scaled",
  Algorithm = "RF",
  TP = mad$TP[1],
  FN = mad$FN[1],
  TN = mad$TN[1],
  FP = mad$FP[1],
  Recall = mad$Recall[1],
  Specificity = mad$Specificity[1],
  Precision = mad$Precision[1],
  F1 = mad$F1[1],
  Balanced_Accuracy = mad$Balanced_Accuracy[1],
  Predicted_Positive_Rate = mad$Predicted_Positive_Rate[1],
  Background_False_Positive_Rate =
    mad$Background_False_Positive_Rate[1],
  Mean_Threshold = mad$Mean_Threshold[1],
  SD_Threshold = mad$SD_Threshold[1],
  Threshold_CV = mad$Threshold_CV[1],
  Minimum_Threshold = mad$Minimum_Threshold[1],
  Maximum_Threshold = mad$Maximum_Threshold[1],
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_sensitivity <- rbind(
  rf_standard,
  mad_row
)
# ------------------------------------------------------------
# ------------------------------------------------------------
ensemble_lookup <- ensemble_thresholds[
  ensemble_thresholds$Threshold_Method %in%
    c(
      "P90",
      "P95",
      "P99",
      "Mean_Plus_2SD",
      "Median_Plus_2MAD_Scaled"
    ),
  ,
  drop = FALSE
]
names(ensemble_lookup)[
  names(ensemble_lookup) == "Threshold"
] <- "Ensemble_Threshold"
names(ensemble_lookup)[
  names(ensemble_lookup) == "Grid_Cell_n"
] <- "Ensemble_Prospective_Cell_n"
names(ensemble_lookup)[
  names(ensemble_lookup) == "Grid_Percent"
] <- "Ensemble_Prospective_Grid_Percent"
rf_sensitivity <- merge(
  rf_sensitivity,
  ensemble_lookup,
  by = "Threshold_Method",
  all.x = TRUE,
  sort = FALSE
)
# ------------------------------------------------------------
# 10. ORDER METHODS
# ------------------------------------------------------------
method_order <- c(
  "P90",
  "P95",
  "P99",
  "Mean_Plus_2SD",
  "Median_Plus_2MAD_Scaled"
)
rf_sensitivity$Threshold_Method <- factor(
  rf_sensitivity$Threshold_Method,
  levels = method_order
)
rf_sensitivity <- rf_sensitivity[
  order(rf_sensitivity$Threshold_Method),
  ,
  drop = FALSE
]
rf_sensitivity$Threshold_Method <- as.character(
  rf_sensitivity$Threshold_Method
)
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_sensitivity$Occurrence_Capture_Percent <-
  100 * rf_sensitivity$Recall
rf_sensitivity$Background_Rejection_Percent <-
  100 * rf_sensitivity$Specificity
rf_sensitivity$Background_False_Positive_Percent <-
  100 * rf_sensitivity$Background_False_Positive_Rate
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_sensitivity <- rf_sensitivity[
  ,
  c(
    "Threshold_Method",
    "Ensemble_Threshold",
    "Ensemble_Prospective_Cell_n",
    "Ensemble_Prospective_Grid_Percent",
    "Mean_Threshold",
    "SD_Threshold",
    "Threshold_CV",
    "Minimum_Threshold",
    "Maximum_Threshold",
    "TP",
    "FN",
    "TN",
    "FP",
    "Recall",
    "Specificity",
    "Precision",
    "F1",
    "Balanced_Accuracy",
    "Occurrence_Capture_Percent",
    "Background_Rejection_Percent",
    "Background_False_Positive_Percent",
    "Predicted_Positive_Rate",
    "Background_False_Positive_Rate"
  ),
  drop = FALSE
]
# ------------------------------------------------------------
# 13. QA/QC
# ------------------------------------------------------------
qaqc <- data.frame(
  Check = c(
    "Exactly five threshold methods",
    "All five methods unique",
    "Each method contains 150 positives",
    "Each method contains 15000 pseudo-background negatives",
    "All ensemble thresholds finite",
    "All ensemble prospective percentages within 0-100",
    "All Recall values within 0-1",
    "All Specificity values within 0-1",
    "All F1 values within 0-1 where defined",
    "All Balanced Accuracy values within 0-1",
    "No unexpected missing final sensitivity values"
  ),
  Result = c(
    nrow(rf_sensitivity) == 5L,
    length(unique(rf_sensitivity$Threshold_Method)) == 5L,
    all(rf_sensitivity$TP + rf_sensitivity$FN == 150L),
    all(rf_sensitivity$TN + rf_sensitivity$FP == 15000L),
    all(is.finite(rf_sensitivity$Ensemble_Threshold)),
    all(
      rf_sensitivity$Ensemble_Prospective_Grid_Percent >= 0 &
        rf_sensitivity$Ensemble_Prospective_Grid_Percent <= 100
    ),
    all(rf_sensitivity$Recall >= 0 & rf_sensitivity$Recall <= 1),
    all(
      rf_sensitivity$Specificity >= 0 &
        rf_sensitivity$Specificity <= 1
    ),
    all(
      is.na(rf_sensitivity$F1) |
        (
          rf_sensitivity$F1 >= 0 &
            rf_sensitivity$F1 <= 1
        )
    ),
    all(
      rf_sensitivity$Balanced_Accuracy >= 0 &
        rf_sensitivity$Balanced_Accuracy <= 1
    ),
    all(
      complete.cases(
        rf_sensitivity[
          ,
          setdiff(
            names(rf_sensitivity),
            c("Precision", "F1")
          ),
          drop = FALSE
        ]
      )
    )
  ),
  stringsAsFactors = FALSE
)
print(qaqc)
if (!all(qaqc$Result)) {
  stop(
    "Final RF threshold-sensitivity QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 14. SAVE
# ------------------------------------------------------------
write.csv(
  rf_sensitivity,
  output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  qaqc,
  qaqc_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# 15. CONSOLE REPORT
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("RF FINAL THRESHOLD SENSITIVITY\n")
cat("============================================\n")
print(
  rf_sensitivity[
    ,
    c(
      "Threshold_Method",
      "Ensemble_Prospective_Grid_Percent",
      "TP",
      "FN",
      "TN",
      "FP",
      "Recall",
      "Specificity",
      "Precision",
      "F1",
      "Balanced_Accuracy",
      "Threshold_CV"
    )
  ],
  row.names = FALSE
)
cat("\n")
cat("============================================\n")
cat("RF THRESHOLD SENSITIVITY QA/QC\n")
cat("============================================\n")
print(
  qaqc,
  row.names = FALSE
)
cat("\n")
cat("RF FINAL THRESHOLD SENSITIVITY: PASS\n")
