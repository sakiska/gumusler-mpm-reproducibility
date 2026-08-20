# ============================================================
# Create_Independent_Validation_Background.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Generate pseudo-background observations used for independent algorithm comparison.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

project_dir <- "Project"
models_dir <- file.path(
  project_dir,
  "Models"
)
tables_dir <- file.path(
  project_dir,
  "Tables"
)
qaqc_dir <- file.path(
  project_dir,
  "QAQC"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
prepare_fold_script <- file.path(
  models_dir,
  "Prepare_Fold_Data.R"
)
fold_repeat_plan_file <- file.path(
  tables_dir,
  "LOOCV_Fold_Repeat_Plan.csv"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
output_validation_file <- file.path(
  tables_dir,
  "Independent_Validation_PseudoBackground.csv"
)
output_summary_file <- file.path(
  tables_dir,
  "Independent_Validation_PseudoBackground_Summary.csv"
)
output_qaqc_file <- file.path(
  qaqc_dir,
  "Independent_Validation_PseudoBackground_QAQC.txt"
)
output_methodology_file <- file.path(
  qaqc_dir,
  "Independent_Validation_PseudoBackground_Methodology.txt"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
# per fold-repeat run.
#
validation_background_n <- 100L
#
validation_seeds <- 1001:1030
# ------------------------------------------------------------
# ------------------------------------------------------------
required_files <- c(
  prepare_fold_script,
  fold_repeat_plan_file
)
missing_files <- required_files[
  !file.exists(required_files)
]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}
if (!dir.exists(qaqc_dir)) {
  dir.create(
    qaqc_dir,
    recursive = TRUE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
source(
  prepare_fold_script
)
if (!exists("prepare_fold_data")) {
  stop(
    "prepare_fold_data() function was not created after sourcing Prepare_Fold_Data.R."
  )
}
# ------------------------------------------------------------
# 7. READ FOLD-REPEAT PLAN
# ------------------------------------------------------------
fold_repeat_plan <- read.csv(
  fold_repeat_plan_file,
  stringsAsFactors = FALSE
)
required_plan_columns <- c(
  "Run_ID",
  "Fold_ID",
  "Fold_Number",
  "Repeat",
  "Background_Seed",
  "Held_Out_CellID"
)
missing_plan_columns <- setdiff(
  required_plan_columns,
  names(fold_repeat_plan)
)
if (length(missing_plan_columns) > 0) {
  stop(
    paste(
      "Missing columns in LOOCV fold-repeat plan:",
      paste(
        missing_plan_columns,
        collapse = ", "
      )
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
repeat_values <- sort(
  unique(
    fold_repeat_plan$Repeat
  )
)
if (length(repeat_values) != 30) {
  stop(
    paste0(
      "Expected 30 repeats, found ",
      length(repeat_values),
      "."
    )
  )
}
if (length(validation_seeds) != length(repeat_values)) {
  stop(
    "Validation seed count does not match repeat count."
  )
}
if (anyDuplicated(validation_seeds) > 0) {
  stop(
    "Validation seed list contains duplicates."
  )
}
if (
  any(
    validation_seeds %in%
      fold_repeat_plan$Background_Seed
  )
) {
  stop(
    paste(
      "Validation seeds overlap",
      "training background seeds."
    )
  )
}
# ------------------------------------------------------------
# 9. STORAGE OBJECTS
# ------------------------------------------------------------
validation_list <- vector(
  "list",
  nrow(fold_repeat_plan)
)
summary_list <- vector(
  "list",
  nrow(fold_repeat_plan)
)
# ------------------------------------------------------------
# ------------------------------------------------------------
for (run_i in seq_len(nrow(fold_repeat_plan))) {
  current_plan <- fold_repeat_plan[
    run_i,
    ,
    drop = FALSE
  ]
  current_run_id <-
    current_plan$Run_ID[1]
  current_repeat <-
    current_plan$Repeat[1]
  current_fold_id <-
    current_plan$Fold_ID[1]
  current_fold_number <-
    current_plan$Fold_Number[1]
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  fold_data <- prepare_fold_data(
    run_id = current_run_id
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  candidate_pool <-
    fold_data$fold_candidate_pool
  training_background <-
    fold_data$selected_background
  training_occurrences <-
    fold_data$training_positive_grid
  held_out_occurrence <-
    fold_data$test_data
  # ----------------------------------------------------------
  # 10.3 QA/QC expected fold structure
  # ----------------------------------------------------------
  if (nrow(training_occurrences) != 4) {
    stop(
      paste(
        current_run_id,
        "does not contain exactly 4 training occurrences."
      )
    )
  }
  if (nrow(training_background) != 100) {
    stop(
      paste(
        current_run_id,
        "does not contain exactly 100 training pseudo-background cells."
      )
    )
  }
  if (nrow(held_out_occurrence) != 1) {
    stop(
      paste(
        current_run_id,
        "does not contain exactly 1 held-out occurrence."
      )
    )
  }
  # ----------------------------------------------------------
  #
  #
  #
  # occurrence here.
  #
  # ----------------------------------------------------------
  validation_candidate_pool <- candidate_pool[
    !candidate_pool$CellID %in%
      training_background$CellID,
    ,
    drop = FALSE
  ]
  rownames(
    validation_candidate_pool
  ) <- NULL
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  candidate_training_background_overlap_n <-
    sum(
      validation_candidate_pool$CellID %in%
        training_background$CellID
    )
  candidate_training_occurrence_overlap_n <-
    sum(
      validation_candidate_pool$CellID %in%
        training_occurrences$CellID
    )
  candidate_heldout_occurrence_overlap_n <-
    sum(
      validation_candidate_pool$CellID %in%
        held_out_occurrence$CellID
    )
  if (
    candidate_training_background_overlap_n != 0
  ) {
    stop(
      paste(
        current_run_id,
        "validation candidate pool overlaps",
        "training background."
      )
    )
  }
  if (
    candidate_training_occurrence_overlap_n != 0
  ) {
    stop(
      paste(
        current_run_id,
        "validation candidate pool overlaps",
        "training occurrences."
      )
    )
  }
  if (
    candidate_heldout_occurrence_overlap_n != 0
  ) {
    stop(
      paste(
        current_run_id,
        "validation candidate pool contains",
        "the held-out occurrence CellID."
      )
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (
    nrow(validation_candidate_pool) <
      validation_background_n
  ) {
    stop(
      paste(
        current_run_id,
        "has only",
        nrow(validation_candidate_pool),
        "eligible validation background cells;",
        validation_background_n,
        "required."
      )
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  validation_seed <-
    validation_seeds[
      match(
        current_repeat,
        repeat_values
      )
    ]
  set.seed(
    validation_seed
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  selected_indices <- sample(
    seq_len(
      nrow(validation_candidate_pool)
    ),
    size = validation_background_n,
    replace = FALSE
  )
  validation_background <-
    validation_candidate_pool[
      selected_indices,
      ,
      drop = FALSE
    ]
  rownames(
    validation_background
  ) <- NULL
  # ----------------------------------------------------------
  # 10.9 Add metadata
  # ----------------------------------------------------------
  validation_background$Class <- 0L
  validation_background$Class_Label <-
    "PseudoBackground"
  validation_background$Data_Role <-
    "Independent_Validation"
  validation_background$Source_Type <-
    "Independent_Validation_PseudoBackground"
  validation_background$Run_ID <-
    current_run_id
  validation_background$Repeat <-
    current_repeat
  validation_background$Fold_ID <-
    current_fold_id
  validation_background$Fold_Number <-
    current_fold_number
  validation_background$Validation_Seed <-
    validation_seed
  validation_background$Held_Out_CellID <-
    current_plan$Held_Out_CellID[1]
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  validation_duplicate_n <-
    sum(
      duplicated(
        validation_background$CellID
      )
    )
  validation_training_background_overlap_n <-
    sum(
      validation_background$CellID %in%
        training_background$CellID
    )
  validation_training_occurrence_overlap_n <-
    sum(
      validation_background$CellID %in%
        training_occurrences$CellID
    )
  validation_heldout_overlap_n <-
    sum(
      validation_background$CellID %in%
        held_out_occurrence$CellID
    )
  minimum_validation_distance_to_training_occurrence <-
    min(
      validation_background$
        Min_Dist_to_Training_Occurrence_m,
      na.rm = TRUE
    )
  run_qaqc_pass <-
    (
      nrow(validation_background) ==
        validation_background_n
    ) &&
    (
      validation_duplicate_n == 0
    ) &&
    (
      validation_training_background_overlap_n == 0
    ) &&
    (
      validation_training_occurrence_overlap_n == 0
    ) &&
    (
      validation_heldout_overlap_n == 0
    ) &&
    (
      minimum_validation_distance_to_training_occurrence >
        50
    )
  if (!run_qaqc_pass) {
    stop(
      paste(
        "Independent validation QA/QC failed for",
        current_run_id
      )
    )
  }
  # ----------------------------------------------------------
  # 10.11 Store outputs
  # ----------------------------------------------------------
  validation_list[[run_i]] <-
    validation_background
  summary_list[[run_i]] <- data.frame(
    Run_ID =
      current_run_id,
    Repeat =
      current_repeat,
    Fold_ID =
      current_fold_id,
    Fold_Number =
      current_fold_number,
    Training_Background_Seed =
      current_plan$Background_Seed[1],
    Validation_Seed =
      validation_seed,
    Fold_Candidate_Pool_n =
      nrow(candidate_pool),
    Validation_Candidate_Pool_n =
      nrow(validation_candidate_pool),
    Validation_Background_n =
      nrow(validation_background),
    Validation_Duplicate_n =
      validation_duplicate_n,
    Validation_Training_Background_Overlap_n =
      validation_training_background_overlap_n,
    Validation_Training_Occurrence_Overlap_n =
      validation_training_occurrence_overlap_n,
    Validation_HeldOut_Occurrence_Overlap_n =
      validation_heldout_overlap_n,
    Minimum_Validation_Distance_to_Training_Occurrence_m =
      minimum_validation_distance_to_training_occurrence,
    QAQC_PASS =
      run_qaqc_pass,
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
validation_background_all <- do.call(
  rbind,
  validation_list
)
rownames(
  validation_background_all
) <- NULL
validation_summary <- do.call(
  rbind,
  summary_list
)
rownames(
  validation_summary
) <- NULL
# ------------------------------------------------------------
# 12. GLOBAL QA/QC
# ------------------------------------------------------------
expected_run_n <-
  nrow(fold_repeat_plan)
expected_validation_row_n <-
  expected_run_n *
  validation_background_n
global_qaqc <- data.frame(
  Check = c(
    "Expected fold-repeat runs",
    "Expected validation rows",
    "Exactly 100 validation backgrounds per run",
    "No duplicate validation CellID within run",
    "No validation-training background overlap",
    "No validation-training occurrence overlap",
    "No validation-held-out occurrence overlap",
    "All validation cells >50 m from training occurrences",
    "All run-level QAQC PASS"
  ),
  Result = c(
    nrow(validation_summary) ==
      expected_run_n,
    nrow(validation_background_all) ==
      expected_validation_row_n,
    all(
      validation_summary$
        Validation_Background_n ==
        validation_background_n
    ),
    all(
      validation_summary$
        Validation_Duplicate_n == 0
    ),
    all(
      validation_summary$
        Validation_Training_Background_Overlap_n == 0
    ),
    all(
      validation_summary$
        Validation_Training_Occurrence_Overlap_n == 0
    ),
    all(
      validation_summary$
        Validation_HeldOut_Occurrence_Overlap_n == 0
    ),
    all(
      validation_summary$
        Minimum_Validation_Distance_to_Training_Occurrence_m >
        50
    ),
    all(
      validation_summary$QAQC_PASS
    )
  ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!all(global_qaqc$Result)) {
  print(
    global_qaqc
  )
  stop(
    paste(
      "Global independent validation",
      "pseudo-background QA/QC failed."
    )
  )
}
# ------------------------------------------------------------
# 14. SAVE OUTPUT TABLES
# ------------------------------------------------------------
write.csv(
  validation_background_all,
  output_validation_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  validation_summary,
  output_summary_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc_lines <- c(
  "INDEPENDENT VALIDATION PSEUDO-BACKGROUND QA/QC",
  "==================================================",
  "",
  paste(
    "Fold-repeat runs:",
    expected_run_n
  ),
  paste(
    "Validation pseudo-background cells per run:",
    validation_background_n
  ),
  paste(
    "Total validation rows:",
    nrow(validation_background_all)
  ),
  "",
  "GLOBAL QA/QC",
  "--------------------------------------------------",
  paste(
    global_qaqc$Check,
    global_qaqc$Result,
    sep = ": "
  ),
  "",
  "LEAKAGE CONTROL",
  "--------------------------------------------------",
  paste(
    "- Validation cells were selected only from the",
    "fold-specific candidate pool constructed from",
    "training-occurrence information."
  ),
  paste(
    "- The held-out occurrence was NOT used to",
    "construct a distance buffer for validation",
    "pseudo-background selection."
  ),
  paste(
    "- Training pseudo-background cells were",
    "explicitly excluded from validation sampling."
  ),
  paste(
    "- All known occurrence CellIDs were excluded",
    "from pseudo-background eligibility by the",
    "underlying fold preparation workflow."
  ),
  paste(
    "- Validation pseudo-background is independent",
    "of model fitting, preprocessing, tuning and",
    "threshold calculation."
  ),
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  output_qaqc_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "INDEPENDENT VALIDATION PSEUDO-BACKGROUND METHODOLOGY",
  "====================================================",
  "",
  paste(
    "For each of the 150 repeated LOOCV runs,",
    validation_background_n,
    "independent pseudo-background cells were selected."
  ),
  "",
  paste(
    "The validation cells were sampled from the same",
    "fold-specific eligibility domain used by the",
    "training workflow, where eligibility was determined",
    "using only the four training occurrences."
  ),
  "",
  paste(
    "The held-out occurrence was not used to define",
    "the 50 m exclusion distance. This prevents test",
    "information from influencing validation-background",
    "selection."
  ),
  "",
  paste(
    "Cells already selected as training pseudo-background",
    "for the corresponding fold-repeat run were removed",
    "before validation sampling."
  ),
  "",
  paste(
    "Validation pseudo-background cells will later be",
    "used together with the single held-out occurrence",
    "to calculate occurrence-versus-pseudo-background",
    "classification metrics."
  ),
  "",
  paste(
    "Pseudo-background cells are not interpreted as",
    "confirmed mineralization absences. Therefore",
    "specificity, precision, F1 and related metrics",
    "must be interpreted as discrimination against",
    "pseudo-background rather than true absence."
  )
)
writeLines(
  methodology_lines,
  output_methodology_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("INDEPENDENT VALIDATION BACKGROUND\n")
cat("============================================\n")
cat(
  "Fold-repeat runs          :",
  expected_run_n,
  "\n"
)
cat(
  "Validation / run          :",
  validation_background_n,
  "\n"
)
cat(
  "Total validation rows     :",
  nrow(validation_background_all),
  "\n"
)
cat("\n")
print(
  global_qaqc
)
cat("\n")
cat(
  "Output:",
  output_validation_file,
  "\n"
)
cat(
  "Summary:",
  output_summary_file,
  "\n"
)
cat(
  "QAQC:",
  output_qaqc_file,
  "\n"
)
cat("\n")
cat(
  "INDEPENDENT VALIDATION PSEUDO-BACKGROUND: PASS\n"
)
