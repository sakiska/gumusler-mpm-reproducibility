# ============================================================
# Create_LOOCV_Folds.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Create the repeated leave-one-occurrence-out fold and repeat plan used by all machine-learning algorithms.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

gc()
# ------------------------------------------------------------
# ------------------------------------------------------------
expected_occurrence_n <- 5
n_repeats <- 30
background_seeds <- 101:130
n_background_per_repeat <- 100
buffer_distance_m <- 50
algorithms <- c(
  "RF",
  "XGBoost",
  "SVM"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrence_file <- paste0(
  "Project/Tables/",
  "Occurrence_Cell_Assignment.csv"
)
candidate_pool_file <- paste0(
  "Project/Data/",
  "PseudoBackground_CandidatePool.csv"
)
output_fold_definition_file <- paste0(
  "Project/Tables/",
  "LOOCV_Fold_Definitions.csv"
)
output_fold_repeat_plan_file <- paste0(
  "Project/Tables/",
  "LOOCV_Fold_Repeat_Plan.csv"
)
output_training_membership_file <- paste0(
  "Project/Tables/",
  "LOOCV_Training_Occurrence_Membership.csv"
)
output_algorithm_plan_file <- paste0(
  "Project/Tables/",
  "LOOCV_Algorithm_Run_Plan.csv"
)
output_qaqc_dir <- "Project/QAQC"
methodology_file <- file.path(
  output_qaqc_dir,
  "LOOCV_Methodology.txt"
)
qaqc_file <- file.path(
  output_qaqc_dir,
  "LOOCV_QAQC.txt"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
dir.create(
  dirname(output_fold_definition_file),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  output_qaqc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!file.exists(occurrence_file)) {
  stop(
    paste(
      "Occurrence file not found:",
      occurrence_file
    )
  )
}
if (!file.exists(candidate_pool_file)) {
  stop(
    paste(
      "Pseudo-background candidate-pool file not found:",
      candidate_pool_file
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrences <- read.csv(
  occurrence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
candidate_pool <- read.csv(
  candidate_pool_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
required_occurrence_columns <- c(
  "Occurrence_Row",
  "X",
  "Y",
  "CellID",
  "Inside_Reference_Mask"
)
required_candidate_columns <- c(
  "CellID",
  "X",
  "Y",
  "Min_Dist_to_Occurrence_m",
  "Background_Eligible"
)
missing_occurrence_columns <- setdiff(
  required_occurrence_columns,
  names(occurrences)
)
missing_candidate_columns <- setdiff(
  required_candidate_columns,
  names(candidate_pool)
)
if (length(missing_occurrence_columns) > 0) {
  stop(
    paste(
      "Missing columns in the occurrence table:",
      paste(
        missing_occurrence_columns,
        collapse = ", "
      )
    )
  )
}
if (length(missing_candidate_columns) > 0) {
  stop(
    paste(
      "Missing columns in the candidate-pool table:",
      paste(
        missing_candidate_columns,
        collapse = ", "
      )
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (nrow(occurrences) != expected_occurrence_n) {
  stop(
    paste0(
      "Expected number of occurrences: ",
      expected_occurrence_n,
      "; observed: ",
      nrow(occurrences),
      "."
    )
  )
}
if (anyDuplicated(occurrences$Occurrence_Row) > 0) {
  stop(
    "Duplicate values were found in the Occurrence_Row column."
  )
}
if (anyDuplicated(occurrences$CellID) > 0) {
  stop(
    "Duplicate CellID values were found in the occurrence table."
  )
}
if (anyDuplicated(candidate_pool$CellID) > 0) {
  stop(
    "Duplicate CellID values were found in the candidate-pool table."
  )
}
if (!all(occurrences$Inside_Reference_Mask)) {
  stop(
    paste(
      "En az bir occurrence referans",
      "falls outside the valid grid mask."
    )
  )
}
if (anyNA(
  occurrences[, c(
    "Occurrence_Row",
    "CellID",
    "X",
    "Y"
  )]
)) {
  stop(
    paste(
      "The occurrence table contains missing Occurrence_Row,",
      "CellID, X, or Y values."
    )
  )
}
if (anyNA(
  candidate_pool[, c(
    "CellID",
    "X",
    "Y"
  )]
)) {
  stop(
    paste(
      "The candidate-pool table contains missing",
      "CellID, X, or Y values."
    )
  )
}
if (length(background_seeds) != n_repeats) {
  stop(
    paste(
      "The number of background seeds and",
      "the number of repeats are not equal."
    )
  )
}
if (anyDuplicated(background_seeds) > 0) {
  stop(
    "Duplicate values were found in the background seed list."
  )
}
if (nrow(candidate_pool) < n_background_per_repeat) {
  stop(
    paste0(
      "The candidate pool contains only ",
      nrow(candidate_pool),
      " cells. In each repeat, ",
      n_background_per_repeat,
      " background cells cannot be sampled."
    )
  )
}
if (any(
  candidate_pool$CellID %in% occurrences$CellID
)) {
  stop(
    paste(
      "An occurrence",
      "cell was found in the candidate pool."
    )
  )
}
if (any(
  candidate_pool$Min_Dist_to_Occurrence_m <=
    buffer_distance_m
)) {
  stop(
    paste(
      "A cell violating the 50 m buffer",
      "rule was found in the candidate pool."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrences <- occurrences[
  order(occurrences$Occurrence_Row),
  ,
  drop = FALSE
]
rownames(occurrences) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
fold_definitions <- data.frame(
  Fold_ID = sprintf(
    "Fold_%02d",
    seq_len(expected_occurrence_n)
  ),
  Fold_Number = seq_len(expected_occurrence_n),
  Held_Out_Occurrence_Row =
    occurrences$Occurrence_Row,
  Held_Out_CellID =
    occurrences$CellID,
  Held_Out_X =
    occurrences$X,
  Held_Out_Y =
    occurrences$Y,
  Training_Occurrence_Count =
    expected_occurrence_n - 1,
  Test_Occurrence_Count = 1L,
  Validation_Method =
    "Leave-One-Occurrence-Out",
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
fold_definitions$Training_Occurrence_Rows <- vapply(
  seq_len(nrow(fold_definitions)),
  function(i) {
    training_rows <- occurrences$Occurrence_Row[
      occurrences$Occurrence_Row !=
        fold_definitions$Held_Out_Occurrence_Row[i]
    ]
    paste(
      training_rows,
      collapse = ";"
    )
  },
  character(1)
)
fold_definitions$Training_CellIDs <- vapply(
  seq_len(nrow(fold_definitions)),
  function(i) {
    training_cellids <- occurrences$CellID[
      occurrences$Occurrence_Row !=
        fold_definitions$Held_Out_Occurrence_Row[i]
    ]
    paste(
      training_cellids,
      collapse = ";"
    )
  },
  character(1)
)
# ------------------------------------------------------------
# ------------------------------------------------------------
training_membership_list <- vector(
  mode = "list",
  length = expected_occurrence_n
)
for (fold_i in seq_len(expected_occurrence_n)) {
  held_out_row <- fold_definitions$
    Held_Out_Occurrence_Row[fold_i]
  fold_membership <- occurrences
  fold_membership$Fold_ID <-
    fold_definitions$Fold_ID[fold_i]
  fold_membership$Fold_Number <-
    fold_definitions$Fold_Number[fold_i]
  fold_membership$Role <- ifelse(
    fold_membership$Occurrence_Row ==
      held_out_row,
    "Test",
    "Training"
  )
  fold_membership$Included_In_Model_Fitting <- (
    fold_membership$Role == "Training"
  )
  fold_membership <- fold_membership[, c(
    "Fold_ID",
    "Fold_Number",
    "Occurrence_Row",
    "CellID",
    "X",
    "Y",
    "Role",
    "Included_In_Model_Fitting"
  )]
  training_membership_list[[fold_i]] <-
    fold_membership
}
training_membership <- do.call(
  rbind,
  training_membership_list
)
rownames(training_membership) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
fold_repeat_plan <- expand.grid(
  Fold_Number = seq_len(expected_occurrence_n),
  Repeat = seq_len(n_repeats),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
fold_repeat_plan <- fold_repeat_plan[
  order(
    fold_repeat_plan$Repeat,
    fold_repeat_plan$Fold_Number
  ),
  ,
  drop = FALSE
]
rownames(fold_repeat_plan) <- NULL
fold_repeat_plan$Fold_ID <- sprintf(
  "Fold_%02d",
  fold_repeat_plan$Fold_Number
)
fold_repeat_plan$Background_Seed <-
  background_seeds[
    fold_repeat_plan$Repeat
  ]
fold_repeat_plan$Background_n <-
  n_background_per_repeat
fold_repeat_plan$Buffer_m <-
  buffer_distance_m
fold_repeat_plan$Candidate_Pool_n <-
  nrow(candidate_pool)
fold_repeat_plan$Held_Out_Occurrence_Row <-
  fold_definitions$Held_Out_Occurrence_Row[
    match(
      fold_repeat_plan$Fold_Number,
      fold_definitions$Fold_Number
    )
  ]
fold_repeat_plan$Held_Out_CellID <-
  fold_definitions$Held_Out_CellID[
    match(
      fold_repeat_plan$Fold_Number,
      fold_definitions$Fold_Number
    )
  ]
fold_repeat_plan$Training_Occurrence_Count <-
  expected_occurrence_n - 1
fold_repeat_plan$Test_Occurrence_Count <- 1L
fold_repeat_plan$Run_ID <- sprintf(
  "Repeat_%02d_%s",
  fold_repeat_plan$Repeat,
  fold_repeat_plan$Fold_ID
)
fold_repeat_plan$Status <- "Planned"
fold_repeat_plan <- fold_repeat_plan[, c(
  "Run_ID",
  "Repeat",
  "Fold_ID",
  "Fold_Number",
  "Background_Seed",
  "Background_n",
  "Buffer_m",
  "Candidate_Pool_n",
  "Training_Occurrence_Count",
  "Test_Occurrence_Count",
  "Held_Out_Occurrence_Row",
  "Held_Out_CellID",
  "Status"
)]
# ------------------------------------------------------------
# ------------------------------------------------------------
algorithm_plan_list <- vector(
  mode = "list",
  length = length(algorithms)
)
for (algorithm_i in seq_along(algorithms)) {
  current_plan <- fold_repeat_plan
  current_plan$Algorithm <-
    algorithms[algorithm_i]
  current_plan$Algorithm_Run_ID <- paste(
    current_plan$Algorithm,
    current_plan$Run_ID,
    sep = "_"
  )
  current_plan <- current_plan[, c(
    "Algorithm_Run_ID",
    "Algorithm",
    "Run_ID",
    "Repeat",
    "Fold_ID",
    "Fold_Number",
    "Background_Seed",
    "Background_n",
    "Buffer_m",
    "Candidate_Pool_n",
    "Training_Occurrence_Count",
    "Test_Occurrence_Count",
    "Held_Out_Occurrence_Row",
    "Held_Out_CellID",
    "Status"
  )]
  algorithm_plan_list[[algorithm_i]] <-
    current_plan
}
algorithm_run_plan <- do.call(
  rbind,
  algorithm_plan_list
)
rownames(algorithm_run_plan) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
expected_fold_n <- expected_occurrence_n
expected_fold_repeat_n <- (
  expected_occurrence_n * n_repeats
)
expected_algorithm_run_n <- (
  expected_occurrence_n *
    n_repeats *
    length(algorithms)
)
fold_count_ok <- (
  nrow(fold_definitions) == expected_fold_n
)
fold_repeat_count_ok <- (
  nrow(fold_repeat_plan) ==
    expected_fold_repeat_n
)
algorithm_run_count_ok <- (
  nrow(algorithm_run_plan) ==
    expected_algorithm_run_n
)
each_occurrence_held_out_once <- all(
  table(
    fold_definitions$Held_Out_Occurrence_Row
  ) == 1
)
each_fold_has_four_training <- all(
  table(
    training_membership$Fold_ID[
      training_membership$Role == "Training"
    ]
  ) == expected_occurrence_n - 1
)
each_fold_has_one_test <- all(
  table(
    training_membership$Fold_ID[
      training_membership$Role == "Test"
    ]
  ) == 1
)
training_test_overlap_n <- 0L
for (fold_i in seq_len(expected_occurrence_n)) {
  fold_data <- training_membership[
    training_membership$Fold_Number == fold_i,
    ,
    drop = FALSE
  ]
  training_ids <- fold_data$CellID[
    fold_data$Role == "Training"
  ]
  test_ids <- fold_data$CellID[
    fold_data$Role == "Test"
  ]
  training_test_overlap_n <-
    training_test_overlap_n +
    sum(training_ids %in% test_ids)
}
duplicate_run_id_n <- sum(
  duplicated(fold_repeat_plan$Run_ID)
)
duplicate_algorithm_run_id_n <- sum(
  duplicated(
    algorithm_run_plan$Algorithm_Run_ID
  )
)
repeat_seed_consistency_ok <- all(
  vapply(
    seq_len(n_repeats),
    function(repeat_i) {
      repeat_seeds <- unique(
        fold_repeat_plan$Background_Seed[
          fold_repeat_plan$Repeat == repeat_i
        ]
      )
      length(repeat_seeds) == 1 &&
        repeat_seeds == background_seeds[repeat_i]
    },
    logical(1)
  )
)
algorithm_plan_consistency_ok <- all(
  table(
    algorithm_run_plan$Algorithm
  ) == expected_fold_repeat_n
)
all_checks_passed <- all(
  fold_count_ok,
  fold_repeat_count_ok,
  algorithm_run_count_ok,
  each_occurrence_held_out_once,
  each_fold_has_four_training,
  each_fold_has_one_test,
  training_test_overlap_n == 0,
  duplicate_run_id_n == 0,
  duplicate_algorithm_run_id_n == 0,
  repeat_seed_consistency_ok,
  algorithm_plan_consistency_ok
)
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  fold_definitions,
  output_fold_definition_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  fold_repeat_plan,
  output_fold_repeat_plan_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  training_membership,
  output_training_membership_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  algorithm_run_plan,
  output_algorithm_plan_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "GUMUSLER MPM - LOOCV METHODOLOGY",
  "==================================================",
  "",
  "Script:",
  "Project/Models/Create_LOOCV_Folds.R",
  "",
  "Purpose:",
  paste(
    "Create a common Leave-One-Occurrence-Out",
    "cross-validation plan for RF, XGBoost and SVM."
  ),
  "",
  "Fixed validation design:",
  paste0(
    "- Total occurrence count: ",
    expected_occurrence_n
  ),
  paste0(
    "- Total LOOCV folds: ",
    expected_fold_n
  ),
  paste0(
    "- Training occurrences per fold: ",
    expected_occurrence_n - 1
  ),
  "- Test occurrences per fold: 1",
  paste0(
    "- Background repeats: ",
    n_repeats
  ),
  paste0(
    "- Pseudo-background cells per repeat: ",
    n_background_per_repeat
  ),
  paste0(
    "- Exclusion-buffer distance: ",
    buffer_distance_m,
    " m"
  ),
  paste0(
    "- Fold-repeat combinations: ",
    expected_fold_repeat_n
  ),
  paste0(
    "- Algorithms: ",
    paste(
      algorithms,
      collapse = ", "
    )
  ),
  paste0(
    "- Total planned algorithm runs: ",
    expected_algorithm_run_n
  ),
  "",
  "LOOCV definition:",
  paste(
    "Each known occurrence is held out exactly once.",
    "The remaining four occurrences form the positive",
    "training class for that fold."
  ),
  "",
  "Shared experimental design:",
  paste(
    "RF, XGBoost and SVM use the same folds, repeat",
    "numbers and pseudo-background seeds. Therefore,",
    "differences among their validation results are not",
    "caused by different occurrence partitions or random",
    "background realizations."
  ),
  "",
  "Background seed rule:",
  paste0(
    "Repeat 1 to ",
    n_repeats,
    " use seeds ",
    min(background_seeds),
    " to ",
    max(background_seeds),
    ", respectively."
  ),
  "",
  "Leakage-control rule:",
  paste(
    "The held-out occurrence is not included in model",
    "fitting. Background sampling, predictor preprocessing,",
    "hyperparameter selection and model fitting must be",
    "executed inside the corresponding training workflow."
  ),
  "",
  "Important:",
  paste(
    "This script defines validation membership and random",
    "seeds only. It does not select pseudo-background cells",
    "and does not fit any model."
  ),
  "",
  "Input files:",
  occurrence_file,
  candidate_pool_file,
  "",
  "Output files:",
  output_fold_definition_file,
  output_fold_repeat_plan_file,
  output_training_membership_file,
  output_algorithm_plan_file
)
writeLines(
  methodology_lines,
  methodology_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc_lines <- c(
  "GUMUSLER MPM - LOOCV QA/QC",
  "==================================================",
  "",
  paste0(
    "Run date: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  paste0(
    "Occurrence count: ",
    nrow(occurrences)
  ),
  paste0(
    "Candidate-pool count: ",
    nrow(candidate_pool)
  ),
  paste0(
    "LOOCV fold count: ",
    nrow(fold_definitions)
  ),
  paste0(
    "Repeat count: ",
    n_repeats
  ),
  paste0(
    "Fold-repeat combinations: ",
    nrow(fold_repeat_plan)
  ),
  paste0(
    "Algorithm count: ",
    length(algorithms)
  ),
  paste0(
    "Total planned algorithm runs: ",
    nrow(algorithm_run_plan)
  ),
  "",
  paste0(
    "Each occurrence held out exactly once: ",
    ifelse(
      each_occurrence_held_out_once,
      "PASS",
      "FAIL"
    )
  ),
  paste0(
    "Each fold has four training occurrences: ",
    ifelse(
      each_fold_has_four_training,
      "PASS",
      "FAIL"
    )
  ),
  paste0(
    "Each fold has one test occurrence: ",
    ifelse(
      each_fold_has_one_test,
      "PASS",
      "FAIL"
    )
  ),
  paste0(
    "Training-test CellID overlaps: ",
    training_test_overlap_n
  ),
  paste0(
    "Duplicate Run_ID values: ",
    duplicate_run_id_n
  ),
  paste0(
    "Duplicate Algorithm_Run_ID values: ",
    duplicate_algorithm_run_id_n
  ),
  paste0(
    "Repeat-seed consistency: ",
    ifelse(
      repeat_seed_consistency_ok,
      "PASS",
      "FAIL"
    )
  ),
  paste0(
    "Algorithm-plan consistency: ",
    ifelse(
      algorithm_plan_consistency_ok,
      "PASS",
      "FAIL"
    )
  ),
  "",
  paste0(
    "FINAL QA/QC STATUS: ",
    ifelse(
      all_checks_passed,
      "PASS",
      "FAIL"
    )
  )
)
writeLines(
  qaqc_lines,
  qaqc_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!all_checks_passed) {
  stop(
    paste(
      "The LOOCV plan was created, but at least one QA/QC",
      "check failed.",
      "Please inspect Project/QAQC/LOOCV_QAQC.txt."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================================\n")
cat("LOOCV FOLD AND REPEAT PLAN CREATED\n")
cat("============================================================\n")
cat(
  "Number of occurrences         :",
  nrow(occurrences),
  "\n"
)
cat(
  "Number of LOOCV folds         :",
  nrow(fold_definitions),
  "\n"
)
cat(
  "Training occurrence / fold    :",
  expected_occurrence_n - 1,
  "\n"
)
cat(
  "Test occurrence / fold        :",
  1,
  "\n"
)
cat(
  "Number of background repeats  :",
  n_repeats,
  "\n"
)
cat(
  "Background / repeat           :",
  n_background_per_repeat,
  "\n"
)
cat(
  "Fold-repeat combination      :",
  nrow(fold_repeat_plan),
  "\n"
)
cat(
  "Number of algorithms          :",
  length(algorithms),
  "\n"
)
cat(
  "Total planned model runs      :",
  nrow(algorithm_run_plan),
  "\n"
)
cat(
  "Candidate pool                :",
  nrow(candidate_pool),
  "cells\n"
)
cat(
  "QA/QC status                  :",
  ifelse(
    all_checks_passed,
    "PASS",
    "FAIL"
  ),
  "\n"
)
cat("\nCREATED FILES:\n")
cat(
  "-",
  output_fold_definition_file,
  "\n"
)
cat(
  "-",
  output_fold_repeat_plan_file,
  "\n"
)
cat(
  "-",
  output_training_membership_file,
  "\n"
)
cat(
  "-",
  output_algorithm_plan_file,
  "\n"
)
cat(
  "-",
  methodology_file,
  "\n"
)
cat(
  "-",
  qaqc_file,
  "\n"
)
cat("============================================================\n")
