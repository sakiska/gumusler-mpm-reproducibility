# ============================================================
# Prepare_Fold_Data.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Prepare leakage-controlled fold-specific training and held-out data for a single outer validation run.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

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
as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x_clean <- toupper(
    trimws(
      as.character(x)
    )
  )
  result <- rep(
    NA,
    length(x_clean)
  )
  result[x_clean %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )] <- TRUE
  result[x_clean %in% c(
    "FALSE",
    "F",
    "0",
    "NO",
    "N"
  )] <- FALSE
  result
}
# ------------------------------------------------------------
# ------------------------------------------------------------
prepare_fold_data <- function(
    run_id,
    full_grid_file = paste0(
      "Project/Data/",
      "MPM_FullGrid_25m.csv"
    ),
    occurrence_file = paste0(
      "Project/Tables/",
      "Occurrence_Cell_Assignment.csv"
    ),
    fold_membership_file = paste0(
      "Project/Tables/",
      "LOOCV_Training_Occurrence_Membership.csv"
    ),
    fold_repeat_plan_file = paste0(
      "Project/Tables/",
      "LOOCV_Fold_Repeat_Plan.csv"
    ),
    buffer_distance_m = 50,
    background_n = 100
) {
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  input_files <- c(
    full_grid_file,
    occurrence_file,
    fold_membership_file,
    fold_repeat_plan_file
  )
  missing_files <- input_files[
    !file.exists(input_files)
  ]
  if (length(missing_files) > 0) {
    stop(
      paste0(
        "The following required input files were not found:\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  full_grid <- read.csv(
    full_grid_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  occurrences <- read.csv(
    occurrence_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fold_membership <- read.csv(
    fold_membership_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fold_repeat_plan <- read.csv(
    fold_repeat_plan_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  assert_required_columns(
    full_grid,
    c(
      "CellID",
      "X",
      "Y"
    ),
    "Full-grid table"
  )
  assert_required_columns(
    occurrences,
    c(
      "Occurrence_Row",
      "CellID",
      "X",
      "Y"
    ),
    "Occurrence table"
  )
  assert_required_columns(
    fold_membership,
    c(
      "Fold_ID",
      "Fold_Number",
      "Occurrence_Row",
      "CellID",
      "X",
      "Y",
      "Role"
    ),
    "LOOCV membership table"
  )
  assert_required_columns(
    fold_repeat_plan,
    c(
      "Run_ID",
      "Repeat",
      "Fold_ID",
      "Fold_Number",
      "Background_Seed",
      "Background_n",
      "Buffer_m",
      "Held_Out_Occurrence_Row",
      "Held_Out_CellID"
    ),
    "LOOCV fold-repeat plan"
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (anyDuplicated(full_grid$CellID) > 0) {
    stop(
      "Duplicate CellID values were found in the full-grid table.",
      call. = FALSE
    )
  }
  if (anyDuplicated(occurrences$Occurrence_Row) > 0) {
    stop(
      paste(
        "Duplicate values were found in the occurrence table for",
        "Occurrence_Row."
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(occurrences$CellID) > 0) {
    stop(
      paste(
        "Duplicate values were found in the occurrence table for",
        "CellID."
      ),
      call. = FALSE
    )
  }
  if (anyNA(
    full_grid[, c(
      "CellID",
      "X",
      "Y"
    )]
  )) {
    stop(
      paste(
        "The full-grid table contains missing",
        "CellID, X, or Y values."
      ),
      call. = FALSE
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
        "The occurrence table contains missing",
        "Occurrence_Row, CellID, X, or Y values."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  selected_run <- fold_repeat_plan[
    fold_repeat_plan$Run_ID == run_id,
    ,
    drop = FALSE
  ]
  if (nrow(selected_run) == 0) {
    stop(
      paste0(
        "Run_ID not found: ",
        run_id
      ),
      call. = FALSE
    )
  }
  if (nrow(selected_run) > 1) {
    stop(
      paste0(
        "Run_ID was found more than once: ",
        run_id
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  current_fold_id <-
    selected_run$Fold_ID[1]
  current_fold_number <-
    selected_run$Fold_Number[1]
  current_repeat <-
    selected_run$Repeat[1]
  current_seed <-
    selected_run$Background_Seed[1]
  planned_background_n <-
    selected_run$Background_n[1]
  planned_buffer_m <-
    selected_run$Buffer_m[1]
  if (planned_background_n != background_n) {
    stop(
      paste0(
        "The function background_n value does not match the plan file. ",
        "Function value: ",
        background_n,
        "; plan: ",
        planned_background_n,
        "."
      ),
      call. = FALSE
    )
  }
  if (planned_buffer_m != buffer_distance_m) {
    stop(
      paste0(
        "The function buffer value does not match the plan file. ",
        "Function value: ",
        buffer_distance_m,
        " m; plan: ",
        planned_buffer_m,
        " m."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  current_membership <- fold_membership[
    fold_membership$Fold_ID == current_fold_id,
    ,
    drop = FALSE
  ]
  if (nrow(current_membership) == 0) {
    stop(
      paste0(
        "Fold not found in the membership table: ",
        current_fold_id
      ),
      call. = FALSE
    )
  }
  training_occurrences <- current_membership[
    current_membership$Role == "Training",
    ,
    drop = FALSE
  ]
  test_occurrence <- current_membership[
    current_membership$Role == "Test",
    ,
    drop = FALSE
  ]
  if (nrow(training_occurrences) != 4) {
    stop(
      paste0(
        current_fold_id,
        " does not have exactly 4 training occurrences: ",
        nrow(training_occurrences)
      ),
      call. = FALSE
    )
  }
  if (nrow(test_occurrence) != 1) {
    stop(
      paste0(
        current_fold_id,
        " does not have exactly 1 test occurrence: ",
        nrow(test_occurrence)
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  if (
    test_occurrence$Occurrence_Row[1] !=
      selected_run$Held_Out_Occurrence_Row[1]
  ) {
    stop(
      paste(
        "The held-out occurrence in the run plan and membership table",
        "do not match."
      ),
      call. = FALSE
    )
  }
  if (
    test_occurrence$CellID[1] !=
      selected_run$Held_Out_CellID[1]
  ) {
    stop(
      paste(
        "The held-out occurrence in the run plan and membership table",
        "The held-out CellID values do not match."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  min_distance_to_training <- rep(
    Inf,
    nrow(full_grid)
  )
  for (
    occurrence_i in
    seq_len(nrow(training_occurrences))
  ) {
    dx <- (
      full_grid$X -
        training_occurrences$X[occurrence_i]
    )
    dy <- (
      full_grid$Y -
        training_occurrences$Y[occurrence_i]
    )
    current_distance <- sqrt(
      dx^2 + dy^2
    )
    min_distance_to_training <- pmin(
      min_distance_to_training,
      current_distance
    )
  }
  full_grid$Min_Dist_to_Training_Occurrence_m <-
    min_distance_to_training
  # ----------------------------------------------------------
  #
  # Rules:
  # ----------------------------------------------------------
  known_occurrence_cellids <-
    occurrences$CellID
  fold_candidate_pool <- full_grid[
    full_grid$Min_Dist_to_Training_Occurrence_m >
      buffer_distance_m &
      !full_grid$CellID %in%
        known_occurrence_cellids,
    ,
    drop = FALSE
  ]
  rownames(fold_candidate_pool) <- NULL
  if (nrow(fold_candidate_pool) < background_n) {
    stop(
      paste0(
        current_fold_id,
        " has only ",
        nrow(fold_candidate_pool),
        " eligible background candidates. Required: ",
        background_n,
        "."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  set.seed(current_seed)
  selected_indices <- sample(
    x = seq_len(
      nrow(fold_candidate_pool)
    ),
    size = background_n,
    replace = FALSE
  )
  selected_background <- fold_candidate_pool[
    selected_indices,
    ,
    drop = FALSE
  ]
  rownames(selected_background) <- NULL
  # ----------------------------------------------------------
  #
  # ----------------------------------------------------------
  training_positive_grid <- full_grid[
    match(
      training_occurrences$CellID,
      full_grid$CellID
    ),
    ,
    drop = FALSE
  ]
  test_positive_grid <- full_grid[
    match(
      test_occurrence$CellID,
      full_grid$CellID
    ),
    ,
    drop = FALSE
  ]
  if (anyNA(training_positive_grid$CellID)) {
    stop(
      paste(
        "At least one training occurrence CellID",
        "was not found in the full grid."
      ),
      call. = FALSE
    )
  }
  if (
    nrow(test_positive_grid) != 1 ||
      is.na(test_positive_grid$CellID[1])
  ) {
    stop(
      paste(
        "Held-out occurrence CellID",
        "was not found in the full grid."
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  training_positive_grid$Class <- 1L
  training_positive_grid$Class_Label <- "Occurrence"
  training_positive_grid$Data_Role <- "Training"
  training_positive_grid$Source_Type <-
    "Training_Occurrence"
  selected_background$Class <- 0L
  selected_background$Class_Label <-
    "PseudoBackground"
  selected_background$Data_Role <- "Training"
  selected_background$Source_Type <-
    "Fold_Specific_PseudoBackground"
  test_positive_grid$Class <- 1L
  test_positive_grid$Class_Label <- "Occurrence"
  test_positive_grid$Data_Role <- "Test"
  test_positive_grid$Source_Type <-
    "Held_Out_Occurrence"
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  metadata_columns <- list(
    Run_ID = run_id,
    Fold_ID = current_fold_id,
    Fold_Number = current_fold_number,
    Repeat = current_repeat,
    Background_Seed = current_seed,
    Buffer_m = buffer_distance_m
  )
  add_metadata <- function(data) {
    for (
      metadata_name in
      names(metadata_columns)
    ) {
      data[[metadata_name]] <-
        metadata_columns[[metadata_name]]
    }
    data
  }
  training_positive_grid <-
    add_metadata(training_positive_grid)
  selected_background <-
    add_metadata(selected_background)
  test_positive_grid <-
    add_metadata(test_positive_grid)
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  common_columns <- intersect(
    names(training_positive_grid),
    names(selected_background)
  )
  training_data <- rbind(
    training_positive_grid[
      ,
      common_columns,
      drop = FALSE
    ],
    selected_background[
      ,
      common_columns,
      drop = FALSE
    ]
  )
  rownames(training_data) <- NULL
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  training_positive_n <-
    sum(training_data$Class == 1)
  training_background_n <-
    sum(training_data$Class == 0)
  selected_background_duplicate_n <-
    sum(
      duplicated(
        selected_background$CellID
      )
    )
  selected_background_occurrence_overlap_n <-
    sum(
      selected_background$CellID %in%
        known_occurrence_cellids
    )
  training_test_overlap_n <-
    sum(
      training_data$CellID %in%
        test_positive_grid$CellID
    )
  minimum_selected_background_distance <-
    min(
      selected_background$
        Min_Dist_to_Training_Occurrence_m
    )
  buffer_rule_ok <- (
    minimum_selected_background_distance >
      buffer_distance_m
  )
  reproducibility_check <- {
    set.seed(current_seed)
    repeated_indices <- sample(
      x = seq_len(
        nrow(fold_candidate_pool)
      ),
      size = background_n,
      replace = FALSE
    )
    identical(
      selected_indices,
      repeated_indices
    )
  }
  all_checks_passed <- all(
    training_positive_n == 4,
    training_background_n == background_n,
    nrow(test_positive_grid) == 1,
    selected_background_duplicate_n == 0,
    selected_background_occurrence_overlap_n == 0,
    training_test_overlap_n == 0,
    buffer_rule_ok,
    reproducibility_check
  )
  if (!all_checks_passed) {
    stop(
      paste0(
        "Fold-data QA/QC failed: ",
        run_id
      ),
      call. = FALSE
    )
  }
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  audit <- data.frame(
    Run_ID = run_id,
    Fold_ID = current_fold_id,
    Fold_Number = current_fold_number,
    Repeat = current_repeat,
    Background_Seed = current_seed,
    Buffer_m = buffer_distance_m,
    Full_Grid_n = nrow(full_grid),
    Training_Occurrence_n =
      nrow(training_occurrences),
    Test_Occurrence_n =
      nrow(test_occurrence),
    Fold_Candidate_Pool_n =
      nrow(fold_candidate_pool),
    Selected_Background_n =
      nrow(selected_background),
    Minimum_Selected_Background_Distance_m =
      minimum_selected_background_distance,
    Background_Duplicate_n =
      selected_background_duplicate_n,
    Background_Known_Occurrence_Overlap_n =
      selected_background_occurrence_overlap_n,
    Training_Test_CellID_Overlap_n =
      training_test_overlap_n,
    Reproducibility_Check =
      reproducibility_check,
    QAQC_Status = ifelse(
      all_checks_passed,
      "PASS",
      "FAIL"
    ),
    stringsAsFactors = FALSE
  )
  # ----------------------------------------------------------
  # ----------------------------------------------------------
  result <- list(
    run_metadata = selected_run,
    training_occurrences =
      training_occurrences,
    test_occurrence =
      test_occurrence,
    fold_candidate_pool =
      fold_candidate_pool,
    selected_background =
      selected_background,
    training_positive_grid =
      training_positive_grid,
    training_data =
      training_data,
    test_data =
      test_positive_grid,
    audit =
      audit
  )
  return(result)
}
# ============================================================
# ============================================================
