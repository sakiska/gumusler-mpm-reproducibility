# ============================================================
# Create_PseudoBackground.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Create the global pseudo-background candidate pool using the occurrence exclusion rule.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

rm(list = ls())
gc()
# ------------------------------------------------------------
# ------------------------------------------------------------
buffer_distance_m <- 50
expected_grid_n       <- 9600
expected_occurrence_n <- 5
predictor_names <- c(
  "Pb_OK",
  "Zn_OK",
  "Cu_OK",
  "Lithology",
  "Dist_Fault",
  "Dist_Silicified",
  "Dist_Brecciated"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
full_grid_file <- "Project/Data/MPM_FullGrid_25m.csv"
occurrence_file <- paste0(
  "Project/Tables/",
  "Occurrence_Cell_Assignment.csv"
)
output_candidate_file <- paste0(
  "Project/Data/",
  "PseudoBackground_CandidatePool.csv"
)
output_buffer_audit_file <- paste0(
  "Project/Tables/",
  "PseudoBackground_Buffer_Audit.csv"
)
output_summary_file <- paste0(
  "Project/Tables/",
  "PseudoBackground_CandidatePool_Summary.csv"
)
output_qaqc_dir <- "Project/QAQC"
methodology_file <- file.path(
  output_qaqc_dir,
  "PseudoBackground_Methodology.txt"
)
qaqc_file <- file.path(
  output_qaqc_dir,
  "PseudoBackground_QAQC.txt"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
dir.create(
  dirname(output_candidate_file),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(output_buffer_audit_file),
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
if (!file.exists(full_grid_file)) {
  stop(
    paste(
      "Full-grid file not found:",
      full_grid_file
    )
  )
}
if (!file.exists(occurrence_file)) {
  stop(
    paste(
      "Occurrence file not found:",
      occurrence_file
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
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
# ------------------------------------------------------------
# ------------------------------------------------------------
required_grid_columns <- c(
  "CellID",
  "X",
  "Y",
  predictor_names,
  "Occurrence"
)
required_occurrence_columns <- c(
  "Occurrence_Row",
  "X",
  "Y",
  "CellID",
  "Inside_Reference_Mask"
)
missing_grid_columns <- setdiff(
  required_grid_columns,
  names(full_grid)
)
missing_occurrence_columns <- setdiff(
  required_occurrence_columns,
  names(occurrences)
)
if (length(missing_grid_columns) > 0) {
  stop(
    paste(
      "Missing columns in the full-grid table:",
      paste(missing_grid_columns, collapse = ", ")
    )
  )
}
if (length(missing_occurrence_columns) > 0) {
  stop(
    paste(
      "Missing columns in the occurrence table:",
      paste(missing_occurrence_columns, collapse = ", ")
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (nrow(full_grid) != expected_grid_n) {
  stop(
    paste0(
      "Expected number of grid cells: ",
      expected_grid_n,
      "; observed: ",
      nrow(full_grid),
      "."
    )
  )
}
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
if (anyDuplicated(full_grid$CellID) > 0) {
  stop("Duplicate CellID values were found in the full-grid table.")
}
if (anyDuplicated(occurrences$CellID) > 0) {
  stop("Duplicate CellID values were found in the occurrence table.")
}
if (!all(occurrences$CellID %in% full_grid$CellID)) {
  stop(
    paste(
      "Some occurrence CellID values",
      "were not found in the full grid."
    )
  )
}
if (!all(occurrences$Inside_Reference_Mask)) {
  stop(
    paste(
      "At least one occurrence is",
      "located outside the mask."
    )
  )
}
if (anyNA(full_grid[, c("CellID", "X", "Y")])) {
  stop(
    "Missing CellID, X, or Y values were found in the full-grid table."
  )
}
if (anyNA(occurrences[, c("CellID", "X", "Y")])) {
  stop(
    "Missing CellID, X, or Y values were found in the occurrence table."
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
distance_matrix <- sapply(
  seq_len(nrow(occurrences)),
  function(i) {
    sqrt(
      (full_grid$X - occurrences$X[i])^2 +
      (full_grid$Y - occurrences$Y[i])^2
    )
  }
)
distance_matrix <- as.matrix(distance_matrix)
full_grid$Min_Dist_to_Occurrence_m <- apply(
  distance_matrix,
  1,
  min
)
# ------------------------------------------------------------
# ------------------------------------------------------------
full_grid$Is_Occurrence_Cell <- (
  full_grid$CellID %in% occurrences$CellID
)
full_grid$Inside_50m_Buffer <- (
  full_grid$Min_Dist_to_Occurrence_m <= buffer_distance_m
)
# ------------------------------------------------------------
# ------------------------------------------------------------
candidate_pool <- full_grid[
  !full_grid$Is_Occurrence_Cell &
    !full_grid$Inside_50m_Buffer,
  ,
  drop = FALSE
]
candidate_pool$Background_Eligible <- TRUE
candidate_pool <- candidate_pool[, c(
  "CellID",
  "X",
  "Y",
  predictor_names,
  "Min_Dist_to_Occurrence_m",
  "Background_Eligible"
)]
rownames(candidate_pool) <- NULL
# ------------------------------------------------------------
# ------------------------------------------------------------
buffer_audit <- full_grid[, c(
  "CellID",
  "X",
  "Y",
  "Occurrence",
  "Is_Occurrence_Cell",
  "Min_Dist_to_Occurrence_m",
  "Inside_50m_Buffer"
)]
buffer_audit$Background_Eligible <- (
  !buffer_audit$Is_Occurrence_Cell &
    !buffer_audit$Inside_50m_Buffer
)
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrence_cells_found <- sum(
  full_grid$Is_Occurrence_Cell
)
buffer_violation_n <- sum(
  candidate_pool$Min_Dist_to_Occurrence_m <=
    buffer_distance_m
)
occurrence_overlap_n <- sum(
  candidate_pool$CellID %in% occurrences$CellID
)
duplicate_candidate_n <- sum(
  duplicated(candidate_pool$CellID)
)
missing_predictor_counts <- sapply(
  candidate_pool[, predictor_names, drop = FALSE],
  function(x) sum(is.na(x))
)
all_checks_passed <- all(
  occurrence_cells_found == expected_occurrence_n,
  buffer_violation_n == 0,
  occurrence_overlap_n == 0,
  duplicate_candidate_n == 0,
  nrow(candidate_pool) > 100
)
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  candidate_pool,
  output_candidate_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  buffer_audit,
  output_buffer_audit_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
summary_table <- data.frame(
  Grid_Cell_Count = nrow(full_grid),
  Occurrence_Count = nrow(occurrences),
  Buffer_Distance_m = buffer_distance_m,
  Cells_Inside_or_Touching_Buffer = sum(
    full_grid$Inside_50m_Buffer
  ),
  Occurrence_Cells = occurrence_cells_found,
  Candidate_Pool_Count = nrow(candidate_pool),
  Buffer_Violations = buffer_violation_n,
  Occurrence_Overlap = occurrence_overlap_n,
  Duplicate_Candidate_CellID = duplicate_candidate_n,
  QAQC_Status = ifelse(
    all_checks_passed,
    "PASS",
    "FAIL"
  )
)
write.csv(
  summary_table,
  output_summary_file,
  row.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
methodology_lines <- c(
  "GUMUSLER MPM - PSEUDO-BACKGROUND METHODOLOGY",
  "==================================================",
  "",
  "Script:",
  "Project/Models/Create_PseudoBackground.R",
  "",
  "Purpose:",
  paste(
    "Create a common pseudo-background candidate pool",
    "for subsequent RF, XGBoost and SVM modeling."
  ),
  "",
  "Fixed methodological decisions:",
  "- Grid resolution: 25 x 25 m",
  paste0(
    "- Known occurrence count: ",
    expected_occurrence_n
  ),
  paste0(
    "- Exclusion-buffer distance: ",
    buffer_distance_m,
    " m"
  ),
  "- Buffer geometry: circular Euclidean distance",
  "- Boundary rule: distances <= 50 m are excluded",
  "",
  "Interpretation:",
  paste(
    "The 50 m buffer does not represent the physical",
    "thickness of mineralized veins. It is an exclusion",
    "zone intended to prevent immediately adjacent and",
    "potentially mineralized cells from being labelled as",
    "pseudo-background."
  ),
  "",
  "Candidate-pool definition:",
  paste(
    "All grid cells located more than 50 m from every",
    "known occurrence and not themselves occurrence cells."
  ),
  "",
  "Important class interpretation:",
  paste(
    "Pseudo-background cells are not confirmed barren or",
    "true absence locations."
  ),
  "",
  "Sampling decision:",
  paste(
    "The selection of 100 pseudo-background cells and the",
    "30 repeated samples will be performed later within",
    "the common model-validation workflow. RF, XGBoost",
    "and SVM will use identical background samples for",
    "each repeat and validation fold."
  ),
  "",
  "Leakage-control decision:",
  paste(
    "During leave-one-occurrence-out validation, background",
    "sampling, preprocessing and model tuning will be",
    "performed inside the training workflow. The held-out",
    "occurrence will never enter model fitting."
  ),
  "",
  paste0(
    "Candidate-pool size produced: ",
    nrow(candidate_pool)
  ),
  "",
  "Predictors:",
  paste(predictor_names, collapse = ", ")
)
writeLines(
  methodology_lines,
  methodology_file
)
# ------------------------------------------------------------
# 16. QA/QC raporu
# ------------------------------------------------------------
qaqc_lines <- c(
  "GUMUSLER MPM - PSEUDO-BACKGROUND QA/QC",
  "==================================================",
  "",
  paste0(
    "Run date: ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  "",
  paste0(
    "Input grid cells: ",
    nrow(full_grid)
  ),
  paste0(
    "Input occurrences: ",
    nrow(occurrences)
  ),
  paste0(
    "Occurrence cells found in grid: ",
    occurrence_cells_found
  ),
  paste0(
    "Exclusion-buffer distance: ",
    buffer_distance_m,
    " m"
  ),
  paste0(
    "Cells inside or touching buffer: ",
    sum(full_grid$Inside_50m_Buffer)
  ),
  paste0(
    "Pseudo-background candidate cells: ",
    nrow(candidate_pool)
  ),
  "",
  paste0(
    "Minimum candidate distance: ",
    round(
      min(candidate_pool$Min_Dist_to_Occurrence_m),
      3
    ),
    " m"
  ),
  paste0(
    "Maximum candidate distance: ",
    round(
      max(candidate_pool$Min_Dist_to_Occurrence_m),
      3
    ),
    " m"
  ),
  "",
  paste0(
    "Buffer violations: ",
    buffer_violation_n
  ),
  paste0(
    "Occurrence-background overlaps: ",
    occurrence_overlap_n
  ),
  paste0(
    "Duplicate candidate CellID values: ",
    duplicate_candidate_n
  ),
  "",
  "Missing predictor values in candidate pool:",
  paste(
    names(missing_predictor_counts),
    missing_predictor_counts,
    sep = " = ",
    collapse = "\n"
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
      "The candidate pool was created, but QA/QC failed.",
      "Project/QAQC/PseudoBackground_QAQC.txt",
      "Please inspect the corresponding QA/QC file."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================================\n")
cat("PSEUDO-BACKGROUND CANDIDATE POOL CREATED\n")
cat("============================================================\n")
cat(
  "Grid cells                    :",
  nrow(full_grid),
  "\n"
)
cat(
  "Number of occurrences         :",
  nrow(occurrences),
  "\n"
)
cat(
  "Buffer distance               :",
  buffer_distance_m,
  "m\n"
)
cat(
  "Cells inside/on buffer boundary:",
  sum(full_grid$Inside_50m_Buffer),
  "\n"
)
cat(
  "Background candidate pool     :",
  nrow(candidate_pool),
  "cell\n"
)
cat(
  "Minimum candidate distance    :",
  round(
    min(candidate_pool$Min_Dist_to_Occurrence_m),
    3
  ),
  "m\n"
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
cat("-", output_candidate_file, "\n")
cat("-", output_buffer_audit_file, "\n")
cat("-", output_summary_file, "\n")
cat("-", methodology_file, "\n")
cat("-", qaqc_file, "\n")
cat("============================================================\n")
