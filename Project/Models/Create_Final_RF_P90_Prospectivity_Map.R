# ============================================================
# Create_Final_RF_P90_Prospectivity_Map.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Create the final Random Forest ensemble prospectivity, predictive-variability, and P90 target rasters.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c(
  "terra",
  "sf"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------
project_dir <- "Project"
rf_summary_file <- file.path(
  project_dir,
  "Models",
  "RF_Results",
  "FullGrid",
  "RF_FullGrid_Ensemble_Summary.csv"
)
threshold_file <- file.path(
  project_dir,
  "Models",
  "Model_Comparison",
  "RF_FullGrid_Candidate_Thresholds.csv"
)
threshold_sensitivity_file <- file.path(
  project_dir,
  "Models",
  "Model_Comparison",
  "RF_Final_Threshold_Sensitivity.csv"
)
reference_grid_file <- file.path(
  project_dir,
  "Raster",
  "reference_grid_25m.tif"
)
vector_dir <- file.path(
  project_dir,
  "Vector"
)
model_output_dir <- file.path(
  project_dir,
  "Models",
  "Final_Prospectivity"
)
raster_output_dir <- file.path(
  project_dir,
  "Raster",
  "Final_Prospectivity"
)
figure_output_dir <- file.path(
  project_dir,
  "Figures",
  "Final_Prospectivity"
)
qaqc_dir <- file.path(
  project_dir,
  "QAQC"
)
dir.create(
  model_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  raster_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  qaqc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
final_csv_file <- file.path(
  model_output_dir,
  "RF_Final_P90_Prospectivity.csv"
)
final_threshold_file <- file.path(
  model_output_dir,
  "RF_Final_P90_Threshold.csv"
)
continuous_raster_file <- file.path(
  raster_output_dir,
  "RF_Final_Continuous_Probability.tif"
)
uncertainty_raster_file <- file.path(
  raster_output_dir,
  "RF_Final_Probability_SD.tif"
)
binary_raster_file <- file.path(
  raster_output_dir,
  "RF_Final_P90_Prospectivity.tif"
)
continuous_figure_file <- file.path(
  figure_output_dir,
  "RF_Final_Continuous_Prospectivity.png"
)
binary_figure_file <- file.path(
  figure_output_dir,
  "RF_Final_P90_Prospectivity.png"
)
qaqc_file <- file.path(
  qaqc_dir,
  "RF_Final_P90_Prospectivity_QAQC.txt"
)
methodology_file <- file.path(
  qaqc_dir,
  "RF_Final_P90_Prospectivity_Methodology.txt"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
required_files <- c(
  rf_summary_file,
  threshold_file,
  threshold_sensitivity_file,
  reference_grid_file
)
missing_files <- required_files[
  !file.exists(required_files)
]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required file(s) missing:\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# 4. READ MODEL OUTPUTS
# ------------------------------------------------------------
rf_summary <- read.csv(
  rf_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
threshold_table <- read.csv(
  threshold_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
threshold_sensitivity <- read.csv(
  threshold_sensitivity_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_rf_columns <- c(
  "CellID",
  "X",
  "Y",
  "RF_Mean_Probability",
  "RF_SD_Probability",
  "RF_Min_Probability",
  "RF_Max_Probability",
  "RF_Prediction_n"
)
missing_rf_columns <- setdiff(
  required_rf_columns,
  names(rf_summary)
)
if (length(missing_rf_columns) > 0) {
  stop(
    paste0(
      "RF summary missing column(s): ",
      paste(missing_rf_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
p90_row <- threshold_table[
  threshold_table$Threshold_Method == "P90",
  ,
  drop = FALSE
]
if (nrow(p90_row) != 1L) {
  stop(
    "Exactly one P90 threshold row was expected.",
    call. = FALSE
  )
}
final_p90_threshold <- as.numeric(
  p90_row$Threshold[1]
)
if (
  !is.finite(final_p90_threshold) ||
    final_p90_threshold < 0 ||
    final_p90_threshold > 1
) {
  stop(
    "Invalid final P90 threshold.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
p90_sensitivity <- threshold_sensitivity[
  threshold_sensitivity$Threshold_Method == "P90",
  ,
  drop = FALSE
]
if (nrow(p90_sensitivity) != 1L) {
  stop(
    "P90 row not found in final threshold-sensitivity table.",
    call. = FALSE
  )
}
if (
  abs(
    as.numeric(p90_sensitivity$Ensemble_Threshold[1]) -
      final_p90_threshold
  ) > 1e-12
) {
  stop(
    "P90 threshold differs between candidate-threshold and sensitivity files.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
rf_summary$Final_Model <- "RF"
rf_summary$Final_Threshold_Method <- "P90"
rf_summary$Final_Threshold <- final_p90_threshold
rf_summary$Prospective_P90 <- ifelse(
  rf_summary$RF_Mean_Probability >= final_p90_threshold,
  1L,
  0L
)
rf_summary$Prospectivity_Class <- ifelse(
  rf_summary$Prospective_P90 == 1L,
  "Prospective_P90",
  "Below_P90"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
if (nrow(rf_summary) != 9600L) {
  stop(
    paste0(
      "Expected 9600 grid cells; found ",
      nrow(rf_summary),
      "."
    ),
    call. = FALSE
  )
}
if (anyDuplicated(rf_summary$CellID) > 0) {
  stop(
    "Duplicated CellID values in RF ensemble summary.",
    call. = FALSE
  )
}
if (
  anyNA(
    rf_summary[
      ,
      c(
        "CellID",
        "X",
        "Y",
        "RF_Mean_Probability",
        "RF_SD_Probability",
        "Prospective_P90"
      )
    ]
  )
) {
  stop(
    "Missing required final-grid values.",
    call. = FALSE
  )
}
if (
  any(
    rf_summary$RF_Mean_Probability < 0 |
      rf_summary$RF_Mean_Probability > 1
  )
) {
  stop(
    "Invalid RF mean probabilities.",
    call. = FALSE
  )
}
if (
  any(
    rf_summary$RF_SD_Probability < 0
  )
) {
  stop(
    "Invalid RF probability SD values.",
    call. = FALSE
  )
}
if (
  !all(
    rf_summary$RF_Prediction_n == 150L
  )
) {
  stop(
    "Not every final grid cell has 150 RF predictions.",
    call. = FALSE
  )
}
prospective_cell_n <- sum(
  rf_summary$Prospective_P90 == 1L
)
prospective_percent <- 100 *
  prospective_cell_n /
  nrow(rf_summary)
if (prospective_cell_n != 960L) {
  warning(
    paste0(
      "P90 produced ",
      prospective_cell_n,
      " prospective cells instead of exactly 960. ",
      "This can occur only if tied probabilities exist at the threshold."
    ),
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_grid <- terra::rast(
  reference_grid_file
)
if (
  terra::ncell(reference_grid) <
    nrow(rf_summary)
) {
  stop(
    "Reference raster has fewer cells than the RF summary.",
    call. = FALSE
  )
}
xy <- as.matrix(
  rf_summary[
    ,
    c("X", "Y")
  ]
)
raster_cells <- terra::cellFromXY(
  reference_grid,
  xy
)
if (
  anyNA(raster_cells)
) {
  stop(
    "At least one RF grid coordinate falls outside the reference raster.",
    call. = FALSE
  )
}
if (
  anyDuplicated(raster_cells) > 0
) {
  stop(
    "Multiple RF rows map to the same reference-raster cell.",
    call. = FALSE
  )
}
continuous_raster <- reference_grid
continuous_raster[] <- NA_real_
continuous_raster[raster_cells] <-
  rf_summary$RF_Mean_Probability
names(continuous_raster) <- "RF_Mean_Probability"
uncertainty_raster <- reference_grid
uncertainty_raster[] <- NA_real_
uncertainty_raster[raster_cells] <-
  rf_summary$RF_SD_Probability
names(uncertainty_raster) <- "RF_Probability_SD"
binary_raster <- reference_grid
binary_raster[] <- NA_real_
binary_raster[raster_cells] <-
  rf_summary$Prospective_P90
names(binary_raster) <- "RF_P90_Prospective"
# ------------------------------------------------------------
# ------------------------------------------------------------
terra::writeRaster(
  continuous_raster,
  continuous_raster_file,
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW")
)
terra::writeRaster(
  uncertainty_raster,
  uncertainty_raster_file,
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW")
)
terra::writeRaster(
  binary_raster,
  binary_raster_file,
  overwrite = TRUE,
  datatype = "INT1U",
  NAflag = 255,
  gdal = c("COMPRESS=LZW")
)
# ------------------------------------------------------------
# ------------------------------------------------------------
find_vector_file <- function(
    vector_dir,
    patterns
) {
  if (!dir.exists(vector_dir)) {
    return(NA_character_)
  }
  shp_files <- list.files(
    vector_dir,
    pattern = "\\.shp$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(shp_files) == 0L) {
    return(NA_character_)
  }
  base_names <- basename(shp_files)
  matches <- rep(FALSE, length(shp_files))
  for (pattern_i in patterns) {
    matches <- matches |
      grepl(
        pattern_i,
        base_names,
        ignore.case = TRUE
      )
  }
  result <- shp_files[matches]
  if (length(result) == 0L) {
    NA_character_
  } else {
    result[1]
  }
}
occurrence_file <- find_vector_file(
  vector_dir,
  c(
    "^occurrences\\.shp$",
    "occurrence",
    "mineral"
  )
)
fault_file <- find_vector_file(
  vector_dir,
  c(
    "fault",
    "fay"
  )
)
silicified_file <- find_vector_file(
  vector_dir,
  c(
    "silic",
    "silis"
  )
)
brecciated_file <- find_vector_file(
  vector_dir,
  c(
    "brecc",
    "bres",
    "bres"
  )
)
read_optional_vector <- function(
    file_name
) {
  if (
    length(file_name) != 1L ||
      is.na(file_name) ||
      !file.exists(file_name)
  ) {
    return(NULL)
  }
  x <- suppressWarnings(
    sf::st_read(
      file_name,
      quiet = TRUE
    )
  )
  if (
    is.na(sf::st_crs(x))
  ) {
    sf::st_crs(x) <- 32635
  }
  x <- sf::st_transform(
    x,
    32635
  )
  terra::vect(x)
}
occurrences <- read_optional_vector(
  occurrence_file
)
faults <- read_optional_vector(
  fault_file
)
silicified <- read_optional_vector(
  silicified_file
)
brecciated <- read_optional_vector(
  brecciated_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
grDevices::png(
  filename = continuous_figure_file,
  width = 2200,
  height = 1800,
  res = 220
)
graphics::par(
  mar = c(4.5, 4.8, 4.2, 5.5)
)
terra::plot(
  continuous_raster,
  main = "Random Forest Mineral Prospectivity",
  xlab = "Easting (m)",
  ylab = "Northing (m)",
  axes = TRUE,
  plg = list(
    title = "RF probability"
  )
)
if (!is.null(silicified)) {
  terra::plot(
    silicified,
    add = TRUE,
    border = "black",
    lwd = 1.3
  )
}
if (!is.null(brecciated)) {
  terra::plot(
    brecciated,
    add = TRUE,
    border = "black",
    lwd = 1.3,
    lty = 2
  )
}
if (!is.null(faults)) {
  terra::plot(
    faults,
    add = TRUE,
    col = "black",
    lwd = 1.2
  )
}
if (!is.null(occurrences)) {
  terra::plot(
    occurrences,
    add = TRUE,
    pch = 8,
    cex = 1.5,
    col = "black"
  )
}
graphics::mtext(
  paste0(
    "Final model: RF ensemble (150 runs); operational threshold: P90 = ",
    format(
      final_p90_threshold,
      digits = 4
    )
  ),
  side = 3,
  line = 0.3,
  cex = 0.8
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
grDevices::png(
  filename = binary_figure_file,
  width = 2200,
  height = 1800,
  res = 220
)
graphics::par(
  mar = c(4.5, 4.8, 4.2, 5.5)
)
terra::plot(
  binary_raster,
  main = "Final RF-P90 Mineral Prospectivity Map",
  xlab = "Easting (m)",
  ylab = "Northing (m)",
  axes = TRUE,
  type = "classes",
  breaks = c(-0.5, 0.5, 1.5),
  col = c("grey90", "grey35"),
  legend = FALSE
)
graphics::legend(
  "topright",
  legend = c(
    "Below P90",
    "Prospective (>= P90)"
  ),
  fill = c(
    "grey90",
    "grey35"
  ),
  title = "Prospectivity",
  bty = "o",
  cex = 0.85
)
if (!is.null(silicified)) {
  terra::plot(
    silicified,
    add = TRUE,
    border = "black",
    lwd = 1.5
  )
}
if (!is.null(brecciated)) {
  terra::plot(
    brecciated,
    add = TRUE,
    border = "black",
    lwd = 1.5,
    lty = 2
  )
}
if (!is.null(faults)) {
  terra::plot(
    faults,
    add = TRUE,
    col = "black",
    lwd = 1.3
  )
}
if (!is.null(occurrences)) {
  terra::plot(
    occurrences,
    add = TRUE,
    pch = 8,
    cex = 1.6,
    col = "black"
  )
}
graphics::mtext(
  paste0(
    "P90 threshold = ",
    format(final_p90_threshold, digits = 4),
    "; prospective grid = ",
    format(prospective_percent, digits = 4),
    "%"
  ),
  side = 3,
  line = 0.3,
  cex = 0.8
)
grDevices::dev.off()
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  rf_summary,
  final_csv_file,
  row.names = FALSE,
  na = ""
)
final_threshold_summary <- data.frame(
  Final_Model = "RF",
  Final_Threshold_Method = "P90",
  Ensemble_P90_Threshold = final_p90_threshold,
  Total_Grid_Cell_n = nrow(rf_summary),
  Prospective_Cell_n = prospective_cell_n,
  Prospective_Grid_Percent = prospective_percent,
  Independent_Validation_TP =
    as.integer(
      p90_sensitivity$TP[1]
    ),
  Independent_Validation_FN =
    as.integer(
      p90_sensitivity$FN[1]
    ),
  Independent_Validation_TN =
    as.integer(
      p90_sensitivity$TN[1]
    ),
  Independent_Validation_FP =
    as.integer(
      p90_sensitivity$FP[1]
    ),
  Independent_Validation_Recall =
    as.numeric(
      p90_sensitivity$Recall[1]
    ),
  Independent_Validation_Specificity =
    as.numeric(
      p90_sensitivity$Specificity[1]
    ),
  Independent_Validation_Precision =
    as.numeric(
      p90_sensitivity$Precision[1]
    ),
  Independent_Validation_F1 =
    as.numeric(
      p90_sensitivity$F1[1]
    ),
  Independent_Validation_Balanced_Accuracy =
    as.numeric(
      p90_sensitivity$Balanced_Accuracy[1]
    ),
  stringsAsFactors = FALSE
)
write.csv(
  final_threshold_summary,
  final_threshold_file,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc <- data.frame(
  Check = c(
    "RF ensemble contains 9600 grid cells",
    "Every grid cell contains 150 RF predictions",
    "P90 threshold is finite and within 0-1",
    "P90 threshold matches final sensitivity table",
    "All RF probabilities are valid",
    "All RF SD values are non-negative",
    "No duplicated final CellID",
    "No duplicated raster-cell assignment",
    "Continuous raster written",
    "Uncertainty raster written",
    "Binary P90 raster written",
    "Continuous figure written",
    "Final P90 figure written",
    "Final CSV written",
    "Final threshold summary written"
  ),
  Result = c(
    nrow(rf_summary) == 9600L,
    all(rf_summary$RF_Prediction_n == 150L),
    is.finite(final_p90_threshold) &&
      final_p90_threshold >= 0 &&
      final_p90_threshold <= 1,
    abs(
      as.numeric(
        p90_sensitivity$Ensemble_Threshold[1]
      ) -
        final_p90_threshold
    ) <= 1e-12,
    all(
      is.finite(
        rf_summary$RF_Mean_Probability
      ) &
        rf_summary$RF_Mean_Probability >= 0 &
        rf_summary$RF_Mean_Probability <= 1
    ),
    all(
      is.finite(
        rf_summary$RF_SD_Probability
      ) &
        rf_summary$RF_SD_Probability >= 0
    ),
    anyDuplicated(
      rf_summary$CellID
    ) == 0L,
    anyDuplicated(
      raster_cells
    ) == 0L,
    file.exists(
      continuous_raster_file
    ),
    file.exists(
      uncertainty_raster_file
    ),
    file.exists(
      binary_raster_file
    ),
    file.exists(
      continuous_figure_file
    ),
    file.exists(
      binary_figure_file
    ),
    file.exists(
      final_csv_file
    ),
    file.exists(
      final_threshold_file
    )
  ),
  stringsAsFactors = FALSE
)
if (!all(qaqc$Result)) {
  print(qaqc)
  stop(
    "Final RF-P90 prospectivity QA/QC failed.",
    call. = FALSE
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
qaqc_lines <- c(
  "FINAL RF-P90 PROSPECTIVITY QA/QC",
  "============================================",
  "",
  paste(
    qaqc$Check,
    qaqc$Result,
    sep = ": "
  ),
  "",
  paste0(
    "Final RF ensemble P90 threshold: ",
    final_p90_threshold
  ),
  paste0(
    "Prospective cells: ",
    prospective_cell_n,
    " / ",
    nrow(rf_summary)
  ),
  paste0(
    "Prospective grid percent: ",
    prospective_percent
  ),
  "",
  "FINAL STATUS: PASS"
)
writeLines(
  qaqc_lines,
  qaqc_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
vector_note <- c(
  paste0(
    "Occurrence vector: ",
    ifelse(
      is.na(occurrence_file),
      "not found",
      occurrence_file
    )
  ),
  paste0(
    "Fault vector: ",
    ifelse(
      is.na(fault_file),
      "not found",
      fault_file
    )
  ),
  paste0(
    "Silicified vector: ",
    ifelse(
      is.na(silicified_file),
      "not found",
      silicified_file
    )
  ),
  paste0(
    "Brecciated vector: ",
    ifelse(
      is.na(brecciated_file),
      "not found",
      brecciated_file
    )
  )
)
methodology_lines <- c(
  "FINAL RF-P90 MINERAL PROSPECTIVITY MAP",
  "============================================",
  "",
  paste(
    "The final prospectivity surface is the cell-wise mean",
    "occurrence probability from 150 accepted leakage-free",
    "Random Forest runs."
  ),
  "",
  paste(
    "Random Forest was selected over XGBoost and SVM using",
    "independent occurrence-versus-pseudo-background validation",
    "and threshold-free and threshold-dependent performance metrics."
  ),
  "",
  paste(
    "The operational threshold was selected as P90 after",
    "sensitivity analysis against P95, P99, Mean+2SD and",
    "scaled Median+2MAD thresholds."
  ),
  "",
  paste0(
    "Final ensemble P90 threshold: ",
    final_p90_threshold,
    "."
  ),
  "",
  paste0(
    "Prospective grid fraction: ",
    prospective_cell_n,
    " of ",
    nrow(rf_summary),
    " cells (",
    prospective_percent,
    "%)."
  ),
  "",
  paste(
    "Known occurrences and available geological vector layers",
    "are displayed only as map overlays and do not alter the",
    "final RF probability or P90 classification."
  ),
  "",
  vector_note
)
writeLines(
  methodology_lines,
  methodology_file
)
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n")
cat("============================================\n")
cat("FINAL RF-P90 PROSPECTIVITY MAP COMPLETE\n")
cat("============================================\n")
cat(
  "Final model                : RF\n"
)
cat(
  "Final threshold method     : P90\n"
)
cat(
  "Final P90 threshold        :",
  final_p90_threshold,
  "\n"
)
cat(
  "Total grid cells           :",
  nrow(rf_summary),
  "\n"
)
cat(
  "Prospective cells          :",
  prospective_cell_n,
  "\n"
)
cat(
  "Prospective grid percent   :",
  prospective_percent,
  "\n"
)
cat("\n")
cat("Detected vector overlays:\n")
cat(
  "  Occurrences :",
  ifelse(is.na(occurrence_file), "NOT FOUND", occurrence_file),
  "\n"
)
cat(
  "  Faults      :",
  ifelse(is.na(fault_file), "NOT FOUND", fault_file),
  "\n"
)
cat(
  "  Silicified  :",
  ifelse(is.na(silicified_file), "NOT FOUND", silicified_file),
  "\n"
)
cat(
  "  Brecciated  :",
  ifelse(is.na(brecciated_file), "NOT FOUND", brecciated_file),
  "\n"
)
cat("\n")
print(
  qaqc,
  row.names = FALSE
)
cat("\n")
cat(
  "Continuous map :",
  continuous_figure_file,
  "\n"
)
cat(
  "Final P90 map  :",
  binary_figure_file,
  "\n"
)
cat(
  "Final GeoTIFF  :",
  binary_raster_file,
  "\n"
)
cat("\n")
cat("FINAL RF-P90 PROSPECTIVITY: PASS\n")
