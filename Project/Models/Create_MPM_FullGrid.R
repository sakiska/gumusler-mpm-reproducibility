# ============================================================
# Create_MPM_FullGrid.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Assemble the seven predictor layers on the common 25 m grid and assign known mineral occurrences to grid cells.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("terra", "sf")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
library(terra)
library(sf)
# ------------------------------------------------------------
# ------------------------------------------------------------
# setwd("C:/.../MPM_Gumusler")
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_file <- "Project/Raster/reference_grid_25m.tif"
raster_files <- c(
  Pb_OK             = "Project/Raster/Pb_OK.tif",
  Zn_OK             = "Project/Raster/Zn_OK.tif",
  Cu_OK             = "Project/Raster/Cu_OK.tif",
  Lithology         = "Project/Raster/lithology_25m.tif",
  Dist_Fault        = "Project/Raster/dist_fault_25m.tif",
  Dist_Silicified   = "Project/Raster/dist_silicified_25m.tif",
  Dist_Brecciated   = "Project/Raster/dist_brecciated_25m.tif"
)
occurrence_file <- "Project/Vector/occurrences.shp"
output_csv <- "Project/Data/MPM_FullGrid_25m.csv"
output_rds <- "Project/Data/MPM_FullGrid_25m.rds"
output_summary <- "Project/Tables/MPM_FullGrid_25m_Summary.csv"
output_occurrences <- "Project/Tables/Occurrence_Cell_Assignment.csv"
# ------------------------------------------------------------
# ------------------------------------------------------------
dir.create("Project/Data", recursive = TRUE, showWarnings = FALSE)
dir.create("Project/Tables", recursive = TRUE, showWarnings = FALSE)
# ------------------------------------------------------------
# ------------------------------------------------------------
all_input_files <- c(
  reference_file,
  unname(raster_files),
  occurrence_file
)
missing_files <- all_input_files[!file.exists(all_input_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following required files were not found:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}
cat("All required input files were found.\n")
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_grid <- rast(reference_file)
if (nlyr(reference_grid) > 1) {
  reference_grid <- reference_grid[[1]]
}
names(reference_grid) <- "Reference"
cat("\nReference grid information:\n")
print(reference_grid)
# ------------------------------------------------------------
# ------------------------------------------------------------
predictor_list <- lapply(
  names(raster_files),
  function(layer_name) {
    raster_object <- rast(raster_files[[layer_name]])
    if (nlyr(raster_object) > 1) {
      stop(
        paste0(
          layer_name,
          " rasteri birden fazla banda sahip. ",
          "Bu asamada her predictor tek bantli olmalidir."
        )
      )
    }
    names(raster_object) <- layer_name
    return(raster_object)
  }
)
names(predictor_list) <- names(raster_files)
cat("\nPredictor rasters were read successfully.\n")
# ------------------------------------------------------------
# ------------------------------------------------------------
check_grid_alignment <- function(reference, candidate, layer_name) {
  same_geometry <- terra::compareGeom(
    reference,
    candidate,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )
  if (!isTRUE(same_geometry)) {
    cat("\n--------------------------------------------\n")
    cat("GRID ALIGNMENT ERROR:", layer_name, "\n")
    cat("--------------------------------------------\n")
    cat("\nReference:\n")
    print(reference)
    cat("\nMisaligned raster:\n")
    print(candidate)
    stop(
      paste0(
        "\n",
        layer_name,
        " raster is not fully aligned with reference_grid_25m.tif.\n",
        "Script guvenlik amaciyla durduruldu.\n",
        "Raster otomatik olarak yeniden orneklenmedi."
      )
    )
  }
  cat(layer_name, ": grid alignment OK\n")
  return(invisible(TRUE))
}
for (layer_name in names(predictor_list)) {
  check_grid_alignment(
    reference = reference_grid,
    candidate = predictor_list[[layer_name]],
    layer_name = layer_name
  )
}
cat("\nAll rasters are fully aligned with the reference grid.\n")
# ------------------------------------------------------------
# ------------------------------------------------------------
predictor_stack <- rast(predictor_list)
cat("\nPredictor stack created:\n")
print(predictor_stack)
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_values <- values(reference_grid, mat = FALSE)
valid_reference_cells <- which(!is.na(reference_values))
if (length(valid_reference_cells) == 0) {
  stop("No valid cells were found within the reference grid.")
}
cat(
  "\nNumber of valid cells within the reference grid:",
  length(valid_reference_cells),
  "\n"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
predictor_values <- values(
  predictor_stack,
  mat = TRUE
)
predictor_values <- predictor_values[
  valid_reference_cells,
  ,
  drop = FALSE
]
coordinates <- xyFromCell(
  reference_grid,
  valid_reference_cells
)
mpm_table <- data.frame(
  CellID = valid_reference_cells,
  X = coordinates[, 1],
  Y = coordinates[, 2],
  predictor_values,
  check.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
na_counts <- colSums(
  is.na(
    mpm_table[
      ,
      names(raster_files),
      drop = FALSE
    ]
  )
)
cat("\nPredictor NA counts:\n")
print(na_counts)
if (any(na_counts > 0)) {
  warning(
    paste0(
      "Some predictor layers contain NA values within the reference grid.\n",
      "Cells containing NA values will be retained in the table.\n",
      "They should be evaluated separately before modeling."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
mpm_table$Lithology <- factor(mpm_table$Lithology)
cat("\nLithology factor levels:\n")
print(levels(mpm_table$Lithology))
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrences_sf <- st_read(
  occurrence_file,
  quiet = TRUE
)
if (nrow(occurrences_sf) == 0) {
  stop("No points were found in the occurrence layer.")
}
cat(
  "\nTotal records in the occurrence file:",
  nrow(occurrences_sf),
  "\n"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
geometry_types <- unique(
  as.character(st_geometry_type(occurrences_sf))
)
allowed_point_types <- c("POINT", "MULTIPOINT")
if (!all(geometry_types %in% allowed_point_types)) {
  stop(
    paste0(
      "The occurrence layer must contain only POINT or MULTIPOINT ",
      "geometrilerinden olusmalidir.\n",
      "Geometry types found: ",
      paste(geometry_types, collapse = ", ")
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_crs <- crs(reference_grid, proj = TRUE)
if (is.na(st_crs(occurrences_sf))) {
  stop(
    "No CRS is defined for the occurrence shapefile. ",
    "QGIS'te dogru CRS atanmalidir."
  )
}
occurrences_sf <- st_transform(
  occurrences_sf,
  crs = reference_crs
)
cat("\nThe occurrence layer was transformed to the reference-grid CRS.\n")
# ------------------------------------------------------------
# ------------------------------------------------------------
if ("MULTIPOINT" %in% geometry_types) {
  occurrences_sf <- st_cast(occurrences_sf, "POINT")
}
# ------------------------------------------------------------
# ------------------------------------------------------------
occurrences_vect <- vect(occurrences_sf)
occurrence_xy <- crds(occurrences_vect)
occurrence_cell_ids <- cellFromXY(
  reference_grid,
  occurrence_xy
)
occurrence_assignment <- data.frame(
  Occurrence_Row = seq_len(nrow(occurrence_xy)),
  X = occurrence_xy[, 1],
  Y = occurrence_xy[, 2],
  CellID = occurrence_cell_ids
)
# ------------------------------------------------------------
# ------------------------------------------------------------
outside_extent <- is.na(occurrence_assignment$CellID)
if (any(outside_extent)) {
  warning(
    paste0(
      sum(outside_extent),
      " occurrence pointsi reference grid extent'i outside kaldi."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
inside_extent_rows <- which(!outside_extent)
occurrence_assignment$Inside_Reference_Mask <- FALSE
if (length(inside_extent_rows) > 0) {
  cells_to_check <- occurrence_assignment$CellID[
    inside_extent_rows
  ]
  mask_values <- reference_values[cells_to_check]
  occurrence_assignment$Inside_Reference_Mask[
    inside_extent_rows
  ] <- !is.na(mask_values)
}
outside_mask <- !occurrence_assignment$Inside_Reference_Mask
if (any(outside_mask)) {
  warning(
    paste0(
      sum(outside_mask),
      " occurrence point(s) fall outside the valid reference-grid mask."
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
valid_occurrence_cells <- unique(
  occurrence_assignment$CellID[
    occurrence_assignment$Inside_Reference_Mask
  ]
)
valid_occurrence_cells <- valid_occurrence_cells[
  !is.na(valid_occurrence_cells)
]
cat(
  "\nNumber of unique cells containing valid occurrence points:",
  length(valid_occurrence_cells),
  "\n"
)
if (length(valid_occurrence_cells) == 0) {
  stop(
    "No occurrence point falls within a valid reference-grid cell."
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
mpm_table$Occurrence <- NA_integer_
mpm_table$Occurrence[
  mpm_table$CellID %in% valid_occurrence_cells
] <- 1L
# ------------------------------------------------------------
# ------------------------------------------------------------
unmatched_cells <- setdiff(
  valid_occurrence_cells,
  mpm_table$CellID
)
if (length(unmatched_cells) > 0) {
  stop(
    paste0(
      "Some occurrence cells were not found in the predictor table:\n",
      paste(unmatched_cells, collapse = ", ")
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
cell_frequency <- table(
  occurrence_assignment$CellID[
    occurrence_assignment$Inside_Reference_Mask
  ]
)
duplicate_occurrence_cells <- cell_frequency[
  cell_frequency > 1
]
if (length(duplicate_occurrence_cells) > 0) {
  cat(
    "\nMultiple occurrences were found within the same 25 m cell:\n"
  )
  print(duplicate_occurrence_cells)
  cat(
    "These cells were retained as a single positive cell in the modeling table.\n"
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
positive_count <- sum(
  mpm_table$Occurrence == 1,
  na.rm = TRUE
)
background_count <- sum(
  is.na(mpm_table$Occurrence)
)
cat("\n--------------------------------------------\n")
cat("FINAL MPM TABLE\n")
cat("--------------------------------------------\n")
cat("Total grid cells      :", nrow(mpm_table), "\n")
cat("Number of predictors  :", length(raster_files), "\n")
cat("Positive cells        :", positive_count, "\n")
cat("Unlabeled cells       :", background_count, "\n")
# ------------------------------------------------------------
# ------------------------------------------------------------
mpm_table <- mpm_table[
  ,
  c(
    "CellID",
    "X",
    "Y",
    "Pb_OK",
    "Zn_OK",
    "Cu_OK",
    "Lithology",
    "Dist_Fault",
    "Dist_Silicified",
    "Dist_Brecciated",
    "Occurrence"
  )
]
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  mpm_table,
  output_csv,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
saveRDS(
  mpm_table,
  output_rds
)
# ------------------------------------------------------------
# ------------------------------------------------------------
write.csv(
  occurrence_assignment,
  output_occurrences,
  row.names = FALSE,
  na = ""
)
# ------------------------------------------------------------
# ------------------------------------------------------------
summary_table <- data.frame(
  Item = c(
    "Total valid grid cells",
    "Number of predictors",
    "Input occurrence records",
    "Valid occurrence records",
    "Unique positive grid cells",
    "Unlabelled grid cells",
    "Cells with Pb NA",
    "Cells with Zn NA",
    "Cells with Cu NA",
    "Cells with Lithology NA",
    "Cells with Dist_Fault NA",
    "Cells with Dist_Silicified NA",
    "Cells with Dist_Brecciated NA"
  ),
  Value = c(
    nrow(mpm_table),
    length(raster_files),
    nrow(occurrences_sf),
    sum(occurrence_assignment$Inside_Reference_Mask),
    positive_count,
    background_count,
    na_counts["Pb_OK"],
    na_counts["Zn_OK"],
    na_counts["Cu_OK"],
    na_counts["Lithology"],
    na_counts["Dist_Fault"],
    na_counts["Dist_Silicified"],
    na_counts["Dist_Brecciated"]
  )
)
write.csv(
  summary_table,
  output_summary,
  row.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
cat("\n--------------------------------------------\n")
cat("PROCESS COMPLETED\n")
cat("--------------------------------------------\n")
cat("\nCreated files:\n")
cat(output_csv, "\n")
cat(output_rds, "\n")
cat(output_occurrences, "\n")
cat(output_summary, "\n")
cat(
  "\nNote: Cells without known occurrences were not assigned class 0; ",
  "they were left as NA.\n",
  sep = ""
)
