# ============================================================
# Run_OK_Rasters_Pb_Zn_Cu.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Generate Pb, Zn, and Cu ordinary-kriging prediction and variance rasters on the 25 m reference grid.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("sp", "gstat", "raster")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    )
  )
}
suppressPackageStartupMessages({
  library(sp)
  library(gstat)
  library(raster)
})
# ------------------------------------------------------------
# ------------------------------------------------------------
data_file <- file.path(
  "Project", "Data", "Suppl_material.csv"
)
template_file <- file.path(
  "Project", "Raster", "reference_grid_25m.tif"
)
model_file <- file.path(
  "Project", "Models",
  "Variogram_OK_Pb_Zn_Cu_All_Results.rds"
)
raster_output_dir <- file.path("Project", "Raster")
figure_output_dir <- file.path("Project", "Figures")
table_output_dir <- file.path("Project", "Tables")
dir.create(raster_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(data_file, template_file, model_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required project file(s) not found:\n",
      paste(paste0("- ", missing_files), collapse = "\n")
    )
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
elements <- c("Pb", "Zn", "Cu")
expected_models <- c(
  Pb = "Exponential",
  Zn = "Exponential",
  Cu = "Gaussian"
)
expected_transformations <- c(
  Pb = "Raw",
  Zn = "Raw",
  Cu = "Raw"
)
minimum_neighbors <- 8
maximum_neighbors <- 20
project_crs <- CRS("+init=epsg:32635")
# ------------------------------------------------------------
# 3. LOAD AND CHECK SAMPLE DATA
# ------------------------------------------------------------
samples <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c("SampleID", "X", "Y", elements)
missing_columns <- setdiff(required_columns, names(samples))
if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing data column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  )
}
if (anyNA(samples[, required_columns])) {
  stop("Missing values were found in the required data columns.")
}
for (element in elements) {
  if (!is.numeric(samples[[element]])) {
    stop(paste0(element, " must be numeric."))
  }
}
coordinates(samples) <- ~ X + Y
proj4string(samples) <- project_crs
# ------------------------------------------------------------
# ------------------------------------------------------------
reference_grid <- raster(template_file)
if (is.na(crs(reference_grid))) {
  crs(reference_grid) <- project_crs
}
grid_resolution <- res(reference_grid)
if (
  abs(grid_resolution[1] - 25) > 0.001 ||
  abs(grid_resolution[2] - 25) > 0.001
) {
  stop(
    paste0(
      "reference_grid_25m.tif is not a 25 m raster. ",
      "Detected resolution: ",
      paste(grid_resolution, collapse = " x "),
      " m."
    )
  )
}
if (!compareCRS(reference_grid, samples)) {
  stop(
    paste0(
      "CRS mismatch between sample data and reference grid.\n",
      "Expected CRS: EPSG:32635."
    )
  )
}
valid_cells <- which(!is.na(getValues(reference_grid)))
if (length(valid_cells) == 0) {
  stop(
    paste0(
      "reference_grid_25m.tif contains no valid cells. ",
      "The study area should have numeric values and the outside area NoData."
    )
  )
}
prediction_coordinates <- xyFromCell(reference_grid, valid_cells)
prediction_grid <- SpatialPointsDataFrame(
  coords = prediction_coordinates,
  data = data.frame(cell_id = valid_cells),
  proj4string = CRS(projection(reference_grid))
)
# ------------------------------------------------------------
# ------------------------------------------------------------
stored_results <- readRDS(model_file)
model_code_lookup <- c(
  Spherical = "Sph",
  Exponential = "Exp",
  Gaussian = "Gau"
)
get_selected_model <- function(element) {
  element_table <- stored_results$results[[element]]
  if (is.null(element_table)) {
    stop(
      paste0(
        "No stored variogram result was found for ",
        element,
        "."
      )
    )
  }
  if (!is.data.frame(element_table)) {
    stop(
      paste0(
        "Stored variogram results for ",
        element,
        " are not in data-frame format."
      )
    )
  }
  required_result_columns <- c(
    "Transformation",
    "Model",
    "Nugget",
    "Partial_Sill",
    "Range_m",
    "Selected"
  )
  missing_result_columns <- setdiff(
    required_result_columns,
    names(element_table)
  )
  if (length(missing_result_columns) > 0) {
    stop(
      paste0(
        element,
        ": required variogram result column(s) missing: ",
        paste(missing_result_columns, collapse = ", ")
      )
    )
  }
  selected_row <- element_table[
    element_table$Selected %in% TRUE,
    ,
    drop = FALSE
  ]
  if (nrow(selected_row) != 1) {
    stop(
      paste0(
        "Exactly one selected variogram model is required for ",
        element,
        "; found ",
        nrow(selected_row),
        "."
      )
    )
  }
  transformation <- as.character(
    selected_row$Transformation[1]
  )
  model_name <- as.character(
    selected_row$Model[1]
  )
  if (transformation != expected_transformations[[element]]) {
    stop(
      paste0(
        element,
        ": expected transformation = ",
        expected_transformations[[element]],
        ", stored selection = ",
        transformation,
        "."
      )
    )
  }
  if (model_name != expected_models[[element]]) {
    stop(
      paste0(
        element,
        ": expected model = ",
        expected_models[[element]],
        ", stored selection = ",
        model_name,
        "."
      )
    )
  }
  model_code <- unname(model_code_lookup[[model_name]])
  nugget_value <- as.numeric(selected_row$Nugget[1])
  partial_sill_value <- as.numeric(
    selected_row$Partial_Sill[1]
  )
  fitted_range <- as.numeric(selected_row$Range_m[1])
  if (
    !is.finite(nugget_value) ||
    !is.finite(partial_sill_value) ||
    !is.finite(fitted_range)
  ) {
    stop(
      paste0(
        "Non-finite variogram parameter detected for ",
        element,
        "."
      )
    )
  }
  if (
    nugget_value < 0 ||
    partial_sill_value <= 0 ||
    fitted_range <= 0
  ) {
    stop(
      paste0(
        "Invalid variogram parameter(s) for ",
        element,
        "."
      )
    )
  }
  fitted_model <- vgm(
    psill = partial_sill_value,
    model = model_code,
    range = fitted_range,
    nugget = nugget_value
  )
  list(
    transformation = transformation,
    model_name = model_name,
    fitted_model = fitted_model,
    fitted_range = fitted_range,
    selected_row = selected_row
  )
}
# ------------------------------------------------------------
# ------------------------------------------------------------
save_map_png <- function(
  raster_object,
  output_file,
  map_title,
  legend_title
) {
  png(
    filename = output_file,
    width = 1800,
    height = 1500,
    res = 200
  )
  par(
    mar = c(4.5, 4.5, 4.0, 6.5)
  )
  plot(
    raster_object,
    main = map_title,
    xlab = "Easting (m)",
    ylab = "Northing (m)",
    colNA = "transparent",
    legend.args = list(
      text = legend_title,
      side = 4,
      line = 3
    )
  )
  points(
    samples,
    pch = 20,
    cex = 0.35
  )
  box()
  dev.off()
}
# ------------------------------------------------------------
# ------------------------------------------------------------
summary_list <- list()
for (element in elements) {
  message("")
  message("============================================")
  message("Ordinary Kriging: ", element)
  message("============================================")
  selected <- get_selected_model(element)
  kriging_formula <- as.formula(
    paste0(element, " ~ 1")
  )
  kriging_output <- krige(
    formula = kriging_formula,
    locations = samples,
    newdata = prediction_grid,
    model = selected$fitted_model,
    nmin = minimum_neighbors,
    nmax = maximum_neighbors,
    debug.level = 0
  )
  predicted_values <- kriging_output$var1.pred
  variance_values <- kriging_output$var1.var
  prediction_raster <- raster(reference_grid)
  variance_raster <- raster(reference_grid)
  prediction_raster[] <- NA_real_
  variance_raster[] <- NA_real_
  prediction_raster[valid_cells] <- predicted_values
  variance_raster[valid_cells] <- variance_values
  prediction_raster_name <- paste0(element, "_OK.tif")
  variance_raster_name <- paste0(
    element,
    "_OK_Variance.tif"
  )
  prediction_raster_file <- file.path(
    raster_output_dir,
    prediction_raster_name
  )
  variance_raster_file <- file.path(
    raster_output_dir,
    variance_raster_name
  )
  writeRaster(
    prediction_raster,
    filename = prediction_raster_file,
    format = "GTiff",
    datatype = "FLT4S",
    overwrite = TRUE,
    options = c("COMPRESS=LZW")
  )
  writeRaster(
    variance_raster,
    filename = variance_raster_file,
    format = "GTiff",
    datatype = "FLT4S",
    overwrite = TRUE,
    options = c("COMPRESS=LZW")
  )
  prediction_png_file <- file.path(
    figure_output_dir,
    paste0(element, "_OK_Concentration_Surface.png")
  )
  variance_png_file <- file.path(
    figure_output_dir,
    paste0(element, "_OK_Kriging_Variance.png")
  )
  save_map_png(
    raster_object = prediction_raster,
    output_file = prediction_png_file,
    map_title = paste0(
      element,
      " Ordinary Kriging Concentration Surface"
    ),
    legend_title = paste0(element, " (ppm)")
  )
  save_map_png(
    raster_object = variance_raster,
    output_file = variance_png_file,
    map_title = paste0(
      element,
      " Ordinary Kriging Variance"
    ),
    legend_title = "Kriging variance"
  )
  finite_predictions <- is.finite(predicted_values)
  finite_variances <- is.finite(variance_values)
  summary_list[[element]] <- data.frame(
    Element = element,
    Transformation = selected$transformation,
    Variogram_Model = selected$model_name,
    Variogram_Range_m = selected$fitted_range,
    Search_Radius_m = NA_real_,
    Minimum_Neighbors = minimum_neighbors,
    Maximum_Neighbors = maximum_neighbors,
    Reference_Grid = basename(template_file),
    Grid_Resolution_m = grid_resolution[1],
    Study_Area_Cells = length(valid_cells),
    Predicted_Cells = sum(finite_predictions),
    Prediction_Coverage_Percent = (
      100 * mean(finite_predictions)
    ),
    NoData_Prediction_Cells = sum(!finite_predictions),
    Minimum_Prediction_ppm = if (
      any(finite_predictions)
    ) {
      min(predicted_values[finite_predictions])
    } else {
      NA_real_
    },
    Maximum_Prediction_ppm = if (
      any(finite_predictions)
    ) {
      max(predicted_values[finite_predictions])
    } else {
      NA_real_
    },
    Mean_Prediction_ppm = if (
      any(finite_predictions)
    ) {
      mean(predicted_values[finite_predictions])
    } else {
      NA_real_
    },
    Negative_Prediction_Cells = sum(
      predicted_values[finite_predictions] < 0
    ),
    Mean_Kriging_Variance = if (
      any(finite_variances)
    ) {
      mean(variance_values[finite_variances])
    } else {
      NA_real_
    },
    Prediction_Raster = prediction_raster_name,
    Variance_Raster = variance_raster_name,
    stringsAsFactors = FALSE
  )
  message(
    element,
    " completed | coverage = ",
    round(100 * mean(finite_predictions), 2),
    "% | nearest-neighbor search without a hard radius"
  )
}
# ------------------------------------------------------------
# 8. SAVE SUMMARY AND METADATA
# ------------------------------------------------------------
summary_table <- do.call(
  rbind,
  summary_list
)
rownames(summary_table) <- NULL
summary_file <- file.path(
  table_output_dir,
  "OK_Raster_Production_Summary.csv"
)
write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)
metadata_file <- file.path(
  table_output_dir,
  "OK_Raster_Production_Metadata.md"
)
metadata_text <- c(
  "# Ordinary Kriging Raster Production",
  "",
  "## Purpose",
  "",
  paste0(
    "Generate complete Pb, Zn and Cu continuous concentration/anomaly ",
    "surfaces and corresponding kriging variance maps."
  ),
  "",
  "## Fixed inputs",
  "",
  "- Sample data: `Project/Data/Suppl_material.csv`",
  paste0(
    "- Reference grid: `Project/Raster/",
    "reference_grid_25m.tif`"
  ),
  "- CRS: EPSG:32635",
  "- Cell size: 25 m",
  "",
  "## Selected models",
  "",
  "- Pb: Raw concentrations, Exponential variogram",
  "- Zn: Raw concentrations, Exponential variogram",
  "- Cu: Raw concentrations, Gaussian variogram",
  "",
  "## Neighborhood",
  "",
  "- Search radius: no hard maximum distance",
  "- Local neighborhood: nearest available samples",
  "- Minimum neighbors: 8",
  "- Maximum neighbors: 20",
  "",
  "## Outputs",
  "",
  "- Three Ordinary Kriging concentration GeoTIFFs",
  "- Three kriging variance GeoTIFFs",
  "- Six base-R quick-look PNG maps",
  "- One production summary CSV",
  "",
  "## Interpretation note",
  "",
  paste0(
    "The OK concentration rasters are continuous geochemical ",
    "surfaces. Formal anomaly classes or thresholds have not ",
    "yet been applied. Those will be produced in the subsequent ",
    "anomaly-extraction step."
  ),
  "",
  "## Quality-control note",
  "",
  paste0(
    "A hard maximum search distance is not imposed because it created ",
    "large internal NoData areas. The nearest 8-20 samples are used at ",
    "each valid reference-grid cell; uncertainty in weakly supported ",
    "areas is represented by the kriging-variance raster. Negative ",
    "predictions, if present, are retained and counted rather than ",
    "silently truncated."
  )
)
writeLines(
  metadata_text,
  metadata_file
)
message("")
message("============================================")
message("ALL Pb-Zn-Cu OK RASTERS COMPLETED")
message("============================================")
print(summary_table)
