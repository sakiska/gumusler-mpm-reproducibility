# ============================================================
# Run_Variogram_OK_Pb_Zn_Cu.R
# Public reproducibility release for the Gumusler MPM study
#
# Purpose: Fit candidate variogram models, perform ordinary-kriging LOOCV, and select the final Pb, Zn, and Cu models.
#
# Notes:
# - Computational logic, model settings, thresholds, and random seeds are unchanged.
# - Human-facing comments and diagnostic messages are provided in English.
# - Run from the repository root so relative Project/... paths resolve correctly.
# ============================================================

required_packages <- c("sp", "gstat", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      "\nRun the following command first:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    )
  )
}
suppressPackageStartupMessages({
  library(sp)
  library(gstat)
  library(openxlsx)
})
# ============================================================
# ============================================================
elements <- c("Pb", "Zn", "Cu")
crs_string <- "+proj=utm +zone=35 +datum=WGS84 +units=m +no_defs"
lag_width_m <- 50
cutoff_m <- 800
initial_range_m <- 300
variogram_models <- c(
  Sph = "Spherical",
  Exp = "Exponential",
  Gau = "Gaussian"
)
# ============================================================
# 3. PATHS
# ============================================================
project_dir <- normalizePath("Project", mustWork = FALSE)
data_file <- file.path(
  project_dir,
  "Data",
  "Suppl_material.csv"
)
fig_dir <- file.path(project_dir, "Figures")
tab_dir <- file.path(project_dir, "Tables")
model_dir <- file.path(project_dir, "Models")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(data_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      data_file,
      "\n\nEnsure that the R working directory is the main MPM_Gumusler folder."
    )
  )
}
# ============================================================
# ============================================================
dat <- read.csv(
  data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c("SampleID", "X", "Y", elements)
if (!all(required_columns %in% names(dat))) {
  missing_columns <- setdiff(required_columns, names(dat))
  stop(
    paste0(
      "Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}
if (anyNA(dat[required_columns])) {
  stop("Missing values were found in required columns.")
}
for (element in elements) {
  if (!is.numeric(dat[[element]])) {
    stop(paste0(element, " column is not numeric."))
  }
  if (any(dat[[element]] <= 0)) {
    stop(
      paste0(
        element,
        " column contains zero or negative values. ",
        "The log10 transformation cannot be applied directly."
      )
    )
  }
  dat[[paste0(element, "_log10")]] <- log10(dat[[element]])
}
coordinates(dat) <- ~ X + Y
proj4string(dat) <- CRS(crs_string)
# ============================================================
# ============================================================
null_or <- function(x, replacement) {
  if (is.null(x) || length(x) == 0) replacement else x
}
fit_variogram_model <- function(empirical_variogram,
                                model_code,
                                observed_values) {
  empirical_max <- max(empirical_variogram$gamma, na.rm = TRUE)
  empirical_min <- min(empirical_variogram$gamma, na.rm = TRUE)
  sample_variance <- var(observed_values, na.rm = TRUE)
  sill_reference <- max(
    empirical_max,
    sample_variance,
    .Machine$double.eps
  )
  starting_ranges <- c(100, 200, 300, 400, 600)
  nugget_fractions <- c(0, 0.10, 0.25, 0.50)
  valid_fits <- list()
  fit_counter <- 0
  for (start_range in starting_ranges) {
    for (nugget_fraction in nugget_fractions) {
      start_nugget <- max(
        0,
        min(
          sill_reference * nugget_fraction,
          empirical_min
        )
      )
      start_partial_sill <- max(
        sill_reference - start_nugget,
        .Machine$double.eps
      )
      starting_model <- vgm(
        psill = start_partial_sill,
        model = model_code,
        range = start_range,
        nugget = start_nugget
      )
      candidate <- try(
        suppressWarnings(
          fit.variogram(
            empirical_variogram,
            starting_model,
            fit.method = 6
          )
        ),
        silent = TRUE
      )
      if (inherits(candidate, "try-error")) {
        next
      }
      structure_rows <- candidate$model != "Nug"
      valid_candidate <- (
        any(structure_rows) &&
        all(is.finite(candidate$psill)) &&
        all(is.finite(candidate$range)) &&
        all(candidate$psill >= 0) &&
        all(candidate$range[structure_rows] > 0)
      )
      if (!valid_candidate) {
        next
      }
      fit_counter <- fit_counter + 1
      valid_fits[[fit_counter]] <- candidate
    }
  }
  if (length(valid_fits) == 0) {
    return(NULL)
  }
  fit_errors <- vapply(
    valid_fits,
    function(x) {
      error_value <- attr(x, "SSErr")
      if (is.null(error_value) || !is.finite(error_value)) {
        Inf
      } else {
        error_value
      }
    },
    numeric(1)
  )
  valid_fits[[which.min(fit_errors)]]
}
extract_model_parameters <- function(fitted_model) {
  nugget_value <- fitted_model$psill[
    fitted_model$model == "Nug"
  ]
  structure_row <- fitted_model[
    fitted_model$model != "Nug",
    ,
    drop = FALSE
  ]
  nugget_value <- if (length(nugget_value) == 0) {
    0
  } else {
    nugget_value[1]
  }
  partial_sill <- structure_row$psill[1]
  range_value <- structure_row$range[1]
  data.frame(
    Nugget = nugget_value,
    Partial_Sill = partial_sill,
    Total_Sill = nugget_value + partial_sill,
    Range_m = range_value,
    Nugget_to_Total_Sill = ifelse(
      nugget_value + partial_sill > 0,
      nugget_value / (nugget_value + partial_sill),
      NA_real_
    ),
    RSS = null_or(attr(fitted_model, "SSErr"), NA_real_)
  )
}
calculate_cv_metrics <- function(formula_object,
                                 fitted_model,
                                 transformation) {
  cv <- krige.cv(
    formula_object,
    dat,
    model = fitted_model,
    nfold = nrow(dat),
    verbose = FALSE
  )
  observed_model_scale <- cv$observed
  predicted_model_scale <- cv$var1.pred
  kriging_variance <- pmax(cv$var1.var, 0)
  zscore <- cv$zscore
  if (transformation == "Raw") {
    observed_original <- observed_model_scale
    predicted_original <- predicted_model_scale
  } else {
    observed_original <- 10 ^ observed_model_scale
    predicted_original <- (
      10 ^ predicted_model_scale
    ) * exp(
      0.5 * (log(10) ^ 2) * kriging_variance
    )
  }
  valid <- is.finite(observed_original) &
    is.finite(predicted_original)
  observed_original <- observed_original[valid]
  predicted_original <- predicted_original[valid]
  residual_original <- predicted_original - observed_original
  r_squared <- if (
    length(observed_original) > 1 &&
    sd(observed_original) > 0 &&
    sd(predicted_original) > 0
  ) {
    cor(observed_original, predicted_original) ^ 2
  } else {
    NA_real_
  }
  metrics <- data.frame(
    ME_Original_ppm = mean(residual_original),
    MAE_Original_ppm = mean(abs(residual_original)),
    RMSE_Original_ppm = sqrt(mean(residual_original ^ 2)),
    R2_Original = r_squared,
    Mean_Standardized_Error = mean(zscore, na.rm = TRUE),
    RMS_Standardized_Error = sqrt(
      mean(zscore ^ 2, na.rm = TRUE)
    ),
    ME_Model_Scale = mean(
      predicted_model_scale - observed_model_scale,
      na.rm = TRUE
    ),
    RMSE_Model_Scale = sqrt(
      mean(
        (predicted_model_scale - observed_model_scale) ^ 2,
        na.rm = TRUE
      )
    )
  )
  list(
    metrics = metrics,
    cv = cv,
    observed_original = observed_original,
    predicted_original = predicted_original
  )
}
save_variogram_plot <- function(empirical_variogram,
                                fitted_models,
                                element,
                                transformation,
                                output_file) {
  png(
    filename = output_file,
    width = 2400,
    height = 1700,
    res = 300
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  plot(
    empirical_variogram$dist,
    empirical_variogram$gamma,
    pch = 16,
    cex = 0.9,
    xlab = "Lag distance (m)",
    ylab = "Semivariance",
    main = paste0(
      element,
      " ",
      transformation,
      ": Experimental Variogram and Model Fits"
    ),
    xlim = c(0, cutoff_m)
  )
  line_types <- c(Sph = 1, Exp = 2, Gau = 3)
  for (model_code in names(fitted_models)) {
    fitted_model <- fitted_models[[model_code]]
    if (is.null(fitted_model)) {
      next
    }
    model_line <- variogramLine(
      fitted_model,
      maxdist = cutoff_m,
      n = 400
    )
    lines(
      model_line$dist,
      model_line$gamma,
      lty = line_types[[model_code]],
      lwd = 2
    )
  }
  legend_labels <- variogram_models[
    names(fitted_models)[
      !vapply(fitted_models, is.null, logical(1))
    ]
  ]
  legend_lty <- line_types[
    names(fitted_models)[
      !vapply(fitted_models, is.null, logical(1))
    ]
  ]
  legend(
    "bottomright",
    legend = legend_labels,
    lty = legend_lty,
    lwd = 2,
    bty = "n"
  )
  grid()
}
save_cv_plot <- function(observed,
                         predicted,
                         element,
                         transformation,
                         model_name,
                         output_file) {
  png(
    filename = output_file,
    width = 1900,
    height = 1900,
    res = 300
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  lims <- range(
    c(observed, predicted),
    finite = TRUE
  )
  plot(
    observed,
    predicted,
    pch = 16,
    cex = 0.8,
    xlab = paste0("Observed ", element, " (ppm)"),
    ylab = paste0("Predicted ", element, " (ppm)"),
    main = paste0(
      element,
      " OK LOOCV: ",
      transformation,
      " + ",
      model_name
    ),
    xlim = lims,
    ylim = lims
  )
  abline(
    a = 0,
    b = 1,
    lty = 2,
    lwd = 2
  )
  grid()
}
run_element_analysis <- function(element) {
  message("\n========================================")
  message("Analysis starting: ", element)
  message("========================================")
  transformation_results <- list()
  for (transformation in c("Raw", "Log10")) {
    response_column <- if (transformation == "Raw") {
      element
    } else {
      paste0(element, "_log10")
    }
    formula_object <- as.formula(
      paste0(response_column, " ~ 1")
    )
    observed_values <- dat[[response_column]]
    empirical_variogram <- variogram(
      formula_object,
      dat,
      width = lag_width_m,
      cutoff = cutoff_m
    )
    fitted_models <- lapply(
      names(variogram_models),
      function(model_code) {
        fit_variogram_model(
          empirical_variogram,
          model_code,
          observed_values
        )
      }
    )
    names(fitted_models) <- names(variogram_models)
    comparison_rows <- list()
    cv_results <- list()
    for (model_code in names(fitted_models)) {
      fitted_model <- fitted_models[[model_code]]
      if (is.null(fitted_model)) {
        warning(
          paste0(
            element,
            " ",
            transformation,
            " ",
            variogram_models[[model_code]],
            " model could not be fitted and was skipped."
          )
        )
        next
      }
      parameters <- extract_model_parameters(
        fitted_model
      )
      cv_result <- calculate_cv_metrics(
        formula_object,
        fitted_model,
        transformation
      )
      cv_results[[model_code]] <- cv_result
      comparison_rows[[model_code]] <- cbind(
        Element = element,
        Transformation = transformation,
        Model = variogram_models[[model_code]],
        parameters,
        cv_result$metrics
      )
    }
    comparison_table <- do.call(
      rbind,
      comparison_rows
    )
    rownames(comparison_table) <- NULL
    write.csv(
      empirical_variogram,
      file.path(
        tab_dir,
        paste0(
          "Experimental_Variogram_",
          element,
          "_",
          transformation,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    save_variogram_plot(
      empirical_variogram,
      fitted_models,
      element,
      transformation,
      file.path(
        fig_dir,
        paste0(
          "Variogram_Model_Fits_",
          element,
          "_",
          transformation,
          ".png"
        )
      )
    )
    transformation_results[[transformation]] <- list(
      empirical = empirical_variogram,
      fits = fitted_models,
      comparison = comparison_table,
      cv = cv_results
    )
  }
  element_comparison <- rbind(
    transformation_results$Raw$comparison,
    transformation_results$Log10$comparison
  )
  element_comparison$Absolute_ME_Original_ppm <- abs(
    element_comparison$ME_Original_ppm
  )
  element_comparison <- element_comparison[
    order(
      element_comparison$RMSE_Original_ppm,
      element_comparison$Absolute_ME_Original_ppm,
      element_comparison$RSS
    ),
  ]
  rownames(element_comparison) <- NULL
  element_comparison$Rank <- seq_len(
    nrow(element_comparison)
  )
  element_comparison$Selected <- (
    element_comparison$Rank == 1
  )
  selected_row <- element_comparison[1, ]
  selected_transformation <- selected_row$Transformation
  selected_model_name <- selected_row$Model
  selected_model_code <- names(variogram_models)[
    variogram_models == selected_model_name
  ][1]
  selected_cv <- transformation_results[[selected_transformation]]$cv[[selected_model_code]]
  save_cv_plot(
    selected_cv$observed_original,
    selected_cv$predicted_original,
    element,
    selected_transformation,
    selected_model_name,
    file.path(
      fig_dir,
      paste0(
        "OK_LOOCV_",
        element,
        "_Selected_Model.png"
      )
    )
  )
  saveRDS(
    list(
      element = element,
      lag_width_m = lag_width_m,
      cutoff_m = cutoff_m,
      raw = transformation_results$Raw,
      log10 = transformation_results$Log10,
      comparison = element_comparison,
      selected = selected_row
    ),
    file.path(
      model_dir,
      paste0(
        "Variogram_OK_",
        element,
        "_Models.rds"
      )
    )
  )
  message(
    "Selected combination: ",
    element,
    " — ",
    selected_transformation,
    " + ",
    selected_model_name
  )
  element_comparison
}
# ============================================================
# 6. RUN Pb, Zn AND Cu
# ============================================================
all_results <- lapply(
  elements,
  run_element_analysis
)
names(all_results) <- elements
combined_comparison <- do.call(
  rbind,
  all_results
)
rownames(combined_comparison) <- NULL
# ============================================================
# 7. EXPORT ONE WORKBOOK
# ============================================================
workbook <- createWorkbook()
header_style <- createStyle(
  fgFill = "#1F4E78",
  fontColour = "#FFFFFF",
  textDecoration = "bold",
  halign = "center",
  border = "Bottom"
)
selected_style <- createStyle(
  textDecoration = "bold",
  border = "TopBottom"
)
addWorksheet(
  workbook,
  "All Model Comparisons"
)
writeData(
  workbook,
  "All Model Comparisons",
  combined_comparison
)
addStyle(
  workbook,
  "All Model Comparisons",
  header_style,
  rows = 1,
  cols = seq_len(ncol(combined_comparison)),
  gridExpand = TRUE
)
selected_rows <- which(
  combined_comparison$Selected
) + 1
if (length(selected_rows) > 0) {
  addStyle(
    workbook,
    "All Model Comparisons",
    selected_style,
    rows = selected_rows,
    cols = seq_len(ncol(combined_comparison)),
    gridExpand = TRUE,
    stack = TRUE
  )
}
freezePane(
  workbook,
  "All Model Comparisons",
  firstRow = TRUE
)
setColWidths(
  workbook,
  "All Model Comparisons",
  cols = seq_len(ncol(combined_comparison)),
  widths = "auto"
)
for (element in elements) {
  element_table <- all_results[[element]]
  sheet_name <- paste0(element, " Comparison")
  addWorksheet(
    workbook,
    sheet_name
  )
  writeData(
    workbook,
    sheet_name,
    element_table
  )
  addStyle(
    workbook,
    sheet_name,
    header_style,
    rows = 1,
    cols = seq_len(ncol(element_table)),
    gridExpand = TRUE
  )
  addStyle(
    workbook,
    sheet_name,
    selected_style,
    rows = 2,
    cols = seq_len(ncol(element_table)),
    gridExpand = TRUE,
    stack = TRUE
  )
  freezePane(
    workbook,
    sheet_name,
    firstRow = TRUE
  )
  setColWidths(
    workbook,
    sheet_name,
    cols = seq_len(ncol(element_table)),
    widths = "auto"
  )
}
saveWorkbook(
  workbook,
  file.path(
    tab_dir,
    "Variogram_and_OK_Model_Comparison_Pb_Zn_Cu.xlsx"
  ),
  overwrite = TRUE
)
# ============================================================
# 8. METADATA
# ============================================================
selected_summary <- do.call(
  rbind,
  lapply(elements, function(element) {
    selected <- all_results[[element]][
      all_results[[element]]$Selected,
      ,
      drop = FALSE
    ]
    data.frame(
      Element = element,
      Transformation = selected$Transformation,
      Model = selected$Model,
      RMSE_ppm = selected$RMSE_Original_ppm,
      ME_ppm = selected$ME_Original_ppm,
      Range_m = selected$Range_m,
      Nugget_to_Total_Sill = selected$Nugget_to_Total_Sill
    )
  })
)
metadata_lines <- c(
  "# Pb–Zn–Cu Variogram and Ordinary Kriging Model Comparison",
  "",
  paste0("- Lag width: ", lag_width_m, " m"),
  paste0("- Cutoff distance: ", cutoff_m, " m"),
  "- Direction: Omnidirectional",
  "- Transformations: Raw and log10",
  "- Candidate models: Spherical, Exponential and Gaussian",
  "- Validation: Leave-one-out Ordinary Kriging cross-validation",
  "- Primary selection criterion: Lowest RMSE on the original ppm scale",
  "- Secondary criteria: Lowest absolute ME and then lowest variogram RSS",
  "- Log10 predictions were back-transformed with an approximate lognormal bias correction.",
  "",
  "## Selected combinations",
  ""
)
for (i in seq_len(nrow(selected_summary))) {
  metadata_lines <- c(
    metadata_lines,
    paste0(
      "- ",
      selected_summary$Element[i],
      ": ",
      selected_summary$Transformation[i],
      " + ",
      selected_summary$Model[i],
      "; RMSE = ",
      round(selected_summary$RMSE_ppm[i], 3),
      " ppm; ME = ",
      round(selected_summary$ME_ppm[i], 3),
      " ppm; range = ",
      round(selected_summary$Range_m[i], 1),
      " m."
    )
  )
}
writeLines(
  metadata_lines,
  file.path(
    tab_dir,
    "Variogram_and_OK_Model_Comparison_Pb_Zn_Cu_Metadata.md"
  )
)
saveRDS(
  list(
    settings = list(
      elements = elements,
      crs = crs_string,
      lag_width_m = lag_width_m,
      cutoff_m = cutoff_m,
      initial_range_m = initial_range_m
    ),
    results = all_results,
    combined_comparison = combined_comparison,
    selected_summary = selected_summary
  ),
  file.path(
    model_dir,
    "Variogram_OK_Pb_Zn_Cu_All_Results.rds"
  )
)
message("\n========================================")
message("ALL ANALYSES COMPLETED")
message("========================================")
print(selected_summary)
