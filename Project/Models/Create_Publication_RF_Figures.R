# ============================================================
# Create_Publication_RF_Figures.R
#
# Publication-quality figures for the Gümüşler MPM reproducibility repository
#
# IMPORTANT CORRECTION
# ------------------------------------------------------------
# This version deliberately avoids coord_sf() for the final figure
# axes. All vector geometries are converted to explicit EPSG:32635
# UTM coordinates and plotted with coord_equal().
#
# Therefore the axes are guaranteed to remain:
#   UTM Easting (m)
#   UTM Northing (m)
#
# Outputs:
#   1) Continuous RF ensemble probability map (fixed 0-1 scale)
#      + P90 boundary + geological overlays
#   2) Binary RF-P90 prospectivity map
#   3) RF ensemble prediction-variability (SD) map
#
# masterArea.shp is optional; if absent, the final raster extent is used.
# Raster values are NOT normalized or rescaled.
# Model grid remains 25 x 25 m.
# P90 threshold and model results are unchanged.
# CRS: EPSG:32635 - WGS 84 / UTM zone 35N
# ============================================================


# ------------------------------------------------------------
# 0. PACKAGE CONTROL
# ------------------------------------------------------------

required_packages <- c(
  "terra",
  "sf",
  "ggplot2",
  "scales"
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
      "\nInstall once with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- "Project"

continuous_file <- file.path(
  project_dir,
  "Raster",
  "Final_Prospectivity",
  "RF_Final_Continuous_Probability.tif"
)

binary_file <- file.path(
  project_dir,
  "Raster",
  "Final_Prospectivity",
  "RF_Final_P90_Prospectivity.tif"
)

uncertainty_file <- file.path(
  project_dir,
  "Raster",
  "Final_Prospectivity",
  "RF_Final_Probability_SD.tif"
)

threshold_file <- file.path(
  project_dir,
  "Models",
  "Final_Prospectivity",
  "RF_Final_P90_Threshold.csv"
)

vector_dir <- file.path(
  project_dir,
  "Vector"
)

output_dir <- file.path(
  project_dir,
  "Figures",
  "Final_Prospectivity",
  "Publication"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. REQUIRED INPUT CHECK
# ------------------------------------------------------------

required_files <- c(
  continuous_file,
  binary_file,
  uncertainty_file,
  threshold_file,
  file.path(vector_dir, "occurrences.shp"),
  file.path(vector_dir, "faults.shp"),
  file.path(vector_dir, "silicifiedZones.shp"),
  file.path(vector_dir, "brecciatedZones.shp")
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
# 3. READ RASTERS
# ------------------------------------------------------------

rf_probability <- terra::rast(
  continuous_file
)

rf_binary <- terra::rast(
  binary_file
)

rf_uncertainty <- terra::rast(
  uncertainty_file
)

threshold_summary <- read.csv(
  threshold_file,
  stringsAsFactors = FALSE
)

p90_threshold <- as.numeric(
  threshold_summary$Ensemble_P90_Threshold[1]
)


# ------------------------------------------------------------
# 4. CRS + RASTER QA/QC
# ------------------------------------------------------------

map_epsg <- 32635

if (
  !terra::compareGeom(
    rf_probability,
    rf_binary,
    rf_uncertainty,
    stopOnError = FALSE
  )
) {
  stop(
    "Final rasters do not have identical geometry.",
    call. = FALSE
  )
}

if (
  !terra::same.crs(
    rf_probability,
    rf_binary
  ) ||
    !terra::same.crs(
      rf_probability,
      rf_uncertainty
    )
) {
  stop(
    "Final rasters do not use the same CRS.",
    call. = FALSE
  )
}

raster_crs <- terra::crs(
  rf_probability,
  proj = TRUE
)

if (
  is.na(raster_crs) ||
    !grepl(
      "32635|utm.*zone=35|zone=35",
      raster_crs,
      ignore.case = TRUE
    )
) {
  stop(
    paste0(
      "Expected EPSG:32635 / UTM zone 35N raster CRS. Current CRS: ",
      raster_crs
    ),
    call. = FALSE
  )
}

if (
  terra::res(rf_probability)[1] != 25 ||
    terra::res(rf_probability)[2] != 25
) {
  stop(
    "Expected 25 x 25 m final raster resolution.",
    call. = FALSE
  )
}

prob_min <- terra::global(
  rf_probability,
  "min",
  na.rm = TRUE
)[1, 1]

prob_max <- terra::global(
  rf_probability,
  "max",
  na.rm = TRUE
)[1, 1]

if (
  prob_min < 0 ||
    prob_max > 1
) {
  stop(
    "RF probability values are outside 0-1.",
    call. = FALSE
  )
}

if (
  !is.finite(p90_threshold) ||
    p90_threshold < 0 ||
    p90_threshold > 1
) {
  stop(
    "Invalid P90 threshold.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 5. READ VECTOR OVERLAYS IN EPSG:32635
# ------------------------------------------------------------

read_sf_utm35 <- function(file_name) {

  x <- suppressWarnings(
    sf::st_read(
      file_name,
      quiet = TRUE
    )
  )

  if (is.na(sf::st_crs(x))) {
    sf::st_crs(x) <- map_epsg
  }

  sf::st_transform(
    x,
    map_epsg
  )
}


occurrences <- read_sf_utm35(
  file.path(vector_dir, "occurrences.shp")
)

faults <- read_sf_utm35(
  file.path(vector_dir, "faults.shp")
)

silicified <- read_sf_utm35(
  file.path(vector_dir, "silicifiedZones.shp")
)

brecciated <- read_sf_utm35(
  file.path(vector_dir, "brecciatedZones.shp")
)

master_area_file <- file.path(
  vector_dir,
  "masterArea.shp"
)

master_area <- NULL

if (file.exists(master_area_file)) {
  master_area <- read_sf_utm35(
    master_area_file
  )
}

lithology_file <- file.path(
  vector_dir,
  "lithology.shp"
)

lithology <- NULL

if (file.exists(lithology_file)) {
  lithology <- read_sf_utm35(
    lithology_file
  )
}


# ------------------------------------------------------------
# 6. RASTERS -> UTM XY DATA FRAMES
# ------------------------------------------------------------

prob_df <- as.data.frame(
  rf_probability,
  xy = TRUE,
  na.rm = TRUE
)

names(prob_df) <- c(
  "x",
  "y",
  "probability"
)


binary_df <- as.data.frame(
  rf_binary,
  xy = TRUE,
  na.rm = TRUE
)

names(binary_df) <- c(
  "x",
  "y",
  "prospective"
)

binary_df$prospective <- factor(
  binary_df$prospective,
  levels = c(0, 1),
  labels = c(
    "Below P90",
    "Prospective (>= P90)"
  )
)


uncertainty_df <- as.data.frame(
  rf_uncertainty,
  xy = TRUE,
  na.rm = TRUE
)

names(uncertainty_df) <- c(
  "x",
  "y",
  "sd_probability"
)


# ------------------------------------------------------------
# 7. P90 BOUNDARY
# ------------------------------------------------------------

p90_polygon <- terra::as.polygons(
  rf_binary,
  dissolve = TRUE,
  na.rm = TRUE
)

p90_polygon <- p90_polygon[
  p90_polygon[[1]] == 1,
]

p90_boundary_sf <- sf::st_as_sf(
  p90_polygon
)

p90_boundary_sf <- sf::st_transform(
  p90_boundary_sf,
  map_epsg
)

p90_boundary_sf <- sf::st_boundary(
  p90_boundary_sf
)


# ------------------------------------------------------------
# 8. SF GEOMETRY -> EXPLICIT UTM PATH DATA
#
# coord_equal() will be used instead of coord_sf(), so no
# geographic degree graticule conversion can occur.
# ------------------------------------------------------------

sf_lines_to_df <- function(
    x,
    polygon_boundary = FALSE
) {

  geom <- sf::st_geometry(x)

  if (polygon_boundary) {
    geom <- sf::st_boundary(geom)
  }

  coords <- sf::st_coordinates(
    geom
  )

  coords <- as.data.frame(
    coords
  )

  names(coords)[
    names(coords) == "X"
  ] <- "x"

  names(coords)[
    names(coords) == "Y"
  ] <- "y"

  level_columns <- grep(
    "^L[0-9]+$",
    names(coords),
    value = TRUE
  )

  if (length(level_columns) > 0) {
    coords$group <- interaction(
      coords[
        ,
        level_columns,
        drop = FALSE
      ],
      drop = TRUE
    )
  } else {
    coords$group <- 1
  }

  coords
}


faults_df <- sf_lines_to_df(
  faults,
  polygon_boundary = FALSE
)

silicified_df <- sf_lines_to_df(
  silicified,
  polygon_boundary = TRUE
)

brecciated_df <- sf_lines_to_df(
  brecciated,
  polygon_boundary = TRUE
)

master_df <- NULL

if (!is.null(master_area)) {
  master_df <- sf_lines_to_df(
    master_area,
    polygon_boundary = TRUE
  )
}

p90_boundary_df <- sf_lines_to_df(
  p90_boundary_sf,
  polygon_boundary = FALSE
)

lithology_df <- NULL

if (!is.null(lithology)) {
  lithology_df <- sf_lines_to_df(
    lithology,
    polygon_boundary = TRUE
  )
}

occurrence_xy <- sf::st_coordinates(
  occurrences
)

occurrences_df <- data.frame(
  x = occurrence_xy[, "X"],
  y = occurrence_xy[, "Y"]
)


# ------------------------------------------------------------
# 9. MAP EXTENT + UTM BREAKS
# ------------------------------------------------------------

if (!is.null(master_area)) {
  study_bbox <- sf::st_bbox(
    master_area
  )

  x_range <- as.numeric(
    c(
      study_bbox["xmin"],
      study_bbox["xmax"]
    )
  )

  y_range <- as.numeric(
    c(
      study_bbox["ymin"],
      study_bbox["ymax"]
    )
  )
} else {
  raster_extent <- terra::ext(rf_probability)
  x_range <- c(raster_extent$xmin, raster_extent$xmax)
  y_range <- c(raster_extent$ymin, raster_extent$ymax)
}

x_break_interval <- 500
y_break_interval <- 500

x_breaks <- seq(
  ceiling(x_range[1] / x_break_interval) *
    x_break_interval,
  floor(x_range[2] / x_break_interval) *
    x_break_interval,
  by = x_break_interval
)

y_breaks <- seq(
  ceiling(y_range[1] / y_break_interval) *
    y_break_interval,
  floor(y_range[2] / y_break_interval) *
    y_break_interval,
  by = y_break_interval
)


# ------------------------------------------------------------
# 10. SCALE BAR + NORTH ARROW GEOMETRY
# ------------------------------------------------------------

x_width <- diff(x_range)
y_height <- diff(y_range)

scale_bar_length <- 500
scale_segment_n <- 4
scale_segment_length <- scale_bar_length /
  scale_segment_n

scale_x0 <- x_range[1] +
  0.035 * x_width

scale_y0 <- y_range[1] +
  0.025 * y_height

scale_height <- 0.012 *
  y_height

scale_bar_df <- data.frame(
  xmin = scale_x0 +
    (0:(scale_segment_n - 1)) *
    scale_segment_length,
  xmax = scale_x0 +
    (1:scale_segment_n) *
    scale_segment_length,
  ymin = scale_y0,
  ymax = scale_y0 + scale_height,
  fill_key = rep(
    c("black", "white"),
    length.out = scale_segment_n
  )
)

north_x <- x_range[2] -
  0.055 * x_width

north_y0 <- y_range[2] -
  0.080 * y_height

north_y1 <- y_range[2] -
  0.025 * y_height


# ------------------------------------------------------------
# 11. COMMON THEME + AXES
# ------------------------------------------------------------

utm_axes <- list(

  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = scales::label_number(
      big.mark = ",",
      accuracy = 1
    ),
    expand = c(0, 0)
  ),

  ggplot2::scale_y_continuous(
    breaks = y_breaks,
    labels = scales::label_number(
      big.mark = ",",
      accuracy = 1
    ),
    expand = c(0, 0)
  ),

  ggplot2::coord_equal(
    xlim = x_range,
    ylim = y_range,
    expand = FALSE
  )
)


map_theme <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(

    panel.grid.major = ggplot2::element_line(
      linewidth = 0.20,
      colour = "grey85"
    ),

    panel.grid.minor = ggplot2::element_blank(),

    axis.title = ggplot2::element_text(
      size = 10
    ),

    axis.text = ggplot2::element_text(
      size = 8
    ),

    plot.title = ggplot2::element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),

    plot.subtitle = ggplot2::element_text(
      size = 9,
      hjust = 0
    ),

    legend.title = ggplot2::element_text(
      size = 9,
      face = "bold"
    ),

    legend.text = ggplot2::element_text(
      size = 8
    ),

    legend.key.width = grid::unit(
      1.2,
      "cm"
    ),

    legend.spacing.y = grid::unit(
      0.08,
      "cm"
    ),

    plot.margin = ggplot2::margin(
      7,
      7,
      7,
      7
    )
  )


# ------------------------------------------------------------
# 12. COMMON GEOLOGICAL OVERLAYS
# ------------------------------------------------------------

add_geological_overlays <- function(p) {

  # Geological overlays are mapped to a dedicated linetype legend.
  # Fixed colours preserve the intended cartographic symbology.
  if (!is.null(lithology_df)) {
    p <- p +
      ggplot2::geom_path(
        data = lithology_df,
        ggplot2::aes(
          x = x,
          y = y,
          group = group,
          linetype = "Lithological boundary"
        ),
        colour = "grey65",
        linewidth = 0.35,
        inherit.aes = FALSE,
        show.legend = TRUE
      )
  }

  p <- p +

    ggplot2::geom_path(
      data = silicified_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = group,
        linetype = "Silicified zone"
      ),
      colour = "#E69F00",
      linewidth = 0.9,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +

    ggplot2::geom_path(
      data = brecciated_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = group,
        linetype = "Brecciated zone"
      ),
      colour = "#CC79A7",
      linewidth = 0.9,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +

    ggplot2::geom_path(
      data = faults_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = group,
        linetype = "Fault"
      ),
      colour = "black",
      linewidth = 0.65,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +

    ggplot2::geom_point(
      data = occurrences_df,
      ggplot2::aes(
        x = x,
        y = y,
        shape = "Mineral occurrence"
      ),
      size = 3.0,
      stroke = 0.8,
      fill = "#D55E00",
      colour = "black",
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +

    ggplot2::scale_linetype_manual(
      name = "Geological overlays",
      breaks = c(
        "Fault",
        "Lithological boundary",
        "Silicified zone",
        "Brecciated zone"
      ),
      values = c(
        "Fault" = "solid",
        "Lithological boundary" = "solid",
        "Silicified zone" = "solid",
        "Brecciated zone" = "22"
      ),
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          colour = c(
            "black",
            "grey65",
            "#E69F00",
            "#CC79A7"
          ),
          linewidth = c(
            0.8,
            0.5,
            1.0,
            1.0
          )
        )
      )
    ) +

    ggplot2::scale_shape_manual(
      name = "Occurrences",
      values = c(
        "Mineral occurrence" = 23
      ),
      guide = ggplot2::guide_legend(
        order = 3,
        override.aes = list(
          fill = "#D55E00",
          colour = "black",
          size = 3
        )
      )
    )

  if (!is.null(master_df)) {
    p <- p + ggplot2::geom_path(
      data = master_df,
      ggplot2::aes(
        x = x,
        y = y,
        group = group
      ),
      colour = "black",
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  p
}




# ------------------------------------------------------------
# 13. FIGURE A - CONTINUOUS RF PROBABILITY
# ------------------------------------------------------------

p_continuous <- ggplot2::ggplot() +

  ggplot2::geom_raster(
    data = prob_df,
    ggplot2::aes(
      x = x,
      y = y,
      fill = probability
    )
  ) +

  ggplot2::scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = scales::label_number(
      accuracy = 0.1
    ),
    oob = scales::squish,
    name = "RF probability"
  ) +

  ggplot2::geom_path(
    data = p90_boundary_df,
    ggplot2::aes(
      x = x,
      y = y,
      group = group,
      colour = "P90 prospective boundary"
    ),
    linewidth = 1.0,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +

  ggplot2::scale_colour_manual(
    name = "Decision boundary",
    values = c(
      "P90 prospective boundary" = "#00D7D7"
    ),
    guide = ggplot2::guide_legend(
      order = 1,
      override.aes = list(
        linewidth = 1.2
      )
    )
  ) +

  ggplot2::labs(
    title = "Random Forest mineral prospectivity",
    subtitle = paste0(
      "Ensemble mean probability (150 runs); ",
      "P90 boundary = ",
      format(
        p90_threshold,
        digits = 4
      )
    ),
    x = "UTM Easting (m)",
    y = "UTM Northing (m)"
  ) +

  utm_axes +

  map_theme


p_continuous <- add_geological_overlays(
  p_continuous
)

# Add scale/north after all fill scales; these cartographic layers use
# explicit fill values without changing the RF probability legend.
p_continuous <- p_continuous +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[1],
    xmax = scale_bar_df$xmax[1],
    ymin = scale_bar_df$ymin[1],
    ymax = scale_bar_df$ymax[1],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[2],
    xmax = scale_bar_df$xmax[2],
    ymin = scale_bar_df$ymin[2],
    ymax = scale_bar_df$ymax[2],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[3],
    xmax = scale_bar_df$xmax[3],
    ymin = scale_bar_df$ymin[3],
    ymax = scale_bar_df$ymax[3],
    fill = "grey80",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[4],
    xmax = scale_bar_df$xmax[4],
    ymin = scale_bar_df$ymin[4],
    ymax = scale_bar_df$ymax[4],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "label",
    x = scale_x0 + scale_bar_length + 25,
    y = scale_y0 + scale_height / 2,
    label = "500 m",
    hjust = 0,
    vjust = 0.5,
    size = 3,
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.08, "lines")
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "black",
    linewidth = 1.5
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "white",
    linewidth = 0.8,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "cm"
      ),
      type = "closed"
    )
  ) +

  ggplot2::annotate(
    "label",
    x = north_x,
    y = north_y1 + 0.012 * y_height,
    label = "N",
    size = 4,
    fontface = "bold",
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.06, "lines")
  )


# ------------------------------------------------------------
# 14. FIGURE B - FINAL BINARY RF-P90 MAP
# ------------------------------------------------------------

p_binary <- ggplot2::ggplot() +

  ggplot2::geom_raster(
    data = binary_df,
    ggplot2::aes(
      x = x,
      y = y,
      fill = prospective
    )
  ) +

  ggplot2::scale_fill_manual(
    values = c(
      "Below P90" = "grey92",
      "Prospective (>= P90)" = "#D55E00"
    ),
    name = "Prospectivity"
  ) +

  ggplot2::labs(
    title = "Final RF-P90 mineral prospectivity map",
    subtitle = paste0(
      "P90 = ",
      format(
        p90_threshold,
        digits = 4
      ),
      "; prospective grid = 10%"
    ),
    x = "UTM Easting (m)",
    y = "UTM Northing (m)"
  ) +

  utm_axes +

  map_theme


p_binary <- add_geological_overlays(
  p_binary
)

p_binary <- p_binary +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[1],
    xmax = scale_bar_df$xmax[1],
    ymin = scale_bar_df$ymin[1],
    ymax = scale_bar_df$ymax[1],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[2],
    xmax = scale_bar_df$xmax[2],
    ymin = scale_bar_df$ymin[2],
    ymax = scale_bar_df$ymax[2],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[3],
    xmax = scale_bar_df$xmax[3],
    ymin = scale_bar_df$ymin[3],
    ymax = scale_bar_df$ymax[3],
    fill = "grey80",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[4],
    xmax = scale_bar_df$xmax[4],
    ymin = scale_bar_df$ymin[4],
    ymax = scale_bar_df$ymax[4],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "label",
    x = scale_x0 + scale_bar_length + 25,
    y = scale_y0 + scale_height / 2,
    label = "500 m",
    hjust = 0,
    vjust = 0.5,
    size = 3,
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.08, "lines")
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "black",
    linewidth = 1.5
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "white",
    linewidth = 0.8,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "cm"
      ),
      type = "closed"
    )
  ) +

  ggplot2::annotate(
    "label",
    x = north_x,
    y = north_y1 + 0.012 * y_height,
    label = "N",
    size = 4,
    fontface = "bold",
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.06, "lines")
  )


# ------------------------------------------------------------
# 15. FIGURE C - RF ENSEMBLE PREDICTION VARIABILITY
# ------------------------------------------------------------

uncertainty_max <- ceiling(
  max(
    uncertainty_df$sd_probability,
    na.rm = TRUE
  ) * 100
) / 100


p_uncertainty <- ggplot2::ggplot() +

  ggplot2::geom_raster(
    data = uncertainty_df,
    ggplot2::aes(
      x = x,
      y = y,
      fill = sd_probability
    )
  ) +

  ggplot2::scale_fill_viridis_c(
    option = "B",
    limits = c(
      0,
      uncertainty_max
    ),
    oob = scales::squish,
    name = "Probability SD"
  ) +

  ggplot2::labs(
    title = "Random Forest ensemble prediction variability",
    subtitle = "Cell-wise SD across 150 accepted RF predictions; higher SD = greater between-run variability",
    x = "UTM Easting (m)",
    y = "UTM Northing (m)"
  ) +

  utm_axes +

  map_theme


p_uncertainty <- add_geological_overlays(
  p_uncertainty
)

p_uncertainty <- p_uncertainty +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[1],
    xmax = scale_bar_df$xmax[1],
    ymin = scale_bar_df$ymin[1],
    ymax = scale_bar_df$ymax[1],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[2],
    xmax = scale_bar_df$xmax[2],
    ymin = scale_bar_df$ymin[2],
    ymax = scale_bar_df$ymax[2],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[3],
    xmax = scale_bar_df$xmax[3],
    ymin = scale_bar_df$ymin[3],
    ymax = scale_bar_df$ymax[3],
    fill = "grey80",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "rect",
    xmin = scale_bar_df$xmin[4],
    xmax = scale_bar_df$xmax[4],
    ymin = scale_bar_df$ymin[4],
    ymax = scale_bar_df$ymax[4],
    fill = "white",
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::annotate(
    "label",
    x = scale_x0 + scale_bar_length + 25,
    y = scale_y0 + scale_height / 2,
    label = "500 m",
    hjust = 0,
    vjust = 0.5,
    size = 3,
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.08, "lines")
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "black",
    linewidth = 1.5
  ) +

  ggplot2::annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    colour = "white",
    linewidth = 0.8,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "cm"
      ),
      type = "closed"
    )
  ) +

  ggplot2::annotate(
    "label",
    x = north_x,
    y = north_y1 + 0.012 * y_height,
    label = "N",
    size = 4,
    fontface = "bold",
    colour = "white",
    fill = "black",
    label.padding = grid::unit(0.06, "lines")
  )


# ------------------------------------------------------------
# 16. SAVE HIGH-RESOLUTION PUBLICATION FIGURES
# ------------------------------------------------------------

save_publication_map <- function(
    plot_object,
    base_name
) {

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        base_name,
        ".png"
      )
    ),
    plot = plot_object,
    width = 180,
    height = 220,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        base_name,
        ".tif"
      )
    ),
    plot = plot_object,
    width = 180,
    height = 220,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        base_name,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = 180,
    height = 220,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}


save_publication_map(
  p_continuous,
  "Figure_RF_Continuous_Probability"
)

save_publication_map(
  p_binary,
  "Figure_RF_P90_Prospectivity"
)

save_publication_map(
  p_uncertainty,
  "Figure_RF_Ensemble_Uncertainty"
)


# ------------------------------------------------------------
# 17. QA/QC REPORT
# ------------------------------------------------------------

output_files <- c(
  file.path(
    output_dir,
    "Figure_RF_Continuous_Probability.png"
  ),
  file.path(
    output_dir,
    "Figure_RF_Continuous_Probability.tif"
  ),
  file.path(
    output_dir,
    "Figure_RF_Continuous_Probability.pdf"
  ),
  file.path(
    output_dir,
    "Figure_RF_P90_Prospectivity.png"
  ),
  file.path(
    output_dir,
    "Figure_RF_P90_Prospectivity.tif"
  ),
  file.path(
    output_dir,
    "Figure_RF_P90_Prospectivity.pdf"
  ),
  file.path(
    output_dir,
    "Figure_RF_Ensemble_Uncertainty.png"
  ),
  file.path(
    output_dir,
    "Figure_RF_Ensemble_Uncertainty.tif"
  ),
  file.path(
    output_dir,
    "Figure_RF_Ensemble_Uncertainty.pdf"
  )
)


qaqc <- data.frame(
  Check = c(
    "Raster CRS is EPSG:32635 / UTM zone 35N",
    "All final rasters use same CRS",
    "Probability raster resolution is 25 m",
    "Probability values remain within 0-1",
    "Probability colour scale fixed at 0-1",
    "Axis system uses explicit projected UTM coordinates",
    "P90 threshold is finite",
    "Continuous map generated",
    "Binary P90 map generated",
    "Prediction-variability map generated",
    "All publication outputs written"
  ),
  Result = c(
    grepl(
      "32635|utm.*zone=35|zone=35",
      raster_crs,
      ignore.case = TRUE
    ),
    terra::same.crs(
      rf_probability,
      rf_binary
    ) &&
      terra::same.crs(
        rf_probability,
        rf_uncertainty
      ),
    all(
      terra::res(
        rf_probability
      ) == 25
    ),
    prob_min >= 0 &&
      prob_max <= 1,
    TRUE,
    all(
      prob_df$x > 500000 &
        prob_df$x < 600000
    ) &&
      all(
        prob_df$y > 4000000 &
          prob_df$y < 5000000
      ),
    is.finite(
      p90_threshold
    ),
    file.exists(
      output_files[1]
    ),
    file.exists(
      output_files[4]
    ),
    file.exists(
      output_files[7]
    ),
    all(
      file.exists(
        output_files
      )
    )
  ),
  stringsAsFactors = FALSE
)


cat("\n")
cat("============================================\n")
cat("PUBLICATION RF PROSPECTIVITY FIGURES COMPLETE\n")
cat("============================================\n")

cat(
  "Map CRS                 : EPSG:32635 - WGS 84 / UTM zone 35N\n"
)

cat(
  "Axis coordinates        : projected UTM metres\n"
)

cat(
  "Easting range           :",
  min(prob_df$x),
  "-",
  max(prob_df$x),
  "m\n"
)

cat(
  "Northing range          :",
  min(prob_df$y),
  "-",
  max(prob_df$y),
  "m\n"
)

cat(
  "Raster resolution       :",
  terra::res(
    rf_probability
  )[1],
  "m\n"
)

cat(
  "RF probability range    :",
  prob_min,
  "-",
  prob_max,
  "\n"
)

cat(
  "Displayed colour scale  : 0 - 1\n"
)

cat(
  "Final P90 threshold     :",
  p90_threshold,
  "\n"
)

cat(
  "Output directory        :",
  output_dir,
  "\n\n"
)

print(
  qaqc,
  row.names = FALSE
)

if (!all(qaqc$Result)) {
  stop(
    "Publication figure QA/QC failed.",
    call. = FALSE
  )
}

cat("\n")
cat("PUBLICATION FIGURES: PASS\n")
