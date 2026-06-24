USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c("terra", "sf", "dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
})

DTM_PATH <- "PRFPD_DTM_leafOff.tiff"
OUTPUT_DIR <- "Indices_SHPs/Enviro_SHPs"

SCALE_CONFIG <- list(
  "10m" = list(
    scale = "10m",
    quadrat_path = "Quad_Scale_SHPs/PR_10m.shp",
    id_col = "sub_id",
    output_stem = "enviro_variables_10m"
  ),
  "20m" = list(
    scale = "20m",
    quadrat_path = "Quad_Scale_SHPs/PR_20m.shp",
    id_col = "Name",
    output_stem = "enviro_variables_20m"
  ),
  "50m" = list(
    scale = "50m",
    quadrat_path = "Quad_Scale_SHPs/PR_50m.shp",
    id_col = "Name",
    output_stem = "enviro_variables_50m"
  )
)

find_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      file.exists(file.path(current_dir, DTM_PATH)) &&
        dir.exists(file.path(current_dir, "Quad_Scale_SHPs")) &&
        dir.exists(file.path(current_dir, "Indices_SHPs"))
    ) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find the Spectral_Diversity project root.", call. = FALSE)
    }
    current_dir <- parent_dir
  }
}

read_quadrats <- function(config, template_raster) {
  quads <- st_read(config$quadrat_path, quiet = TRUE) %>%
    dplyr::select(-matches("^Dscrptn")) %>%
    st_transform(st_crs(terra::crs(template_raster))) %>%
    mutate(
      quad_id = as.character(.data[[config$id_col]]),
      scale = config$scale
    )

  quads %>%
    relocate(quad_id, scale)
}

calculate_riley_tri_values <- function(values) {
  center_index <- ceiling(length(values) / 2)
  center_value <- values[[center_index]]

  if (!is.finite(center_value)) {
    return(NA_real_)
  }

  neighbor_values <- values[-center_index]
  differences <- neighbor_values - center_value
  differences <- differences[is.finite(differences)]

  if (length(differences) == 0) {
    return(NA_real_)
  }

  sqrt(sum(differences^2))
}

window_sum_matrix <- function(values, window_size) {
  n_rows <- nrow(values)
  n_cols <- ncol(values)
  half_window <- floor(window_size / 2)

  values[!is.finite(values)] <- 0

  integral <- matrix(0, nrow = n_rows + 1, ncol = n_cols + 1)
  col_cumulative <- apply(values, 2, cumsum)
  integral[-1, -1] <- t(apply(col_cumulative, 1, cumsum))

  row_min <- pmax(seq_len(n_rows) - half_window, 1)
  row_max <- pmin(seq_len(n_rows) + half_window, n_rows)
  col_min <- pmax(seq_len(n_cols) - half_window, 1)
  col_max <- pmin(seq_len(n_cols) + half_window, n_cols)

  result <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)
  for (row_index in seq_len(n_rows)) {
    result[row_index, ] <- integral[row_max[[row_index]] + 1, col_max + 1] -
      integral[row_min[[row_index]], col_max + 1] -
      integral[row_max[[row_index]] + 1, col_min] +
      integral[row_min[[row_index]], col_min]
  }

  result
}

matrix_to_raster <- function(values, template_raster) {
  output_raster <- template_raster
  terra::values(output_raster) <- as.vector(t(values))
  output_raster
}

calculate_riley_tri_raster <- function(dtm_raster, window_size) {
  if (window_size %% 2 != 1 || window_size < 3) {
    stop("Riley TRI window size must be an odd integer of at least 3.", call. = FALSE)
  }

  dtm_values <- terra::as.matrix(dtm_raster, wide = TRUE)
  finite_center <- is.finite(dtm_values)
  center_values <- dtm_values
  center_values[!finite_center] <- 0

  window_sum <- window_sum_matrix(dtm_values, window_size)
  window_sum_squares <- window_sum_matrix(dtm_values^2, window_size)
  window_count <- window_sum_matrix(matrix(as.numeric(finite_center), nrow(dtm_values)), window_size)

  neighbor_sum <- window_sum - center_values
  neighbor_sum_squares <- window_sum_squares - center_values^2
  neighbor_count <- window_count - as.numeric(finite_center)

  tri_squares <- neighbor_sum_squares -
    (2 * center_values * neighbor_sum) +
    (neighbor_count * center_values^2)
  tri_squares[tri_squares < 0] <- 0

  tri_values <- sqrt(tri_squares)
  tri_values[!finite_center | neighbor_count <= 0] <- NA_real_

  matrix_to_raster(tri_values, dtm_raster)
}

extract_mean_values <- function(raster_layer, quads, output_name) {
  extracted <- terra::extract(
    x = raster_layer,
    y = terra::vect(quads),
    fun = mean,
    na.rm = TRUE
  )

  tibble::tibble(
    quad_id = quads$quad_id,
    !!output_name := extracted[[2]]
  )
}

write_scale_outputs <- function(enviro_sf, config) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(OUTPUT_DIR, paste0(config$output_stem, ".csv"))
  shp_path <- file.path(OUTPUT_DIR, paste0(config$output_stem, ".shp"))

  readr::write_csv(st_drop_geometry(enviro_sf), csv_path)
  st_write(enviro_sf, shp_path, delete_layer = TRUE, quiet = TRUE)

  list(csv = csv_path, shapefile = shp_path)
}

process_scale <- function(config, dtm_raster, tri_rasters) {
  message("Processing environmental variables for ", config$scale, " quadrats")

  quads <- read_quadrats(config, dtm_raster)

  elevation_values <- extract_mean_values(
    raster_layer = dtm_raster,
    quads = quads,
    output_name = "elev_mean"
  )

  tri_values <- lapply(names(tri_rasters), function(metric_name) {
    extract_mean_values(
      raster_layer = tri_rasters[[metric_name]],
      quads = quads,
      output_name = metric_name
    )
  })

  enviro_values <- Reduce(
    function(left, right) left_join(left, right, by = "quad_id"),
    c(list(elevation_values), tri_values)
  )

  enviro_sf <- quads %>%
    left_join(enviro_values, by = "quad_id") %>%
    relocate(quad_id, scale, elev_mean, tri5_mean, tri11_mean, tri21_mean)

  outputs <- write_scale_outputs(enviro_sf, config)

  message(
    "Wrote ", config$scale, " outputs: ",
    outputs$csv, " and ", outputs$shapefile
  )

  enviro_sf
}

run_environmental_variables_workflow <- function(scales = names(SCALE_CONFIG)) {
  project_root <- find_project_root()
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  dtm_raster <- terra::rast(DTM_PATH)

  message("Calculating Riley TRI rasters")
  tri_rasters <- list(
    tri5_mean = calculate_riley_tri_raster(dtm_raster, 5),
    tri11_mean = calculate_riley_tri_raster(dtm_raster, 11),
    tri21_mean = calculate_riley_tri_raster(dtm_raster, 21)
  )

  results <- lapply(scales, function(scale_name) {
    process_scale(
      config = SCALE_CONFIG[[scale_name]],
      dtm_raster = dtm_raster,
      tri_rasters = tri_rasters
    )
  })
  names(results) <- scales

  invisible(results)
}

if (sys.nframe() == 0) {
  run_environmental_variables_workflow()
}
