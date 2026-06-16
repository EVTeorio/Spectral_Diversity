library(testthat)

r_library_paths <- c(
  "C:/Users/PaintRock/AppData/Local/R/win-library/4.2",
  "C:/Program Files/R/R-4.2.3/library"
)
.libPaths(r_library_paths[dir.exists(r_library_paths)])

library(terra)
library(sf)
library(dplyr)

terraOptions(memfrac = 0.45, progress = 0)

start_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
repo_root <- if (basename(start_dir) == "tests") {
  normalizePath(file.path(start_dir, ".."), winslash = "/", mustWork = TRUE)
} else {
  start_dir
}
hsi_dir <- file.path(repo_root, "HSI_NA_trimmed")

get_hsi_files <- function() {
  all_files <- list.files(hsi_dir, full.names = TRUE)
  hsi_files <- all_files[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", all_files)]
  readable <- vapply(hsi_files, function(f) {
    !inherits(try(rast(f), silent = TRUE), "try-error")
  }, logical(1))

  skipped <- basename(hsi_files[!readable])
  if (length(skipped) > 0) {
    warning(
      "Skipping unreadable source HSI tiles, likely cloud-only OneDrive files: ",
      paste(skipped, collapse = ", "),
      call. = FALSE
    )
  }

  hsi_files[readable]
}

read_quads <- function(scale) {
  shp_path <- file.path(repo_root, "Quad_Scale_SHPs", paste0("PR_", scale, ".shp"))
  quads <- st_read(shp_path, quiet = TRUE)
  quads <- quads[, !grepl("^Dscrptn", names(quads))]

  if (scale == "10m") {
    quads <- quads %>%
      mutate(Name = sub_id) %>%
      select(-sub_id)
  }

  quads
}

pick_existing_quads <- function(scale, quads, hsi_files) {
  scale_dir <- file.path(repo_root, "Quad_Spectra", scale)
  quad_names <- as.character(st_drop_geometry(quads)$Name)
  existing_paths <- file.path(scale_dir, quad_names)
  existing <- quad_names[file.exists(existing_paths)]
  existing_paths <- file.path(scale_dir, existing)
  existing_readable <- vapply(existing_paths, function(p) {
    !inherits(try(rast(p), silent = TRUE), "try-error")
  }, logical(1))
  existing <- existing[existing_readable]

  if (length(existing) == 0) {
    stop("No readable existing quad spectra found for ", scale)
  }

  hsi_template <- rast(hsi_files[1])
  quads_vect <- project(vect(quads), hsi_template)
  covered <- rep(FALSE, length(quad_names))

  for (f in hsi_files) {
    r <- rast(f)
    covered <- covered | relate(quads_vect, as.polygons(ext(r), crs = crs(r)), "intersects")
  }

  covered_existing <- quad_names[covered & quad_names %in% existing]
  if (length(covered_existing) == 0) {
    stop("No existing ", scale, " quad intersects the readable HSI source rasters")
  }

  covered_existing
}

recreate_quad_crop <- function(scale, quad_id, out_dir, hsi_files) {
  quads <- read_quads(scale)
  quad_names <- as.character(st_drop_geometry(quads)$Name)
  quad_sf <- quads[quad_names == quad_id, ][1, ]
  quad_vect <- vect(quad_sf)

  hsi_template <- rast(hsi_files[1])
  quad_vect <- project(quad_vect, hsi_template)
  quad_template <- rast(
    ext(quad_vect),
    resolution = res(hsi_template),
    crs = crs(hsi_template)
  )

  spectral_stack <- list()
  distance_stack <- list()

  for (f in hsi_files) {
    r <- rast(f)
    e <- ext(r)

    if (!relate(quad_vect, as.polygons(e, crs = crs(r)), "intersects")) {
      next
    }

    r_quad <- crop(r, quad_vect)
    r_quad <- mask(r_quad, quad_vect)
    r_quad <- resample(r_quad, quad_template, method = "near")

    if (all(is.na(values(r_quad[[1]])))) {
      next
    }

    x_center <- (e$xmin + e$xmax) / 2
    y_center <- (e$ymin + e$ymax) / 2
    x_raster <- init(quad_template, "x")
    y_raster <- init(quad_template, "y")
    dist_r <- sqrt((x_raster - x_center)^2 + (y_raster - y_center)^2)

    no_data <- app(r_quad, function(x) ifelse(all(is.na(x)), 1, NA))
    dist_r[no_data == 1] <- Inf

    spectral_stack[[length(spectral_stack) + 1]] <- r_quad
    distance_stack[[length(distance_stack) + 1]] <- dist_r
  }

  if (length(spectral_stack) == 0) {
    stop("No intersecting HSI raster produced data for ", scale, " quad ", quad_id)
  }

  distances <- rast(distance_stack)
  winner <- app(distances, which.min)

  r_out <- spectral_stack[[1]]
  vals_out <- matrix(NA_real_, nrow = ncell(r_out), ncol = nlyr(r_out))

  for (k in seq_along(spectral_stack)) {
    idx <- which(values(winner) == k)
    if (length(idx) == 0) {
      next
    }
    vals_out[idx, ] <- values(spectral_stack[[k]])[idx, ]
  }

  values(r_out) <- vals_out

  out_path <- file.path(out_dir, paste0(scale, "_", quad_id))
  writeRaster(r_out, out_path, filetype = "ENVI", overwrite = TRUE)
  rast(out_path)
}

compare_quad <- function(scale, quad_id, generated) {
  existing_path <- file.path(repo_root, "Quad_Spectra", scale, quad_id)
  existing <- rast(existing_path)

  expect_true(
    compareGeom(existing, generated, stopOnError = FALSE, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE),
    info = paste("Geometry mismatch for", scale, quad_id)
  )
  expect_equal(nlyr(existing), nlyr(generated), info = paste("Band count mismatch for", scale, quad_id))

  layer_ids <- unique(round(seq(1, nlyr(existing), length.out = min(16, nlyr(existing)))))
  existing_layers <- existing[[layer_ids]]
  generated_layers <- generated[[layer_ids]]

  na_mismatch <- global(
    ifel(is.na(existing_layers) != is.na(generated_layers), 1, 0),
    "sum",
    na.rm = TRUE
  )
  max_na_mismatch <- max(na_mismatch[, 1], na.rm = TRUE)
  expect_equal(max_na_mismatch, 0, info = paste("NA pattern mismatch for", scale, quad_id))

  diff_raster <- abs(existing_layers - generated_layers)
  max_diff <- global(diff_raster, "max", na.rm = TRUE)
  max_abs_diff <- max(max_diff[, 1], na.rm = TRUE)
  expect_lt(max_abs_diff, 1e-6, info = paste("Value mismatch for", scale, quad_id))

  data.frame(
    scale = scale,
    quad_id = quad_id,
    layers_checked = length(layer_ids),
    max_abs_diff = max_abs_diff,
    max_na_mismatch = max_na_mismatch
  )
}

test_that("stored quadrat spectra match HSI_quad_crop_refined crop logic", {
  scales <- c("10m", "20m", "50m")
  out_dir <- file.path(repo_root, "tests", "_tmp_quad_crop_repro")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  hsi_files <- get_hsi_files()

  results <- lapply(scales, function(scale) {
    quads <- read_quads(scale)
    quad_ids <- pick_existing_quads(scale, quads, hsi_files)

    for (quad_id in quad_ids) {
      generated <- try(recreate_quad_crop(scale, quad_id, out_dir, hsi_files), silent = TRUE)
      if (!inherits(generated, "try-error")) {
        return(compare_quad(scale, quad_id, generated))
      }
    }

    stop("No readable ", scale, " candidate produced data from readable HSI source rasters")
  })

  results <- bind_rows(results)
  write.csv(
    results,
    file.path(repo_root, "reports", "validation", "quad_crop_reproducibility_results.csv"),
    row.names = FALSE
  )
})
