USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c("testthat", "terra")
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

library(testthat)
library(terra)

find_test_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    script_path <- file.path(
      current_dir,
      "scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R"
    )
    if (file.exists(script_path)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find the Spectral_Diversity project root.", call. = FALSE)
    }
    current_dir <- parent_dir
  }
}

project_root <- find_test_project_root()
source(file.path(
  project_root,
  "scripts/2_Indices Creation/Enviro_Variables/environmental_variables_all_scales.R"
))

test_that("Riley TRI helper calculates square root of summed squared differences", {
  values <- c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  )

  expected <- sqrt(sum((values[-5] - 5)^2))

  expect_equal(calculate_riley_tri_values(values), expected)
})

test_that("Riley TRI helper ignores missing neighbors but requires finite center", {
  values <- c(
    1, NA, 3,
    4, 5, 6,
    7, 8, 9
  )

  expected <- sqrt(sum((values[c(1, 3, 4, 6, 7, 8, 9)] - 5)^2))

  expect_equal(calculate_riley_tri_values(values), expected)
  expect_true(is.na(calculate_riley_tri_values(rep(NA_real_, 9))))
})

test_that("Riley TRI raster rejects invalid window sizes", {
  raster_layer <- rast(nrows = 3, ncols = 3, vals = 1:9)

  expect_error(
    calculate_riley_tri_raster(raster_layer, 4),
    "odd integer"
  )
  expect_error(
    calculate_riley_tri_raster(raster_layer, 1),
    "at least 3"
  )
})

test_that("Riley TRI raster matches helper for a center cell", {
  raster_layer <- rast(nrows = 3, ncols = 3, vals = 1:9)
  tri_layer <- calculate_riley_tri_raster(raster_layer, 3)

  expected <- calculate_riley_tri_values(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  ))

  expect_equal(terra::values(tri_layer)[[5]], expected)
})

test_that("rolling window sums match brute-force matrix sums", {
  values <- matrix(seq_len(25), nrow = 5, byrow = TRUE)
  window_sums <- window_sum_matrix(values, 3)

  brute_force_sum <- function(row_index, col_index) {
    rows <- max(1, row_index - 1):min(nrow(values), row_index + 1)
    cols <- max(1, col_index - 1):min(ncol(values), col_index + 1)
    sum(values[rows, cols])
  }

  expected <- matrix(NA_real_, nrow = 5, ncol = 5)
  for (row_index in seq_len(5)) {
    for (col_index in seq_len(5)) {
      expected[row_index, col_index] <- brute_force_sum(row_index, col_index)
    }
  }

  expect_equal(window_sums, expected)
})
