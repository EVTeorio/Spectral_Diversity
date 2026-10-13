user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
OUTPUT_DIR <- "reports/tables/alpha_hull_all_pixels_values"
RESULTS_CSV <- file.path(OUTPUT_DIR, "alpha_hull_all_pixels_values.csv")
SUMMARY_CSV <- file.path(OUTPUT_DIR, "alpha_hull_all_pixels_values_summary.csv")
VALIDATION_CSV <- file.path(OUTPUT_DIR, "alpha_hull_all_pixels_values_validation.csv")
ANALYSIS_REPORT <- "reports/analysis/20261011_alpha_hull_all_pixels_values_50m_alpha1.md"
TASK_REPORT <- "reports/tasks/20261011_alpha_hull_all_pixels_values_50m_alpha1.md"
VALIDATION_REPORT <- "reports/validation/20261011_alpha_hull_all_pixels_values_50m_alpha1_validation.md"

PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds"
STANDARDIZED_PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds"

BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
PCA_N_COMPONENTS <- 3L

ALPHA_VALUES <- as.numeric(strsplit(Sys.getenv("ALPHA_VALUES", "1"), ",")[[1]])
SCALE_VALUES <- strsplit(Sys.getenv("ALPHA_SCALES", "50m"), ",")[[1]]
BASIS_VALUES <- strsplit(Sys.getenv("ALPHA_BASES", "vector_normalized_pca"), ",")[[1]]
ALPHA_TIMEOUT <- as.numeric(Sys.getenv("ALPHA_TIMEOUT", "60"))

alpha_label <- paste(gsub("\\.", "_", format(ALPHA_VALUES, trim = TRUE, scientific = FALSE)), collapse = "_")
scale_label <- paste(SCALE_VALUES, collapse = "_")
ANALYSIS_REPORT <- file.path(
  "reports/analysis",
  paste0("20261011_alpha_hull_all_pixels_values_", scale_label, "_alpha", alpha_label, ".md")
)
TASK_REPORT <- file.path(
  "reports/tasks",
  paste0("20261011_alpha_hull_all_pixels_values_", scale_label, "_alpha", alpha_label, ".md")
)
VALIDATION_REPORT <- file.path(
  "reports/validation",
  paste0("20261011_alpha_hull_all_pixels_values_", scale_label, "_alpha", alpha_label, "_validation.md")
)

SCALE_CONFIG <- data.frame(
  scale = c("10m", "20m", "50m"),
  spec_dir = c(
    "Quad_Spectra/10m_smooth_5nm",
    "Quad_Spectra/20m_smooth_5nm",
    "Quad_Spectra/50m_smooth_5nm"
  ),
  quadrat_csv = c(
    "quadrat_analysis_10m.csv",
    "quadrat_analysis_20m.csv",
    "quadrat_analysis_50m.csv"
  ),
  spectral_csv = c(
    "Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv",
    "Quad_Values/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv",
    "Quad_Values/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv"
  ),
  stringsAsFactors = FALSE
)

required_packages <- c("terra", "alphahull")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is unavailable: ", package, call. = FALSE)
  }
}

setwd(PROJECT_DIR)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

list_raster_files <- function(spec_dir) {
  files <- list.files(spec_dir, full.names = TRUE)
  files[!grepl("\\.(hdr|aux|xml|enp|sta)$", files, ignore.case = TRUE)]
}

get_quad_id <- function(file, scale) {
  quad_id <- basename(file)
  if (scale == "20m") {
    return(regmatches(quad_id, regexpr("^\\d+", quad_id)))
  }
  quad_id
}

read_masked_spectra <- function(file) {
  raster <- terra::rast(file)
  band_mask <- raster[[BEST_BAND]]
  sunlit_mask <- band_mask > SHADOW_THRESHOLD
  raster_masked <- terra::mask(raster, sunlit_mask, maskvalues = 0)
  x <- terra::values(raster_masked, mat = TRUE)
  if (nrow(x) == 0) {
    return(x)
  }
  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]
}

vector_normalize <- function(x) {
  norms <- sqrt(rowSums(x^2))
  keep <- is.finite(norms) & norms > 0
  x <- x[keep, , drop = FALSE]
  norms <- norms[keep]
  sweep(x, 1, norms, "/")
}

project_scores <- function(x, pca_object, preprocess) {
  if (preprocess == "vector_normalized") {
    x <- vector_normalize(x)
  }
  if (nrow(x) == 0) {
    return(matrix(numeric(), ncol = PCA_N_COMPONENTS))
  }
  x_scaled <- sweep(x, 2, pca_object$center, "-")
  x_scaled <- sweep(x_scaled, 2, pca_object$scale, "/")
  x_scaled %*% pca_object$pca$rotation[, seq_len(PCA_N_COMPONENTS), drop = FALSE]
}

calculate_alpha_area <- function(points, alpha_value) {
  points_unique <- unique(points)
  if (nrow(points_unique) < 4) {
    return(list(
      area = NA_real_,
      status = "insufficient_unique_points",
      n_unique_points = nrow(points_unique),
      arc_count = NA_integer_,
      elapsed_seconds = NA_real_
    ))
  }

  start_time <- proc.time()[["elapsed"]]
  tryCatch({
    hull <- alphahull::ahull(points_unique[, 1], points_unique[, 2], alpha = alpha_value)
    area <- alphahull::areaahull(hull, timeout = ALPHA_TIMEOUT)
    elapsed <- proc.time()[["elapsed"]] - start_time
    list(
      area = area,
      status = if (is.finite(area)) "ok" else "nonfinite_area",
      n_unique_points = nrow(points_unique),
      arc_count = if (!is.null(hull$arcs)) nrow(hull$arcs) else NA_integer_,
      elapsed_seconds = elapsed
    )
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - start_time
    list(
      area = NA_real_,
      status = paste0("error: ", conditionMessage(e)),
      n_unique_points = nrow(points_unique),
      arc_count = NA_integer_,
      elapsed_seconds = elapsed
    )
  })
}

read_existing_results <- function() {
  if (!file.exists(RESULTS_CSV)) {
    return(data.frame())
  }
  read.csv(RESULTS_CSV, stringsAsFactors = FALSE)
}

is_completed_combo <- function(scale, quad_id, basis, alpha_value) {
  existing <- read_existing_results()
  if (nrow(existing) == 0) {
    return(FALSE)
  }
  any(
    existing$scale == scale &
      existing$quad_id == quad_id &
      existing$basis == basis &
      abs(existing$alpha_value - alpha_value) < 1e-12
  )
}

append_result_row <- function(result_row) {
  write.table(
    result_row,
    RESULTS_CSV,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(RESULTS_CSV),
    append = file.exists(RESULTS_CSV)
  )
}

format_number <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

write_markdown_table <- function(data) {
  if (nrow(data) == 0) {
    return("_No rows._")
  }
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  rows <- apply(data, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  paste(c(header, separator, rows), collapse = "\n")
}

load_target_rows <- function(config) {
  spectral_df <- read.csv(config$spectral_csv, stringsAsFactors = FALSE)
  quadrat_df <- read.csv(config$quadrat_csv, stringsAsFactors = FALSE)
  raster_files <- list_raster_files(config$spec_dir)
  raster_index <- data.frame(
    quad_id = vapply(raster_files, get_quad_id, character(1), scale = config$scale),
    raster_file = raster_files,
    stringsAsFactors = FALSE
  )

  targets <- merge(
    spectral_df[, c(
      "quad_id",
      "scale",
      "n_pixels",
      "manual_excluded",
      "metric_method",
      "alpha_hull_area",
      "alpha_hull_method",
      "alpha_hull_n_points",
      "alpha_hull_alpha",
      "standardized_PCA_alpha_hull_area",
      "standardized_PCA_alpha_hull_method",
      "standardized_PCA_alpha_hull_n_points",
      "standardized_PCA_alpha_hull_alpha"
    )],
    raster_index,
    by = "quad_id",
    all.x = TRUE
  )

  if ("phy_afaith" %in% names(quadrat_df)) {
    targets <- merge(
      targets,
      quadrat_df[, c("quad_id", "phy_afaith")],
      by = "quad_id",
      all.x = TRUE
    )
  } else {
    targets$phy_afaith <- NA_real_
  }

  targets[order(targets$quad_id), , drop = FALSE]
}

pca_objects <- list(
  original_pca = list(
    pca_object = readRDS(PCA_RDS),
    preprocess = "raw",
    production_column = "alpha_hull_area"
  ),
  vector_normalized_pca = list(
    pca_object = readRDS(STANDARDIZED_PCA_RDS),
    preprocess = "vector_normalized",
    production_column = "standardized_PCA_alpha_hull_area"
  )
)

run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
target_configs <- SCALE_CONFIG[SCALE_CONFIG$scale %in% SCALE_VALUES, , drop = FALSE]
if (nrow(target_configs) == 0) {
  stop("No matching scale configuration for ALPHA_SCALES.", call. = FALSE)
}

for (config_i in seq_len(nrow(target_configs))) {
  config <- target_configs[config_i, ]
  targets <- load_target_rows(config)

  for (row_i in seq_len(nrow(targets))) {
    row <- targets[row_i, ]
    message("Processing ", row$scale, " ", row$quad_id, " (", row_i, "/", nrow(targets), ")")

    for (basis_name in intersect(names(pca_objects), BASIS_VALUES)) {
      for (alpha_value in ALPHA_VALUES) {
        if (is_completed_combo(row$scale, row$quad_id, basis_name, alpha_value)) {
          message("Skipping completed ", row$scale, " ", row$quad_id, " ", basis_name, " alpha=", alpha_value)
          next
        }

        if (basis_name == "original_pca") {
          previous_sampled_alpha_area <- row$alpha_hull_area
          previous_sampled_method <- row$alpha_hull_method
          previous_sampled_n_points <- row$alpha_hull_n_points
          previous_sampled_alpha <- row$alpha_hull_alpha
        } else {
          previous_sampled_alpha_area <- row$standardized_PCA_alpha_hull_area
          previous_sampled_method <- row$standardized_PCA_alpha_hull_method
          previous_sampled_n_points <- row$standardized_PCA_alpha_hull_n_points
          previous_sampled_alpha <- row$standardized_PCA_alpha_hull_alpha
        }

        base_row <- data.frame(
          scale = row$scale,
          quad_id = row$quad_id,
          basis = basis_name,
          alpha_value = alpha_value,
          manual_excluded = as.logical(row$manual_excluded),
          metric_method = row$metric_method,
          phy_afaith = row$phy_afaith,
          n_retained_pixels = NA_integer_,
          n_unique_pc12_points = NA_integer_,
          alpha_area_all_pixels = NA_real_,
          status = NA_character_,
          arc_count = NA_integer_,
          elapsed_seconds = NA_real_,
          raster_file = row$raster_file,
          previous_sampled_alpha_area = previous_sampled_alpha_area,
          previous_sampled_method = previous_sampled_method,
          previous_sampled_n_points = previous_sampled_n_points,
          previous_sampled_alpha = previous_sampled_alpha,
          run_timestamp = run_timestamp,
          stringsAsFactors = FALSE
        )

        if (isTRUE(as.logical(row$manual_excluded))) {
          base_row$status <- "skipped_manual_excluded"
          append_result_row(base_row)
          next
        }

        if (!is.character(row$raster_file) || is.na(row$raster_file) || !file.exists(row$raster_file)) {
          base_row$status <- "missing_raster"
          append_result_row(base_row)
          next
        }

        x <- read_masked_spectra(row$raster_file)
        scores <- project_scores(x, pca_objects[[basis_name]]$pca_object, pca_objects[[basis_name]]$preprocess)
        points <- scores[, 1:2, drop = FALSE]
        alpha_result <- calculate_alpha_area(points, alpha_value)

        base_row$n_retained_pixels <- nrow(x)
        base_row$n_unique_pc12_points <- alpha_result$n_unique_points
        base_row$alpha_area_all_pixels <- alpha_result$area
        base_row$status <- alpha_result$status
        base_row$arc_count <- alpha_result$arc_count
        base_row$elapsed_seconds <- alpha_result$elapsed_seconds

        append_result_row(base_row)
      }
    }
  }
}

results <- read_existing_results()
current_results <- results[
  results$scale %in% SCALE_VALUES &
    results$basis %in% BASIS_VALUES &
    results$alpha_value %in% ALPHA_VALUES,
  ,
  drop = FALSE
]

summary_rows <- aggregate(
  elapsed_seconds ~ scale + basis + alpha_value + status,
  data = current_results,
  FUN = function(x) round(sum(x, na.rm = TRUE), 3)
)
names(summary_rows)[names(summary_rows) == "elapsed_seconds"] <- "total_elapsed_seconds"
count_rows <- aggregate(
  quad_id ~ scale + basis + alpha_value + status,
  data = current_results,
  FUN = length
)
names(count_rows)[names(count_rows) == "quad_id"] <- "rows"
summary_rows <- merge(summary_rows, count_rows, by = c("scale", "basis", "alpha_value", "status"), all = TRUE)
summary_rows <- summary_rows[order(summary_rows$scale, summary_rows$basis, summary_rows$alpha_value, summary_rows$status), ]
write.csv(summary_rows, SUMMARY_CSV, row.names = FALSE)

validation <- data.frame(
  check = c(
    "requested_scale",
    "requested_alpha",
    "requested_basis",
    "current_rows",
    "ok_rows",
    "skipped_manual_excluded_rows",
    "non_ok_non_skipped_rows"
  ),
  value = c(
    paste(SCALE_VALUES, collapse = ","),
    paste(ALPHA_VALUES, collapse = ","),
    paste(BASIS_VALUES, collapse = ","),
    nrow(current_results),
    sum(current_results$status == "ok"),
    sum(current_results$status == "skipped_manual_excluded"),
    sum(!(current_results$status %in% c("ok", "skipped_manual_excluded")))
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, VALIDATION_CSV, row.names = FALSE)

report_summary <- summary_rows
report_summary$total_elapsed_seconds <- format_number(report_summary$total_elapsed_seconds)

analysis_lines <- c(
  "# All-Pixel Alpha-Hull Values: 50 m Alpha 1",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Purpose",
  "",
  "Calculate all-retained-pixel 2D alpha-hull spectral heterogeneity values for the current 50 m quadrats using alpha = 1.",
  "The output table is designed to accumulate future alpha values and spatial scales in the same CSV.",
  "",
  "## Scope",
  "",
  paste("- Scale(s):", paste(SCALE_VALUES, collapse = ", ")),
  paste("- Alpha value(s):", paste(ALPHA_VALUES, collapse = ", ")),
  paste("- PCA basis:", paste(BASIS_VALUES, collapse = ", ")),
  "- Metric space: vector-normalized PCA PC1-PC2",
  "- Manual atmospheric/cloud exclusions are retained as explicit skipped rows.",
  "",
  "## Summary",
  "",
  write_markdown_table(report_summary),
  "",
  "## Output Tables",
  "",
  paste0("- `", RESULTS_CSV, "`"),
  paste0("- `", SUMMARY_CSV, "`"),
  paste0("- `", VALIDATION_CSV, "`"),
  "",
  "## Notes",
  "",
  "- `alpha_area_all_pixels` is the all-retained-pixel PC1-PC2 alpha-hull area.",
  "- Duplicate PC1-PC2 coordinates are collapsed before alpha-hull geometry because duplicate points do not change the boundary.",
  "- `previous_sampled_alpha_area` carries the earlier 10,000-point sampled vector-normalized PCA alpha-hull value for comparison."
)
writeLines(analysis_lines, ANALYSIS_REPORT)

task_lines <- c(
  "# Task Report: All-Pixel Alpha-Hull Values",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Actions",
  "",
  "- Added `scripts/3_Analysis/alpha_hull_all_pixels_values.R`.",
  "- Created an append/resume CSV for all-pixel alpha-hull values across alpha values, scales, and PCA bases.",
  "- Ran the workflow for 50 m, alpha = 1, vector-normalized PCA.",
  "- Preserved manual exclusions as skipped rows rather than silently dropping them."
)
writeLines(task_lines, TASK_REPORT)

validation_lines <- c(
  "# Validation Report: All-Pixel Alpha-Hull Values",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Checks",
  "",
  write_markdown_table(validation),
  "",
  "## Result",
  "",
  if (validation$value[validation$check == "non_ok_non_skipped_rows"] == "0") {
    "Pass. Requested all-pixel alpha-hull rows are complete or explicitly skipped for manual exclusions."
  } else {
    "Review required. One or more requested rows has a non-ok status."
  }
)
writeLines(validation_lines, VALIDATION_REPORT)

message("Completed all-pixel alpha-hull value workflow.")
