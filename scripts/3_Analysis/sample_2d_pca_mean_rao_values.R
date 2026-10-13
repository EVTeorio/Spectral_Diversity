user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
OUTPUT_DIR <- "reports/tables/sample_2d_pca_mean_rao_values"
RESULTS_CSV <- file.path(OUTPUT_DIR, "sample_2d_pca_mean_rao_values.csv")
SUMMARY_CSV <- file.path(OUTPUT_DIR, "sample_2d_pca_mean_rao_values_summary.csv")
CORRELATION_CSV <- file.path(OUTPUT_DIR, "sample_2d_pca_mean_rao_phy_afaith_correlations.csv")
VALIDATION_CSV <- file.path(OUTPUT_DIR, "sample_2d_pca_mean_rao_validation.csv")
ANALYSIS_REPORT <- "reports/analysis/20261011_sample_2d_pca_mean_rao_50m.md"
TASK_REPORT <- "reports/tasks/20261011_sample_2d_pca_mean_rao_50m.md"
VALIDATION_REPORT <- "reports/validation/20261011_sample_2d_pca_mean_rao_50m_validation.md"

STANDARDIZED_PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds"

BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
PCA_N_COMPONENTS <- 3L
SAMPLE_MAX_POINTS <- as.integer(Sys.getenv("SAMPLED_2D_MAX_POINTS", "10000"))
RANDOM_SEED <- as.integer(Sys.getenv("SAMPLED_2D_RANDOM_SEED", "42"))
SCALE_VALUES <- strsplit(Sys.getenv("SAMPLED_2D_SCALES", "50m"), ",")[[1]]
BASIS_VALUES <- strsplit(Sys.getenv("SAMPLED_2D_BASES", "vector_normalized_pca"), ",")[[1]]

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

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Required R package is unavailable: terra", call. = FALSE)
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

seed_from_values <- function(scale, quad_id, offset = 0) {
  RANDOM_SEED + offset + sum(utf8ToInt(paste(scale, quad_id, sep = "_")))
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

sample_scores_2d <- function(scores, scale, quad_id) {
  scores_2d <- scores[, 1:2, drop = FALSE]
  n_available <- nrow(scores_2d)
  if (n_available > SAMPLE_MAX_POINTS) {
    set.seed(seed_from_values(scale, quad_id, offset = 2000))
    sample_index <- sample.int(n_available, SAMPLE_MAX_POINTS)
    scores_2d <- scores_2d[sample_index, , drop = FALSE]
    method <- "sampled_pixels"
  } else {
    method <- "all_pixels"
  }

  list(points = scores_2d, method = method, n_available = n_available, n_used = nrow(scores_2d))
}

calculate_2d_metrics <- function(points) {
  if (nrow(points) < 2) {
    return(list(
      mean_distance_2d = NA_real_,
      rao_q_2d = NA_real_,
      status = "insufficient_points"
    ))
  }

  center <- colMeans(points)
  centered <- sweep(points, 2, center, "-")
  squared_radius <- rowSums(centered^2)
  distances <- sqrt(squared_radius)

  list(
    mean_distance_2d = mean(distances),
    rao_q_2d = 2 * mean(squared_radius),
    status = "ok"
  )
}

read_existing_results <- function() {
  if (!file.exists(RESULTS_CSV)) {
    return(data.frame())
  }
  read.csv(RESULTS_CSV, stringsAsFactors = FALSE)
}

is_completed_combo <- function(scale, quad_id, basis, sample_max_points) {
  existing <- read_existing_results()
  if (nrow(existing) == 0) {
    return(FALSE)
  }
  any(
    existing$scale == scale &
      existing$quad_id == quad_id &
      existing$basis == basis &
      existing$sample_max_points == sample_max_points
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

write_markdown_table <- function(data) {
  if (nrow(data) == 0) {
    return("_No rows._")
  }
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  rows <- apply(data, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  paste(c(header, separator, rows), collapse = "\n")
}

format_number <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), digits = digits, format = "f"))
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
      "standardized_PCA_pca_euclidean_mean",
      "standardized_PCA_rao_q_pca"
    )],
    raster_index,
    by = "quad_id",
    all.x = TRUE
  )
  targets <- merge(targets, quadrat_df[, c("quad_id", "phy_afaith")], by = "quad_id", all.x = TRUE)
  targets[order(targets$quad_id), , drop = FALSE]
}

correlation_rows <- function(results) {
  ok <- results[
    results$status == "ok" &
      is.finite(results$mean_distance_2d) &
      is.finite(results$rao_q_2d) &
      is.finite(results$phy_afaith),
    ,
    drop = FALSE
  ]
  metrics <- c("mean_distance_2d", "rao_q_2d")
  labels <- c("sampled_2d_mean_distance", "sampled_2d_rao_q")

  do.call(rbind, lapply(seq_along(metrics), function(i) {
    values <- ok[[metrics[i]]]
    n <- length(values)
    pearson_r <- NA_real_
    pearson_p <- NA_real_
    spearman_r <- NA_real_
    spearman_p <- NA_real_
    if (n >= 3 && stats::sd(values) > 0 && stats::sd(ok$phy_afaith) > 0) {
      pearson_test <- stats::cor.test(values, ok$phy_afaith, method = "pearson")
      spearman_test <- suppressWarnings(stats::cor.test(values, ok$phy_afaith, method = "spearman", exact = FALSE))
      pearson_r <- unname(pearson_test$estimate)
      pearson_p <- pearson_test$p.value
      spearman_r <- unname(spearman_test$estimate)
      spearman_p <- spearman_test$p.value
    }
    data.frame(
      scale = paste(unique(ok$scale), collapse = ","),
      basis = paste(unique(ok$basis), collapse = ","),
      sample_max_points = SAMPLE_MAX_POINTS,
      metric_id = labels[i],
      n = n,
      pearson_r = pearson_r,
      r_squared = pearson_r^2,
      pearson_p = pearson_p,
      spearman_r = spearman_r,
      spearman_p = spearman_p,
      stringsAsFactors = FALSE
    )
  }))
}

pca_objects <- list(
  vector_normalized_pca = list(
    pca_object = readRDS(STANDARDIZED_PCA_RDS),
    preprocess = "vector_normalized"
  )
)

run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
target_configs <- SCALE_CONFIG[SCALE_CONFIG$scale %in% SCALE_VALUES, , drop = FALSE]
if (nrow(target_configs) == 0) {
  stop("No matching scale configuration for SAMPLED_2D_SCALES.", call. = FALSE)
}

for (config_i in seq_len(nrow(target_configs))) {
  config <- target_configs[config_i, ]
  targets <- load_target_rows(config)

  for (row_i in seq_len(nrow(targets))) {
    row <- targets[row_i, ]
    message("Processing ", row$scale, " ", row$quad_id, " (", row_i, "/", nrow(targets), ")")

    for (basis_name in intersect(names(pca_objects), BASIS_VALUES)) {
      if (is_completed_combo(row$scale, row$quad_id, basis_name, SAMPLE_MAX_POINTS)) {
        message("Skipping completed ", row$scale, " ", row$quad_id, " ", basis_name)
        next
      }

      result_row <- data.frame(
        scale = row$scale,
        quad_id = row$quad_id,
        basis = basis_name,
        sample_max_points = SAMPLE_MAX_POINTS,
        manual_excluded = as.logical(row$manual_excluded),
        metric_method = row$metric_method,
        phy_afaith = row$phy_afaith,
        n_retained_pixels = NA_integer_,
        n_sampled_pixels = NA_integer_,
        sample_method = NA_character_,
        mean_distance_2d = NA_real_,
        rao_q_2d = NA_real_,
        status = NA_character_,
        elapsed_seconds = NA_real_,
        raster_file = row$raster_file,
        previous_all_pixel_pc1_pc3_mean_distance = row$standardized_PCA_pca_euclidean_mean,
        previous_all_pixel_pc1_pc3_rao_q = row$standardized_PCA_rao_q_pca,
        run_timestamp = run_timestamp,
        stringsAsFactors = FALSE
      )

      if (isTRUE(as.logical(row$manual_excluded))) {
        result_row$status <- "skipped_manual_excluded"
        append_result_row(result_row)
        next
      }

      if (!is.character(row$raster_file) || is.na(row$raster_file) || !file.exists(row$raster_file)) {
        result_row$status <- "missing_raster"
        append_result_row(result_row)
        next
      }

      start_time <- proc.time()[["elapsed"]]
      x <- read_masked_spectra(row$raster_file)
      scores <- project_scores(x, pca_objects[[basis_name]]$pca_object, pca_objects[[basis_name]]$preprocess)
      sampled <- sample_scores_2d(scores, row$scale, row$quad_id)
      metrics <- calculate_2d_metrics(sampled$points)
      elapsed <- proc.time()[["elapsed"]] - start_time

      result_row$n_retained_pixels <- sampled$n_available
      result_row$n_sampled_pixels <- sampled$n_used
      result_row$sample_method <- sampled$method
      result_row$mean_distance_2d <- metrics$mean_distance_2d
      result_row$rao_q_2d <- metrics$rao_q_2d
      result_row$status <- metrics$status
      result_row$elapsed_seconds <- elapsed
      append_result_row(result_row)
    }
  }
}

results <- read_existing_results()
current_results <- results[
  results$scale %in% SCALE_VALUES &
    results$basis %in% BASIS_VALUES &
    results$sample_max_points == SAMPLE_MAX_POINTS,
  ,
  drop = FALSE
]

summary_counts <- as.data.frame(table(current_results$status), stringsAsFactors = FALSE)
names(summary_counts) <- c("status", "rows")
summary_counts$total_elapsed_seconds <- vapply(
  summary_counts$status,
  function(status) {
    round(sum(current_results$elapsed_seconds[current_results$status == status], na.rm = TRUE), 3)
  },
  numeric(1)
)
write.csv(summary_counts, SUMMARY_CSV, row.names = FALSE)

correlations <- correlation_rows(current_results)
write.csv(correlations, CORRELATION_CSV, row.names = FALSE)

validation <- data.frame(
  check = c(
    "requested_scale",
    "requested_basis",
    "sample_max_points",
    "current_rows",
    "ok_rows",
    "skipped_manual_excluded_rows",
    "non_ok_non_skipped_rows"
  ),
  value = c(
    paste(SCALE_VALUES, collapse = ","),
    paste(BASIS_VALUES, collapse = ","),
    SAMPLE_MAX_POINTS,
    nrow(current_results),
    sum(current_results$status == "ok"),
    sum(current_results$status == "skipped_manual_excluded"),
    sum(!(current_results$status %in% c("ok", "skipped_manual_excluded")))
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, VALIDATION_CSV, row.names = FALSE)

report_correlations <- correlations
numeric_cols <- c("pearson_r", "r_squared", "pearson_p", "spearman_r", "spearman_p")
for (col in numeric_cols) {
  report_correlations[[col]] <- format_number(report_correlations[[col]], digits = 3)
}

analysis_lines <- c(
  "# Sampled 2D PCA Mean Distance and Rao's Q",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Purpose",
  "",
  "Calculate vector-normalized PCA PC1-PC2 mean Euclidean distance and Rao's Q using the same 10,000-pixel sampling cap used by the sampled alpha-hull workflow.",
  "",
  "## Scope",
  "",
  paste("- Scale(s):", paste(SCALE_VALUES, collapse = ", ")),
  paste("- PCA basis:", paste(BASIS_VALUES, collapse = ", ")),
  paste("- Maximum sampled pixels per quadrat:", SAMPLE_MAX_POINTS),
  "- Metric space: vector-normalized PCA PC1-PC2",
  "- Rao's Q uses equal pixel weights and squared Euclidean distance, computed as `2 * mean(squared radius)` around the sampled centroid.",
  "",
  "## Status Summary",
  "",
  write_markdown_table(summary_counts),
  "",
  "## Correlations With Abundance-Weighted Faith's PD",
  "",
  write_markdown_table(report_correlations),
  "",
  "## Output Tables",
  "",
  paste0("- `", RESULTS_CSV, "`"),
  paste0("- `", SUMMARY_CSV, "`"),
  paste0("- `", CORRELATION_CSV, "`"),
  paste0("- `", VALIDATION_CSV, "`")
)
writeLines(analysis_lines, ANALYSIS_REPORT)

task_lines <- c(
  "# Task Report: Sampled 2D PCA Mean Distance and Rao's Q",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Actions",
  "",
  "- Added `scripts/3_Analysis/sample_2d_pca_mean_rao_values.R`.",
  "- Calculated 10,000-pixel sampled vector-normalized PCA PC1-PC2 mean distance and Rao's Q for 50 m quadrats.",
  "- Wrote values, summary, validation, and abundance-weighted Faith's PD correlation tables."
)
writeLines(task_lines, TASK_REPORT)

validation_lines <- c(
  "# Validation Report: Sampled 2D PCA Mean Distance and Rao's Q",
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
    "Pass. Requested sampled 2D PCA metric rows are complete or explicitly skipped for manual exclusions."
  } else {
    "Review required. One or more requested rows has a non-ok status."
  }
)
writeLines(validation_lines, VALIDATION_REPORT)

message("Completed sampled 2D PCA mean/Rao workflow.")
