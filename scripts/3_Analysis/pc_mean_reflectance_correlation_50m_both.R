USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Missing required package: terra", call. = FALSE)
}

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

PCA_CONFIGS <- data.frame(
  pca_basis = c("regular_PCA", "standardized_PCA"),
  pca_label = c("Regular PCA", "Standardized PCA"),
  pca_rds = c(
    "Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds",
    "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds"
  ),
  stringsAsFactors = FALSE
)

SPEC_DIR <- "Quad_Spectra/50m_smooth_5nm"
TABLE_DIR <- "reports/tables/pc1_mean_reflectance_correlation"
FIGURE_DIR <- "reports/figures/pc1_mean_reflectance_correlation"
REPORT_PATH <- "reports/analysis/20260707_50m_pc_mean_reflectance_correlation.md"
TASK_REPORT <- "reports/tasks/20260707_50m_pc_mean_reflectance_correlation.md"
VALIDATION_REPORT <- "reports/validation/20260707_50m_pc_mean_reflectance_correlation_validation.md"

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

SIDECAR_PATTERN <- "\\.(hdr|aux|xml|enp|sta)$"
PIXELS_PER_QUADRAT <- 250L
RANDOM_SEED <- 20260706L
PC_AXES <- paste0("pc", 1:3)

BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
DIRECTION <- "<"

MANUAL_EXCLUDED_50M <- c(
  "sub50_80", "sub50_79", "sub50_71", "sub50_70", "sub50_62", "sub50_53"
)

list_raster_files <- function(spec_dir) {
  files <- list.files(spec_dir, full.names = TRUE)
  files[!grepl(SIDECAR_PATTERN, files, ignore.case = TRUE)]
}

apply_shadow_mask <- function(raster) {
  band_mask <- raster[[BEST_BAND]]
  sunlit_mask <- if (DIRECTION == ">") {
    band_mask < SHADOW_THRESHOLD
  } else {
    band_mask > SHADOW_THRESHOLD
  }
  terra::mask(raster, sunlit_mask, maskvalues = 0)
}

clean_spectra <- function(x) {
  if (nrow(x) == 0) return(x)
  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]
}

vector_normalize_spectra <- function(x) {
  norms <- sqrt(rowSums(x^2))
  keep <- is.finite(norms) & norms > 0
  x <- x[keep, , drop = FALSE]
  norms <- norms[keep]
  sweep(x, 1, norms, "/")
}

preprocess_pca_spectra <- function(x, preprocess) {
  if (identical(preprocess, "vector_normalized")) {
    return(vector_normalize_spectra(x))
  }
  x
}

seed_from_quad_id <- function(quad_id) {
  RANDOM_SEED + sum(utf8ToInt(as.character(quad_id)))
}

read_masked_spectra <- function(file) {
  raster <- terra::rast(file)
  raster_masked <- apply_shadow_mask(raster)
  clean_spectra(terra::values(raster_masked, mat = TRUE))
}

project_pcs <- function(x, pca_object, pc_axes = PC_AXES) {
  x <- preprocess_pca_spectra(x, pca_object$preprocess)
  x_scaled <- sweep(x, 2, pca_object$center, "-")
  x_scaled <- sweep(x_scaled, 2, pca_object$scale, "/")
  pc_numbers <- as.integer(sub("^pc", "", pc_axes))
  pc_scores <- x_scaled %*% pca_object$pca$rotation[, pc_numbers, drop = FALSE]
  colnames(pc_scores) <- pc_axes
  as.data.frame(pc_scores)
}

sample_quadrat_pixels <- function(file) {
  quad_id <- basename(file)
  spectra <- read_masked_spectra(file)
  n_available <- nrow(spectra)
  if (n_available < PIXELS_PER_QUADRAT) {
    stop("Quadrat ", quad_id, " has only ", n_available, " valid sunlit pixels", call. = FALSE)
  }
  set.seed(seed_from_quad_id(quad_id))
  sample_index <- sample.int(n_available, PIXELS_PER_QUADRAT, replace = FALSE)
  sampled_spectra <- spectra[sample_index, , drop = FALSE]
  data.frame(
    quad_id = quad_id,
    scale = "50m",
    manual_excluded_50m = quad_id %in% MANUAL_EXCLUDED_50M,
    n_available_pixels = n_available,
    sample_index = sample_index,
    mean_reflectance = rowMeans(sampled_spectra),
    sampled_spectra,
    check.names = FALSE
  )
}

correlation_row <- function(data, pca_basis, pca_label, level, pc_axis, n_quadrats = NA_integer_) {
  test <- stats::cor.test(data$mean_reflectance, data[[pc_axis]], method = "pearson")
  r <- unname(test$estimate)
  data.frame(
    pca_basis = pca_basis,
    pca_label = pca_label,
    analysis_level = level,
    pc_axis = toupper(pc_axis),
    n_rows = nrow(data),
    n_quadrats = n_quadrats,
    pearson_r = r,
    r_squared = r^2,
    p_value = test$p.value,
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    stringsAsFactors = FALSE
  )
}

plot_scatter <- function(data, x_col, y_col, title, subtitle, xlab, ylab, file) {
  png(file, width = 2400, height = 1700, res = 300)
  par(mar = c(4.5, 4.5, 3.5, 1), las = 1)
  status_col <- ifelse(data$manual_excluded_50m, "#b35c44", "#3b6f6d")
  plot(data[[x_col]], data[[y_col]], pch = 16, cex = 0.45, col = adjustcolor(status_col, alpha.f = 0.35), xlab = xlab, ylab = ylab, main = title)
  mtext(subtitle, side = 3, line = 0.3, cex = 0.75)
  fit <- stats::lm(data[[y_col]] ~ data[[x_col]])
  abline(fit, col = "#2f4858", lwd = 2)
  legend("topleft", legend = c("Included", "Manual PCA exclusion"), col = c("#3b6f6d", "#b35c44"), pch = 16, bty = "n")
  dev.off()
}

markdown_table <- function(data, digits = 4) {
  out <- data
  for (col in names(out)) {
    if (is.numeric(out[[col]])) {
      out[[col]] <- ifelse(is.na(out[[col]]), "NA", formatC(out[[col]], format = "f", digits = digits))
    }
  }
  header <- paste0("| ", paste(names(out), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(out)), collapse = " | "), " |")
  rows <- apply(out, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, divider, rows), collapse = "\n")
}

raster_files <- list_raster_files(SPEC_DIR)
raster_files <- raster_files[order(basename(raster_files))]

sample_with_spectra <- do.call(rbind, lapply(raster_files, sample_quadrat_pixels))
metadata_cols <- c("quad_id", "scale", "manual_excluded_50m", "n_available_pixels", "sample_index", "mean_reflectance")
spectra_cols <- setdiff(names(sample_with_spectra), metadata_cols)
sampled_spectra <- as.matrix(sample_with_spectra[, spectra_cols, drop = FALSE])
pixel_metadata <- sample_with_spectra[, metadata_cols, drop = FALSE]

all_pixel_scores <- list()
all_quadrat_means <- list()
all_correlations <- list()

for (i in seq_len(nrow(PCA_CONFIGS))) {
  config <- PCA_CONFIGS[i, ]
  pca_object <- readRDS(config$pca_rds)
  pc_scores <- project_pcs(sampled_spectra, pca_object)
  pixel_sample <- cbind(
    data.frame(pca_basis = config$pca_basis, pca_label = config$pca_label, pixel_metadata, stringsAsFactors = FALSE),
    pc_scores
  )

  quadrat_means <- aggregate(
    pixel_sample[, c("mean_reflectance", PC_AXES)],
    by = pixel_sample[, c("pca_basis", "pca_label", "scale", "quad_id", "manual_excluded_50m", "n_available_pixels")],
    FUN = mean
  )
  names(quadrat_means)[names(quadrat_means) == "mean_reflectance"] <- "mean_reflectance"
  quadrat_counts <- aggregate(sample_index ~ quad_id, data = pixel_sample, FUN = length)
  names(quadrat_counts)[2] <- "sampled_pixels"
  quadrat_means <- merge(quadrat_means, quadrat_counts, by = "quad_id", all.x = TRUE)
  quadrat_means <- quadrat_means[order(quadrat_means$quad_id), ]

  correlation_summary <- do.call(rbind, lapply(PC_AXES, function(pc_axis) {
    rbind(
      correlation_row(pixel_sample, config$pca_basis, config$pca_label, "pixel_level", pc_axis, length(unique(pixel_sample$quad_id))),
      correlation_row(quadrat_means, config$pca_basis, config$pca_label, "quadrat_mean_level", pc_axis, nrow(quadrat_means))
    )
  }))

  for (pc_axis in PC_AXES) {
    pc_label <- toupper(pc_axis)
    pixel_result <- correlation_summary[correlation_summary$pc_axis == pc_label & correlation_summary$analysis_level == "pixel_level", ]
    quadrat_result <- correlation_summary[correlation_summary$pc_axis == pc_label & correlation_summary$analysis_level == "quadrat_mean_level", ]
    plot_scatter(
      pixel_sample,
      "mean_reflectance",
      pc_axis,
      paste(config$pca_label, "50 m Sampled Pixels: Mean Reflectance vs", pc_label),
      paste0("n = ", nrow(pixel_sample), " pixels; r = ", round(pixel_result$pearson_r, 4), ", R2 = ", round(pixel_result$r_squared, 4)),
      "Pixel mean reflectance across 5 nm bands",
      paste(pc_label, "score"),
      file.path(FIGURE_DIR, paste0("50m_", config$pca_basis, "_", pc_axis, "_mean_reflectance_pixel_scatter.png"))
    )
    plot_scatter(
      quadrat_means,
      "mean_reflectance",
      pc_axis,
      paste(config$pca_label, "50 m Quadrat Means: Mean Reflectance vs Mean", pc_label),
      paste0("n = ", nrow(quadrat_means), " quadrats; r = ", round(quadrat_result$pearson_r, 4), ", R2 = ", round(quadrat_result$r_squared, 4)),
      "Quadrat mean reflectance across sampled pixels",
      paste("Quadrat mean", pc_label, "score"),
      file.path(FIGURE_DIR, paste0("50m_", config$pca_basis, "_", pc_axis, "_mean_reflectance_quadrat_mean_scatter.png"))
    )
  }

  write.csv(pixel_sample, file.path(TABLE_DIR, paste0("50m_", config$pca_basis, "_pc_mean_reflectance_pixel_sample.csv")), row.names = FALSE)
  write.csv(quadrat_means, file.path(TABLE_DIR, paste0("50m_", config$pca_basis, "_pc_mean_reflectance_quadrat_means.csv")), row.names = FALSE)
  write.csv(correlation_summary, file.path(TABLE_DIR, paste0("50m_", config$pca_basis, "_pc_mean_reflectance_correlation_summary.csv")), row.names = FALSE)

  all_pixel_scores[[config$pca_basis]] <- pixel_sample
  all_quadrat_means[[config$pca_basis]] <- quadrat_means
  all_correlations[[config$pca_basis]] <- correlation_summary
}

combined_pixel_sample <- do.call(rbind, all_pixel_scores)
combined_quadrat_means <- do.call(rbind, all_quadrat_means)
combined_correlation_summary <- do.call(rbind, all_correlations)

validation_summary <- data.frame(
  check = c("50m raster files found", "pixels per quadrat", "total sampled pixels per PCA basis", "quadrat mean rows per PCA basis", "manual excluded quadrats included"),
  expected = c(80, PIXELS_PER_QUADRAT, 80 * PIXELS_PER_QUADRAT, 80, length(MANUAL_EXCLUDED_50M)),
  observed = c(
    length(raster_files),
    paste(unique(table(combined_pixel_sample$quad_id, combined_pixel_sample$pca_basis)), collapse = ", "),
    paste(tapply(combined_pixel_sample$sample_index, combined_pixel_sample$pca_basis, length), collapse = ", "),
    paste(tapply(combined_quadrat_means$quad_id, combined_quadrat_means$pca_basis, length), collapse = ", "),
    sum(unique(combined_quadrat_means[, c("quad_id", "manual_excluded_50m")])$manual_excluded_50m)
  ),
  stringsAsFactors = FALSE
)

write.csv(combined_pixel_sample[, c("pca_basis", "pca_label", metadata_cols, PC_AXES)], file.path(TABLE_DIR, "50m_pc_mean_reflectance_pixel_sample.csv"), row.names = FALSE)
write.csv(combined_quadrat_means, file.path(TABLE_DIR, "50m_pc_mean_reflectance_quadrat_means.csv"), row.names = FALSE)
write.csv(combined_correlation_summary, file.path(TABLE_DIR, "50m_pc_mean_reflectance_correlation_summary.csv"), row.names = FALSE)
write.csv(validation_summary, file.path(TABLE_DIR, "50m_pc_mean_reflectance_validation_summary.csv"), row.names = FALSE)

standardized_summary <- combined_correlation_summary[combined_correlation_summary$pca_basis == "standardized_PCA", ]
regular_summary <- combined_correlation_summary[combined_correlation_summary$pca_basis == "regular_PCA", ]

report_lines <- c(
  "# 50 m PC1-PC3 and Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Purpose",
  "",
  "Evaluate how strongly PC1, PC2, and PC3 are associated with simple mean spectral reflectance for both the regular PCA and vector-normalized standardized PCA. The key interpretation is based on standardized PCA because it is the brightness-reduced PCA basis.",
  "",
  "## Design",
  "",
  paste0("- Sampled pixels: ", PIXELS_PER_QUADRAT, " valid sunlit pixels per 50 m quadrat."),
  paste0("- Total sample per PCA basis: ", format(80 * PIXELS_PER_QUADRAT, big.mark = ","), " pixels from 80 quadrats."),
  "- Spectra source: `Quad_Spectra/50m_smooth_5nm/`.",
  "- PCA sources: regular and standardized PCA RDS files in `Quad_Values/Spectral_diversitySHPs/`.",
  paste0("- Shadow mask: retain pixels with `", BEST_BAND, "` reflectance greater than ", SHADOW_THRESHOLD, "."),
  "- Manual PCA-excluded 50 m quadrats were included in this diagnostic sample and are flagged in the outputs.",
  "",
  "## Key Results: Standardized PCA",
  "",
  markdown_table(standardized_summary, digits = 6),
  "",
  "## Supporting Regular PCA Results",
  "",
  markdown_table(regular_summary, digits = 6),
  "",
  "## Interpretation Notes",
  "",
  "- The pixel-level result is the requested 20,000-pixel diagnostic per PCA basis.",
  "- The quadrat-mean result is a sensitivity check where each quadrat contributes one mean point.",
  "- PC axis signs are arbitrary, so the sign of the correlation reflects the orientation of each saved PCA object.",
  "",
  "## Output Tables",
  "",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_pixel_sample.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_quadrat_means.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_regular_PCA_pc_mean_reflectance_correlation_summary.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_standardized_PCA_pc_mean_reflectance_correlation_summary.csv`"
)
writeLines(report_lines, REPORT_PATH)

task_lines <- c(
  "# Task Report: 50 m PC1-PC3 Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Objective",
  "",
  "Compare mean spectral reflectance to PC1, PC2, and PC3 for regular and standardized PCA using 250 sampled sunlit pixels from each current 50 m quadrat.",
  "",
  "## Key Output",
  "",
  "The combined summary table is `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`; basis-specific summaries are also written for regular and standardized PCA.",
  "",
  "## Standardized PCA Results",
  "",
  markdown_table(standardized_summary, digits = 6)
)
writeLines(task_lines, TASK_REPORT)

validation_lines <- c(
  "# Validation: 50 m PC1-PC3 Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Checks",
  "",
  markdown_table(validation_summary),
  "",
  "## Result",
  "",
  "Pass when the script completes. It samples 250 valid sunlit pixels from each of 80 current 50 m quadrat rasters and projects the same sampled spectra into both regular and standardized PCA space."
)
writeLines(validation_lines, VALIDATION_REPORT)

cat("50 m PC1-PC3 mean reflectance correlation complete\n")
cat("Report:", REPORT_PATH, "\n")
