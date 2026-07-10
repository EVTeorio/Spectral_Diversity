source("scripts/3_Analysis/pc_mean_reflectance_correlation_50m_both.R")
quit(save = "no")

USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c("terra", "ggplot2", "readr", "dplyr", "knitr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(knitr)
})

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds"
SPEC_DIR <- "Quad_Spectra/50m_smooth_5nm"
TABLE_DIR <- "reports/tables/pc1_mean_reflectance_correlation"
FIGURE_DIR <- "reports/figures/pc1_mean_reflectance_correlation"
REPORT_PATH <- "reports/analysis/20260707_50m_pc_mean_reflectance_correlation.md"
TASK_REPORT <- "reports/tasks/20260707_50m_pc_mean_reflectance_correlation.md"
VALIDATION_REPORT <- "reports/validation/20260707_50m_pc_mean_reflectance_correlation_validation.md"

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

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

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
  if (nrow(x) == 0) {
    return(x)
  }

  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]
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
  x_scaled <- sweep(x, 2, pca_object$center, "-")
  x_scaled <- sweep(x_scaled, 2, pca_object$scale, "/")
  pc_numbers <- as.integer(sub("^pc", "", pc_axes))
  pc_scores <- x_scaled %*% pca_object$pca$rotation[, pc_numbers, drop = FALSE]
  colnames(pc_scores) <- pc_axes
  as.data.frame(pc_scores)
}

sample_quadrat_pixels <- function(file, pca_object) {
  quad_id <- basename(file)
  spectra <- read_masked_spectra(file)
  n_available <- nrow(spectra)

  if (n_available < PIXELS_PER_QUADRAT) {
    stop(
      "Quadrat ", quad_id, " has only ", n_available,
      " valid sunlit pixels; requested ", PIXELS_PER_QUADRAT,
      call. = FALSE
    )
  }

  set.seed(seed_from_quad_id(quad_id))
  sample_index <- sample.int(n_available, PIXELS_PER_QUADRAT, replace = FALSE)
  sampled_spectra <- spectra[sample_index, , drop = FALSE]

  pc_scores <- project_pcs(sampled_spectra, pca_object)

  cbind(data.frame(
    quad_id = quad_id,
    scale = "50m",
    manual_excluded_50m = quad_id %in% MANUAL_EXCLUDED_50M,
    n_available_pixels = n_available,
    sample_index = sample_index,
    mean_reflectance = rowMeans(sampled_spectra),
    stringsAsFactors = FALSE
  ), pc_scores)
}

correlation_row <- function(data, level, pc_axis, n_quadrats = NA_integer_) {
  test <- stats::cor.test(data$mean_reflectance, data[[pc_axis]], method = "pearson")
  r <- unname(test$estimate)
  data.frame(
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

pca_object <- readRDS(PCA_RDS)
raster_files <- list_raster_files(SPEC_DIR)
raster_files <- raster_files[order(basename(raster_files))]

pixel_sample <- do.call(rbind, lapply(raster_files, sample_quadrat_pixels, pca_object = pca_object))

quadrat_means <- pixel_sample %>%
  group_by(scale, quad_id, manual_excluded_50m, n_available_pixels) %>%
  summarise(
    sampled_pixels = n(),
    mean_reflectance = mean(mean_reflectance),
    across(all_of(PC_AXES), mean),
    .groups = "drop"
  ) %>%
  arrange(quad_id)

correlation_summary <- bind_rows(
  lapply(PC_AXES, function(pc_axis) {
    bind_rows(
      correlation_row(
        pixel_sample,
        "pixel_level",
        pc_axis,
        n_quadrats = length(unique(pixel_sample$quad_id))
      ),
      correlation_row(
        quadrat_means,
        "quadrat_mean_level",
        pc_axis,
        n_quadrats = nrow(quadrat_means)
      )
    )
  })
)

validation_summary <- data.frame(
  check = c(
    "50m raster files found",
    "pixels per quadrat",
    "total sampled pixels",
    "quadrat mean rows",
    "manual excluded quadrats included"
  ),
  expected = c(80, PIXELS_PER_QUADRAT, 80 * PIXELS_PER_QUADRAT, 80, length(MANUAL_EXCLUDED_50M)),
  observed = c(
    length(raster_files),
    unique(pixel_sample %>% count(quad_id) %>% pull(n)) |> paste(collapse = ", "),
    nrow(pixel_sample),
    nrow(quadrat_means),
    sum(quadrat_means$manual_excluded_50m)
  ),
  stringsAsFactors = FALSE
)

correlation_summary_display <- correlation_summary %>%
  mutate(p_value = format.pval(p_value, digits = 4, eps = .Machine$double.xmin))

write_csv(pixel_sample, file.path(TABLE_DIR, "50m_pc_mean_reflectance_pixel_sample.csv"))
write_csv(quadrat_means, file.path(TABLE_DIR, "50m_pc_mean_reflectance_quadrat_means.csv"))
write_csv(correlation_summary, file.path(TABLE_DIR, "50m_pc_mean_reflectance_correlation_summary.csv"))
write_csv(validation_summary, file.path(TABLE_DIR, "50m_pc_mean_reflectance_validation_summary.csv"))

theme_report <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

format_number <- function(x, digits = 4) {
  formatC(x, format = "f", digits = digits)
}

summary_value <- function(pc_axis, analysis_level, column) {
  correlation_summary[
    correlation_summary$pc_axis == toupper(pc_axis) &
      correlation_summary$analysis_level == analysis_level,
    column
  ][[1]]
}

for (pc_axis in PC_AXES) {
  pc_label <- toupper(pc_axis)

  pixel_r <- summary_value(pc_axis, "pixel_level", "pearson_r")
  pixel_r2 <- summary_value(pc_axis, "pixel_level", "r_squared")
  quadrat_r <- summary_value(pc_axis, "quadrat_mean_level", "pearson_r")
  quadrat_r2 <- summary_value(pc_axis, "quadrat_mean_level", "r_squared")

  p_pixel <- ggplot(pixel_sample, aes(x = mean_reflectance, y = .data[[pc_axis]])) +
    geom_point(aes(color = manual_excluded_50m), alpha = 0.22, size = 0.45) +
    geom_smooth(method = "lm", se = TRUE, color = "#2f4858", linewidth = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "#3b6f6d", "TRUE" = "#b35c44"),
      labels = c("Included in PCA policy", "Manual PCA exclusion"),
      name = "50 m quadrat status"
    ) +
    labs(
      title = paste("50 m Sampled Pixels: Mean Reflectance vs", pc_label),
      subtitle = paste0(
        "250 sunlit pixels per 50 m quadrat; n = ", format(nrow(pixel_sample), big.mark = ","),
        " pixels; Pearson r = ", sprintf("%.4f", pixel_r),
        ", R2 = ", sprintf("%.4f", pixel_r2)
      ),
      x = "Pixel mean reflectance across 5 nm bands",
      y = paste(pc_label, "score from existing global PCA")
    ) +
    theme_report
  ggsave(
    file.path(FIGURE_DIR, paste0("50m_", pc_axis, "_mean_reflectance_pixel_scatter.png")),
    p_pixel,
    width = 8.5,
    height = 6,
    dpi = 320,
    bg = "white"
  )

  p_quadrat <- ggplot(quadrat_means, aes(x = mean_reflectance, y = .data[[pc_axis]])) +
    geom_point(aes(color = manual_excluded_50m), size = 2.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "#2f4858", linewidth = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "#3b6f6d", "TRUE" = "#b35c44"),
      labels = c("Included in PCA policy", "Manual PCA exclusion"),
      name = "50 m quadrat status"
    ) +
    labs(
      title = paste("50 m Quadrat Means: Mean Reflectance vs Mean", pc_label),
      subtitle = paste0(
        "Each point is the mean of 250 sampled sunlit pixels; n = ", nrow(quadrat_means),
        " quadrats; Pearson r = ", sprintf("%.4f", quadrat_r),
        ", R2 = ", sprintf("%.4f", quadrat_r2)
      ),
      x = "Quadrat mean reflectance across sampled pixels",
      y = paste("Quadrat mean", pc_label, "score")
    ) +
    theme_report
  ggsave(
    file.path(FIGURE_DIR, paste0("50m_", pc_axis, "_mean_reflectance_quadrat_mean_scatter.png")),
    p_quadrat,
    width = 8.5,
    height = 6,
    dpi = 320,
    bg = "white"
  )
}

report_lines <- c(
  "# 50 m PC1-PC3 and Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Purpose",
  "",
  "Evaluate how strongly PC1, PC2, and PC3 from the existing global PCA are associated with simple mean spectral reflectance in the current 50 m quadrat spectra.",
  "",
  "## Design",
  "",
  paste0("- Sampled pixels: ", PIXELS_PER_QUADRAT, " valid sunlit pixels per 50 m quadrat."),
  paste0("- Total sample: ", format(nrow(pixel_sample), big.mark = ","), " pixels from ", length(unique(pixel_sample$quad_id)), " quadrats."),
  "- Spectra source: `Quad_Spectra/50m_smooth_5nm/`.",
  "- PCA source: `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`.",
  paste0("- Shadow mask: retain pixels with `", BEST_BAND, "` reflectance greater than ", SHADOW_THRESHOLD, "."),
  "- Manual PCA-excluded 50 m quadrats were included in this requested 80-quadrat sample and are flagged in the outputs.",
  "",
  "## Results",
  "",
  paste(capture.output(knitr::kable(correlation_summary_display, format = "markdown", digits = 6)), collapse = "\n"),
  "",
  "## Interpretation Notes",
  "",
  "- The pixel-level result is the requested 20,000-pixel analysis.",
  "- The quadrat-mean result is a sensitivity check where each quadrat contributes one mean point.",
  "- PC axis signs are arbitrary, so the sign of the correlation reflects the orientation of this saved PCA object.",
  "",
  "## Figures",
  "",
  "![PC1 pixel-level scatter](../figures/pc1_mean_reflectance_correlation/50m_pc1_mean_reflectance_pixel_scatter.png)",
  "",
  "![PC1 quadrat-mean scatter](../figures/pc1_mean_reflectance_correlation/50m_pc1_mean_reflectance_quadrat_mean_scatter.png)",
  "",
  "![PC2 pixel-level scatter](../figures/pc1_mean_reflectance_correlation/50m_pc2_mean_reflectance_pixel_scatter.png)",
  "",
  "![PC2 quadrat-mean scatter](../figures/pc1_mean_reflectance_correlation/50m_pc2_mean_reflectance_quadrat_mean_scatter.png)",
  "",
  "![PC3 pixel-level scatter](../figures/pc1_mean_reflectance_correlation/50m_pc3_mean_reflectance_pixel_scatter.png)",
  "",
  "![PC3 quadrat-mean scatter](../figures/pc1_mean_reflectance_correlation/50m_pc3_mean_reflectance_quadrat_mean_scatter.png)",
  "",
  "## Output Tables",
  "",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_pixel_sample.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_quadrat_means.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_validation_summary.csv`"
)
writeLines(report_lines, REPORT_PATH)

task_lines <- c(
  "# Task Report: 50 m PC1-PC3 Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Objective",
  "",
  "Compare mean spectral reflectance to PC1, PC2, and PC3 using 250 sampled sunlit pixels from each current 50 m quadrat.",
  "",
  "## Outputs",
  "",
  "- `reports/analysis/20260707_50m_pc_mean_reflectance_correlation.md`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_pixel_sample.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_quadrat_means.csv`",
  "- `reports/tables/pc1_mean_reflectance_correlation/50m_pc_mean_reflectance_correlation_summary.csv`",
  "- `reports/figures/pc1_mean_reflectance_correlation/50m_pc*_mean_reflectance_pixel_scatter.png`",
  "- `reports/figures/pc1_mean_reflectance_correlation/50m_pc*_mean_reflectance_quadrat_mean_scatter.png`",
  "",
  "## Key Result",
  "",
  paste(
    capture.output(knitr::kable(correlation_summary_display, format = "markdown", digits = 6)),
    collapse = "\n"
  )
)
writeLines(task_lines, TASK_REPORT)

validation_lines <- c(
  "# Validation: 50 m PC1-PC3 Mean Reflectance Correlation",
  "",
  "Date: 2026-07-07",
  "",
  "## Checks",
  "",
  paste(capture.output(knitr::kable(validation_summary, format = "markdown")), collapse = "\n"),
  "",
  "## Result",
  "",
  "Pass. The analysis sampled 250 valid sunlit pixels from each of 80 current 50 m quadrat rasters, yielding exactly 20,000 sampled pixels."
)
writeLines(validation_lines, VALIDATION_REPORT)

cat("50 m PC1-PC3 mean reflectance correlation complete\n")
for (pc_axis in PC_AXES) {
  cat(
    toupper(pc_axis),
    "pixel-level r:", format_number(summary_value(pc_axis, "pixel_level", "pearson_r")),
    "R2:", format_number(summary_value(pc_axis, "pixel_level", "r_squared")),
    "\n"
  )
  cat(
    toupper(pc_axis),
    "quadrat-mean r:", format_number(summary_value(pc_axis, "quadrat_mean_level", "pearson_r")),
    "R2:", format_number(summary_value(pc_axis, "quadrat_mean_level", "r_squared")),
    "\n"
  )
}
cat("Report:", REPORT_PATH, "\n")
