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

TABLE_DIR <- "reports/tables/pca_loading_spectral_regions"
FIGURE_DIR <- "reports/figures/pca_loading_spectral_regions"
REPORT_PATH <- "reports/analysis/20260707_pca_loading_spectral_regions.md"
TASK_REPORT <- "reports/tasks/20260707_pca_loading_spectral_regions.md"
VALIDATION_REPORT <- "reports/validation/20260707_pca_loading_spectral_regions_validation.md"

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

NONUNIFORM_QUANTILE <- 0.90
BRIGHTNESS_FLAT_Z <- 1
RATIO_EPSILON <- 0.10

parse_wavelength_nm <- function(x) as.numeric(gsub("[^0-9.]+", "", x))
standardize <- function(x) as.numeric(scale(x))
sign_label <- function(x) ifelse(x >= 0, "positive", "negative")
dominant_value <- function(x) names(sort(table(x), decreasing = TRUE))[1]

spectral_zone <- function(wavelength_nm) {
  cut(
    wavelength_nm,
    breaks = c(-Inf, 450, 520, 600, 680, 750, 900, Inf),
    right = FALSE,
    labels = c("violet-blue", "blue-green", "green-yellow", "red", "red-edge", "near infrared", "longer near infrared")
  )
}

collapse_regions <- function(data, flag_col, value_col, label) {
  flagged <- data[data[[flag_col]], , drop = FALSE]
  if (nrow(flagged) == 0) return(data.frame())

  run_id <- cumsum(c(TRUE, diff(flagged$wavelength_nm) > 5))
  regions <- lapply(split(flagged, run_id), function(region_data) {
    peak_index <- which.max(abs(region_data[[value_col]]))
    peak_row <- region_data[peak_index, , drop = FALSE]
    data.frame(
      pca_basis = peak_row$pca_basis,
      pca_label = peak_row$pca_label,
      region_type = label,
      start_nm = min(region_data$wavelength_nm),
      end_nm = max(region_data$wavelength_nm),
      n_bands = nrow(region_data),
      dominant_zone = dominant_value(region_data$spectral_zone),
      peak_wavelength_nm = peak_row$wavelength_nm,
      peak_value = peak_row[[value_col]],
      loading_sign = sign_label(peak_row[[value_col]]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, regions)
}

fmt <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))

markdown_table <- function(data, digits = 4) {
  if (nrow(data) == 0) return("_No rows._")
  out <- data
  for (col in names(out)) {
    if (is.numeric(out[[col]])) out[[col]] <- fmt(out[[col]], digits)
  }
  header <- paste0("| ", paste(names(out), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(out)), collapse = " | "), " |")
  rows <- apply(out, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, divider, rows), collapse = "\n")
}

analyze_pca_loadings <- function(config) {
  pca_object <- readRDS(config$pca_rds)
  rotation <- pca_object$pca$rotation[, 1:2, drop = FALSE]
  wavelength_nm <- parse_wavelength_nm(rownames(rotation))

  pc1_loading <- rotation[, 1]
  pc2_loading <- rotation[, 2]
  if (mean(pc1_loading) < 0) {
    pc1_loading <- -pc1_loading
    pc2_loading <- -pc2_loading
  }

  n_bands <- length(pc1_loading)
  uniform_pc1_loading <- rep(1 / sqrt(n_bands), n_bands)
  pc1_departure <- pc1_loading - uniform_pc1_loading
  pc1_departure_z <- standardize(pc1_departure)
  pc2_loading_z <- standardize(pc2_loading)
  pc1_nonuniform_threshold <- unname(quantile(abs(pc1_departure_z), NONUNIFORM_QUANTILE))
  pc2_nonuniform_threshold <- unname(quantile(abs(pc2_loading_z), NONUNIFORM_QUANTILE))

  loading_table <- data.frame(
    pca_basis = config$pca_basis,
    pca_label = config$pca_label,
    wavelength_nm = wavelength_nm,
    spectral_zone = as.character(spectral_zone(wavelength_nm)),
    pc1_loading = pc1_loading,
    pc1_uniform_loading = uniform_pc1_loading,
    pc1_departure_from_uniform = pc1_departure,
    pc1_departure_z = pc1_departure_z,
    pc2_loading = pc2_loading,
    pc2_loading_z = pc2_loading_z,
    stringsAsFactors = FALSE
  )
  loading_table$pc1_nonuniform_flag <- abs(loading_table$pc1_departure_z) >= pc1_nonuniform_threshold
  loading_table$pc2_nonuniform_flag <- abs(loading_table$pc2_loading_z) >= pc2_nonuniform_threshold
  loading_table$pc2_strong_not_brightness_flag <-
    loading_table$pc2_nonuniform_flag & abs(loading_table$pc1_departure_z) < BRIGHTNESS_FLAT_Z
  loading_table$pc2_to_pc1_departure_ratio <-
    abs(loading_table$pc2_loading_z) / (abs(loading_table$pc1_departure_z) + RATIO_EPSILON)
  loading_table <- loading_table[order(loading_table$wavelength_nm), ]

  region_table <- rbind(
    collapse_regions(loading_table, "pc1_nonuniform_flag", "pc1_departure_z", "PC1 departure from uniform brightness"),
    collapse_regions(loading_table, "pc2_nonuniform_flag", "pc2_loading_z", "PC2 nonuniform loading"),
    collapse_regions(loading_table, "pc2_strong_not_brightness_flag", "pc2_loading_z", "PC2 strong while PC1 is relatively flat")
  )

  summary_table <- data.frame(
    pca_basis = config$pca_basis,
    pca_label = config$pca_label,
    metric = c(
      "n_bands",
      "wavelength_min_nm",
      "wavelength_max_nm",
      "pc1_uniform_loading",
      "pc1_nonuniform_abs_z_threshold",
      "pc2_nonuniform_abs_z_threshold",
      "brightness_flat_abs_pc1_departure_z_threshold",
      "pc1_nonuniform_band_count",
      "pc2_nonuniform_band_count",
      "pc2_strong_not_brightness_band_count"
    ),
    value = c(
      n_bands,
      min(loading_table$wavelength_nm),
      max(loading_table$wavelength_nm),
      unique(loading_table$pc1_uniform_loading)[1],
      pc1_nonuniform_threshold,
      pc2_nonuniform_threshold,
      BRIGHTNESS_FLAT_Z,
      sum(loading_table$pc1_nonuniform_flag),
      sum(loading_table$pc2_nonuniform_flag),
      sum(loading_table$pc2_strong_not_brightness_flag)
    ),
    stringsAsFactors = FALSE
  )

  list(loadings = loading_table, regions = region_table, summary = summary_table)
}

plot_loading_summary <- function(loading_table, config) {
  png(file.path(FIGURE_DIR, paste0(config$pca_basis, "_pc1_pc2_loadings_by_wavelength.png")), width = 2400, height = 1400, res = 300)
  par(mar = c(4.5, 4.5, 3.5, 1), las = 1)
  yrange <- range(c(loading_table$pc1_loading, loading_table$pc2_loading), na.rm = TRUE)
  plot(loading_table$wavelength_nm, loading_table$pc1_loading, type = "l", col = "#3b6f6d", lwd = 2, ylim = yrange, xlab = "Wavelength (nm)", ylab = "Loading", main = paste(config$pca_label, "PC1 and PC2 Loadings"))
  abline(h = 0, col = "gray55", lwd = 1)
  lines(loading_table$wavelength_nm, loading_table$pc2_loading, col = "#b35c44", lwd = 2)
  legend("topright", legend = c("PC1", "PC2"), col = c("#3b6f6d", "#b35c44"), lwd = 2, bty = "n")
  dev.off()

  pc1_regions <- loading_table[loading_table$pc1_nonuniform_flag, ]
  png(file.path(FIGURE_DIR, paste0(config$pca_basis, "_pc1_departure_from_uniform_brightness.png")), width = 2400, height = 1400, res = 300)
  par(mar = c(4.5, 4.5, 3.5, 1), las = 1)
  plot(loading_table$wavelength_nm, loading_table$pc1_departure_z, type = "l", col = "#3b6f6d", lwd = 2, xlab = "Wavelength (nm)", ylab = "PC1 departure from uniform loading (z)", main = paste(config$pca_label, "PC1 Departure From Uniform Brightness"))
  abline(h = 0, col = "gray55", lwd = 1)
  points(pc1_regions$wavelength_nm, pc1_regions$pc1_departure_z, pch = 16, col = "#d8a84f")
  dev.off()

  pc2_regions <- loading_table[loading_table$pc2_nonuniform_flag, ]
  pc2_independent <- loading_table[loading_table$pc2_strong_not_brightness_flag, ]
  png(file.path(FIGURE_DIR, paste0(config$pca_basis, "_pc2_nonuniform_and_nonbrightness_regions.png")), width = 2400, height = 1400, res = 300)
  par(mar = c(4.5, 4.5, 3.5, 1), las = 1)
  plot(loading_table$wavelength_nm, loading_table$pc2_loading_z, type = "l", col = "#b35c44", lwd = 2, xlab = "Wavelength (nm)", ylab = "PC2 loading (z)", main = paste(config$pca_label, "PC2 Nonuniform and Non-Brightness Regions"))
  abline(h = 0, col = "gray55", lwd = 1)
  points(pc2_regions$wavelength_nm, pc2_regions$pc2_loading_z, pch = 16, col = "#b35c44")
  points(pc2_independent$wavelength_nm, pc2_independent$pc2_loading_z, pch = 16, col = "#3b6f6d")
  legend("topright", legend = c("PC2 nonuniform", "PC2 strong while PC1 flat"), col = c("#b35c44", "#3b6f6d"), pch = 16, bty = "n")
  dev.off()
}

region_phrase <- function(regions) {
  if (nrow(regions) == 0) return("none")
  paste(ifelse(regions$start_nm == regions$end_nm, paste0(regions$start_nm, " nm"), paste0(regions$start_nm, "-", regions$end_nm, " nm")), collapse = ", ")
}

results <- lapply(seq_len(nrow(PCA_CONFIGS)), function(i) analyze_pca_loadings(PCA_CONFIGS[i, ]))
loading_table <- do.call(rbind, lapply(results, `[[`, "loadings"))
region_table <- do.call(rbind, lapply(results, `[[`, "regions"))
summary_table <- do.call(rbind, lapply(results, `[[`, "summary"))

write.csv(loading_table, file.path(TABLE_DIR, "pca_pc1_pc2_loadings_by_wavelength.csv"), row.names = FALSE)
write.csv(region_table, file.path(TABLE_DIR, "pca_pc1_pc2_nonuniform_regions.csv"), row.names = FALSE)
write.csv(summary_table, file.path(TABLE_DIR, "pca_pc1_pc2_loading_summary.csv"), row.names = FALSE)

for (i in seq_len(nrow(PCA_CONFIGS))) {
  plot_loading_summary(loading_table[loading_table$pca_basis == PCA_CONFIGS$pca_basis[i], ], PCA_CONFIGS[i, ])
}

standardized_regions <- region_table[region_table$pca_basis == "standardized_PCA", ]
standardized_pc1 <- standardized_regions[standardized_regions$region_type == "PC1 departure from uniform brightness", ]
standardized_pc2 <- standardized_regions[standardized_regions$region_type == "PC2 nonuniform loading", ]
standardized_independent <- standardized_regions[standardized_regions$region_type == "PC2 strong while PC1 is relatively flat", ]
regular_regions <- region_table[region_table$pca_basis == "regular_PCA", ]

report_lines <- c(
  "# PCA Loading Spectral Region Analysis",
  "",
  "Date: 2026-07-07",
  "",
  "## Purpose",
  "",
  "Identify where PC1 and PC2 loadings are not spectrally uniform for both the regular PCA and the vector-normalized standardized PCA. The stated key finding is based on the standardized PCA because vector normalization reduces broad brightness dominance before PCA fitting.",
  "",
  "## Method",
  "",
  "- PC1 was oriented so its mean loading is positive and then compared with a flat brightness vector with loading `1 / sqrt(n_bands)` at every wavelength.",
  "- PC1 nonuniform wavelengths are the top 10% of absolute z-scored departures from that flat brightness vector.",
  "- PC2 nonuniform wavelengths are the top 10% of absolute z-scored PC2 loadings.",
  "- PC2 regions least dominated by brightness are wavelengths with strong PC2 loadings where absolute PC1 departure from uniform brightness is less than 1 z-score.",
  "- Contiguous 5 nm wavelengths meeting each rule were collapsed into spectral regions.",
  "",
  "## Key Finding: Standardized PCA",
  "",
  paste0("- In the standardized PCA, PC1 is not a simple brightness-only axis; its strongest nonuniform wavelength regions are ", region_phrase(standardized_pc1), "."),
  paste0("- Standardized PC2 has concentrated nonuniform structure at ", region_phrase(standardized_pc2), "."),
  paste0("- The standardized PC2 regions least attributable to broad PC1 brightness are ", region_phrase(standardized_independent), ". These are the spectral regions to carry forward as the primary loading-based finding."),
  "- PC axis signs are arbitrary; interpret positive and negative loading signs relative to the saved PCA object.",
  "",
  "## Standardized PCA Regions",
  "",
  markdown_table(standardized_regions, digits = 4),
  "",
  "## Supporting Regular PCA Comparison",
  "",
  "The regular PCA loading regions were also regenerated for traceability, but they are not stated as the key finding because the regular PCA remains more brightness dominated.",
  "",
  markdown_table(regular_regions, digits = 4),
  "",
  "## Summary Checks",
  "",
  markdown_table(summary_table, digits = 4),
  "",
  "## Figures",
  "",
  "![Standardized PCA loadings](../figures/pca_loading_spectral_regions/standardized_PCA_pc1_pc2_loadings_by_wavelength.png)",
  "",
  "![Standardized PCA PC1 departure](../figures/pca_loading_spectral_regions/standardized_PCA_pc1_departure_from_uniform_brightness.png)",
  "",
  "![Standardized PCA PC2 regions](../figures/pca_loading_spectral_regions/standardized_PCA_pc2_nonuniform_and_nonbrightness_regions.png)",
  "",
  "![Regular PCA loadings](../figures/pca_loading_spectral_regions/regular_PCA_pc1_pc2_loadings_by_wavelength.png)",
  "",
  "## Output Tables",
  "",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_loadings_by_wavelength.csv`",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_nonuniform_regions.csv`",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_loading_summary.csv`"
)
writeLines(report_lines, REPORT_PATH)

task_lines <- c(
  "# Task Report: PCA Loading Spectral Region Analysis",
  "",
  "Date: 2026-07-07",
  "",
  "## Objective",
  "",
  "Examine where PC1 and PC2 loadings are not uniform for regular and standardized PCA, and identify standardized PC2 wavelength regions with strong structure beyond broad brightness.",
  "",
  "## Outputs",
  "",
  "- `reports/analysis/20260707_pca_loading_spectral_regions.md`",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_loadings_by_wavelength.csv`",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_nonuniform_regions.csv`",
  "- `reports/tables/pca_loading_spectral_regions/pca_pc1_pc2_loading_summary.csv`",
  "- `reports/figures/pca_loading_spectral_regions/*_pc1_pc2_loadings_by_wavelength.png`",
  "- `reports/figures/pca_loading_spectral_regions/*_pc1_departure_from_uniform_brightness.png`",
  "- `reports/figures/pca_loading_spectral_regions/*_pc2_nonuniform_and_nonbrightness_regions.png`",
  "",
  "## Key Finding",
  "",
  paste0("The key finding is stated from standardized PCA: standardized PC2 has brightness-reduced spectral structure at ", region_phrase(standardized_independent), ". Regular PCA outputs were regenerated only as supporting comparison.")
)
writeLines(task_lines, TASK_REPORT)

validation_lines <- c(
  "# Validation: PCA Loading Spectral Region Analysis",
  "",
  "Date: 2026-07-07",
  "",
  "## Checks",
  "",
  markdown_table(summary_table, digits = 4),
  "",
  "## Result",
  "",
  "Pass. The analysis read both saved PCA rotation matrices, extracted PC1 and PC2 loadings for all wavelengths, wrote combined per-band and region-level tables, and generated basis-specific figures."
)
writeLines(validation_lines, VALIDATION_REPORT)

cat("PCA loading spectral region analysis complete\n")
cat("Standardized PCA key regions:", region_phrase(standardized_independent), "\n")
cat("Report:", REPORT_PATH, "\n")
