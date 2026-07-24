PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/multiscale_spectral_biodiversity")
TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/multiscale_spectral_biodiversity")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260710_sv_diversity_pairwise_correlations.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md")

SV_MEASURES <- c("spec_spca_mean", "spec_sa")
DIVERSITY_MEASURES <- c(
  "phy_faith", "phy_rao", "phy_afaith",
  "sp_rich", "sp_shannon", "sp_simpson", "sp_even"
)

DISPLAY_NAMES <- c(
  spec_spca_mean = "Standardized PCA mean Euclidean distance",
  spec_sa = "SA entropy",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_rich = "Species richness",
  sp_shannon = "Shannon diversity",
  sp_simpson = "Simpson diversity",
  sp_even = "Species evenness"
)

DIVERSITY_GROUPS <- c(
  phy_faith = "Phylogenetic diversity",
  phy_rao = "Phylogenetic diversity",
  phy_afaith = "Phylogenetic diversity",
  sp_rich = "Species diversity",
  sp_shannon = "Species diversity",
  sp_simpson = "Species diversity",
  sp_even = "Species diversity"
)

HEATMAP_NAMES <- c(
  spec_spca_mean = "Std PCA distance",
  spec_sa = "SA entropy",
  phy_faith = "Faith PD",
  phy_rao = "Phylo Rao Q",
  phy_afaith = "Abund. Faith PD",
  sp_rich = "Richness",
  sp_shannon = "Shannon",
  sp_simpson = "Simpson",
  sp_even = "Evenness"
)

EDGE_20M <- c(
  "0", "1", "100", "2", "200", "3", "300", "4", "400", "5", "500", "6", "600",
  "7", "700", "8", "800", "9", "900", "10", "1000", "11", "1100", "12", "1200",
  "13", "14", "1300", "15", "1400", "16", "1500", "17", "1600", "18", "1700",
  "19", "1800", "20", "1900", "21", "22", "1901", "23", "24", "1903", "124",
  "1905", "224", "1906", "324", "1907", "424", "1908", "524", "1909", "624",
  "1910", "1911", "724", "1912", "824", "1913", "924", "1914", "1024", "1915",
  "1124", "1916", "1224", "1324", "1917", "1424", "1919", "1524", "1920",
  "1624", "1921", "1724", "1922", "1824", "1923", "1924", "1904", "1902", "1918"
)

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

display_name <- function(x) {
  unname(DISPLAY_NAMES[x])
}

heatmap_name <- function(x) {
  unname(HEATMAP_NAMES[x])
}

fmt_num <- function(x, digits = 3) {
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(as.numeric(x[ok]), digits = digits, format = "f")
  out
}

fmt_p <- function(x) {
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- ifelse(
    as.numeric(x[ok]) < 0.001,
    formatC(as.numeric(x[ok]), digits = 2, format = "e"),
    formatC(as.numeric(x[ok]), digits = 3, format = "f")
  )
  out
}

read_scale_table <- function(scale_name) {
  file_path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale_name, ".csv"))
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- as.character(data$scale)
  data
}

read_sampling_metadata <- function(scale_name) {
  file_path <- file.path(
    PROJECT_DIR,
    "Quad_Values",
    paste0(scale_name, "_spectral_heterogeneity_smooth_masked_5nm_summary.csv")
  )
  metadata <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  metadata$quad_id <- as.character(metadata$quad_id)
  metadata$scale <- scale_name

  sa_all_pixels_sampled <- rep(NA, nrow(metadata))
  has_sa_metadata <- !is.na(metadata$sa_n_pixels) & !is.na(metadata$sa_method)
  sa_all_pixels_sampled[has_sa_metadata] <- metadata$sa_method[has_sa_metadata] == "exact_all_pixels" |
    (metadata$sa_method[has_sa_metadata] == "bootstrap_mean" & metadata$sa_n_pixels[has_sa_metadata] <= 5000)

  spca_euclidean_all_pixels_sampled <- rep(NA, nrow(metadata))
  has_spca_metadata <- !is.na(metadata$standardized_PCA_metric_method)
  spca_euclidean_all_pixels_sampled[has_spca_metadata] <-
    metadata$standardized_PCA_metric_method[has_spca_metadata] == "all_pixels"

  data.frame(
    quad_id = metadata$quad_id,
    scale = metadata$scale,
    sa_n_pixels = metadata$sa_n_pixels,
    sa_method = metadata$sa_method,
    sa_all_pixels_sampled = sa_all_pixels_sampled,
    spca_n_pixels = metadata$standardized_PCA_n_pixels,
    spca_metric_method = metadata$standardized_PCA_metric_method,
    spca_euclidean_all_pixels_sampled = spca_euclidean_all_pixels_sampled,
    stringsAsFactors = FALSE
  )
}

add_edge_flags <- function(data) {
  parent_20m <- data$quad_id
  is_10m <- data$scale == "10m"
  parent_20m[is_10m] <- sub("_[a-d]$", "", parent_20m[is_10m])

  edge_flag <- rep(FALSE, nrow(data))
  edge_flag[data$scale == "20m" & data$quad_id %in% EDGE_20M] <- TRUE
  edge_flag[data$scale == "10m" & parent_20m %in% EDGE_20M] <- TRUE

  data$parent_20m <- parent_20m
  data$edge_flag <- edge_flag
  data$primary_analysis <- ifelse(data$scale %in% c("10m", "20m"), !edge_flag, TRUE)
  data
}

load_analysis_data <- function() {
  data <- do.call(
    rbind,
    list(read_scale_table("10m"), read_scale_table("20m"), read_scale_table("50m"))
  )
  metadata <- do.call(
    rbind,
    list(read_sampling_metadata("10m"), read_sampling_metadata("20m"), read_sampling_metadata("50m"))
  )
  data$row_order <- seq_len(nrow(data))
  data <- merge(data, metadata, by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  data <- data[order(data$row_order), ]
  data$row_order <- NULL
  add_edge_flags(data)
}

correlate_pair <- function(data, scale_name, response, predictor) {
  keep <- data$scale == scale_name &
    data$primary_analysis &
    !is.na(data[[response]]) &
    !is.na(data[[predictor]])
  pair_df <- data[keep, c("quad_id", "scale", response, predictor)]
  names(pair_df) <- c("quad_id", "scale", "response_value", "predictor_value")

  empty_result <- data.frame(
    scale = scale_name,
    sv_measure = response,
    sv_label = display_name(response),
    diversity_measure = predictor,
    diversity_label = display_name(predictor),
    diversity_group = unname(DIVERSITY_GROUPS[predictor]),
    n = nrow(pair_df),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_df1 = NA_real_,
    f_df2 = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    slope_p_value = NA_real_,
    spearman_r = NA_real_,
    spearman_p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(pair_df) < 3) {
    return(empty_result)
  }

  fit <- lm(response_value ~ predictor_value, data = pair_df)
  fit_summary <- summary(fit)
  f_values <- fit_summary$fstatistic
  pearson <- cor.test(pair_df$predictor_value, pair_df$response_value, method = "pearson")
  spearman <- suppressWarnings(cor.test(
    pair_df$predictor_value,
    pair_df$response_value,
    method = "spearman",
    exact = FALSE
  ))
  f_p <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)

  empty_result$pearson_r <- unname(pearson$estimate)
  empty_result$r_squared <- unname(pearson$estimate)^2
  empty_result$f_statistic <- unname(f_values[["value"]])
  empty_result$f_df1 <- unname(f_values[["numdf"]])
  empty_result$f_df2 <- unname(f_values[["dendf"]])
  empty_result$f_p_value <- f_p
  empty_result$slope <- unname(coef(fit)[["predictor_value"]])
  empty_result$intercept <- unname(coef(fit)[["(Intercept)"]])
  empty_result$slope_p_value <- coef(fit_summary)["predictor_value", "Pr(>|t|)"]
  empty_result$spearman_r <- unname(spearman$estimate)
  empty_result$spearman_p_value <- spearman$p.value
  empty_result
}

run_correlations <- function(data) {
  results <- list()
  index <- 1
  for (scale_name in c("10m", "20m", "50m")) {
    for (response in SV_MEASURES) {
      for (predictor in DIVERSITY_MEASURES) {
        results[[index]] <- correlate_pair(data, scale_name, response, predictor)
        index <- index + 1
      }
    }
  }
  do.call(rbind, results)
}

make_coverage_summary <- function(data) {
  scales <- c("10m", "20m", "50m")
  rows <- lapply(scales, function(scale_name) {
    in_scale <- data$scale == scale_name
    primary <- in_scale & data$primary_analysis
    data.frame(
      scale = scale_name,
      total_quadrats = sum(in_scale),
      primary_quadrats = sum(primary),
      edge_flagged = sum(in_scale & data$edge_flag),
      complete_sa_entropy = sum(primary & !is.na(data$spec_sa)),
      complete_standardized_pca_distance = sum(primary & !is.na(data$spec_spca_mean)),
      sa_all_pixels_sampled = sum(primary & data$sa_all_pixels_sampled %in% TRUE),
      spca_euclidean_all_pixels_sampled = sum(primary & data$spca_euclidean_all_pixels_sampled %in% TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

top_results_by_sv <- function(correlations) {
  rows <- list()
  index <- 1
  for (scale_name in c("10m", "20m", "50m")) {
    for (response in SV_MEASURES) {
      subset_rows <- correlations[correlations$scale == scale_name & correlations$sv_measure == response, ]
      best <- subset_rows[which.max(abs(subset_rows$pearson_r)), ]
      rows[[index]] <- best
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

color_for_r <- function(r) {
  palette <- grDevices::colorRampPalette(c("#4B77BE", "#F7F7F2", "#B03A2E"))(201)
  idx <- round((pmax(pmin(r, 1), -1) + 1) * 100) + 1
  palette[idx]
}

save_heatmap <- function(correlations) {
  file_path <- file.path(FIGURE_DIR, "01_sv_diversity_pairwise_correlation_heatmap.png")
  grDevices::png(file_path, width = 4200, height = 1900, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(1, 3), mar = c(7.5, 5.2, 4, 1.5), oma = c(0, 0, 2, 0))
  for (scale_name in c("10m", "20m", "50m")) {
    scale_rows <- correlations[correlations$scale == scale_name, ]
    plot(
      NA,
      xlim = c(0.5, length(DIVERSITY_MEASURES) + 0.5),
      ylim = c(0.5, length(SV_MEASURES) + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = scale_name,
      bty = "n"
    )
    axis(1, at = seq_along(DIVERSITY_MEASURES), labels = heatmap_name(DIVERSITY_MEASURES), las = 2, cex.axis = 0.8)
    axis(2, at = seq_along(SV_MEASURES), labels = heatmap_name(SV_MEASURES), las = 2, cex.axis = 0.85)
    for (i in seq_along(DIVERSITY_MEASURES)) {
      for (j in seq_along(SV_MEASURES)) {
        row <- scale_rows[
          scale_rows$diversity_measure == DIVERSITY_MEASURES[i] &
            scale_rows$sv_measure == SV_MEASURES[j],
        ]
        rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col = color_for_r(row$pearson_r), border = "white", lwd = 2)
        text(
          i,
          j,
          labels = paste0("r=", fmt_num(row$pearson_r, 2), "\nR2=", fmt_num(row$r_squared, 2)),
          cex = 0.72
        )
      }
    }
  }
  mtext("Pairwise SV-diversity correlations", outer = TRUE, cex = 1.3, font = 2)
  file_path
}

save_scatter <- function(data, response) {
  file_name <- if (response == "spec_sa") {
    "02_sa_entropy_diversity_scatterplots.png"
  } else {
    "03_standardized_pca_distance_diversity_scatterplots.png"
  }
  file_path <- file.path(FIGURE_DIR, file_name)
  colors <- c("10m" = "#2F6F73", "20m" = "#A35D2D", "50m" = "#6B5CA5")

  grDevices::png(file_path, width = 3900, height = 2400, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(2, 4), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
  for (predictor in DIVERSITY_MEASURES) {
    keep <- data$primary_analysis & !is.na(data[[response]]) & !is.na(data[[predictor]])
    x <- data[[predictor]][keep]
    y <- data[[response]][keep]
    scale_values <- data$scale[keep]
    plot(
      x,
      y,
      col = adjustcolor(colors[scale_values], alpha.f = 0.45),
      pch = 16,
      cex = 0.55,
      xlab = display_name(predictor),
      ylab = display_name(response),
      main = display_name(predictor)
    )
    for (scale_name in c("10m", "20m", "50m")) {
      scale_keep <- keep & data$scale == scale_name
      if (sum(scale_keep) >= 3) {
        fit <- lm(data[[response]][scale_keep] ~ data[[predictor]][scale_keep])
        abline(fit, col = colors[scale_name], lwd = 2)
      }
    }
  }
  plot.new()
  legend("center", legend = names(colors), col = colors, pch = 16, lwd = 2, bty = "n", title = "Scale")
  mtext(paste0(display_name(response), " versus diversity measures"), outer = TRUE, cex = 1.3, font = 2)
  file_path
}

save_figures <- function(data, correlations) {
  c(
    save_heatmap(correlations),
    save_scatter(data, "spec_sa"),
    save_scatter(data, "spec_spca_mean")
  )
}

markdown_table <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  rows <- apply(data, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, separator, rows)
}

relative_path <- function(paths) {
  normalized_project <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = FALSE)
  normalized_paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalized_project), "/"), "", normalized_paths)
}

write_analysis_report <- function(correlations, coverage, figure_files) {
  top_results <- top_results_by_sv(correlations)
  top_table <- data.frame(
    Scale = top_results$scale,
    `SV measure` = top_results$sv_label,
    `Strongest diversity pairing` = top_results$diversity_label,
    Group = top_results$diversity_group,
    n = top_results$n,
    r = fmt_num(top_results$pearson_r, 3),
    R2 = fmt_num(top_results$r_squared, 3),
    F = fmt_num(top_results$f_statistic, 2),
    `F p-value` = fmt_p(top_results$f_p_value),
    check.names = FALSE
  )

  full_table <- data.frame(
    Scale = correlations$scale,
    `SV measure` = correlations$sv_label,
    `Diversity measure` = correlations$diversity_label,
    Group = correlations$diversity_group,
    n = correlations$n,
    r = fmt_num(correlations$pearson_r, 3),
    R2 = fmt_num(correlations$r_squared, 3),
    F = fmt_num(correlations$f_statistic, 2),
    `F p-value` = fmt_p(correlations$f_p_value),
    check.names = FALSE
  )

  coverage_table <- data.frame(
    Scale = coverage$scale,
    `Total quadrats` = coverage$total_quadrats,
    `Primary quadrats` = coverage$primary_quadrats,
    `Edge flagged` = coverage$edge_flagged,
    `Complete SA entropy` = coverage$complete_sa_entropy,
    `Complete standardized PCA distance` = coverage$complete_standardized_pca_distance,
    `SA all pixels sampled` = coverage$sa_all_pixels_sampled,
    `Std PCA all pixels sampled` = coverage$spca_euclidean_all_pixels_sampled,
    check.names = FALSE
  )

  report_lines <- c(
    "# SV-Diversity Pairwise Correlation Analysis",
    "",
    "Date: 2026-07-10",
    "",
    "## Research Question",
    "",
    "What is the relationship between spectral variation and phylogenetic/species diversity across 10 m, 20 m, and 50 m quadrats?",
    "",
    "This analysis is intentionally direct: each primary spectral variation measure is independently paired with each phylogenetic and species diversity measure. No environmental covariates, interaction terms, or multivariable model-selection steps are included in this first relationship layer.",
    "",
    "## Primary Spectral Variation Measures",
    "",
    "- `spec_spca_mean`: standardized PCA mean Euclidean distance in PC1-PC3 space after vector-normalizing spectra.",
    "- `spec_sa`: spectral angle entropy mean from sunlit, shadow-masked smoothed 5 nm spectra.",
    "",
    "## Diversity Measures",
    "",
    "- Phylogenetic diversity: `phy_faith`, `phy_rao`, `phy_afaith`.",
    "- Species diversity: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.",
    "",
    "## Model Form",
    "",
    "For each scale, spectral variation measure, and diversity measure, the workflow fits:",
    "",
    "`SV_measure ~ diversity_measure`",
    "",
    "Reported values include Pearson `r`, `R2`, simple-regression `F`, F-test p-value, slope, intercept, and Spearman rank correlation. In these one-predictor models, `R2` is the squared Pearson correlation.",
    "",
    "`sv_diversity_analysis_dataset.csv` also includes sampling/provenance flags from the upstream spectral heterogeneity summaries. `sa_all_pixels_sampled` is TRUE when the SA entropy workflow sampled all retained pixels for that quadrat and FALSE when the 5,000-pixel cap was used. `spca_euclidean_all_pixels_sampled` is TRUE when the standardized PCA Euclidean-distance metric used all retained pixels for that quadrat.",
    "",
    "Documented 10 m and 20 m edge quadrats are excluded from primary analysis; all 50 m quadrats are retained because no separate 50 m edge rule is documented.",
    "",
    "## Coverage",
    "",
    markdown_table(coverage_table),
    "",
    "## Strongest Pairing Per SV Measure And Scale",
    "",
    markdown_table(top_table),
    "",
    "## Full Pairwise Results",
    "",
    markdown_table(full_table),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Output Tables",
    "",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`",
    "",
    "## Superseded Analysis Direction",
    "",
    "Earlier candidate-model tables that ranked environmental and interaction models are retained only as historical context if present in the output folder. The active direction for this stage is the direct pairwise relationship between the two primary SV measures and each diversity measure."
  )

  writeLines(report_lines, ANALYSIS_REPORT)
}

write_task_report <- function(correlations, coverage) {
  task_lines <- c(
    "# Multiscale Spectral-Biodiversity Analysis",
    "",
    "Last updated: 2026-07-10",
    "",
    "## Current Direction",
    "",
    "The active analysis has been revamped around the first research question: what is the relationship between spectral variation and phylogenetic/species diversity?",
    "",
    "## Inputs",
    "",
    "- `quadrat_analysis_10m.csv`",
    "- `quadrat_analysis_20m.csv`",
    "- `quadrat_analysis_50m.csv`",
    "- `RESEARCH_OBJECTIVES.md` for scientific framing.",
    "",
    "## Primary Spectral Variation Measures",
    "",
    "- `spec_spca_mean`: standardized PCA mean Euclidean distance.",
    "- `spec_sa`: spectral angle entropy.",
    "",
    "## Diversity Predictors",
    "",
    "- Phylogenetic: `phy_faith`, `phy_rao`, `phy_afaith`.",
    "- Species: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.",
    "",
    "## Methods",
    "",
    "- Each SV measure is independently paired with each diversity measure within each spatial scale.",
    "- Each pair is summarized with Pearson `r`, `R2`, simple-regression `F`, F-test p-value, slope, intercept, and Spearman `r`.",
    "- `sv_diversity_analysis_dataset.csv` includes `sa_all_pixels_sampled` and `spca_euclidean_all_pixels_sampled`, which flag whether each primary SV metric used all retained pixels for each quadrat.",
    "- No environmental covariates or multivariable candidate-model ranking are included in this first relationship layer.",
    "- Documented 10 m and 20 m edge quadrats are excluded from primary analysis; 50 m uses all quadrats.",
    "",
    "## Outputs",
    "",
    "- `reports/analysis/20260710_sv_diversity_pairwise_correlations.md`",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_top_pairings.csv`",
    "- `reports/figures/multiscale_spectral_biodiversity/01_sv_diversity_pairwise_correlation_heatmap.png`",
    "- `reports/figures/multiscale_spectral_biodiversity/02_sa_entropy_diversity_scatterplots.png`",
    "- `reports/figures/multiscale_spectral_biodiversity/03_standardized_pca_distance_diversity_scatterplots.png`",
    "",
    "## Result Size",
    "",
    paste0("- Pairwise correlation rows: ", nrow(correlations)),
    paste0("- Analysis dataset rows: ", sum(coverage$total_quadrats)),
    "",
    "## Superseded Material",
    "",
    "The previous environment-adjusted AIC model-ranking workflow is superseded for the current stage. It can be revisited later as a second layer after the direct SV-diversity relationships are interpreted."
  )
  writeLines(task_lines, TASK_REPORT)
}

write_validation_report <- function(correlations, coverage, figure_files) {
  expected_rows <- 3 * length(SV_MEASURES) * length(DIVERSITY_MEASURES)
  table_files <- c(
    file.path(TABLE_DIR, "sv_diversity_pairwise_correlations.csv"),
    file.path(TABLE_DIR, "sv_diversity_analysis_dataset.csv"),
    file.path(TABLE_DIR, "sv_diversity_top_pairings.csv")
  )

  coverage_table <- data.frame(
    Scale = coverage$scale,
    `Total quadrats` = coverage$total_quadrats,
    `Primary quadrats` = coverage$primary_quadrats,
    `Complete SA entropy` = coverage$complete_sa_entropy,
    `Complete standardized PCA distance` = coverage$complete_standardized_pca_distance,
    `SA all pixels sampled` = coverage$sa_all_pixels_sampled,
    `Std PCA all pixels sampled` = coverage$spca_euclidean_all_pixels_sampled,
    check.names = FALSE
  )

  validation_lines <- c(
    "# Multiscale Spectral-Biodiversity Analysis Validation",
    "",
    "Last updated: 2026-07-10",
    "",
    "## Checks",
    "",
    paste0("- Expected pairwise rows: ", expected_rows),
    paste0("- Observed pairwise rows: ", nrow(correlations)),
    paste0("- Pairwise row count matches expectation: ", nrow(correlations) == expected_rows),
    paste0("- Missing Pearson r values: ", sum(is.na(correlations$pearson_r))),
    paste0("- Missing F statistics: ", sum(is.na(correlations$f_statistic))),
    paste0("- Table CSV count: ", sum(file.exists(table_files))),
    paste0("- Figure PNG count: ", sum(file.exists(figure_files))),
    "",
    "## Coverage Summary",
    "",
    markdown_table(coverage_table),
    "",
    "## Figure Files",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Table Files",
    "",
    paste0("- `", relative_path(table_files), "`"),
    "",
    "## Notes",
    "",
    "- This validation checks that the requested pairwise correlation layer is complete.",
    "- Spatial and environmental sensitivity models remain future second-layer analyses, not part of this direct correlation output."
  )

  writeLines(validation_lines, VALIDATION_REPORT)
}

run_multiscale_spectral_biodiversity_analysis <- function() {
  data <- load_analysis_data()
  missing_columns <- setdiff(c(SV_MEASURES, DIVERSITY_MEASURES), names(data))
  if (length(missing_columns) > 0) {
    stop("Missing expected columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  correlations <- run_correlations(data)
  coverage <- make_coverage_summary(data)
  top_pairings <- top_results_by_sv(correlations)
  figure_files <- save_figures(data, correlations)

  write.csv(data, file.path(TABLE_DIR, "sv_diversity_analysis_dataset.csv"), row.names = FALSE)
  write.csv(correlations, file.path(TABLE_DIR, "sv_diversity_pairwise_correlations.csv"), row.names = FALSE)
  write.csv(top_pairings, file.path(TABLE_DIR, "sv_diversity_top_pairings.csv"), row.names = FALSE)

  write_analysis_report(correlations, coverage, figure_files)
  write_task_report(correlations, coverage)
  write_validation_report(correlations, coverage, figure_files)

  message("SV-diversity pairwise analysis complete.")
  invisible(list(data = data, correlations = correlations, coverage = coverage))
}

if (sys.nframe() == 0) {
  run_multiscale_spectral_biodiversity_analysis()
}
