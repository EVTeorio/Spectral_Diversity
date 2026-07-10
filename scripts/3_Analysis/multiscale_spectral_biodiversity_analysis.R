USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c(
  "dplyr", "ggplot2", "readr", "tibble", "tidyr", "stringr",
  "broom", "gridExtra", "grid", "spdep", "scales", "png", "beepr"
)
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
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(broom)
  library(gridExtra)
  library(grid)
  library(spdep)
  library(scales)
})

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/multiscale_spectral_biodiversity")
TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/multiscale_spectral_biodiversity")
PDF_DIR <- file.path(PROJECT_DIR, "Documents/PDFs")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260624_multiscale_spectral_biodiversity_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260624_multiscale_spectral_biodiversity_analysis_validation.md")

PRIMARY_RESPONSE <- "spec_sa"
SECONDARY_RESPONSES <- c(
  "spec_pca_mean", "spec_rao_q", "spec_alpha",
  "spec_spca_mean", "spec_spca_rao", "spec_spca_alpha"
)
BIODIVERSITY_PREDICTORS <- c("phy_faith", "phy_afaith", "sp_shannon")
ENVIRONMENT_PREDICTORS <- c("env_elev", "env_tri11")
DISPLAY_NAMES <- c(
  spec_sa = "SA entropy mean",
  spec_pca_mean = "PCA mean distance",
  spec_rao_q = "Spectral Rao's Q",
  spec_alpha = "Alpha-hull area",
  spec_spca_mean = "standardized_PCA mean distance",
  spec_spca_rao = "standardized_PCA Rao's Q",
  spec_spca_alpha = "standardized_PCA alpha-hull area",
  phy_faith = "Faith's PD",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_shannon = "Shannon diversity",
  env_elev = "Elevation",
  env_tri5 = "TRI 5x5",
  env_tri11 = "TRI 11x11",
  env_tri21 = "TRI 21x21"
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

theme_report <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "#263238", size = base_size + 3),
      plot.subtitle = element_text(color = "#546E7A"),
      axis.title = element_text(color = "#263238"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", color = "#263238"),
      plot.margin = margin(9, 12, 9, 12)
    )
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PDF_DIR, recursive = TRUE, showWarnings = FALSE)

read_scale_table <- function(scale_name) {
  file_path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale_name, ".csv"))
  read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      quad_id = as.character(quad_id),
      scale = factor(scale, levels = c("10m", "20m", "50m"))
    )
}

add_edge_flags <- function(data) {
  data %>%
    mutate(
      parent_20m = if_else(
        scale == "10m",
        str_replace(quad_id, "_[a-d]$", ""),
        quad_id
      ),
      edge_flag = case_when(
        scale == "20m" & quad_id %in% EDGE_20M ~ TRUE,
        scale == "10m" & parent_20m %in% EDGE_20M ~ TRUE,
        TRUE ~ FALSE
      ),
      primary_analysis = case_when(
        scale %in% c("10m", "20m") ~ !edge_flag,
        TRUE ~ TRUE
      )
    )
}

z_score <- function(x) {
  if (all(is.na(x)) || isTRUE(sd(x, na.rm = TRUE) == 0)) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

prep_model_data <- function(data) {
  scale_vars <- unique(c(PRIMARY_RESPONSE, SECONDARY_RESPONSES, BIODIVERSITY_PREDICTORS, ENVIRONMENT_PREDICTORS))
  data %>%
    mutate(across(all_of(scale_vars), z_score, .names = "{.col}_z"))
}

model_formula <- function(response, model_type, predictor = NULL) {
  response_z <- paste0(response, "_z")
  if (model_type == "Null") {
    return(as.formula(paste(response_z, "~ 1")))
  }
  if (model_type == "Biodiversity") {
    return(as.formula(paste(response_z, "~", paste0(predictor, "_z"))))
  }
  if (model_type == "Environment") {
    return(as.formula(paste(response_z, "~ env_elev_z + env_tri11_z")))
  }
  if (model_type == "Biodiversity + environment") {
    return(as.formula(paste(response_z, "~", paste0(predictor, "_z"), "+ env_elev_z + env_tri11_z")))
  }
  if (model_type == "Biodiversity x elevation") {
    return(as.formula(paste(response_z, "~", paste0(predictor, "_z"), "* env_elev_z + env_tri11_z")))
  }
  stop("Unknown model type: ", model_type, call. = FALSE)
}

fit_one_model <- function(data, scale_name, response, model_type, predictor = NA_character_) {
  formula <- model_formula(response, model_type, predictor)
  model_vars <- all.vars(formula)
  model_df <- data %>%
    filter(scale == scale_name, primary_analysis) %>%
    select(quad_id, center_x, center_y, all_of(model_vars)) %>%
    tidyr::drop_na()

  if (nrow(model_df) < length(model_vars) + 8) {
    return(NULL)
  }

  fit <- lm(formula, data = model_df)
  glance_fit <- broom::glance(fit)
  coefficients <- broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      scale = scale_name,
      response = response,
      model_type = model_type,
      predictor = predictor,
      n = nrow(model_df),
      r_squared = glance_fit$r.squared,
      adj_r_squared = glance_fit$adj.r.squared,
      aic = AIC(fit),
      bic = BIC(fit),
      residual_sd = glance_fit$sigma
    )

  list(
    fit = fit,
    model_df = model_df,
    coefficients = coefficients,
    summary = tibble(
      scale = scale_name,
      response = response,
      model_type = model_type,
      predictor = predictor,
      n = nrow(model_df),
      r_squared = glance_fit$r.squared,
      adj_r_squared = glance_fit$adj.r.squared,
      aic = AIC(fit),
      bic = BIC(fit),
      residual_sd = glance_fit$sigma
    )
  )
}

fit_model_grid <- function(data) {
  responses <- c(PRIMARY_RESPONSE, SECONDARY_RESPONSES)
  scales <- levels(data$scale)
  fits <- list()
  index <- 1

  for (scale_name in scales) {
    for (response in responses) {
      fits[[index]] <- fit_one_model(data, scale_name, response, "Null")
      index <- index + 1
      fits[[index]] <- fit_one_model(data, scale_name, response, "Environment")
      index <- index + 1
      for (predictor in BIODIVERSITY_PREDICTORS) {
        fits[[index]] <- fit_one_model(data, scale_name, response, "Biodiversity", predictor)
        index <- index + 1
        fits[[index]] <- fit_one_model(data, scale_name, response, "Biodiversity + environment", predictor)
        index <- index + 1
        if (response == PRIMARY_RESPONSE) {
          fits[[index]] <- fit_one_model(data, scale_name, response, "Biodiversity x elevation", predictor)
          index <- index + 1
        }
      }
    }
  }

  fits <- fits[!vapply(fits, is.null, logical(1))]
  summaries <- bind_rows(lapply(fits, `[[`, "summary")) %>%
    group_by(scale, response) %>%
    mutate(delta_aic = aic - min(aic, na.rm = TRUE)) %>%
    ungroup()
  coefficients <- bind_rows(lapply(fits, `[[`, "coefficients"))

  list(fits = fits, summaries = summaries, coefficients = coefficients)
}

correlation_tables <- function(data) {
  combos <- expand_grid(
    scale = levels(data$scale),
    response = c(PRIMARY_RESPONSE, SECONDARY_RESPONSES),
    predictor = c(BIODIVERSITY_PREDICTORS, ENVIRONMENT_PREDICTORS)
  )

  bind_rows(lapply(seq_len(nrow(combos)), function(i) {
    scale_i <- combos$scale[[i]]
    response_i <- combos$response[[i]]
    predictor_i <- combos$predictor[[i]]
    tmp <- data %>%
      filter(.data$scale == scale_i, primary_analysis) %>%
      select(all_of(c(response_i, predictor_i))) %>%
      tidyr::drop_na()

    tibble(
      scale = scale_i,
      response = response_i,
      predictor = predictor_i,
      n = nrow(tmp),
      pearson_r = if (nrow(tmp) > 3) cor(tmp[[response_i]], tmp[[predictor_i]], method = "pearson") else NA_real_,
      spearman_r = if (nrow(tmp) > 3) cor(tmp[[response_i]], tmp[[predictor_i]], method = "spearman") else NA_real_,
      pearson_p = if (nrow(tmp) > 3) cor.test(tmp[[response_i]], tmp[[predictor_i]], method = "pearson")$p.value else NA_real_
    )
  }))
}

moran_diagnostic <- function(model_item, k = 8) {
  model_df <- model_item$model_df
  if (nrow(model_df) <= k + 3) {
    return(NULL)
  }
  coords <- as.matrix(model_df[, c("center_x", "center_y")])
  knn <- spdep::knearneigh(coords, k = min(k, nrow(model_df) - 1))
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  test <- spdep::moran.test(residuals(model_item$fit), lw, zero.policy = TRUE)
  tibble(
    scale = model_item$summary$scale,
    response = model_item$summary$response,
    model_type = model_item$summary$model_type,
    predictor = model_item$summary$predictor,
    n = nrow(model_df),
    moran_i = unname(test$estimate[["Moran I statistic"]]),
    expected_i = unname(test$estimate[["Expectation"]]),
    p_value = test$p.value
  )
}

choose_primary_models <- function(model_results) {
  model_results$summaries %>%
    filter(response == PRIMARY_RESPONSE, model_type != "Null") %>%
    group_by(scale) %>%
    arrange(aic, desc(adj_r_squared), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}

save_plot <- function(plot, filename, width = 11, height = 7, dpi = 320) {
  path <- file.path(FIGURE_DIR, filename)
  ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  path
}

make_figures <- function(data, correlations, model_results, moran_results, primary_models) {
  figures <- list()

  coverage <- data %>%
    group_by(scale) %>%
    summarize(
      total_quadrats = n(),
      primary_quadrats = sum(primary_analysis),
      complete_spec_sa = sum(primary_analysis & !is.na(spec_sa)),
      complete_pca_mean = sum(primary_analysis & !is.na(spec_pca_mean)),
      edge_flagged = sum(edge_flag),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(primary_quadrats, complete_spec_sa, complete_pca_mean, edge_flagged),
      names_to = "metric",
      values_to = "count"
    ) %>%
    mutate(metric = recode(
      metric,
      primary_quadrats = "Primary analysis quadrats",
      complete_spec_sa = "Complete SA entropy mean",
      complete_pca_mean = "Complete PCA mean distance",
      edge_flagged = "Edge flagged"
    ))

  p_coverage <- ggplot(coverage, aes(x = scale, y = count, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = c("#1B4D89", "#4F9D69", "#D95F02", "#8D99AE")) +
    labs(
      title = "Analysis coverage by scale",
      subtitle = "Primary analysis excludes documented 10 m and 20 m edge quadrats; 50 m has no separate edge rule.",
      x = "Quadrat scale",
      y = "Quadrat count",
      fill = NULL
    ) +
    theme_report()
  figures$coverage <- save_plot(p_coverage, "01_analysis_coverage_by_scale.png")

  dist_data <- data %>%
    filter(primary_analysis) %>%
    select(scale, all_of(c(PRIMARY_RESPONSE, SECONDARY_RESPONSES))) %>%
    pivot_longer(-scale, names_to = "metric_id", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(metric = recode(metric_id, !!!as.list(DISPLAY_NAMES)))

  p_dist_combined <- ggplot(dist_data, aes(x = value, fill = scale)) +
    geom_histogram(bins = 28, color = "white", alpha = 0.85) +
    facet_grid(metric ~ scale, scales = "free") +
    scale_fill_manual(values = c("10m" = "#1B4D89", "20m" = "#4F9D69", "50m" = "#D95F02")) +
    labs(
      title = "Spectral heterogeneity distributions",
      subtitle = "SA entropy is usually the mean of 70 bootstrap iterations; secondary metrics are shown for primary-analysis quadrats.",
      x = "Metric value",
      y = "Count",
      fill = "Scale"
    ) +
    theme_report(10) +
    theme(legend.position = "none")
  figures$combined_distribution <- save_plot(p_dist_combined, "02_spectral_metric_distributions.png", height = 9)

  distribution_files <- list()
  distribution_order <- c(PRIMARY_RESPONSE, SECONDARY_RESPONSES)
  distribution_names <- c(
    spec_sa = "02a_distribution_sa_entropy_mean.png",
    spec_pca_mean = "02b_distribution_pca_mean_distance.png",
    spec_rao_q = "02c_distribution_spectral_rao_q.png",
    spec_alpha = "02d_distribution_alpha_hull_area.png",
    spec_spca_mean = "02e_distribution_standardized_pca_mean_distance.png",
    spec_spca_rao = "02f_distribution_standardized_pca_rao_q.png",
    spec_spca_alpha = "02g_distribution_standardized_pca_alpha_hull_area.png"
  )
  for (metric_id in distribution_order) {
    metric_label <- DISPLAY_NAMES[[metric_id]]
    metric_note <- if (metric_id == PRIMARY_RESPONSE) {
      "Most quadrat values are means of 70 bootstrap iterations using up to 5,000 retained pixels; the small exact subset used all retained pixel pairs."
    } else {
      "Distributions are shown for primary-analysis quadrats with non-missing values."
    }
    p_dist <- ggplot(dist_data %>% filter(metric_id == !!metric_id), aes(x = value, fill = scale)) +
      geom_histogram(bins = 28, color = "white", alpha = 0.85) +
      facet_wrap(~ scale, scales = "free") +
      scale_fill_manual(values = c("10m" = "#1B4D89", "20m" = "#4F9D69", "50m" = "#D95F02")) +
      labs(
        title = paste(metric_label, "distribution by scale"),
        subtitle = metric_note,
        x = metric_label,
        y = "Count",
        fill = "Scale"
      ) +
      theme_report(10) +
      theme(legend.position = "none")
    distribution_files[[metric_id]] <- save_plot(
      p_dist,
      distribution_names[[metric_id]],
      height = 5.8
    )
  }
  figures$distributions <- distribution_files

  p_map <- ggplot(data %>% filter(primary_analysis, !is.na(spec_sa)), aes(center_x, center_y, color = spec_sa)) +
    geom_point(size = 1.55, alpha = 0.88) +
    facet_wrap(~ scale) +
    coord_equal() +
    scale_color_viridis_c(option = "C", name = "SA entropy mean") +
    labs(
      title = "Spatial pattern of SA entropy means",
      subtitle = "Most SA values are means of 70 bootstrap iterations; centroid maps use plant-diversity quadrat centroids in NAD83 / UTM zone 16N.",
      x = "UTM easting",
      y = "UTM northing"
    ) +
    theme_report()
  figures$spatial_sa <- save_plot(p_map, "03_spatial_pattern_spec_sa.png", height = 7.5)

  scatter_for <- function(pred, file, title) {
    ggplot(data %>% filter(primary_analysis), aes(.data[[pred]], spec_sa)) +
      geom_point(aes(color = env_elev), alpha = 0.58, size = 1.7, na.rm = TRUE) +
      geom_smooth(method = "lm", se = TRUE, color = "#263238", linewidth = 0.9, na.rm = TRUE) +
      facet_wrap(~ scale, scales = "free") +
      scale_color_viridis_c(option = "D", name = "Elevation") +
      labs(
        title = title,
        subtitle = "Lines show simple linear fits within scale; most SA values are bootstrap means, with exact entropy only for small quadrats.",
        x = DISPLAY_NAMES[[pred]],
        y = "SA entropy mean"
      ) +
      theme_report()
  }
  figures$faith_scatter <- save_plot(
    scatter_for("phy_faith", "04_spec_sa_vs_faith_pd.png", "SA entropy mean versus Faith's phylogenetic diversity"),
    "04_spec_sa_vs_faith_pd.png"
  )
  figures$afaith_scatter <- save_plot(
    scatter_for("phy_afaith", "05_spec_sa_vs_abundance_weighted_faith_pd.png", "SA entropy mean versus abundance-weighted Faith's PD"),
    "05_spec_sa_vs_abundance_weighted_faith_pd.png"
  )
  figures$shannon_scatter <- save_plot(
    scatter_for("sp_shannon", "06_spec_sa_vs_shannon_diversity.png", "SA entropy mean versus Shannon diversity"),
    "06_spec_sa_vs_shannon_diversity.png"
  )

  p_corr <- correlations %>%
    filter(response == PRIMARY_RESPONSE, predictor %in% c(BIODIVERSITY_PREDICTORS, ENVIRONMENT_PREDICTORS)) %>%
    mutate(
      predictor_label = recode(predictor, !!!as.list(DISPLAY_NAMES)),
      sig = case_when(
        pearson_p < 0.001 ~ "***",
        pearson_p < 0.01 ~ "**",
        pearson_p < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    ggplot(aes(x = predictor_label, y = scale, fill = pearson_r)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = paste0(sprintf("%.2f", pearson_r), sig)), size = 3.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
    labs(
      title = "Primary correlations with SA entropy means",
      subtitle = "Pearson r values for primary-analysis quadrats. Asterisks mark p < 0.05, 0.01, and 0.001.",
      x = NULL,
      y = "Scale",
      fill = "Pearson r"
    ) +
    theme_report() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  figures$correlation_heatmap <- save_plot(p_corr, "07_primary_correlation_heatmap.png", height = 6.2)

  p_r2 <- model_results$summaries %>%
    filter(response == PRIMARY_RESPONSE, model_type != "Null") %>%
    mutate(
      model_label = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      model_label = str_replace(model_label, "Biodiversity: ", ""),
      model_label = str_replace(model_label, "Biodiversity \\+ environment: ", "+ env: "),
      model_label = str_replace(model_label, "Biodiversity x elevation: ", "x elev: ")
    ) %>%
    ggplot(aes(x = reorder(model_label, adj_r_squared), y = adj_r_squared, fill = model_type)) +
    geom_col(width = 0.72) +
    coord_flip() +
    facet_wrap(~ scale, scales = "free_y") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c(
      "Biodiversity" = "#1B4D89",
      "Environment" = "#7A9E3F",
      "Biodiversity + environment" = "#D95F02",
      "Biodiversity x elevation" = "#6A4C93"
    )) +
    labs(
      title = "Model explanatory strength for SA entropy means",
      subtitle = "Adjusted R2 compares standardized OLS candidate models within each scale.",
      x = NULL,
      y = "Adjusted R2",
      fill = "Model family"
    ) +
    theme_report(9)
  figures$model_r2 <- save_plot(p_r2, "08_model_adjusted_r2_comparison.png", height = 8.5)

  p_coef <- model_results$coefficients %>%
    filter(
      response == PRIMARY_RESPONSE,
      model_type == "Biodiversity + environment",
      term %in% paste0(BIODIVERSITY_PREDICTORS, "_z")
    ) %>%
    mutate(
      predictor_label = recode(str_remove(term, "_z$"), !!!as.list(DISPLAY_NAMES)),
      significant_label = if_else(p.value < 0.05, "Yes", "No")
    ) %>%
    ggplot(aes(x = estimate, y = predictor_label, color = significant_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_wrap(~ scale) +
    scale_color_manual(values = c("No" = "#455A64", "Yes" = "#B2182B"), drop = FALSE) +
    labs(
      title = "Standardized biodiversity effects after environmental controls",
      subtitle = "Points are standardized coefficients; bars are 95 percent confidence intervals.",
      x = "Standardized coefficient",
      y = NULL,
      color = "p < 0.05"
    ) +
    theme_report()
  figures$coefficient_forest <- save_plot(p_coef, "09_primary_coefficient_forest.png", height = 6.5)

  p_moran <- moran_results %>%
    filter(response == PRIMARY_RESPONSE) %>%
    mutate(
      label = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      significant_label = if_else(p_value < 0.05, "Yes", "No")
    ) %>%
    ggplot(aes(x = moran_i, y = reorder(label, moran_i), color = significant_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 2.2) +
    facet_wrap(~ scale, scales = "free_y") +
    scale_color_manual(values = c("No" = "#455A64", "Yes" = "#B2182B"), drop = FALSE) +
    labs(
      title = "Spatial autocorrelation in model residuals",
      subtitle = "Moran's I is calculated from 8-nearest-neighbor residuals. Significant residual structure cautions against over-interpreting p-values.",
      x = "Residual Moran's I",
      y = NULL,
      color = "p < 0.05"
    ) +
    theme_report(9)
  figures$residual_moran <- save_plot(p_moran, "10_residual_moran_diagnostics.png", height = 8.5)

  p_scale <- data %>%
    filter(primary_analysis) %>%
    select(scale, all_of(c(PRIMARY_RESPONSE, BIODIVERSITY_PREDICTORS, ENVIRONMENT_PREDICTORS))) %>%
    pivot_longer(-scale, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>%
    mutate(value_z = z_score(value)) %>%
    ungroup() %>%
    mutate(metric = recode(metric, !!!as.list(DISPLAY_NAMES))) %>%
    ggplot(aes(scale, value_z, fill = scale)) +
    geom_boxplot(outlier.alpha = 0.22, width = 0.65) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("10m" = "#1B4D89", "20m" = "#4F9D69", "50m" = "#D95F02")) +
    labs(
      title = "Scale comparison of standardized variables",
      subtitle = "Values are standardized within each variable to make scale shifts visually comparable.",
      x = "Quadrat scale",
      y = "Standardized value"
    ) +
    theme_report(9) +
    theme(legend.position = "none")
  figures$scale_boxplots <- save_plot(p_scale, "11_scale_comparison_standardized_variables.png", height = 8.5)

  figures
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

build_findings <- function(data, correlations, model_results, primary_models, moran_results) {
  best_rows <- primary_models %>%
    mutate(
      model_readable = if_else(
        is.na(predictor),
        model_type,
        paste0(model_type, " using ", recode(predictor, !!!as.list(DISPLAY_NAMES)))
      )
    )

  corr_faith <- correlations %>%
    filter(response == PRIMARY_RESPONSE, predictor == "phy_faith") %>%
    mutate(direction = if_else(pearson_r >= 0, "positive", "negative"))

  list(
    title = "Paint Rock spectral-biodiversity relationships",
    summary = c(
      "The analysis evaluates whether quadrat-level hyperspectral heterogeneity tracks biodiversity and phylogenetic diversity across 10 m, 20 m, and 50 m grains.",
      "SA entropy mean is treated as the primary response. For most quadrats it is the mean of 70 bootstrap iterations using up to 5,000 retained pixels; the small exact subset uses all retained pixel pairs. PCA mean distance, spectral Rao's Q, alpha-hull area, and their standardized_PCA analogs are secondary sensitivity metrics.",
      "Models use standardized predictors so coefficients are comparable within scale. Documented 10 m and 20 m edge quadrats are excluded from primary inference; 50 m uses all quadrats because no separate 50 m edge rule is documented."
    ),
    best_model_table = best_rows,
    faith_correlation = corr_faith,
    residual_caution = moran_results %>%
      filter(response == PRIMARY_RESPONSE, p_value < 0.05) %>%
      count(scale, name = "significant_residual_models")
  )
}

text_page <- function(title, body_lines, footer = NULL, title_size = 20, body_size = 11) {
  grid.newpage()
  pushViewport(viewport(width = 0.88, height = 0.86))
  grid.text(title, x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = title_size, fontface = "bold", col = "#263238"))
  grid.text(
    paste(body_lines, collapse = "\n\n"),
    x = 0,
    y = 0.88,
    just = c("left", "top"),
    gp = gpar(fontsize = body_size, col = "#263238", lineheight = 1.1)
  )
  if (!is.null(footer)) {
    grid.text(footer, x = 1, y = 0, just = c("right", "bottom"), gp = gpar(fontsize = 8.5, col = "#607D8B"))
  }
  popViewport()
}

table_page <- function(title, table_data, footer = NULL, rows = 18) {
  grid.newpage()
  pushViewport(viewport(width = 0.92, height = 0.88))
  grid.text(title, x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 17, fontface = "bold", col = "#263238"))
  if (nrow(table_data) == 0) {
    grid.text("No rows available.", x = 0, y = 0.88, just = c("left", "top"))
  } else {
    table_data <- head(table_data, rows)
    tg <- tableGrob(table_data, rows = NULL, theme = ttheme_minimal(
      core = list(fg_params = list(fontsize = 8.2), padding = unit(c(3, 3), "mm")),
      colhead = list(fg_params = list(fontsize = 8.5, fontface = "bold", col = "white"),
                     bg_params = list(fill = "#263238", col = NA))
    ))
    grid.draw(editGrob(tg, vp = viewport(x = 0.5, y = 0.45, width = 1, height = 0.78)))
  }
  if (!is.null(footer)) {
    grid.text(footer, x = 1, y = 0, just = c("right", "bottom"), gp = gpar(fontsize = 8.5, col = "#607D8B"))
  }
  popViewport()
}

image_page <- function(title, image_path, footer = NULL) {
  grid.newpage()
  pushViewport(viewport(width = 0.92, height = 0.9))
  grid.text(title, x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 17, fontface = "bold", col = "#263238"))
  img <- png::readPNG(image_path)
  grid.raster(img, x = 0.5, y = 0.48, width = 0.96, height = 0.78, interpolate = TRUE)
  if (!is.null(footer)) {
    grid.text(footer, x = 1, y = 0, just = c("right", "bottom"), gp = gpar(fontsize = 8.5, col = "#607D8B"))
  }
  popViewport()
}

write_main_pdf <- function(findings, figures, coverage_table, primary_models, correlations, moran_results) {
  pdf_path <- file.path(PDF_DIR, "spectral_biodiversity_multiscale_findings.pdf")
  pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  text_page(
    findings$title,
    c(
      "Main question: can hyperspectral spectral heterogeneity serve as a proxy for biodiversity and phylogenetic diversity in the Paint Rock Forest Dynamics Plot?",
      findings$summary,
      "Primary response: SA entropy mean (`spec_sa`). For most quadrats it is the mean of 70 bootstrap iterations using up to 5,000 retained pixels; the small exact subset uses all retained pixel pairs. Primary predictors: Faith's PD, abundance-weighted Faith's PD, Shannon diversity, elevation, and 11x11 Riley TRI.",
      "Outputs in this report are descriptive and inferential OLS summaries with spatial residual checks. When residual Moran's I remains significant, p-values should be treated cautiously."
    ),
    footer = "Generated from current quadrat_analysis_*m.csv tables"
  )

  coverage_display <- coverage_table %>%
    mutate(across(where(is.numeric), as.integer)) %>%
    rename(
      Scale = scale,
      `All quadrats` = total_quadrats,
      `Primary n` = primary_quadrats,
      `Complete SA` = complete_spec_sa,
      `Complete PCA mean` = complete_pca_mean,
      `Edge flagged` = edge_flagged
    )
  table_page("Analysis sample sizes", coverage_display)

  image_page("Coverage and data completeness", figures$coverage)
  image_page("Spatial pattern of primary spectral heterogeneity", figures$spatial_sa)
  image_page("SA entropy mean versus Faith's PD", figures$faith_scatter)
  image_page("SA entropy mean versus abundance-weighted Faith's PD", figures$afaith_scatter)
  image_page("SA entropy mean versus Shannon diversity", figures$shannon_scatter)
  image_page("Primary correlations", figures$correlation_heatmap)
  image_page("Model comparison", figures$model_r2)
  image_page("Biodiversity effects after environmental controls", figures$coefficient_forest)
  image_page("Residual spatial autocorrelation", figures$residual_moran)

  best_display <- primary_models %>%
    transmute(
      Scale = scale,
      Model = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      n = n,
      `Adj. R2` = sprintf("%.3f", adj_r_squared),
      AIC = sprintf("%.1f", aic),
      `Delta AIC` = sprintf("%.1f", delta_aic)
    )
  table_page("Best-supported primary-response model by scale", best_display)

  primary_corr_display <- correlations %>%
    filter(response == PRIMARY_RESPONSE, predictor %in% BIODIVERSITY_PREDICTORS) %>%
    transmute(
      Scale = scale,
      Predictor = recode(predictor, !!!as.list(DISPLAY_NAMES)),
      n = n,
      `Pearson r` = sprintf("%.3f", pearson_r),
      `Spearman r` = sprintf("%.3f", spearman_r),
      p = format_p(pearson_p)
    )
  table_page("Primary biodiversity correlations", primary_corr_display)

  moran_display <- moran_results %>%
    filter(response == PRIMARY_RESPONSE) %>%
    arrange(scale, p_value) %>%
    transmute(
      Scale = scale,
      Model = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      n = n,
      `Moran I` = sprintf("%.3f", moran_i),
      p = format_p(p_value)
    )
  table_page("Spatial residual diagnostic summary", moran_display, rows = 24)

  text_page(
    "Interpretation guardrails",
    c(
      "1. A positive biodiversity coefficient supports the spectral-variation hypothesis only when its confidence interval is mostly above zero and residual spatial structure is weak or acknowledged.",
      "2. Environmental controls are included because topography can affect both canopy composition and observed reflectance. Models that retain biodiversity signal after elevation/TRI controls are stronger evidence for biodiversity-spectral linkage.",
      "3. The 50 m scale is useful for scale dependence, but sample size is small. Large coefficients at 50 m should be read alongside confidence intervals and residual diagnostics.",
      "4. Known spectral missingness and PCA exclusion rules are preserved. No spectral values were imputed.",
      "5. `spec_sa` is usually the mean of 70 bootstrap iterations using up to 5,000 retained pixels; only the small exact subset used all retained pixel pairs. Bootstrap QC fields are not in the combined tables, so sensitivity checks using `boot_cv` and bootstrap CI width should be added before final manuscript inference."
    )
  )

  invisible(pdf_path)
}

write_appendix_pdf <- function(figures, model_results, correlations, moran_results) {
  pdf_path <- file.path(PDF_DIR, "spectral_biodiversity_model_appendix.pdf")
  pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  text_page(
    "Model appendix",
    c(
      "This appendix contains supplemental figures and detailed model tables for the multiscale spectral-biodiversity analysis.",
      "All models use standardized responses and predictors. Candidate sets include null, single biodiversity, environment-only, biodiversity plus environment, and biodiversity-by-elevation models for the primary SA entropy mean response.",
      "Tables are written in full to `reports/tables/multiscale_spectral_biodiversity/`."
    )
  )

  for (metric_id in names(figures$distributions)) {
    image_page(paste(DISPLAY_NAMES[[metric_id]], "distribution"), figures$distributions[[metric_id]])
  }
  image_page("Scale comparison of standardized variables", figures$scale_boxplots)
  image_page("Primary correlations", figures$correlation_heatmap)
  image_page("Adjusted R2 comparison", figures$model_r2)
  image_page("Coefficient forest", figures$coefficient_forest)
  image_page("Residual Moran's I diagnostics", figures$residual_moran)

  top_models <- model_results$summaries %>%
    arrange(response, scale, delta_aic) %>%
    group_by(response, scale) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    transmute(
      Scale = scale,
      Response = recode(response, !!!as.list(DISPLAY_NAMES)),
      Model = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      n = n,
      `Adj. R2` = sprintf("%.3f", adj_r_squared),
      AIC = sprintf("%.1f", aic),
      `Delta AIC` = sprintf("%.1f", delta_aic)
    )
  table_page("Top candidate models by response and scale", top_models, rows = 26)

  coef_display <- model_results$coefficients %>%
    filter(term != "(Intercept)", response == PRIMARY_RESPONSE) %>%
    arrange(scale, model_type, predictor, term) %>%
    transmute(
      Scale = scale,
      Model = if_else(is.na(predictor), model_type, paste(model_type, recode(predictor, !!!as.list(DISPLAY_NAMES)), sep = ": ")),
      Term = str_replace_all(term, "_z", ""),
      Estimate = sprintf("%.3f", estimate),
      SE = sprintf("%.3f", std.error),
      `CI low` = sprintf("%.3f", conf.low),
      `CI high` = sprintf("%.3f", conf.high),
      p = format_p(p.value)
    )
  table_page("Primary-response coefficients", coef_display, rows = 24)

  corr_display <- correlations %>%
    transmute(
      Scale = scale,
      Response = recode(response, !!!as.list(DISPLAY_NAMES)),
      Predictor = recode(predictor, !!!as.list(DISPLAY_NAMES)),
      n = n,
      `Pearson r` = sprintf("%.3f", pearson_r),
      `Spearman r` = sprintf("%.3f", spearman_r),
      p = format_p(pearson_p)
    )
  table_page("Correlation table excerpt", corr_display, rows = 24)

  invisible(pdf_path)
}

write_reports <- function(outputs) {
  coverage_table <- outputs$data %>%
    group_by(scale) %>%
    summarize(
      total_quadrats = n(),
      primary_quadrats = sum(primary_analysis),
      complete_spec_sa = sum(primary_analysis & !is.na(spec_sa)),
      complete_pca_mean = sum(primary_analysis & !is.na(spec_pca_mean)),
      edge_flagged = sum(edge_flag),
      .groups = "drop"
    )

  findings <- build_findings(
    outputs$data,
    outputs$correlations,
    outputs$model_results,
    outputs$primary_models,
    outputs$moran_results
  )

  main_pdf <- write_main_pdf(
    findings,
    outputs$figures,
    coverage_table,
    outputs$primary_models,
    outputs$correlations,
    outputs$moran_results
  )
  appendix_pdf <- write_appendix_pdf(
    outputs$figures,
    outputs$model_results,
    outputs$correlations,
    outputs$moran_results
  )

  c(main_pdf = main_pdf, appendix_pdf = appendix_pdf)
}

write_markdown_reports <- function(outputs, pdf_paths) {
  coverage_table <- outputs$data %>%
    group_by(scale) %>%
    summarize(
      total_quadrats = n(),
      primary_quadrats = sum(primary_analysis),
      complete_spec_sa = sum(primary_analysis & !is.na(spec_sa)),
      complete_pca_mean = sum(primary_analysis & !is.na(spec_pca_mean)),
      edge_flagged = sum(edge_flag),
      .groups = "drop"
    )

  task_lines <- c(
    "# Multiscale Spectral-Biodiversity Analysis",
    "",
    "Last updated: 2026-06-24",
    "",
    "## Task",
    "",
    "Created a reproducible downstream analysis and PDF reporting workflow for the current combined quadrat analysis tables.",
    "",
    "## Inputs",
    "",
    "- `quadrat_analysis_10m.csv`",
    "- `quadrat_analysis_20m.csv`",
    "- `quadrat_analysis_50m.csv`",
    "- `scripts/3_Analysis/LLM.R` and `scripts/3_Analysis/Analysis_PDF.R` were used as style and method references.",
    "- `RESEARCH_OBJECTIVES.md` supplied the scientific questions and hypothesis framing.",
    "",
    "## Outputs",
    "",
    paste0("- `", normalizePath(pdf_paths[["main_pdf"]], winslash = "/", mustWork = FALSE), "`"),
    paste0("- `", normalizePath(pdf_paths[["appendix_pdf"]], winslash = "/", mustWork = FALSE), "`"),
    "- `reports/figures/multiscale_spectral_biodiversity/`",
    "- `reports/tables/multiscale_spectral_biodiversity/`",
    "",
    "## Methods",
    "",
    "- Primary response: `spec_sa`, the SA entropy mean. For most quadrats it is the mean of 70 bootstrap iterations using up to 5,000 retained pixels; the small exact subset uses all retained pixel pairs.",
    "- Secondary responses: `spec_pca_mean`, `spec_rao_q`, `spec_alpha`, `spec_spca_mean`, `spec_spca_rao`, and `spec_spca_alpha`.",
    "- Primary biodiversity predictors: `phy_faith`, `phy_afaith`, and `sp_shannon`.",
    "- Environmental controls: `env_elev` and `env_tri11`.",
    "- Predictors and responses were standardized within model datasets.",
    "- Candidate OLS models were compared using adjusted R2 and AIC.",
    "- Residual spatial autocorrelation was evaluated with 8-nearest-neighbor Moran's I.",
    "- Documented 10 m and 20 m edge quadrats were excluded from primary analysis; 50 m used all quadrats because no separate 50 m edge rule is documented.",
    "",
    "## Coverage Summary",
    "",
    "| Scale | Total quadrats | Primary quadrats | Complete SA mean | Complete PCA mean | Edge flagged |",
    "|---|---:|---:|---:|---:|---:|",
    apply(coverage_table, 1, function(row) {
      paste0(
        "| ", row[["scale"]], " | ",
        row[["total_quadrats"]], " | ",
        row[["primary_quadrats"]], " | ",
        row[["complete_spec_sa"]], " | ",
        row[["complete_pca_mean"]], " | ",
        row[["edge_flagged"]], " |"
      )
    })
  )
  writeLines(task_lines, TASK_REPORT)

  fig_files <- list.files(FIGURE_DIR, pattern = "\\.png$", full.names = FALSE)
  table_files <- list.files(TABLE_DIR, pattern = "\\.csv$", full.names = FALSE)
  validation_lines <- c(
    "# Multiscale Spectral-Biodiversity Analysis Validation",
    "",
    "Last updated: 2026-06-24",
    "",
    "## Checks",
    "",
    paste0("- Main PDF exists: ", file.exists(pdf_paths[["main_pdf"]])),
    paste0("- Appendix PDF exists: ", file.exists(pdf_paths[["appendix_pdf"]])),
    paste0("- Figure PNG count: ", length(fig_files)),
    paste0("- Table CSV count: ", length(table_files)),
    paste0("- Model summary rows: ", nrow(outputs$model_results$summaries)),
    paste0("- Coefficient rows: ", nrow(outputs$model_results$coefficients)),
    paste0("- Correlation rows: ", nrow(outputs$correlations)),
    paste0("- Moran diagnostic rows: ", nrow(outputs$moran_results)),
    "",
    "## Figure Files",
    "",
    paste0("- `reports/figures/multiscale_spectral_biodiversity/", fig_files, "`"),
    "",
    "## Table Files",
    "",
    paste0("- `reports/tables/multiscale_spectral_biodiversity/", table_files, "`"),
    "",
    "## Notes",
    "",
    "- PDF rendering should be checked with Poppler page previews after script execution.",
    "- Bootstrap quality-control sensitivity checks remain recommended before final manuscript inference."
  )
  writeLines(validation_lines, VALIDATION_REPORT)
}

run_multiscale_spectral_biodiversity_analysis <- function() {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(PROJECT_DIR)

  data <- bind_rows(
    read_scale_table("10m"),
    read_scale_table("20m"),
    read_scale_table("50m")
  ) %>%
    add_edge_flags() %>%
    prep_model_data()

  model_results <- fit_model_grid(data)
  primary_models <- choose_primary_models(model_results)
  correlations <- correlation_tables(data)

  primary_fit_keys <- primary_models %>%
    select(scale, response, model_type, predictor)
  primary_fit_indices <- vapply(model_results$fits, function(item) {
    any(
      primary_fit_keys$scale == item$summary$scale &
        primary_fit_keys$response == item$summary$response &
        primary_fit_keys$model_type == item$summary$model_type &
        (
          (is.na(primary_fit_keys$predictor) & is.na(item$summary$predictor)) |
            primary_fit_keys$predictor == item$summary$predictor
        )
    )
  }, logical(1))
  diagnostic_fits <- c(
    model_results$fits[primary_fit_indices],
    model_results$fits[vapply(model_results$fits, function(item) {
      item$summary$response == PRIMARY_RESPONSE && item$summary$model_type %in% c("Null", "Environment")
    }, logical(1))]
  )
  moran_results <- bind_rows(lapply(diagnostic_fits, moran_diagnostic))

  figures <- make_figures(data, correlations, model_results, moran_results, primary_models)

  readr::write_csv(data, file.path(TABLE_DIR, "analysis_dataset_with_flags.csv"))
  readr::write_csv(correlations, file.path(TABLE_DIR, "correlation_results.csv"))
  readr::write_csv(model_results$summaries, file.path(TABLE_DIR, "model_summary_results.csv"))
  readr::write_csv(model_results$coefficients, file.path(TABLE_DIR, "model_coefficient_results.csv"))
  readr::write_csv(primary_models, file.path(TABLE_DIR, "primary_best_models_by_scale.csv"))
  readr::write_csv(moran_results, file.path(TABLE_DIR, "residual_moran_diagnostics.csv"))

  outputs <- list(
    data = data,
    correlations = correlations,
    model_results = model_results,
    primary_models = primary_models,
    moran_results = moran_results,
    figures = figures
  )

  pdf_paths <- write_reports(outputs)
  write_markdown_reports(outputs, pdf_paths)

  if (requireNamespace("beepr", quietly = TRUE)) {
    beepr::beep(3)
  }

  print(pdf_paths)
  invisible(c(outputs, list(pdf_paths = pdf_paths)))
}

if (sys.nframe() == 0) {
  run_multiscale_spectral_biodiversity_analysis()
}
