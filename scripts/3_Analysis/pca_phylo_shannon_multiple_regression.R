PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-10-09"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/pca_phylo_shannon_multiple_regression")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20261009_pca_phylo_shannon_multiple_regression.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20261009_pca_phylo_shannon_multiple_regression.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20261009_pca_phylo_shannon_multiple_regression_validation.md")

SCALES <- c("10m", "20m", "50m")
SPECTRAL_METRICS <- c(
  "spec_alpha",
  "spec_pca_mean",
  "spec_rao_q",
  "spec_spca_alpha",
  "spec_spca_mean",
  "spec_spca_rao"
)
PHYLO_PREDICTORS <- c("phy_faith", "phy_rao", "phy_afaith")
BIODIVERSITY_PREDICTORS <- c(PHYLO_PREDICTORS, "sp_shannon")

DISPLAY_NAMES <- c(
  spec_alpha = "Original PCA alpha hull",
  spec_pca_mean = "Original PCA mean distance",
  spec_rao_q = "Original PCA Rao's Q",
  spec_spca_alpha = "Vector-normalized PCA alpha hull",
  spec_spca_mean = "Vector-normalized PCA mean distance",
  spec_spca_rao = "Vector-normalized PCA Rao's Q",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_shannon = "Shannon diversity"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

display_name <- function(metric) {
  ifelse(metric %in% names(DISPLAY_NAMES), unname(DISPLAY_NAMES[metric]), metric)
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

read_scale_table <- function(scale_name) {
  file_path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale_name, ".csv"))
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- as.character(data$scale)
  data
}

read_analysis_data <- function() {
  data <- do.call(rbind, lapply(SCALES, read_scale_table))
  required <- c("quad_id", "scale", SPECTRAL_METRICS, BIODIVERSITY_PREDICTORS)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data
}

complete_model_data <- function(data, scale_name, response, predictors) {
  scale_data <- data[data$scale == scale_name, c("quad_id", "scale", response, predictors), drop = FALSE]
  for (column in c(response, predictors)) {
    scale_data[[column]] <- as.numeric(scale_data[[column]])
  }
  complete <- stats::complete.cases(scale_data[, c(response, predictors), drop = FALSE])
  scale_data[complete, , drop = FALSE]
}

fit_model <- function(data, response, predictors) {
  if (nrow(data) <= length(predictors) + 2) {
    return(NULL)
  }
  if (stats::sd(data[[response]]) == 0) {
    return(NULL)
  }
  variable_sds <- vapply(data[, predictors, drop = FALSE], stats::sd, numeric(1))
  if (any(variable_sds == 0)) {
    return(NULL)
  }
  stats::lm(stats::as.formula(paste(response, "~", paste(predictors, collapse = " + "))), data = data)
}

model_summary_row <- function(scale_name, response, model_name, predictors, fit, n_total) {
  row <- data.frame(
    scale = scale_name,
    spectral_metric = response,
    spectral_label = display_name(response),
    model_name = model_name,
    predictors = paste(predictors, collapse = " + "),
    n = ifelse(is.null(fit), 0L, stats::nobs(fit)),
    n_total_scale_rows = n_total,
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    aic = NA_real_,
    bic = NA_real_,
    rmse = NA_real_,
    residual_abs_cor_fitted = NA_real_,
    residual_shapiro_p = NA_real_,
    max_cooks_distance = NA_real_,
    stringsAsFactors = FALSE
  )
  if (is.null(fit)) {
    return(row)
  }

  fit_summary <- summary(fit)
  f_values <- fit_summary$fstatistic
  residuals <- stats::residuals(fit)
  fitted_values <- stats::fitted(fit)
  row$r_squared <- fit_summary$r.squared
  row$adj_r_squared <- fit_summary$adj.r.squared
  row$f_statistic <- unname(f_values[["value"]])
  row$f_p_value <- stats::pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
  row$aic <- stats::AIC(fit)
  row$bic <- stats::BIC(fit)
  row$rmse <- sqrt(mean(residuals^2))
  if (length(residuals) >= 3 && stats::sd(abs(residuals)) > 0 && stats::sd(fitted_values) > 0) {
    row$residual_abs_cor_fitted <- unname(stats::cor(abs(residuals), fitted_values))
  }
  if (length(residuals) >= 3 && length(residuals) <= 5000) {
    row$residual_shapiro_p <- tryCatch(stats::shapiro.test(residuals)$p.value, error = function(e) NA_real_)
  }
  cooks <- stats::cooks.distance(fit)
  row$max_cooks_distance <- max(cooks, na.rm = TRUE)
  row
}

standardized_fit <- function(data, response, predictors) {
  scaled <- data[, c(response, predictors), drop = FALSE]
  scaled[] <- lapply(scaled, function(x) as.numeric(scale(x)))
  fit_model(scaled, response, predictors)
}

coefficient_rows <- function(scale_name, response, fit, std_fit) {
  if (is.null(fit)) {
    return(data.frame())
  }
  coef_table <- as.data.frame(summary(fit)$coefficients)
  names(coef_table) <- c("estimate", "std_error", "t_value", "p_value")
  coef_table$term <- row.names(coef_table)
  row.names(coef_table) <- NULL

  std_coef <- rep(NA_real_, nrow(coef_table))
  if (!is.null(std_fit)) {
    std_values <- stats::coef(std_fit)
    std_coef <- unname(std_values[coef_table$term])
  }

  data.frame(
    scale = scale_name,
    spectral_metric = response,
    spectral_label = display_name(response),
    term = coef_table$term,
    term_label = ifelse(coef_table$term %in% names(DISPLAY_NAMES), display_name(coef_table$term), coef_table$term),
    estimate = coef_table$estimate,
    std_error = coef_table$std_error,
    t_value = coef_table$t_value,
    p_value = coef_table$p_value,
    standardized_estimate = std_coef,
    stringsAsFactors = FALSE
  )
}

safe_r2 <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }
  summary(fit)$r.squared
}

predictor_increment_rows <- function(scale_name, response, data, full_fit, phylo_fit, shannon_fit) {
  full_r2 <- safe_r2(full_fit)
  rows <- data.frame(
    scale = scale_name,
    spectral_metric = response,
    spectral_label = display_name(response),
    increment_type = c("phylogenetic_block_after_shannon", "shannon_after_phylogenetic_block"),
    added_predictor = c("phy_faith + phy_rao + phy_afaith", "sp_shannon"),
    added_label = c("Phylogenetic block", "Shannon diversity"),
    base_r_squared = c(safe_r2(shannon_fit), safe_r2(phylo_fit)),
    full_r_squared = full_r2,
    incremental_r_squared = c(full_r2 - safe_r2(shannon_fit), full_r2 - safe_r2(phylo_fit)),
    partial_f = NA_real_,
    partial_p_value = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(full_fit) && !is.null(shannon_fit)) {
    test <- stats::anova(shannon_fit, full_fit)
    rows$partial_f[1] <- test$F[2]
    rows$partial_p_value[1] <- test$`Pr(>F)`[2]
  }
  if (!is.null(full_fit) && !is.null(phylo_fit)) {
    test <- stats::anova(phylo_fit, full_fit)
    rows$partial_f[2] <- test$F[2]
    rows$partial_p_value[2] <- test$`Pr(>F)`[2]
  }

  individual_rows <- lapply(BIODIVERSITY_PREDICTORS, function(predictor) {
    reduced_predictors <- setdiff(BIODIVERSITY_PREDICTORS, predictor)
    reduced_fit <- fit_model(data, response, reduced_predictors)
    reduced_r2 <- safe_r2(reduced_fit)
    partial_f <- NA_real_
    partial_p <- NA_real_
    if (!is.null(full_fit) && !is.null(reduced_fit)) {
      test <- stats::anova(reduced_fit, full_fit)
      partial_f <- test$F[2]
      partial_p <- test$`Pr(>F)`[2]
    }
    data.frame(
      scale = scale_name,
      spectral_metric = response,
      spectral_label = display_name(response),
      increment_type = "single_predictor_after_other_three",
      added_predictor = predictor,
      added_label = display_name(predictor),
      base_r_squared = reduced_r2,
      full_r_squared = full_r2,
      incremental_r_squared = full_r2 - reduced_r2,
      partial_f = partial_f,
      partial_p_value = partial_p,
      stringsAsFactors = FALSE
    )
  })

  rbind(rows, do.call(rbind, individual_rows))
}

vif_rows <- function(scale_name, response, data) {
  rows <- lapply(BIODIVERSITY_PREDICTORS, function(predictor) {
    other_predictors <- setdiff(BIODIVERSITY_PREDICTORS, predictor)
    fit <- fit_model(data, predictor, other_predictors)
    r2 <- safe_r2(fit)
    data.frame(
      scale = scale_name,
      spectral_metric = response,
      spectral_label = display_name(response),
      predictor = predictor,
      predictor_label = display_name(predictor),
      predictor_r_squared_with_other_predictors = r2,
      vif = ifelse(is.na(r2), NA_real_, 1 / (1 - r2)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

single_predictor_rows <- function(scale_name, response, data) {
  rows <- lapply(BIODIVERSITY_PREDICTORS, function(predictor) {
    fit <- fit_model(data, response, predictor)
    summary <- model_summary_row(scale_name, response, paste0("single_", predictor), predictor, fit, nrow(data))
    summary$predictor <- predictor
    summary$predictor_label <- display_name(predictor)
    summary
  })
  do.call(rbind, rows)
}

run_models <- function(data) {
  model_rows <- list()
  coef_rows <- list()
  increment_rows <- list()
  vif_output_rows <- list()
  single_rows <- list()
  row_index <- 1

  for (scale_name in SCALES) {
    for (response in SPECTRAL_METRICS) {
      model_data <- complete_model_data(data, scale_name, response, BIODIVERSITY_PREDICTORS)
      n_total <- nrow(data[data$scale == scale_name, , drop = FALSE])

      full_fit <- fit_model(model_data, response, BIODIVERSITY_PREDICTORS)
      phylo_fit <- fit_model(model_data, response, PHYLO_PREDICTORS)
      shannon_fit <- fit_model(model_data, response, "sp_shannon")
      std_fit <- standardized_fit(model_data, response, BIODIVERSITY_PREDICTORS)

      model_rows[[row_index]] <- model_summary_row(
        scale_name, response, "full_phylogenetic_plus_shannon", BIODIVERSITY_PREDICTORS, full_fit, n_total
      )
      coef_rows[[row_index]] <- coefficient_rows(scale_name, response, full_fit, std_fit)
      increment_rows[[row_index]] <- predictor_increment_rows(scale_name, response, model_data, full_fit, phylo_fit, shannon_fit)
      vif_output_rows[[row_index]] <- vif_rows(scale_name, response, model_data)
      single_rows[[row_index]] <- single_predictor_rows(scale_name, response, model_data)
      row_index <- row_index + 1
    }
  }

  list(
    model_summary = do.call(rbind, model_rows),
    coefficient_summary = do.call(rbind, coef_rows),
    incremental_r2 = do.call(rbind, increment_rows),
    predictor_vif = do.call(rbind, vif_output_rows),
    single_predictor_summary = do.call(rbind, single_rows)
  )
}

make_top_full_models <- function(model_summary) {
  full <- model_summary[model_summary$model_name == "full_phylogenetic_plus_shannon", , drop = FALSE]
  full <- full[order(-full$r_squared), , drop = FALSE]
  full
}

make_strongest_partial_predictors <- function(incremental_r2) {
  single <- incremental_r2[incremental_r2$increment_type == "single_predictor_after_other_three", , drop = FALSE]
  single <- single[order(-single$incremental_r_squared), , drop = FALSE]
  single
}

write_outputs <- function(outputs) {
  model_file <- file.path(TABLE_DIR, "pca_phylo_shannon_full_model_summary.csv")
  coef_file <- file.path(TABLE_DIR, "pca_phylo_shannon_full_model_coefficients.csv")
  increment_file <- file.path(TABLE_DIR, "pca_phylo_shannon_incremental_r2.csv")
  vif_file <- file.path(TABLE_DIR, "pca_phylo_shannon_predictor_vif.csv")
  single_file <- file.path(TABLE_DIR, "pca_phylo_shannon_single_predictor_models.csv")
  top_file <- file.path(TABLE_DIR, "pca_phylo_shannon_top_full_models.csv")
  partial_file <- file.path(TABLE_DIR, "pca_phylo_shannon_strongest_partial_predictors.csv")

  top_models <- make_top_full_models(outputs$model_summary)
  partial_predictors <- make_strongest_partial_predictors(outputs$incremental_r2)

  write.csv(outputs$model_summary, model_file, row.names = FALSE)
  write.csv(outputs$coefficient_summary, coef_file, row.names = FALSE)
  write.csv(outputs$incremental_r2, increment_file, row.names = FALSE)
  write.csv(outputs$predictor_vif, vif_file, row.names = FALSE)
  write.csv(outputs$single_predictor_summary, single_file, row.names = FALSE)
  write.csv(top_models, top_file, row.names = FALSE)
  write.csv(partial_predictors, partial_file, row.names = FALSE)

  list(
    model_file = model_file,
    coef_file = coef_file,
    increment_file = increment_file,
    vif_file = vif_file,
    single_file = single_file,
    top_file = top_file,
    partial_file = partial_file,
    top_models = top_models,
    partial_predictors = partial_predictors
  )
}

write_analysis_report <- function(outputs, files) {
  top_table <- head(files$top_models, 12)
  top_table <- data.frame(
    Scale = top_table$scale,
    Spectral = top_table$spectral_label,
    n = top_table$n,
    R2 = fmt_num(top_table$r_squared, 3),
    Adj_R2 = fmt_num(top_table$adj_r_squared, 3),
    p = fmt_p(top_table$f_p_value),
    stringsAsFactors = FALSE
  )

  partial_table <- head(files$partial_predictors, 12)
  partial_table <- data.frame(
    Scale = partial_table$scale,
    Spectral = partial_table$spectral_label,
    Added = partial_table$added_label,
    Base_R2 = fmt_num(partial_table$base_r_squared, 3),
    Full_R2 = fmt_num(partial_table$full_r_squared, 3),
    Incremental_R2 = fmt_num(partial_table$incremental_r_squared, 3),
    p = fmt_p(partial_table$partial_p_value),
    stringsAsFactors = FALSE
  )

  block <- outputs$incremental_r2[outputs$incremental_r2$increment_type %in% c(
    "phylogenetic_block_after_shannon",
    "shannon_after_phylogenetic_block"
  ), , drop = FALSE]
  block <- block[order(block$scale, block$spectral_metric, block$increment_type), , drop = FALSE]
  block_table <- data.frame(
    Scale = block$scale,
    Spectral = block$spectral_label,
    Added = block$added_label,
    Base_R2 = fmt_num(block$base_r_squared, 3),
    Full_R2 = fmt_num(block$full_r_squared, 3),
    Incremental_R2 = fmt_num(block$incremental_r_squared, 3),
    p = fmt_p(block$partial_p_value),
    stringsAsFactors = FALSE
  )

  lines <- c(
    "# PCA Phylogenetic-Shannon Multiple Regression",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis fits scale-specific multiple regressions for the six PCA-based spectral heterogeneity metrics. Each model uses Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity as biodiversity predictors.",
    "",
    "## Model Form",
    "",
    "`spectral_metric ~ phy_faith + phy_rao + phy_afaith + sp_shannon`",
    "",
    "Models are fit separately for each quadrat scale and each spectral metric. The analysis uses complete cases from the current combined quadrat tables and does not add edge, environmental, brightness, or spatial covariates.",
    "",
    "## Top Full Models",
    "",
    markdown_table(top_table),
    "",
    "## Strongest Single-Predictor Partial Contributions",
    "",
    "These rows show the additional R2 from adding one biodiversity predictor after the other three predictors are already in the model.",
    "",
    markdown_table(partial_table),
    "",
    "## Biodiversity Block Contributions",
    "",
    "The phylogenetic block is Faith's PD, phylogenetic Rao's Q, and abundance-weighted Faith's PD fitted after Shannon diversity. The Shannon row is Shannon diversity fitted after the three phylogenetic predictors.",
    "",
    markdown_table(block_table),
    "",
    "## Output Tables",
    "",
    paste0("- `", relative_path(files$model_file), "`"),
    paste0("- `", relative_path(files$coef_file), "`"),
    paste0("- `", relative_path(files$increment_file), "`"),
    paste0("- `", relative_path(files$vif_file), "`"),
    paste0("- `", relative_path(files$single_file), "`"),
    paste0("- `", relative_path(files$top_file), "`"),
    paste0("- `", relative_path(files$partial_file), "`"),
    "",
    "## Interpretation Notes",
    "",
    "- The vector-normalized PCA metrics remain the primary PCA evidence layer; original PCA metrics are retained as brightness-sensitive comparisons.",
    "- These are ordinary least-squares models and should be interpreted cautiously because prior diagnostics found spatial autocorrelation in several spectral metrics, biodiversity metrics, and residuals.",
    "- Predictor collinearity is expected because the biodiversity metrics are not independent. Use the standardized coefficients, incremental R2 rows, and VIF table together rather than interpreting raw coefficients alone.",
    "- The 50 m models have the strongest R2 values but smaller complete-case sample sizes."
  )
  writeLines(lines, ANALYSIS_REPORT)
}

write_task_report <- function(files) {
  lines <- c(
    "# PCA Phylogenetic-Shannon Multiple Regression Task Report",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Completed",
    "",
    "- Added a focused multiple-regression analysis for six PCA-based spectral heterogeneity metrics.",
    "- Modeled each spectral metric as a response to the three phylogenetic diversity metrics plus Shannon diversity.",
    "- Wrote full-model summaries, coefficients, incremental R2 diagnostics, VIF diagnostics, single-predictor model summaries, and ranked top-model tables.",
    "- Wrote analysis and validation Markdown notes.",
    "",
    "## Outputs",
    "",
    paste0("- `", relative_path(files$model_file), "`"),
    paste0("- `", relative_path(files$coef_file), "`"),
    paste0("- `", relative_path(files$increment_file), "`"),
    paste0("- `", relative_path(files$vif_file), "`"),
    paste0("- `", relative_path(files$single_file), "`"),
    paste0("- `", relative_path(files$top_file), "`"),
    paste0("- `", relative_path(files$partial_file), "`"),
    paste0("- `", relative_path(ANALYSIS_REPORT), "`"),
    paste0("- `", relative_path(VALIDATION_REPORT), "`")
  )
  writeLines(lines, TASK_REPORT)
}

write_validation_report <- function(outputs, data) {
  expected_models <- length(SCALES) * length(SPECTRAL_METRICS)
  expected_coefficients <- expected_models * (length(BIODIVERSITY_PREDICTORS) + 1)
  expected_increment_rows <- expected_models * (2 + length(BIODIVERSITY_PREDICTORS))
  expected_vif_rows <- expected_models * length(BIODIVERSITY_PREDICTORS)

  checks <- data.frame(
    check = c(
      "required_columns_present",
      "full_model_row_count",
      "coefficient_row_count",
      "incremental_r2_row_count",
      "vif_row_count",
      "all_full_models_estimable",
      "finite_full_model_r2",
      "complete_case_n_positive"
    ),
    expected = c(
      "all required fields",
      expected_models,
      expected_coefficients,
      expected_increment_rows,
      expected_vif_rows,
      expected_models,
      expected_models,
      expected_models
    ),
    observed = c(
      "all required fields",
      nrow(outputs$model_summary),
      nrow(outputs$coefficient_summary),
      nrow(outputs$incremental_r2),
      nrow(outputs$predictor_vif),
      sum(outputs$model_summary$n > length(BIODIVERSITY_PREDICTORS) + 2),
      sum(is.finite(outputs$model_summary$r_squared)),
      sum(outputs$model_summary$n > 0)
    ),
    passed = c(
      TRUE,
      nrow(outputs$model_summary) == expected_models,
      nrow(outputs$coefficient_summary) == expected_coefficients,
      nrow(outputs$incremental_r2) == expected_increment_rows,
      nrow(outputs$predictor_vif) == expected_vif_rows,
      sum(outputs$model_summary$n > length(BIODIVERSITY_PREDICTORS) + 2) == expected_models,
      sum(is.finite(outputs$model_summary$r_squared)) == expected_models,
      sum(outputs$model_summary$n > 0) == expected_models
    ),
    stringsAsFactors = FALSE
  )

  validation_file <- file.path(TABLE_DIR, "pca_phylo_shannon_validation_checks.csv")
  write.csv(checks, validation_file, row.names = FALSE)

  complete_n <- outputs$model_summary[, c("scale", "spectral_metric", "spectral_label", "n", "n_total_scale_rows")]
  complete_file <- file.path(TABLE_DIR, "pca_phylo_shannon_complete_case_counts.csv")
  write.csv(complete_n, complete_file, row.names = FALSE)

  lines <- c(
    "# PCA Phylogenetic-Shannon Multiple Regression Validation",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Validation Checks",
    "",
    markdown_table(checks),
    "",
    "## Complete-Case Counts",
    "",
    markdown_table(data.frame(
      Scale = complete_n$scale,
      Spectral = complete_n$spectral_label,
      n = complete_n$n,
      Total_rows = complete_n$n_total_scale_rows,
      stringsAsFactors = FALSE
    )),
    "",
    "## Validation Outputs",
    "",
    paste0("- `", relative_path(validation_file), "`"),
    paste0("- `", relative_path(complete_file), "`"),
    "",
    "## Notes",
    "",
    "- All models use complete cases for the response and the four biodiversity predictors.",
    "- The combined quadrat tables already exclude unavailable PCA values with `NA`, so manually excluded atmospheric/cloud quadrats are absent from complete-case fits for PCA-dependent metrics."
  )
  writeLines(lines, VALIDATION_REPORT)

  if (!all(checks$passed)) {
    stop("Validation failed; see ", VALIDATION_REPORT, call. = FALSE)
  }
}

data <- read_analysis_data()
outputs <- run_models(data)
files <- write_outputs(outputs)
write_analysis_report(outputs, files)
write_task_report(files)
write_validation_report(outputs, data)

message("Wrote PCA phylogenetic-Shannon multiple regression outputs to: ", TABLE_DIR)
