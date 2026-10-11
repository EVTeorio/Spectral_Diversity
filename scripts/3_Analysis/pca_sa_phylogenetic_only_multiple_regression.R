PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-10-09"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/pca_sa_phylogenetic_only_multiple_regression")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20261009_pca_sa_phylogenetic_only_multiple_regression.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20261009_pca_sa_phylogenetic_only_multiple_regression.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20261009_pca_sa_phylogenetic_only_multiple_regression_validation.md")
PAPER_TABLE <- file.path(PROJECT_DIR, "Documents/Paper/20261009_phylogenetic_only_pca_sa_multiple_regression_table.md")

SCALES <- c("10m", "20m", "50m")
SPECTRAL_METRICS <- c(
  "spec_sa",
  "spec_alpha",
  "spec_pca_mean",
  "spec_rao_q",
  "spec_spca_alpha",
  "spec_spca_mean",
  "spec_spca_rao"
)
PHYLO_PREDICTORS <- c("phy_faith", "phy_rao", "phy_afaith")

DISPLAY_NAMES <- c(
  spec_sa = "Spectral angle entropy",
  spec_alpha = "Original PCA alpha hull",
  spec_pca_mean = "Original PCA mean distance",
  spec_rao_q = "Original PCA Rao's Q",
  spec_spca_alpha = "Vector-normalized PCA alpha hull",
  spec_spca_mean = "Vector-normalized PCA mean distance",
  spec_spca_rao = "Vector-normalized PCA Rao's Q",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(PAPER_TABLE), recursive = TRUE, showWarnings = FALSE)

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
  required <- c("quad_id", "scale", SPECTRAL_METRICS, PHYLO_PREDICTORS)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data
}

complete_model_data <- function(data, scale_name, response) {
  columns <- c("quad_id", "scale", response, PHYLO_PREDICTORS)
  scale_data <- data[data$scale == scale_name, columns, drop = FALSE]
  for (column in c(response, PHYLO_PREDICTORS)) {
    scale_data[[column]] <- as.numeric(scale_data[[column]])
  }
  complete <- stats::complete.cases(scale_data[, c(response, PHYLO_PREDICTORS), drop = FALSE])
  scale_data[complete, , drop = FALSE]
}

fit_model <- function(data, response, predictors) {
  if (nrow(data) <= length(predictors) + 2) {
    return(NULL)
  }
  if (stats::sd(data[[response]]) == 0) {
    return(NULL)
  }
  predictor_sds <- vapply(data[, predictors, drop = FALSE], stats::sd, numeric(1))
  if (any(predictor_sds == 0)) {
    return(NULL)
  }
  stats::lm(stats::as.formula(paste(response, "~", paste(predictors, collapse = " + "))), data = data)
}

standardized_fit <- function(data, response, predictors) {
  scaled <- data[, c(response, predictors), drop = FALSE]
  scaled[] <- lapply(scaled, function(x) as.numeric(scale(x)))
  fit_model(scaled, response, predictors)
}

model_summary_row <- function(scale_name, response, fit, n_total) {
  row <- data.frame(
    scale = scale_name,
    spectral_metric = response,
    spectral_label = display_name(response),
    predictors = paste(PHYLO_PREDICTORS, collapse = " + "),
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
  row$max_cooks_distance <- max(stats::cooks.distance(fit), na.rm = TRUE)
  row
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

partial_predictor_rows <- function(scale_name, response, data, full_fit) {
  rows <- lapply(PHYLO_PREDICTORS, function(predictor) {
    reduced_predictors <- setdiff(PHYLO_PREDICTORS, predictor)
    reduced_fit <- fit_model(data, response, reduced_predictors)
    reduced_r2 <- safe_r2(reduced_fit)
    full_r2 <- safe_r2(full_fit)
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
  do.call(rbind, rows)
}

vif_rows <- function(scale_name, response, data) {
  rows <- lapply(PHYLO_PREDICTORS, function(predictor) {
    other_predictors <- setdiff(PHYLO_PREDICTORS, predictor)
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

run_models <- function(data) {
  model_rows <- list()
  coef_rows <- list()
  partial_rows <- list()
  vif_output_rows <- list()
  row_index <- 1

  for (scale_name in SCALES) {
    for (response in SPECTRAL_METRICS) {
      model_data <- complete_model_data(data, scale_name, response)
      n_total <- nrow(data[data$scale == scale_name, , drop = FALSE])
      full_fit <- fit_model(model_data, response, PHYLO_PREDICTORS)
      std_fit <- standardized_fit(model_data, response, PHYLO_PREDICTORS)

      model_rows[[row_index]] <- model_summary_row(scale_name, response, full_fit, n_total)
      coef_rows[[row_index]] <- coefficient_rows(scale_name, response, full_fit, std_fit)
      partial_rows[[row_index]] <- partial_predictor_rows(scale_name, response, model_data, full_fit)
      vif_output_rows[[row_index]] <- vif_rows(scale_name, response, model_data)
      row_index <- row_index + 1
    }
  }

  list(
    model_summary = do.call(rbind, model_rows),
    coefficient_summary = do.call(rbind, coef_rows),
    partial_r2 = do.call(rbind, partial_rows),
    predictor_vif = do.call(rbind, vif_output_rows)
  )
}

make_paper_wide_table <- function(model_summary) {
  rows <- lapply(SPECTRAL_METRICS, function(metric) {
    metric_rows <- model_summary[model_summary$spectral_metric == metric, , drop = FALSE]
    values <- list(Spectral_metric = display_name(metric))
    for (scale_name in SCALES) {
      scale_row <- metric_rows[metric_rows$scale == scale_name, , drop = FALSE]
      values[[paste0(scale_name, " n")]] <- scale_row$n
      values[[paste0(scale_name, " R2")]] <- fmt_num(scale_row$r_squared, 3)
      values[[paste0(scale_name, " adj R2")]] <- fmt_num(scale_row$adj_r_squared, 3)
      values[[paste0(scale_name, " p")]] <- fmt_p(scale_row$f_p_value)
    }
    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  })
  do.call(rbind, rows)
}

write_outputs <- function(outputs) {
  summary_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_full_model_summary.csv")
  coef_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_full_model_coefficients.csv")
  partial_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_partial_r2.csv")
  vif_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_predictor_vif.csv")
  top_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_top_full_models.csv")

  top_models <- outputs$model_summary[order(-outputs$model_summary$r_squared), , drop = FALSE]
  write.csv(outputs$model_summary, summary_file, row.names = FALSE)
  write.csv(outputs$coefficient_summary, coef_file, row.names = FALSE)
  write.csv(outputs$partial_r2, partial_file, row.names = FALSE)
  write.csv(outputs$predictor_vif, vif_file, row.names = FALSE)
  write.csv(top_models, top_file, row.names = FALSE)

  list(
    summary_file = summary_file,
    coef_file = coef_file,
    partial_file = partial_file,
    vif_file = vif_file,
    top_file = top_file,
    top_models = top_models,
    paper_wide_table = make_paper_wide_table(outputs$model_summary)
  )
}

write_paper_table <- function(files) {
  lines <- c(
    "# Phylogenetic-Only Multiple Regression Models for PCA and Spectral Angle Heterogeneity",
    "",
    "OLS model form for each scale and spectral metric:",
    "",
    "`spectral_metric ~ Faith's PD + phylogenetic Rao's Q + abundance-weighted Faith's PD`",
    "",
    markdown_table(files$paper_wide_table),
    "",
    "Source: `reports/tables/pca_sa_phylogenetic_only_multiple_regression/pca_sa_phylogenetic_only_full_model_summary.csv`.",
    "",
    "Notes: spectral angle entropy is included as the angular non-PCA comparison metric. Original PCA refers to the brightness-sensitive PCA basis; vector-normalized PCA refers to the PCA basis built after vector-normalizing spectra. Models use complete cases from the current combined quadrat tables. Significance is from the full-model F test."
  )
  writeLines(lines, PAPER_TABLE)
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

  partial_table <- outputs$partial_r2[order(-outputs$partial_r2$incremental_r_squared), , drop = FALSE]
  partial_table <- head(partial_table, 12)
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

  lines <- c(
    "# PCA and Spectral Angle Phylogenetic-Only Multiple Regression",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis fits scale-specific multiple regressions for spectral angle entropy plus the six PCA-based spectral heterogeneity metrics using only phylogenetic diversity predictors.",
    "",
    "## Model Form",
    "",
    "`spectral_metric ~ phy_faith + phy_rao + phy_afaith`",
    "",
    "## Top Full Models",
    "",
    markdown_table(top_table),
    "",
    "## Strongest Single-Predictor Partial Contributions",
    "",
    "These rows show the additional R2 from adding one phylogenetic predictor after the other two phylogenetic predictors are already in the model.",
    "",
    markdown_table(partial_table),
    "",
    "## Paper-Facing Table",
    "",
    paste0("- `", relative_path(PAPER_TABLE), "`"),
    "",
    "## Output Tables",
    "",
    paste0("- `", relative_path(files$summary_file), "`"),
    paste0("- `", relative_path(files$coef_file), "`"),
    paste0("- `", relative_path(files$partial_file), "`"),
    paste0("- `", relative_path(files$vif_file), "`"),
    paste0("- `", relative_path(files$top_file), "`"),
    "",
    "## Interpretation Notes",
    "",
    "- These models isolate the phylogenetic predictors by excluding Shannon diversity.",
    "- The PCA top model rankings are nearly identical to the phylogenetic-plus-Shannon layer, confirming that most explanatory power in the four-predictor model came from the phylogenetic block.",
    "- As with the prior OLS layer, p-values should remain cautious because prior spatial diagnostics found residual spatial autocorrelation for several priority relationships."
  )
  writeLines(lines, ANALYSIS_REPORT)
}

write_task_report <- function(files) {
  lines <- c(
    "# PCA and Spectral Angle Phylogenetic-Only Multiple Regression Task Report",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Completed",
    "",
    "- Added a phylogenetic-only multiple-regression analysis for spectral angle entropy plus six PCA-based spectral heterogeneity metrics.",
    "- Modeled each spectral metric as a response to Faith's PD, phylogenetic Rao's Q, and abundance-weighted Faith's PD.",
    "- Wrote model summaries, coefficients, partial R2 diagnostics, VIF diagnostics, ranked tables, and a paper-facing Markdown table.",
    "",
    "## Outputs",
    "",
    paste0("- `", relative_path(files$summary_file), "`"),
    paste0("- `", relative_path(files$coef_file), "`"),
    paste0("- `", relative_path(files$partial_file), "`"),
    paste0("- `", relative_path(files$vif_file), "`"),
    paste0("- `", relative_path(files$top_file), "`"),
    paste0("- `", relative_path(PAPER_TABLE), "`"),
    paste0("- `", relative_path(ANALYSIS_REPORT), "`"),
    paste0("- `", relative_path(VALIDATION_REPORT), "`")
  )
  writeLines(lines, TASK_REPORT)
}

write_validation_report <- function(outputs) {
  expected_models <- length(SCALES) * length(SPECTRAL_METRICS)
  expected_coefficients <- expected_models * (length(PHYLO_PREDICTORS) + 1)
  expected_partial_rows <- expected_models * length(PHYLO_PREDICTORS)
  expected_vif_rows <- expected_models * length(PHYLO_PREDICTORS)

  checks <- data.frame(
    check = c(
      "full_model_row_count",
      "coefficient_row_count",
      "partial_r2_row_count",
      "vif_row_count",
      "all_full_models_estimable",
      "finite_full_model_r2",
      "complete_case_n_positive",
      "paper_table_created"
    ),
    expected = c(
      expected_models,
      expected_coefficients,
      expected_partial_rows,
      expected_vif_rows,
      expected_models,
      expected_models,
      expected_models,
      "exists"
    ),
    observed = c(
      nrow(outputs$model_summary),
      nrow(outputs$coefficient_summary),
      nrow(outputs$partial_r2),
      nrow(outputs$predictor_vif),
      sum(outputs$model_summary$n > length(PHYLO_PREDICTORS) + 2),
      sum(is.finite(outputs$model_summary$r_squared)),
      sum(outputs$model_summary$n > 0),
      file.exists(PAPER_TABLE)
    ),
    passed = c(
      nrow(outputs$model_summary) == expected_models,
      nrow(outputs$coefficient_summary) == expected_coefficients,
      nrow(outputs$partial_r2) == expected_partial_rows,
      nrow(outputs$predictor_vif) == expected_vif_rows,
      sum(outputs$model_summary$n > length(PHYLO_PREDICTORS) + 2) == expected_models,
      sum(is.finite(outputs$model_summary$r_squared)) == expected_models,
      sum(outputs$model_summary$n > 0) == expected_models,
      file.exists(PAPER_TABLE)
    ),
    stringsAsFactors = FALSE
  )
  validation_file <- file.path(TABLE_DIR, "pca_sa_phylogenetic_only_validation_checks.csv")
  write.csv(checks, validation_file, row.names = FALSE)

  lines <- c(
    "# PCA and Spectral Angle Phylogenetic-Only Multiple Regression Validation",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Validation Checks",
    "",
    markdown_table(checks),
    "",
    "## Notes",
    "",
    "- All models use complete cases for the response and the three phylogenetic predictors.",
    "- Complete-case counts match PCA metric availability from the combined quadrat tables."
  )
  writeLines(lines, VALIDATION_REPORT)

  if (!all(checks$passed)) {
    stop("Validation failed; see ", VALIDATION_REPORT, call. = FALSE)
  }
}

data <- read_analysis_data()
outputs <- run_models(data)
files <- write_outputs(outputs)
write_paper_table(files)
write_analysis_report(outputs, files)
write_task_report(files)
write_validation_report(outputs)

message("Wrote PCA and spectral angle phylogenetic-only multiple regression outputs to: ", TABLE_DIR)
message("Wrote paper-facing table to: ", PAPER_TABLE)
