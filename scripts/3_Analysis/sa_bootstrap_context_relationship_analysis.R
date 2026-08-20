PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-17"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/bootstrap_metric_relationships")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/bootstrap_metric_relationships")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260817_sa_bootstrap_context_relationship_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260817_sa_bootstrap_context_relationship_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260817_sa_bootstrap_context_relationship_analysis_validation.md")

SCALES <- c("10m", "20m", "50m")

BOOT_STATS <- c("boot_sd", "boot_cv", "boot_se")
BOOT_STAT_LABELS <- c(
  boot_sd = "SA bootstrap SD",
  boot_cv = "SA bootstrap CV",
  boot_se = "SA bootstrap SE"
)

CONTINUOUS_PREDICTORS <- data.frame(
  predictor = c(
    "env_elev",
    "present_species_count",
    "crown_equivalent_individuals",
    "mean_pixel_brightness",
    "mean_blue_pixel_brightness",
    "mean_green_pixel_brightness",
    "mean_red_pixel_brightness",
    "mean_nir_pixel_brightness"
  ),
  predictor_label = c(
    "Elevation",
    "Species count",
    "Individual count",
    "Overall brightness",
    "Blue brightness",
    "Green brightness",
    "Red brightness",
    "NIR brightness"
  ),
  color = c(
    "#5A5A5A",
    "#2E7D32",
    "#7B4D2B",
    "#4B4B4B",
    "#2166AC",
    "#1B9E77",
    "#D73027",
    "#542788"
  ),
  stringsAsFactors = FALSE
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

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

read_scale_analysis_table <- function(scale) {
  path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale, ".csv"))
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data[, c("scale", "quad_id", "env_elev"), drop = FALSE]
}

read_context_dataset <- function() {
  final_sa <- read.csv(
    file.path(TABLE_DIR, "final_sa_entropy_bootstrap_metric_stats.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  final_sa <- final_sa[, c("scale", "quad_id", "metric_value", "metric_sd", "metric_cv", "metric_se"), drop = FALSE]
  names(final_sa) <- c("scale", "quad_id", "sa_entropy", "boot_sd", "boot_cv", "boot_se")

  elevation <- do.call(rbind, lapply(SCALES, read_scale_analysis_table))

  counts <- read.csv(
    file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  clusters <- read.csv(
    file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  clusters <- clusters[, c("scale", "quad_id", "composition_cluster", "composition_type"), drop = FALSE]

  brightness <- read.csv(
    file.path(PROJECT_DIR, "reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  data <- merge(final_sa, elevation, by = c("scale", "quad_id"), all.x = TRUE)
  data <- merge(data, counts, by = c("scale", "quad_id"), all.x = TRUE)
  data <- merge(data, clusters, by = c("scale", "quad_id"), all.x = TRUE)
  data <- merge(data, brightness, by = c("scale", "quad_id"), all.x = TRUE)
  data
}

continuous_relationship_row <- function(data, predictor, boot_stat) {
  complete <- !is.na(data[[predictor]]) & !is.na(data[[boot_stat]])
  row <- data.frame(
    n = sum(complete),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    stringsAsFactors = FALSE
  )
  if (sum(complete) >= 3) {
    model_data <- data[complete, c(predictor, boot_stat)]
    names(model_data) <- c("x", "y")
    fit <- lm(y ~ x, data = model_data)
    fit_summary <- summary(fit)
    f_values <- fit_summary$fstatistic
    pearson <- cor.test(model_data$x, model_data$y, method = "pearson")
    row$pearson_r <- unname(pearson$estimate)
    row$r_squared <- fit_summary$r.squared
    row$f_statistic <- unname(f_values[["value"]])
    row$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
    row$slope <- unname(coef(fit)[["x"]])
    row$intercept <- unname(coef(fit)[["(Intercept)"]])
  }
  row
}

composition_relationship_row <- function(data, boot_stat) {
  complete <- !is.na(data$composition_type) & nzchar(data$composition_type) & !is.na(data[[boot_stat]])
  row <- data.frame(
    n = sum(complete),
    n_groups = length(unique(data$composition_type[complete])),
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    eta_squared = NA_real_,
    stringsAsFactors = FALSE
  )
  if (sum(complete) >= 3 && row$n_groups >= 2) {
    model_data <- data[complete, c("composition_type", boot_stat)]
    names(model_data) <- c("composition_type", "y")
    fit <- lm(y ~ composition_type, data = model_data)
    anova_fit <- anova(fit)
    row$f_statistic <- anova_fit$`F value`[1]
    row$f_p_value <- anova_fit$`Pr(>F)`[1]
    row$eta_squared <- anova_fit$`Sum Sq`[1] / sum(anova_fit$`Sum Sq`)
  }
  row
}

summarize_relationships <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    for (boot_stat in BOOT_STATS) {
      for (i in seq_len(nrow(CONTINUOUS_PREDICTORS))) {
        predictor <- CONTINUOUS_PREDICTORS$predictor[i]
        row <- continuous_relationship_row(scale_data, predictor, boot_stat)
        rows[[index]] <- cbind(
          data.frame(
            relationship_type = "continuous",
            scale = scale_name,
            boot_stat = boot_stat,
            boot_stat_label = BOOT_STAT_LABELS[[boot_stat]],
            predictor = predictor,
            predictor_label = CONTINUOUS_PREDICTORS$predictor_label[i],
            stringsAsFactors = FALSE
          ),
          row,
          data.frame(n_groups = NA_integer_, eta_squared = NA_real_, stringsAsFactors = FALSE)
        )
        index <- index + 1
      }
      comp <- composition_relationship_row(scale_data, boot_stat)
      rows[[index]] <- data.frame(
        relationship_type = "composition",
        scale = scale_name,
        boot_stat = boot_stat,
        boot_stat_label = BOOT_STAT_LABELS[[boot_stat]],
        predictor = "composition_type",
        predictor_label = "Species composition type",
        n = comp$n,
        pearson_r = NA_real_,
        r_squared = NA_real_,
        f_statistic = comp$f_statistic,
        f_p_value = comp$f_p_value,
        slope = NA_real_,
        intercept = NA_real_,
        n_groups = comp$n_groups,
        eta_squared = comp$eta_squared,
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

plot_continuous_panel <- function(scale_data, predictor, predictor_label, boot_stat, color, relationships) {
  complete <- !is.na(scale_data[[predictor]]) & !is.na(scale_data[[boot_stat]])
  plot(
    scale_data[[predictor]],
    scale_data[[boot_stat]],
    pch = 16,
    cex = 0.45,
    col = adjustcolor(color, alpha.f = 0.58),
    xlab = predictor_label,
    ylab = BOOT_STAT_LABELS[[boot_stat]],
    main = paste(BOOT_STAT_LABELS[[boot_stat]], "vs", predictor_label)
  )
  row <- relationships[
    relationships$relationship_type == "continuous" &
      relationships$scale == unique(scale_data$scale)[1] &
      relationships$boot_stat == boot_stat &
      relationships$predictor == predictor,
  ]
  if (nrow(row) == 1 && !is.na(row$pearson_r) && sum(complete) >= 3) {
    fit <- lm(scale_data[[boot_stat]][complete] ~ scale_data[[predictor]][complete])
    abline(fit, col = "#111111", lwd = 2)
    legend(
      "topleft",
      legend = c(
        paste0("n=", row$n),
        paste0("r=", fmt_num(row$pearson_r, 2)),
        paste0("R2=", fmt_num(row$r_squared, 2)),
        paste0("p=", fmt_p(row$f_p_value))
      ),
      bty = "n",
      cex = 0.70
    )
  }
}

save_continuous_figures <- function(data, relationships) {
  files <- character(length(SCALES))
  for (i in seq_along(SCALES)) {
    scale_name <- SCALES[i]
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    file_path <- file.path(FIGURE_DIR, paste0("03_sa_bootstrap_stats_vs_continuous_context_", scale_name, ".png"))
    png(file_path, width = 6400, height = 2700, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(3, 8), mar = c(4.1, 4.2, 2.3, 0.6), oma = c(0, 0, 2.1, 0))
    for (boot_stat in BOOT_STATS) {
      for (j in seq_len(nrow(CONTINUOUS_PREDICTORS))) {
        plot_continuous_panel(
          scale_data,
          CONTINUOUS_PREDICTORS$predictor[j],
          CONTINUOUS_PREDICTORS$predictor_label[j],
          boot_stat,
          CONTINUOUS_PREDICTORS$color[j],
          relationships
        )
      }
    }
    mtext(paste0(scale_name, " SA bootstrap statistics vs continuous context variables"), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    dev.off()
    files[i] <- file_path
  }
  files
}

save_composition_figures <- function(data, relationships) {
  files <- character(length(SCALES))
  palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#1F78B4")
  for (i in seq_along(SCALES)) {
    scale_name <- SCALES[i]
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    type_levels <- sort(unique(scale_data$composition_type[!is.na(scale_data$composition_type) & nzchar(scale_data$composition_type)]))
    scale_data$composition_type <- factor(scale_data$composition_type, levels = type_levels)
    type_colors <- setNames(rep(palette, length.out = length(type_levels)), type_levels)
    file_path <- file.path(FIGURE_DIR, paste0("04_sa_bootstrap_stats_vs_species_composition_", scale_name, ".png"))
    png(file_path, width = 3600, height = 1300, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(1, 3), mar = c(4.9, 4.4, 2.5, 0.8), oma = c(0, 0, 2.1, 0))
    for (boot_stat in BOOT_STATS) {
      complete <- !is.na(scale_data$composition_type) & !is.na(scale_data[[boot_stat]])
      y <- scale_data[[boot_stat]][complete]
      group <- droplevels(scale_data$composition_type[complete])
      x <- as.numeric(group)
      plot(
        jitter(x, amount = 0.16),
        y,
        pch = 16,
        cex = 0.52,
        col = adjustcolor(type_colors[as.character(group)], alpha.f = 0.58),
        xaxt = "n",
        xlab = "Composition type",
        ylab = BOOT_STAT_LABELS[[boot_stat]],
        main = paste(BOOT_STAT_LABELS[[boot_stat]], "by composition")
      )
      axis(1, at = seq_along(type_levels), labels = seq_along(type_levels), cex.axis = 0.78)
      boxplot(y ~ group, add = TRUE, outline = FALSE, border = "#111111", col = adjustcolor("#FFFFFF", alpha.f = 0), xaxt = "n")
      row <- relationships[
        relationships$relationship_type == "composition" &
          relationships$scale == scale_name &
          relationships$boot_stat == boot_stat,
      ]
      if (nrow(row) == 1) {
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("groups=", row$n_groups),
            paste0("eta2=", fmt_num(row$eta_squared, 2)),
            paste0("p=", fmt_p(row$f_p_value))
          ),
          bty = "n",
          cex = 0.72
        )
      }
      if (boot_stat == BOOT_STATS[1]) {
        legend(
          "topright",
          legend = paste0(seq_along(type_levels), " = ", type_levels),
          fill = type_colors[type_levels],
          border = NA,
          bty = "n",
          cex = 0.58
        )
      }
    }
    mtext(paste0(scale_name, " SA bootstrap statistics vs species composition type"), outer = TRUE, cex = 1.1, font = 2)
    par(old_par)
    dev.off()
    files[i] <- file_path
  }
  files
}

write_reports <- function(data, relationships, figure_files) {
  top_cont <- relationships[relationships$relationship_type == "continuous", , drop = FALSE]
  top_cont <- top_cont[order(-abs(top_cont$pearson_r)), ]
  top_comp <- relationships[relationships$relationship_type == "composition", , drop = FALSE]
  top_comp <- top_comp[order(-top_comp$eta_squared), ]
  top_table <- rbind(
    data.frame(
      Type = "continuous",
      Scale = head(top_cont$scale, 10),
      Statistic = head(top_cont$boot_stat_label, 10),
      Predictor = head(top_cont$predictor_label, 10),
      n = head(top_cont$n, 10),
      Effect = fmt_num(head(top_cont$pearson_r, 10), 3),
      R2_or_eta2 = fmt_num(head(top_cont$r_squared, 10), 3),
      `p-value` = fmt_p(head(top_cont$f_p_value, 10)),
      check.names = FALSE
    ),
    data.frame(
      Type = "composition",
      Scale = head(top_comp$scale, 6),
      Statistic = head(top_comp$boot_stat_label, 6),
      Predictor = head(top_comp$predictor_label, 6),
      n = head(top_comp$n, 6),
      Effect = "ANOVA",
      R2_or_eta2 = fmt_num(head(top_comp$eta_squared, 6), 3),
      `p-value` = fmt_p(head(top_comp$f_p_value, 6)),
      check.names = FALSE
    )
  )

  analysis_lines <- c(
    "# SA Bootstrap Context Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis plots the three final SA entropy bootstrap statistics against elevation, species count, crown-equivalent individual count, species composition type, overall brightness, and blue, green, red, and near-infrared brightness.",
    "",
    "## Data",
    "",
    "- SA bootstrap statistics come from `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`.",
    "- Elevation comes from the current `quadrat_analysis_*m.csv` tables.",
    "- Species count and crown-equivalent individual count come from `quadrat_crown_equivalent_individual_totals.csv`.",
    "- Species composition types come from the presence/absence cluster table `species_composition_clusters.csv`.",
    "- Brightness variables come from `quadrat_pixel_brightness_summary.csv`.",
    "",
    "## Strongest Context Relationships",
    "",
    markdown_table(top_table),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Tables",
    "",
    "- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_dataset.csv`",
    "- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_relationships.csv`"
  )
  writeLines(analysis_lines, ANALYSIS_REPORT)

  task_lines <- c(
    "# Task Report: SA Bootstrap Context Relationship Figures",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Create scatterplot figures for final SA entropy bootstrap SD, CV, and SE against elevation, species count, individual count, species composition, and brightness variables at each quadrat scale.",
    "",
    "## Outputs",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_dataset.csv`",
    "- `reports/tables/bootstrap_metric_relationships/sa_bootstrap_context_relationships.csv`"
  )
  writeLines(task_lines, TASK_REPORT)

  validation_lines <- c(
    "# Validation: SA Bootstrap Context Relationship Figures",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- Context dataset rows: ", nrow(data)),
    paste0("- Relationship rows: ", nrow(relationships)),
    paste0("- Continuous relationship rows: ", sum(relationships$relationship_type == "continuous")),
    paste0("- Composition relationship rows: ", sum(relationships$relationship_type == "composition")),
    paste0("- Output figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    paste0("- Missing elevation values: ", sum(is.na(data$env_elev))),
    paste0("- Missing species count values: ", sum(is.na(data$present_species_count))),
    paste0("- Missing individual count values: ", sum(is.na(data$crown_equivalent_individuals))),
    paste0("- Missing composition values: ", sum(is.na(data$composition_type) | !nzchar(data$composition_type))),
    paste0("- Missing overall brightness values: ", sum(is.na(data$mean_pixel_brightness))),
    "",
    "## Result",
    "",
    "The requested SA bootstrap context figures and tables were produced."
  )
  writeLines(validation_lines, VALIDATION_REPORT)
}

run_sa_bootstrap_context_relationship_analysis <- function() {
  data <- read_context_dataset()
  relationships <- summarize_relationships(data)
  figure_files <- c(
    save_continuous_figures(data, relationships),
    save_composition_figures(data, relationships)
  )

  write.csv(data, file.path(TABLE_DIR, "sa_bootstrap_context_dataset.csv"), row.names = FALSE)
  write.csv(relationships, file.path(TABLE_DIR, "sa_bootstrap_context_relationships.csv"), row.names = FALSE)
  write_reports(data, relationships, figure_files)

  message("SA bootstrap context relationship analysis complete.")
  invisible(list(data = data, relationships = relationships, figure_files = figure_files))
}

if (sys.nframe() == 0) {
  run_sa_bootstrap_context_relationship_analysis()
}
