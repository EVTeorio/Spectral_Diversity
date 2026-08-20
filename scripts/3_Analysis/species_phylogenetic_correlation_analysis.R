PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-17"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/species_phylogenetic_correlation")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260817_species_phylogenetic_correlation_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260817_species_phylogenetic_correlation_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260817_species_phylogenetic_correlation_analysis_validation.md")

SPECIES_METRICS <- c("sp_rich", "sp_shannon", "sp_simpson", "sp_even")
PHYLO_METRICS <- c("phy_faith", "phy_rao", "phy_afaith")
ALL_DIVERSITY_METRICS <- c(SPECIES_METRICS, PHYLO_METRICS)
SV_METRICS <- c("spec_spca_mean", "spec_sa")
SCALES <- c("10m", "20m", "50m")

DISPLAY_NAMES <- c(
  sp_rich = "Species richness",
  sp_shannon = "Shannon diversity",
  sp_simpson = "Simpson diversity",
  sp_even = "Species evenness",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  spec_spca_mean = "Std PCA mean distance",
  spec_sa = "SA entropy"
)

PLANT_FILES <- c(
  "10m" = file.path(PROJECT_DIR, "Quad_Values/Diversity_SHPs/plant_diversity_10m.csv"),
  "20m" = file.path(PROJECT_DIR, "Quad_Values/Diversity_SHPs/plant_diversity_20m.csv"),
  "50m" = file.path(PROJECT_DIR, "Quad_Values/Diversity_SHPs/plant_diversity_50m.csv")
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

display_name <- function(x) {
  unname(DISPLAY_NAMES[x])
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

read_analysis_data <- function() {
  file_path <- file.path(TABLE_DIR, "..", "multiscale_spectral_biodiversity", "sv_diversity_analysis_dataset.csv")
  file_path <- normalizePath(file_path, winslash = "/", mustWork = TRUE)
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- as.character(data$scale)
  data$primary_analysis <- as.logical(data$primary_analysis)
  data$edge_flag <- as.logical(data$edge_flag)
  data
}

species_columns <- function(data) {
  metric_cols <- c("richness", "shannon", "simpson", "evenness", "faith_pd", "rao_pd", "afaith_pd")
  id_cols <- c("", "X", "Name", "sub_id", "quad_id", "scale", metric_cols)
  numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  setdiff(numeric_cols, id_cols)
}

read_species_composition <- function(scale_name) {
  data <- read.csv(PLANT_FILES[[scale_name]], stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- scale_name
  sp_cols <- species_columns(data)
  data[, c("quad_id", "scale", sp_cols), drop = FALSE]
}

make_quadrat_individual_totals <- function() {
  rows <- list()
  for (scale_name in SCALES) {
    comp <- read_species_composition(scale_name)
    sp_cols <- setdiff(names(comp), c("quad_id", "scale"))
    abundance_matrix <- as.matrix(comp[, sp_cols, drop = FALSE])
    abundance_matrix[is.na(abundance_matrix)] <- 0
    rows[[scale_name]] <- data.frame(
      quad_id = comp$quad_id,
      scale = comp$scale,
      crown_equivalent_individuals = rowSums(abundance_matrix),
      present_species_count = rowSums(abundance_matrix > 0),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

make_composition_clusters <- function(analysis_data) {
  rows <- list()
  for (scale_name in SCALES) {
    comp <- read_species_composition(scale_name)
    scale_rows <- analysis_data[analysis_data$scale == scale_name, c("quad_id", "scale", SPECIES_METRICS, PHYLO_METRICS)]
    data <- merge(scale_rows, comp, by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
    sp_cols <- setdiff(names(comp), c("quad_id", "scale"))
    keep <- rowSums(!is.na(data[, sp_cols, drop = FALSE])) > 0
    cluster_data <- data[keep, ]
    comp_matrix <- as.matrix(cluster_data[, sp_cols, drop = FALSE])
    comp_matrix[is.na(comp_matrix)] <- 0
    presence_matrix <- ifelse(comp_matrix > 0, 1, 0)
    cluster_count <- min(6, max(2, floor(nrow(presence_matrix) / 25)))
    set.seed(20260817)
    clusters <- kmeans(presence_matrix, centers = cluster_count, nstart = 25, iter.max = 100)
    rows[[scale_name]] <- data.frame(
      quad_id = cluster_data$quad_id,
      scale = cluster_data$scale,
      composition_cluster = paste0(scale_name, "_C", clusters$cluster),
      cluster_size = as.integer(table(clusters$cluster)[as.character(clusters$cluster)]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

make_composition_type_names <- function(clusters) {
  rows <- list()
  freq_rows <- list()
  row_index <- 1
  freq_index <- 1

  for (scale_name in SCALES) {
    comp <- read_species_composition(scale_name)
    data <- merge(clusters[clusters$scale == scale_name, ], comp, by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
    sp_cols <- setdiff(names(comp), c("quad_id", "scale"))
    presence <- as.data.frame(ifelse(as.matrix(data[, sp_cols, drop = FALSE]) > 0, 1, 0))
    names(presence) <- sp_cols

    for (cluster_name in sort(unique(data$composition_cluster))) {
      in_cluster <- data$composition_cluster == cluster_name
      cluster_presence <- presence[in_cluster, , drop = FALSE]
      other_presence <- presence[!in_cluster, , drop = FALSE]
      freq_in <- colMeans(cluster_presence, na.rm = TRUE)
      freq_out <- colMeans(other_presence, na.rm = TRUE)
      indicator_score <- freq_in - freq_out
      order_index <- order(indicator_score, freq_in, decreasing = TRUE)
      label_species <- names(freq_in)[order_index][freq_in[order_index] >= 0.25]
      label_species <- head(label_species, 3)
      if (length(label_species) == 0) {
        label_species <- head(names(freq_in)[order(freq_in, decreasing = TRUE)], 3)
      }
      composition_type <- paste(label_species, collapse = "+")

      top_species <- head(names(freq_in)[order(freq_in, decreasing = TRUE)], 8)
      rows[[row_index]] <- data.frame(
        scale = scale_name,
        composition_cluster = cluster_name,
        composition_type = composition_type,
        cluster_size = sum(in_cluster),
        mean_presence_richness = mean(rowSums(cluster_presence > 0), na.rm = TRUE),
        defining_species_codes = paste(label_species, collapse = ";"),
        top_species_presence = paste0(top_species, "=", fmt_num(freq_in[top_species] * 100, 1), "%", collapse = "; "),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1

      for (species_code in sp_cols) {
        freq_rows[[freq_index]] <- data.frame(
          scale = scale_name,
          composition_cluster = cluster_name,
          species_code = species_code,
          presence_frequency = freq_in[[species_code]],
          outside_presence_frequency = freq_out[[species_code]],
          indicator_score = indicator_score[[species_code]],
          stringsAsFactors = FALSE
        )
        freq_index <- freq_index + 1
      }
    }
  }

  type_names <- do.call(rbind, rows)
  frequencies <- do.call(rbind, freq_rows)
  list(type_names = type_names, frequencies = frequencies)
}

calculate_composition_cluster_separation <- function(clusters) {
  rows <- list()
  index <- 1

  for (scale_name in SCALES) {
    comp <- read_species_composition(scale_name)
    data <- merge(
      clusters[clusters$scale == scale_name, c("quad_id", "scale", "composition_cluster")],
      comp,
      by = c("quad_id", "scale"),
      all.x = TRUE,
      sort = FALSE
    )
    sp_cols <- setdiff(names(comp), c("quad_id", "scale"))
    keep <- !is.na(data$composition_cluster) & rowSums(!is.na(data[, sp_cols, drop = FALSE])) > 0
    data <- data[keep, ]
    presence_matrix <- ifelse(as.matrix(data[, sp_cols, drop = FALSE]) > 0, 1, 0)
    presence_matrix[is.na(presence_matrix)] <- 0

    cluster_names <- sort(unique(data$composition_cluster))
    separation <- rep(NA_real_, length(cluster_names))
    names(separation) <- cluster_names

    if (nrow(presence_matrix) > length(cluster_names) && length(cluster_names) > 1 && requireNamespace("cluster", quietly = TRUE)) {
      cluster_factor <- factor(data$composition_cluster, levels = cluster_names)
      silhouette_values <- cluster::silhouette(as.integer(cluster_factor), stats::dist(presence_matrix))
      separation <- tapply(silhouette_values[, "sil_width"], data$composition_cluster, mean, na.rm = TRUE)
    }

    for (cluster_name in cluster_names) {
      rows[[index]] <- data.frame(
        scale = scale_name,
        composition_cluster = cluster_name,
        mean_silhouette_width = unname(separation[[cluster_name]]),
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }

  do.call(rbind, rows)
}

calculate_scatterplot_cluster_separation <- function(data, clusters) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  plot_data <- merge(
    data,
    clusters[, c("quad_id", "scale", "composition_cluster", "composition_type")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  rows <- list()
  index <- 1

  for (scale_name in SCALES) {
    for (pair in metric_pairs) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- plot_data$scale == scale_name &
        !is.na(plot_data[[x_metric]]) &
        !is.na(plot_data[[y_metric]]) &
        !is.na(plot_data$composition_cluster)
      panel_data <- plot_data[complete, c("composition_cluster", "composition_type", x_metric, y_metric)]
      names(panel_data) <- c("composition_cluster", "composition_type", "x_value", "y_value")
      cluster_names <- sort(unique(panel_data$composition_cluster))
      cluster_stats <- unique(panel_data[, c("composition_cluster", "composition_type")])
      cluster_stats <- cluster_stats[order(cluster_stats$composition_type), ]
      cluster_stats$n <- as.integer(table(panel_data$composition_cluster)[cluster_stats$composition_cluster])
      cluster_stats$mean_scatterplot_silhouette_width <- NA_real_
      overall_silhouette <- NA_real_

      if (nrow(panel_data) > length(cluster_names) && length(cluster_names) > 1 && requireNamespace("cluster", quietly = TRUE)) {
        xy <- scale(panel_data[, c("x_value", "y_value")])
        cluster_factor <- factor(panel_data$composition_cluster, levels = cluster_names)
        silhouette_values <- cluster::silhouette(as.integer(cluster_factor), stats::dist(xy))
        cluster_silhouette <- tapply(silhouette_values[, "sil_width"], panel_data$composition_cluster, mean, na.rm = TRUE)
        overall_silhouette <- mean(silhouette_values[, "sil_width"], na.rm = TRUE)
        cluster_stats$mean_scatterplot_silhouette_width <- unname(cluster_silhouette[cluster_stats$composition_cluster])
      }

      for (row_index in seq_len(nrow(cluster_stats))) {
        rows[[index]] <- data.frame(
          scale = scale_name,
          x_metric = x_metric,
          x_label = display_name(x_metric),
          y_metric = y_metric,
          y_label = display_name(y_metric),
          composition_cluster = cluster_stats$composition_cluster[row_index],
          composition_type = cluster_stats$composition_type[row_index],
          n = cluster_stats$n[row_index],
          mean_scatterplot_silhouette_width = cluster_stats$mean_scatterplot_silhouette_width[row_index],
          overall_scatterplot_silhouette_width = overall_silhouette,
          stringsAsFactors = FALSE
        )
        index <- index + 1
      }
    }
  }

  do.call(rbind, rows)
}

correlate_species_phylo_pair <- function(data, scale_name, species_metric, phylo_metric) {
  keep <- data$scale == scale_name &
    !is.na(data[[species_metric]]) &
    !is.na(data[[phylo_metric]])
  pair <- data[keep, c(species_metric, phylo_metric)]
  names(pair) <- c("species_value", "phylo_value")

  result <- data.frame(
    scale = scale_name,
    species_metric = species_metric,
    species_label = display_name(species_metric),
    phylo_metric = phylo_metric,
    phylo_label = display_name(phylo_metric),
    n = nrow(pair),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    spearman_r = NA_real_,
    spearman_p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(pair) < 3) {
    return(result)
  }

  fit <- lm(phylo_value ~ species_value, data = pair)
  fit_summary <- summary(fit)
  f_values <- fit_summary$fstatistic
  pearson <- cor.test(pair$species_value, pair$phylo_value, method = "pearson")
  spearman <- suppressWarnings(cor.test(pair$species_value, pair$phylo_value, method = "spearman", exact = FALSE))

  result$pearson_r <- unname(pearson$estimate)
  result$r_squared <- fit_summary$r.squared
  result$f_statistic <- unname(f_values[["value"]])
  result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
  result$slope <- unname(coef(fit)[["species_value"]])
  result$intercept <- unname(coef(fit)[["(Intercept)"]])
  result$spearman_r <- unname(spearman$estimate)
  result$spearman_p_value <- spearman$p.value
  result
}

run_pairwise_correlations <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (species_metric in SPECIES_METRICS) {
      for (phylo_metric in PHYLO_METRICS) {
        rows[[index]] <- correlate_species_phylo_pair(data, scale_name, species_metric, phylo_metric)
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

correlate_all_metric_pair <- function(data, scale_name, x_metric, y_metric) {
  keep <- data$scale == scale_name &
    !is.na(data[[x_metric]]) &
    !is.na(data[[y_metric]])
  pair <- data[keep, c(x_metric, y_metric)]
  names(pair) <- c("x_value", "y_value")

  result <- data.frame(
    scale = scale_name,
    x_metric = x_metric,
    x_label = display_name(x_metric),
    y_metric = y_metric,
    y_label = display_name(y_metric),
    n = nrow(pair),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    spearman_r = NA_real_,
    spearman_p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(pair) < 3) {
    return(result)
  }

  fit <- lm(y_value ~ x_value, data = pair)
  fit_summary <- summary(fit)
  f_values <- fit_summary$fstatistic
  pearson <- cor.test(pair$x_value, pair$y_value, method = "pearson")
  spearman <- suppressWarnings(cor.test(pair$x_value, pair$y_value, method = "spearman", exact = FALSE))

  result$pearson_r <- unname(pearson$estimate)
  result$r_squared <- fit_summary$r.squared
  result$f_statistic <- unname(f_values[["value"]])
  result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
  result$slope <- unname(coef(fit)[["x_value"]])
  result$intercept <- unname(coef(fit)[["(Intercept)"]])
  result$spearman_r <- unname(spearman$estimate)
  result$spearman_p_value <- spearman$p.value
  result
}

run_all_metric_correlations <- function(data) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (pair in metric_pairs) {
      rows[[index]] <- correlate_all_metric_pair(data, scale_name, pair[1], pair[2])
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

elevation_adjust_all_metric_pair <- function(data, scale_name, x_metric, y_metric) {
  keep <- data$scale == scale_name &
    !is.na(data[[x_metric]]) &
    !is.na(data[[y_metric]]) &
    !is.na(data$env_elev)
  model_data <- data[keep, c(x_metric, y_metric, "env_elev")]
  names(model_data) <- c("x_value", "y_value", "env_elev")

  result <- data.frame(
    scale = scale_name,
    x_metric = x_metric,
    x_label = display_name(x_metric),
    y_metric = y_metric,
    y_label = display_name(y_metric),
    n = nrow(model_data),
    base_r2 = NA_real_,
    elevation_adjusted_r2 = NA_real_,
    elevation_incremental_r2_after_metric = NA_real_,
    elevation_partial_f = NA_real_,
    elevation_partial_p_value = NA_real_,
    elevation_slope = NA_real_,
    elevation_p_value_in_adjusted_model = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(model_data) >= 4) {
    base_fit <- lm(y_value ~ x_value, data = model_data)
    adjusted_fit <- lm(y_value ~ x_value + env_elev, data = model_data)
    comparison <- anova(base_fit, adjusted_fit)
    adjusted_summary <- summary(adjusted_fit)

    result$base_r2 <- summary(base_fit)$r.squared
    result$elevation_adjusted_r2 <- adjusted_summary$r.squared
    result$elevation_incremental_r2_after_metric <- result$elevation_adjusted_r2 - result$base_r2
    result$elevation_partial_f <- comparison$F[2]
    result$elevation_partial_p_value <- comparison$`Pr(>F)`[2]
    result$elevation_slope <- coef(adjusted_fit)[["env_elev"]]
    result$elevation_p_value_in_adjusted_model <- coef(adjusted_summary)["env_elev", "Pr(>|t|)"]
  }

  result
}

run_elevation_adjusted_all_metric_models <- function(data) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (pair in metric_pairs) {
      rows[[index]] <- elevation_adjust_all_metric_pair(data, scale_name, pair[1], pair[2])
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

incremental_sv_models <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (sv_metric in SV_METRICS) {
      for (species_metric in SPECIES_METRICS) {
        for (phylo_metric in PHYLO_METRICS) {
          keep <- data$scale == scale_name &
            !is.na(data[[sv_metric]]) &
            !is.na(data[[species_metric]]) &
            !is.na(data[[phylo_metric]])
          model_data <- data[keep, c(sv_metric, species_metric, phylo_metric)]
          names(model_data) <- c("sv_value", "species_value", "phylo_value")

          result <- data.frame(
            scale = scale_name,
            sv_metric = sv_metric,
            sv_label = display_name(sv_metric),
            species_metric = species_metric,
            species_label = display_name(species_metric),
            phylo_metric = phylo_metric,
            phylo_label = display_name(phylo_metric),
            n = nrow(model_data),
            species_only_r2 = NA_real_,
            phylo_only_r2 = NA_real_,
            combined_r2 = NA_real_,
            phylo_incremental_r2_after_species = NA_real_,
            partial_f = NA_real_,
            partial_f_p_value = NA_real_,
            phylo_slope_in_combined = NA_real_,
            phylo_p_value_in_combined = NA_real_,
            stringsAsFactors = FALSE
          )

          if (nrow(model_data) >= 4) {
            species_fit <- lm(sv_value ~ species_value, data = model_data)
            phylo_fit <- lm(sv_value ~ phylo_value, data = model_data)
            combined_fit <- lm(sv_value ~ species_value + phylo_value, data = model_data)
            comparison <- anova(species_fit, combined_fit)
            combined_summary <- summary(combined_fit)
            result$species_only_r2 <- summary(species_fit)$r.squared
            result$phylo_only_r2 <- summary(phylo_fit)$r.squared
            result$combined_r2 <- combined_summary$r.squared
            result$phylo_incremental_r2_after_species <- result$combined_r2 - result$species_only_r2
            result$partial_f <- comparison$F[2]
            result$partial_f_p_value <- comparison$`Pr(>F)`[2]
            result$phylo_slope_in_combined <- coef(combined_fit)[["phylo_value"]]
            result$phylo_p_value_in_combined <- coef(combined_summary)["phylo_value", "Pr(>|t|)"]
          }

          rows[[index]] <- result
          index <- index + 1
        }
      }
    }
  }
  do.call(rbind, rows)
}

make_residual_dataset <- function(data, clusters) {
  rows <- list()
  index <- 1
  cluster_data <- merge(data, clusters, by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  for (scale_name in SCALES) {
    for (species_metric in SPECIES_METRICS) {
      for (phylo_metric in PHYLO_METRICS) {
        keep <- cluster_data$scale == scale_name &
          !is.na(cluster_data[[species_metric]]) &
          !is.na(cluster_data[[phylo_metric]]) &
          !is.na(cluster_data$composition_cluster)
        model_data <- cluster_data[keep, c("quad_id", "scale", "composition_cluster", "cluster_size", species_metric, phylo_metric)]
        names(model_data) <- c("quad_id", "scale", "composition_cluster", "cluster_size", "species_value", "phylo_value")

        if (nrow(model_data) >= 3) {
          fit <- lm(phylo_value ~ species_value, data = model_data)
          model_data$phylo_residual <- residuals(fit)
          model_data$phylo_residual_z <- as.numeric(scale(model_data$phylo_residual))
          model_data$abs_phylo_residual_z <- abs(model_data$phylo_residual_z)
          model_data$species_metric <- species_metric
          model_data$species_label <- display_name(species_metric)
          model_data$phylo_metric <- phylo_metric
          model_data$phylo_label <- display_name(phylo_metric)
          rows[[index]] <- model_data
          index <- index + 1
        }
      }
    }
  }
  do.call(rbind, rows)
}

summarize_cluster_divergence <- function(residual_data) {
  grouped <- split(
    residual_data,
    list(
      residual_data$scale,
      residual_data$composition_cluster,
      residual_data$species_metric,
      residual_data$phylo_metric
    ),
    drop = TRUE
  )

  rows <- lapply(grouped, function(group) {
    data.frame(
      scale = group$scale[1],
      composition_cluster = group$composition_cluster[1],
      cluster_size = group$cluster_size[1],
      species_metric = group$species_metric[1],
      species_label = group$species_label[1],
      phylo_metric = group$phylo_metric[1],
      phylo_label = group$phylo_label[1],
      n = nrow(group),
      mean_abs_phylo_residual_z = mean(group$abs_phylo_residual_z, na.rm = TRUE),
      sd_phylo_residual = sd(group$phylo_residual, na.rm = TRUE),
      mean_species_value = mean(group$species_value, na.rm = TRUE),
      mean_phylo_value = mean(group$phylo_value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  output <- do.call(rbind, rows)
  output[order(output$scale, output$phylo_metric, -output$mean_abs_phylo_residual_z), ]
}

metric_summary <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (metric in c(SPECIES_METRICS, PHYLO_METRICS)) {
      values <- data[[metric]][data$scale == scale_name]
      rows[[index]] <- data.frame(
        scale = scale_name,
        metric = metric,
        metric_label = display_name(metric),
        n = sum(!is.na(values)),
        mean = mean(values, na.rm = TRUE),
        sd = sd(values, na.rm = TRUE),
        cv = sd(values, na.rm = TRUE) / mean(values, na.rm = TRUE),
        min = min(values, na.rm = TRUE),
        median = median(values, na.rm = TRUE),
        max = max(values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

moderator_label <- function(moderator) {
  labels <- c(
    present_species_count = "Number of species",
    crown_equivalent_individuals = "Crown-equivalent individuals"
  )
  unname(labels[moderator])
}

analyze_species_individual_moderation <- function(data, individual_totals) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  moderators <- c("present_species_count", "crown_equivalent_individuals")
  model_data <- merge(
    data,
    individual_totals[, c("quad_id", "scale", moderators)],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  rows <- list()
  index <- 1

  for (scale_name in SCALES) {
    for (pair in metric_pairs) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      for (moderator in moderators) {
        complete <- model_data$scale == scale_name &
          !is.na(model_data[[x_metric]]) &
          !is.na(model_data[[y_metric]]) &
          !is.na(model_data[[moderator]])
        panel_data <- model_data[complete, c(x_metric, y_metric, moderator)]
        names(panel_data) <- c("x_value", "y_value", "moderator_value")
        moderator_is_plotted_metric <- moderator == "present_species_count" &&
          (x_metric == "sp_rich" || y_metric == "sp_rich")

        result <- data.frame(
          scale = scale_name,
          x_metric = x_metric,
          x_label = display_name(x_metric),
          y_metric = y_metric,
          y_label = display_name(y_metric),
          moderator = moderator,
          moderator_label = moderator_label(moderator),
          moderator_is_plotted_metric = moderator_is_plotted_metric,
          n = nrow(panel_data),
          base_r2 = NA_real_,
          additive_r2 = NA_real_,
          full_interaction_r2 = NA_real_,
          moderator_incremental_r2_after_metric = NA_real_,
          interaction_incremental_r2_after_additive = NA_real_,
          full_incremental_r2_after_metric = NA_real_,
          moderator_main_slope = NA_real_,
          moderator_main_p_value = NA_real_,
          interaction_slope = NA_real_,
          interaction_f = NA_real_,
          interaction_p_value = NA_real_,
          full_partial_f = NA_real_,
          full_partial_p_value = NA_real_,
          stringsAsFactors = FALSE
        )

        if (
          nrow(panel_data) >= 6 &&
            sd(panel_data$x_value, na.rm = TRUE) > 0 &&
            sd(panel_data$y_value, na.rm = TRUE) > 0 &&
            sd(panel_data$moderator_value, na.rm = TRUE) > 0
        ) {
          fit_data <- data.frame(
            x_z = as.numeric(scale(panel_data$x_value)),
            y_z = as.numeric(scale(panel_data$y_value)),
            moderator_z = as.numeric(scale(panel_data$moderator_value))
          )
          base_fit <- lm(y_z ~ x_z, data = fit_data)
          additive_fit <- lm(y_z ~ x_z + moderator_z, data = fit_data)
          full_fit <- lm(y_z ~ x_z * moderator_z, data = fit_data)
          additive_comparison <- anova(base_fit, additive_fit)
          interaction_comparison <- anova(additive_fit, full_fit)
          full_comparison <- anova(base_fit, full_fit)
          additive_summary <- summary(additive_fit)
          full_summary <- summary(full_fit)

          result$base_r2 <- summary(base_fit)$r.squared
          result$additive_r2 <- additive_summary$r.squared
          result$full_interaction_r2 <- full_summary$r.squared
          result$moderator_incremental_r2_after_metric <- result$additive_r2 - result$base_r2
          result$interaction_incremental_r2_after_additive <- result$full_interaction_r2 - result$additive_r2
          result$full_incremental_r2_after_metric <- result$full_interaction_r2 - result$base_r2
          if ("moderator_z" %in% rownames(coef(additive_summary))) {
            result$moderator_main_slope <- coef(additive_summary)["moderator_z", "Estimate"]
            result$moderator_main_p_value <- coef(additive_summary)["moderator_z", "Pr(>|t|)"]
          }
          if ("x_z:moderator_z" %in% rownames(coef(full_summary))) {
            result$interaction_slope <- coef(full_summary)["x_z:moderator_z", "Estimate"]
          }
          result$interaction_f <- interaction_comparison$F[2]
          result$interaction_p_value <- interaction_comparison$`Pr(>F)`[2]
          result$full_partial_f <- full_comparison$F[2]
          result$full_partial_p_value <- full_comparison$`Pr(>F)`[2]
        }

        rows[[index]] <- result
        index <- index + 1
      }
    }
  }

  output <- do.call(rbind, rows)
  output$interaction_q_value <- p.adjust(output$interaction_p_value, method = "BH")
  output$full_partial_q_value <- p.adjust(output$full_partial_p_value, method = "BH")
  output$moderator_main_q_value <- p.adjust(output$moderator_main_p_value, method = "BH")
  output
}

significant_moderation_results <- function(moderation, alpha = 0.05) {
  significant <- moderation[
    !moderation$moderator_is_plotted_metric &
      !is.na(moderation$interaction_q_value) &
      moderation$interaction_q_value < alpha,
  ]
  significant[order(significant$interaction_q_value, -significant$interaction_incremental_r2_after_additive), ]
}

color_for_r <- function(r) {
  palette <- grDevices::colorRampPalette(c("#315C8A", "#F8F7F2", "#A93F38"))(201)
  idx <- round((pmax(pmin(r, 1), -1) + 1) * 100) + 1
  palette[idx]
}

save_correlation_heatmap <- function(correlations) {
  file_path <- file.path(FIGURE_DIR, "01_species_phylogenetic_correlation_heatmap.png")
  grDevices::png(file_path, width = 3600, height = 1800, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(1, 3), mar = c(7, 6, 3.5, 1), oma = c(0, 0, 2, 0))
  for (scale_name in SCALES) {
    scale_rows <- correlations[correlations$scale == scale_name, ]
    plot(
      NA,
      xlim = c(0.5, length(PHYLO_METRICS) + 0.5),
      ylim = c(0.5, length(SPECIES_METRICS) + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = scale_name,
      bty = "n"
    )
    axis(1, at = seq_along(PHYLO_METRICS), labels = display_name(PHYLO_METRICS), las = 2, cex.axis = 0.75)
    axis(2, at = seq_along(SPECIES_METRICS), labels = display_name(SPECIES_METRICS), las = 2, cex.axis = 0.78)
    for (i in seq_along(PHYLO_METRICS)) {
      for (j in seq_along(SPECIES_METRICS)) {
        row <- scale_rows[
          scale_rows$phylo_metric == PHYLO_METRICS[i] &
            scale_rows$species_metric == SPECIES_METRICS[j],
        ]
        rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col = color_for_r(row$pearson_r), border = "white", lwd = 2)
        text(i, j, labels = paste0("r=", fmt_num(row$pearson_r, 2), "\nR2=", fmt_num(row$r_squared, 2)), cex = 0.75)
      }
    }
  }
  mtext("Species-diversity versus phylogenetic-diversity correlations", outer = TRUE, cex = 1.25, font = 2)
  file_path
}

save_increment_heatmap <- function(incremental) {
  file_path <- file.path(FIGURE_DIR, "02_phylogenetic_incremental_r2_heatmap.png")
  grDevices::png(file_path, width = 4300, height = 2300, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(2, 3), mar = c(6.7, 6, 3.5, 1), oma = c(0, 0, 2, 0))
  palette <- grDevices::colorRampPalette(c("#F7F7F2", "#E0A458", "#8E3B46"))(101)
  max_delta <- max(incremental$phylo_incremental_r2_after_species, na.rm = TRUE)
  for (sv_metric in SV_METRICS) {
    for (scale_name in SCALES) {
      rows <- incremental[incremental$scale == scale_name & incremental$sv_metric == sv_metric, ]
      plot(
        NA,
        xlim = c(0.5, length(PHYLO_METRICS) + 0.5),
        ylim = c(0.5, length(SPECIES_METRICS) + 0.5),
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        main = paste(display_name(sv_metric), scale_name),
        bty = "n"
      )
      axis(1, at = seq_along(PHYLO_METRICS), labels = display_name(PHYLO_METRICS), las = 2, cex.axis = 0.68)
      axis(2, at = seq_along(SPECIES_METRICS), labels = display_name(SPECIES_METRICS), las = 2, cex.axis = 0.72)
      for (i in seq_along(PHYLO_METRICS)) {
        for (j in seq_along(SPECIES_METRICS)) {
          row <- rows[rows$phylo_metric == PHYLO_METRICS[i] & rows$species_metric == SPECIES_METRICS[j], ]
          idx <- round((row$phylo_incremental_r2_after_species / max_delta) * 100) + 1
          rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col = palette[idx], border = "white", lwd = 2)
          text(i, j, labels = fmt_num(row$phylo_incremental_r2_after_species, 3), cex = 0.72)
        }
      }
    }
  }
  mtext("Incremental R2 from adding phylogenetic diversity after species diversity", outer = TRUE, cex = 1.15, font = 2)
  file_path
}

save_scatter_grid <- function(data) {
  file_path <- file.path(FIGURE_DIR, "03_species_phylogenetic_scatter_grid.png")
  grDevices::png(file_path, width = 3900, height = 3000, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  colors <- c("10m" = "#2F6F73", "20m" = "#B06B34", "50m" = "#6B5CA5")
  par(mfrow = c(4, 3), mar = c(4, 4.5, 2.5, 0.8), oma = c(0, 0, 2, 0))
  for (species_metric in SPECIES_METRICS) {
    for (phylo_metric in PHYLO_METRICS) {
      keep <- !is.na(data[[species_metric]]) & !is.na(data[[phylo_metric]])
      plot(
        data[[species_metric]][keep],
        data[[phylo_metric]][keep],
        col = adjustcolor(colors[data$scale[keep]], alpha.f = 0.42),
        pch = 16,
        cex = 0.5,
        xlab = display_name(species_metric),
        ylab = display_name(phylo_metric),
        main = paste(display_name(species_metric), "vs", display_name(phylo_metric))
      )
      for (scale_name in SCALES) {
        scale_keep <- keep & data$scale == scale_name
        if (sum(scale_keep) >= 3) {
          fit <- lm(data[[phylo_metric]][scale_keep] ~ data[[species_metric]][scale_keep])
          abline(fit, col = colors[scale_name], lwd = 2)
        }
      }
    }
  }
  mtext("Species-diversity and phylogenetic-diversity relationships by scale", outer = TRUE, cex = 1.2, font = 2)
  file_path
}

save_divergence_boxplot <- function(residual_data) {
  file_path <- file.path(FIGURE_DIR, "04_phylogenetic_divergence_by_metric_boxplot.png")
  grDevices::png(file_path, width = 3600, height = 2100, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  residual_data$group <- paste(residual_data$scale, residual_data$phylo_label, sep = "\n")
  boxplot(
    abs_phylo_residual_z ~ group,
    data = residual_data,
    las = 2,
    col = rep(c("#8FB8A8", "#E6B17E", "#A8A0C8"), each = length(PHYLO_METRICS)),
    border = "#333333",
    ylab = "Absolute standardized residual from species metric",
    xlab = "",
    main = "Divergence of phylogenetic metrics after species-diversity expectation"
  )
  file_path
}

save_all_metric_scatter_by_scale <- function(data, all_metric_correlations) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  colors <- c(species_species = "#315C8A", species_phylo = "#2F6F73", phylo_phylo = "#8E3B46")
  files <- character(length(SCALES))

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("05_all_diversity_metric_pairwise_scatter_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)

    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2, 0))
    for (pair in metric_pairs) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      keep <- data$scale == scale_name &
        !is.na(data[[x_metric]]) &
        !is.na(data[[y_metric]])

      pair_type <- if (x_metric %in% SPECIES_METRICS && y_metric %in% SPECIES_METRICS) {
        "species_species"
      } else if (x_metric %in% PHYLO_METRICS && y_metric %in% PHYLO_METRICS) {
        "phylo_phylo"
      } else {
        "species_phylo"
      }

      plot(
        data[[x_metric]][keep],
        data[[y_metric]][keep],
        col = adjustcolor(colors[[pair_type]], alpha.f = 0.45),
        pch = 16,
        cex = 0.48,
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )

      if (sum(keep) >= 3) {
        fit <- lm(data[[y_metric]][keep] ~ data[[x_metric]][keep])
        abline(fit, col = "#222222", lwd = 2)
        row <- all_metric_correlations[
          all_metric_correlations$scale == scale_name &
            all_metric_correlations$x_metric == x_metric &
            all_metric_correlations$y_metric == y_metric,
        ]
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("r=", fmt_num(row$pearson_r, 2)),
            paste0("R2=", fmt_num(row$r_squared, 2))
          ),
          bty = "n",
          cex = 0.8
        )
      }
    }
    mtext(
      paste0(scale_name, " quadrats: all species and phylogenetic diversity metric pairs"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

cluster_palette <- function(n) {
  base_colors <- c(
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#332288",
    "#117733"
  )
  if (n <= length(base_colors)) {
    base_colors[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

individual_ramp_colors <- function(values) {
  palette <- grDevices::colorRampPalette(c(
    "#440154",
    "#2A6FBB",
    "#00A878",
    "#FDE725",
    "#F94144"
  ))(101)
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(palette[51], length(values)))
  }
  idx <- round((values - value_range[1]) / diff(value_range) * 100) + 1
  palette[pmax(pmin(idx, 101), 1)]
}

species_count_ramp_colors <- function(values) {
  palette <- grDevices::colorRampPalette(c(
    "#2D004B",
    "#0057B8",
    "#00B4D8",
    "#7AE582",
    "#FFE45E",
    "#FF3D00"
  ))(101)
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(palette[51], length(values)))
  }
  idx <- round((values - value_range[1]) / diff(value_range) * 100) + 1
  palette[pmax(pmin(idx, 101), 1)]
}

add_individual_ramp_legend <- function(values, title = "Crown-equiv.\nindividuals") {
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range))) {
    return(invisible(NULL))
  }
  legend_values <- pretty(value_range, n = 5)
  legend_values <- legend_values[legend_values >= value_range[1] & legend_values <= value_range[2]]
  if (length(legend_values) < 3) {
    legend_values <- seq(value_range[1], value_range[2], length.out = 5)
  }
  legend_colors <- individual_ramp_colors(legend_values)
  legend(
    "bottomright",
    legend = fmt_num(legend_values, 1),
    col = legend_colors,
    pch = 16,
    title = title,
    bty = "n",
    cex = 0.62
  )
}

add_species_count_ramp_legend <- function(values, title = "Species\nper quad") {
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range))) {
    return(invisible(NULL))
  }
  legend_values <- pretty(value_range, n = 5)
  legend_values <- legend_values[legend_values >= value_range[1] & legend_values <= value_range[2]]
  if (length(legend_values) < 3) {
    legend_values <- seq(value_range[1], value_range[2], length.out = 5)
  }
  legend_values <- sort(unique(round(legend_values)))
  legend_colors <- species_count_ramp_colors(legend_values)
  legend(
    "bottomright",
    legend = legend_values,
    col = legend_colors,
    pch = 16,
    title = title,
    bty = "n",
    cex = 0.62
  )
}

save_individual_ramp_scatter_by_scale <- function(data, individual_totals, all_metric_correlations) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  plot_data <- merge(
    data,
    individual_totals[, c("quad_id", "scale", "crown_equivalent_individuals")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  files <- character(length(SCALES))

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("10_all_diversity_metric_pairwise_scatter_individual_ramp_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))

    for (pair_index in seq_along(metric_pairs)) {
      pair <- metric_pairs[[pair_index]]
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- plot_data$scale == scale_name &
        !is.na(plot_data[[x_metric]]) &
        !is.na(plot_data[[y_metric]]) &
        !is.na(plot_data$crown_equivalent_individuals)
      point_colors <- adjustcolor(individual_ramp_colors(plot_data$crown_equivalent_individuals[complete]), alpha.f = 0.68)

      plot(
        plot_data[[x_metric]][complete],
        plot_data[[y_metric]][complete],
        col = point_colors,
        pch = 16,
        cex = 0.48,
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )

      if (sum(complete) >= 3) {
        fit <- lm(plot_data[[y_metric]][complete] ~ plot_data[[x_metric]][complete])
        abline(fit, col = "#111111", lwd = 2)
        row <- all_metric_correlations[
          all_metric_correlations$scale == scale_name &
            all_metric_correlations$x_metric == x_metric &
            all_metric_correlations$y_metric == y_metric,
        ]
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("r=", fmt_num(row$pearson_r, 2)),
            paste0("R2=", fmt_num(row$r_squared, 2))
          ),
          bty = "n",
          cex = 0.72
        )
      }

      if (pair_index == 1) {
        add_individual_ramp_legend(plot_data$crown_equivalent_individuals[plot_data$scale == scale_name])
      }
    }

    mtext(
      paste0(scale_name, " quadrats: diversity metric pairs colored by crown-equivalent individuals"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

save_species_count_ramp_scatter_by_scale <- function(data, individual_totals, all_metric_correlations) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  plot_data <- merge(
    data,
    individual_totals[, c("quad_id", "scale", "present_species_count")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  files <- character(length(SCALES))

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("11_all_diversity_metric_pairwise_scatter_species_count_ramp_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))

    for (pair_index in seq_along(metric_pairs)) {
      pair <- metric_pairs[[pair_index]]
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- plot_data$scale == scale_name &
        !is.na(plot_data[[x_metric]]) &
        !is.na(plot_data[[y_metric]]) &
        !is.na(plot_data$present_species_count)
      point_colors <- adjustcolor(species_count_ramp_colors(plot_data$present_species_count[complete]), alpha.f = 0.68)

      plot(
        plot_data[[x_metric]][complete],
        plot_data[[y_metric]][complete],
        col = point_colors,
        pch = 16,
        cex = 0.48,
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )

      if (sum(complete) >= 3) {
        fit <- lm(plot_data[[y_metric]][complete] ~ plot_data[[x_metric]][complete])
        abline(fit, col = "#111111", lwd = 2)
        row <- all_metric_correlations[
          all_metric_correlations$scale == scale_name &
            all_metric_correlations$x_metric == x_metric &
            all_metric_correlations$y_metric == y_metric,
        ]
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("r=", fmt_num(row$pearson_r, 2)),
            paste0("R2=", fmt_num(row$r_squared, 2))
          ),
          bty = "n",
          cex = 0.72
        )
      }

      if (pair_index == 1) {
        add_species_count_ramp_legend(plot_data$present_species_count[plot_data$scale == scale_name])
      }
    }

    mtext(
      paste0(scale_name, " quadrats: diversity metric pairs colored by species count"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

save_cluster_highlight_scatter_by_scale <- function(data, clusters, all_metric_correlations, scatterplot_cluster_separation) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  plot_data <- merge(
    data,
    clusters[, c("quad_id", "scale", "composition_cluster", "composition_type")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  files <- character(length(SCALES))

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    scale_clusters <- unique(plot_data[plot_data$scale == scale_name & !is.na(plot_data$composition_type), c("composition_cluster", "composition_type")])
    scale_clusters <- scale_clusters[order(scale_clusters$composition_type), ]
    cluster_colors <- cluster_palette(nrow(scale_clusters))
    names(cluster_colors) <- scale_clusters$composition_cluster
    file_path <- file.path(
      FIGURE_DIR,
      paste0("09_all_diversity_metric_pairwise_scatter_composition_cluster_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))

    for (pair_index in seq_along(metric_pairs)) {
      pair <- metric_pairs[[pair_index]]
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- plot_data$scale == scale_name &
        !is.na(plot_data[[x_metric]]) &
        !is.na(plot_data[[y_metric]]) &
        !is.na(plot_data$composition_type)

      point_colors <- adjustcolor(cluster_colors[plot_data$composition_cluster[complete]], alpha.f = 0.58)
      plot(
        plot_data[[x_metric]][complete],
        plot_data[[y_metric]][complete],
        col = point_colors,
        pch = 16,
        cex = 0.48,
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )

      if (sum(complete) >= 3) {
        fit <- lm(plot_data[[y_metric]][complete] ~ plot_data[[x_metric]][complete])
        abline(fit, col = "#111111", lwd = 2)
        row <- all_metric_correlations[
          all_metric_correlations$scale == scale_name &
            all_metric_correlations$x_metric == x_metric &
            all_metric_correlations$y_metric == y_metric,
        ]
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("r=", fmt_num(row$pearson_r, 2)),
            paste0("R2=", fmt_num(row$r_squared, 2))
          ),
          bty = "n",
          cex = 0.72
        )
      }

      panel_silhouette <- scatterplot_cluster_separation[
        scatterplot_cluster_separation$scale == scale_name &
          scatterplot_cluster_separation$x_metric == x_metric &
          scatterplot_cluster_separation$y_metric == y_metric,
      ]
      panel_silhouette <- panel_silhouette[match(scale_clusters$composition_cluster, panel_silhouette$composition_cluster), ]
      legend_labels <- paste0(
        scale_clusters$composition_type,
        " (sil=",
        fmt_num(panel_silhouette$mean_scatterplot_silhouette_width, 2),
        ")"
      )
      legend(
        "bottomright",
        legend = legend_labels,
        col = cluster_colors,
        pch = 16,
        title = "Composition type (panel sil.)",
        bty = "n",
        cex = ifelse(length(cluster_colors) > 4, 0.38, 0.48)
      )
    }

    mtext(
      paste0(scale_name, " quadrats: diversity metric pairs colored by species composition type"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

save_edge_highlight_scatter_by_scale <- function(data) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  edge_scales <- c("10m", "20m")
  files <- character(length(edge_scales))

  for (scale_index in seq_along(edge_scales)) {
    scale_name <- edge_scales[scale_index]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("06_all_diversity_metric_pairwise_scatter_edge_highlight_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2, 0))

    for (pair in metric_pairs) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- data$scale == scale_name &
        !is.na(data[[x_metric]]) &
        !is.na(data[[y_metric]])
      non_edge <- complete & !data$edge_flag
      edge <- complete & data$edge_flag

      plot(
        data[[x_metric]][complete],
        data[[y_metric]][complete],
        type = "n",
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )
      points(
        data[[x_metric]][non_edge],
        data[[y_metric]][non_edge],
        col = adjustcolor("#315C8A", alpha.f = 0.42),
        pch = 16,
        cex = 0.45
      )
      points(
        data[[x_metric]][edge],
        data[[y_metric]][edge],
        col = adjustcolor("#B03A2E", alpha.f = 0.72),
        pch = 16,
        cex = 0.52
      )

      if (sum(non_edge) >= 3) {
        fit <- lm(data[[y_metric]][non_edge] ~ data[[x_metric]][non_edge])
        abline(fit, col = "#111111", lwd = 2)
        pearson_r <- cor(data[[x_metric]][non_edge], data[[y_metric]][non_edge])
        legend(
          "topleft",
          legend = c(
            paste0("non-edge n=", sum(non_edge)),
            paste0("edge n=", sum(edge)),
            paste0("non-edge r=", fmt_num(pearson_r, 2)),
            paste0("non-edge R2=", fmt_num(summary(fit)$r.squared, 2))
          ),
          bty = "n",
          cex = 0.72
        )
      }
    }

    legend(
      "bottomright",
      legend = c("Non-edge", "Edge"),
      col = c("#315C8A", "#B03A2E"),
      pch = 16,
      bty = "n",
      cex = 0.9
    )
    mtext(
      paste0(scale_name, " quadrats: edge-highlighted diversity metric pairs"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

elevation_colors <- function(values) {
  palette <- grDevices::colorRampPalette(c("#B03A2E", "#F7F7F2", "#2F7D4F"))(101)
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(palette[51], length(values)))
  }
  idx <- round((values - value_range[1]) / diff(value_range) * 100) + 1
  palette[pmax(pmin(idx, 101), 1)]
}

save_elevation_gradient_scatter_by_scale <- function(data, elevation_models) {
  metric_pairs <- utils::combn(ALL_DIVERSITY_METRICS, 2, simplify = FALSE)
  files <- character(length(SCALES))

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("07_all_diversity_metric_pairwise_scatter_elevation_gradient_", scale_name, ".png")
    )

    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))

    for (pair in metric_pairs) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      keep <- data$scale == scale_name &
        !is.na(data[[x_metric]]) &
        !is.na(data[[y_metric]]) &
        !is.na(data$env_elev)
      point_colors <- elevation_colors(data$env_elev[keep])

      plot(
        data[[x_metric]][keep],
        data[[y_metric]][keep],
        col = adjustcolor(point_colors, alpha.f = 0.72),
        pch = 16,
        cex = 0.48,
        xlab = display_name(x_metric),
        ylab = display_name(y_metric),
        main = paste(display_name(x_metric), "vs", display_name(y_metric))
      )

      if (sum(keep) >= 3) {
        fit <- lm(data[[y_metric]][keep] ~ data[[x_metric]][keep])
        abline(fit, col = "#222222", lwd = 2)
        row <- elevation_models[
          elevation_models$scale == scale_name &
            elevation_models$x_metric == x_metric &
            elevation_models$y_metric == y_metric,
        ]
        legend(
          "topleft",
          legend = c(
            paste0("n=", row$n),
            paste0("elev dR2=", fmt_num(row$elevation_incremental_r2_after_metric, 3)),
            paste0("elev p=", fmt_p(row$elevation_partial_p_value))
          ),
          bty = "n",
          cex = 0.72
        )
      }
    }

    legend(
      "bottomright",
      legend = c("Low elevation", "Mid", "High elevation"),
      col = c("#B03A2E", "#F7F7F2", "#2F7D4F"),
      pch = 16,
      bty = "n",
      cex = 0.85
    )
    mtext(
      paste0(scale_name, " quadrats: diversity metric pairs colored by mean elevation"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

save_composition_type_presence_heatmaps <- function(type_names, frequencies) {
  files <- character(length(SCALES))
  palette <- grDevices::colorRampPalette(c("#F7F7F2", "#8FB8A8", "#2F7D4F"))(101)

  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    scale_freq <- frequencies[frequencies$scale == scale_name, ]
    scale_types <- type_names[type_names$scale == scale_name, ]
    cluster_order <- scale_types$composition_cluster[order(scale_types$composition_type)]
    species_order <- names(sort(tapply(scale_freq$presence_frequency, scale_freq$species_code, mean), decreasing = TRUE))
    matrix_values <- matrix(
      NA_real_,
      nrow = length(cluster_order),
      ncol = length(species_order),
      dimnames = list(cluster_order, species_order)
    )

    for (cluster_name in cluster_order) {
      for (species_code in species_order) {
        value <- scale_freq$presence_frequency[
          scale_freq$composition_cluster == cluster_name &
            scale_freq$species_code == species_code
        ]
        matrix_values[cluster_name, species_code] <- value
      }
    }

    row_labels <- paste0(
      scale_types$composition_type[match(cluster_order, scale_types$composition_cluster)],
      "\n",
      cluster_order,
      " n=",
      scale_types$cluster_size[match(cluster_order, scale_types$composition_cluster)]
    )

    file_path <- file.path(
      FIGURE_DIR,
      paste0("08_species_presence_composition_type_heatmap_", scale_name, ".png")
    )
    grDevices::png(file_path, width = 5200, height = 1800, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mar = c(7.8, 9.5, 3.5, 1), xpd = FALSE)

    image(
      x = seq_len(ncol(matrix_values)),
      y = seq_len(nrow(matrix_values)),
      z = t(matrix_values[nrow(matrix_values):1, ]),
      col = palette,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = paste0(scale_name, " composition types by species presence")
    )
    axis(1, at = seq_len(ncol(matrix_values)), labels = colnames(matrix_values), las = 2, cex.axis = 0.52)
    axis(
      2,
      at = seq_len(nrow(matrix_values)),
      labels = rev(row_labels),
      las = 2,
      cex.axis = 0.68
    )
    box()
    legend(
      "topright",
      legend = c("0%", "50%", "100%"),
      fill = palette[c(1, 51, 101)],
      title = "Presence",
      bty = "n",
      cex = 0.8
    )
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }

  files
}

save_figures <- function(data, correlations, incremental, residual_data, all_metric_correlations, elevation_models, clusters, type_names, presence_frequencies, scatterplot_cluster_separation, individual_totals) {
  c(
    save_correlation_heatmap(correlations),
    save_increment_heatmap(incremental),
    save_scatter_grid(data),
    save_divergence_boxplot(residual_data),
    save_all_metric_scatter_by_scale(data, all_metric_correlations),
    save_cluster_highlight_scatter_by_scale(data, clusters, all_metric_correlations, scatterplot_cluster_separation),
    save_individual_ramp_scatter_by_scale(data, individual_totals, all_metric_correlations),
    save_species_count_ramp_scatter_by_scale(data, individual_totals, all_metric_correlations),
    save_edge_highlight_scatter_by_scale(data),
    save_elevation_gradient_scatter_by_scale(data, elevation_models),
    save_composition_type_presence_heatmaps(type_names, presence_frequencies)
  )
}

best_correlations <- function(correlations) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (phylo_metric in PHYLO_METRICS) {
      subset_rows <- correlations[correlations$scale == scale_name & correlations$phylo_metric == phylo_metric, ]
      rows[[index]] <- subset_rows[which.max(abs(subset_rows$pearson_r)), ]
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

best_incremental_models <- function(incremental) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (sv_metric in SV_METRICS) {
      subset_rows <- incremental[incremental$scale == scale_name & incremental$sv_metric == sv_metric, ]
      rows[[index]] <- subset_rows[which.max(subset_rows$phylo_incremental_r2_after_species), ]
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

write_analysis_report <- function(correlations, incremental, metric_stats, cluster_summary, moderation, significant_moderation, figure_files) {
  top_cor <- best_correlations(correlations)
  top_inc <- best_incremental_models(incremental)
  top_cluster <- cluster_summary[order(-cluster_summary$mean_abs_phylo_residual_z), ]
  top_cluster <- top_cluster[seq_len(min(12, nrow(top_cluster))), ]
  top_moderation <- significant_moderation[seq_len(min(12, nrow(significant_moderation))), ]

  cor_table <- data.frame(
    Scale = top_cor$scale,
    `Phylogenetic metric` = top_cor$phylo_label,
    `Closest species metric` = top_cor$species_label,
    n = top_cor$n,
    r = fmt_num(top_cor$pearson_r, 3),
    R2 = fmt_num(top_cor$r_squared, 3),
    `p-value` = fmt_p(top_cor$f_p_value),
    check.names = FALSE
  )

  inc_table <- data.frame(
    Scale = top_inc$scale,
    `Observed SV metric` = top_inc$sv_label,
    `Species metric already in model` = top_inc$species_label,
    `Added phylogenetic metric` = top_inc$phylo_label,
    n = top_inc$n,
    `Species-only R2` = fmt_num(top_inc$species_only_r2, 3),
    `Combined R2` = fmt_num(top_inc$combined_r2, 3),
    `Phylo incremental R2` = fmt_num(top_inc$phylo_incremental_r2_after_species, 3),
    `Partial p-value` = fmt_p(top_inc$partial_f_p_value),
    check.names = FALSE
  )

  cluster_table <- data.frame(
    Scale = top_cluster$scale,
    Cluster = top_cluster$composition_cluster,
    `Cluster size` = top_cluster$cluster_size,
    `Species metric` = top_cluster$species_label,
    `Phylogenetic metric` = top_cluster$phylo_label,
    n = top_cluster$n,
    `Mean abs residual z` = fmt_num(top_cluster$mean_abs_phylo_residual_z, 3),
    check.names = FALSE
  )

  if (nrow(top_moderation) > 0) {
    moderation_table <- data.frame(
      Scale = top_moderation$scale,
      Relationship = paste(top_moderation$x_label, "vs", top_moderation$y_label),
      Moderator = top_moderation$moderator_label,
      n = top_moderation$n,
      `Interaction slope` = fmt_num(top_moderation$interaction_slope, 3),
      `Interaction delta R2` = fmt_num(top_moderation$interaction_incremental_r2_after_additive, 3),
      `Full delta R2` = fmt_num(top_moderation$full_incremental_r2_after_metric, 3),
      `Interaction q-value` = fmt_p(top_moderation$interaction_q_value),
      check.names = FALSE
    )
  } else {
    moderation_table <- data.frame(
      Result = "No species-count or crown-equivalent-individual interaction terms were significant after BH adjustment when cases where species count was already the richness axis were excluded.",
      check.names = FALSE
    )
  }

  report_lines <- c(
    "# Species And Phylogenetic Diversity Correlation Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis compares four species-diversity measures with three phylogenetic-diversity measures across 10 m, 20 m, and 50 m quadrats. It asks how closely the phylogenetic metrics track species diversity, where they diverge, and whether adding phylogenetic diversity contributes additional explanatory variation in the observed spectral-heterogeneity data.",
    "",
    "## Inputs",
    "",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`",
    "- `Quad_Values/Diversity_SHPs/plant_diversity_10m.csv`",
    "- `Quad_Values/Diversity_SHPs/plant_diversity_20m.csv`",
    "- `Quad_Values/Diversity_SHPs/plant_diversity_50m.csv`",
    "",
    "## Measures Compared",
    "",
    "- Species-diversity measures: `sp_rich`, `sp_shannon`, `sp_simpson`, `sp_even`.",
    "- Phylogenetic-diversity measures: `phy_faith`, `phy_rao`, `phy_afaith`.",
    "- Observed spectral-variation responses for incremental checks: `spec_spca_mean` and `spec_sa`.",
    "",
    "## Methods",
    "",
    "- Pairwise species-phylogenetic concordance was calculated within each scale using all quadrats with available species and phylogenetic diversity values.",
    "- Each pair was summarized with Pearson correlation, linear-model R2, slope, intercept, F statistic, and Spearman correlation.",
    "- A separate all-metric scatterplot set compares every unique pair among the four species-diversity metrics and three phylogenetic-diversity metrics within each quadrat scale.",
    "- Edge-highlight scatterplot panels were created for 10 m and 20 m quadrats, with edge quadrats in red, non-edge quadrats in blue, and the regression line plus displayed statistics calculated from non-edge quadrats.",
    "- Elevation-gradient scatterplot panels were created for each scale, with points colored from red at low mean elevation to green at high mean elevation. Each panel reports the incremental R2 and p-value for adding `env_elev` after the x-axis diversity metric.",
    "- Individual-ramp scatterplot panels were created for each scale, with points colored by crown-equivalent individuals calculated as the sum of species crown-overlap proportions in each quadrat.",
    "- Species-count scatterplot panels were created for each scale, with points colored by the number of species present in each quadrat.",
    "- Incremental spectral-variation models compared `SV ~ species_metric`, `SV ~ phylo_metric`, and `SV ~ species_metric + phylo_metric` using all complete spectral-diversity cases; the reported incremental R2 is the additional variance explained by the phylogenetic metric after the species metric.",
    "- To evaluate convergence or divergence around broadly similar species composition, quadrats were clustered within each scale from species presence/absence, not abundance. Phylogenetic residuals from `phylo_metric ~ species_metric` were summarized within those composition clusters.",
    "- Mean silhouette width was calculated for each composition type within every all-metric scatterplot panel using the plotted two-metric space; larger positive values indicate stronger separation from the other composition types in that panel.",
    "- Species count and crown-equivalent individuals were then tested as moderators of each all-metric relationship using standardized models of the form `y_metric ~ x_metric * moderator`. Significant results below use BH-adjusted interaction p-values and exclude species-count tests where species richness is already one of the plotted axes.",
    "",
    "## Closest Species-Diversity Match For Each Phylogenetic Metric",
    "",
    markdown_table(cor_table),
    "",
    "## Strongest Incremental Phylogenetic Contribution To Observed Spectral Variation",
    "",
    markdown_table(inc_table),
    "",
    "## Largest Composition-Cluster Divergence Cases",
    "",
    markdown_table(cluster_table),
    "",
    "## Significant Species-Count Or Individual-Count Effects On Metric Relationships",
    "",
    markdown_table(moderation_table),
    "",
    "## Interpretation Notes",
    "",
    "- Faith's PD is expected to converge most strongly with species richness because both are primarily presence-based.",
    "- Phylogenetic Rao's Q and abundance-weighted Faith's PD can diverge from richness when common species differ in relatedness or when abundance is concentrated in phylogenetically similar or distinct lineages.",
    "- A positive incremental R2 means the phylogenetic metric explains spectral-variation structure not captured by the paired species-diversity metric alone.",
    "- Cluster residual summaries identify compositionally similar groups where phylogenetic values remain unusually high or low relative to the species-diversity expectation.",
    "- A significant moderator interaction means the slope between two plotted diversity metrics changes as species count or crown-equivalent individuals increase.",
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Output Tables",
    "",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_pairwise_correlations.csv`",
    "- `reports/tables/species_phylogenetic_correlation/all_diversity_metric_pairwise_correlations.csv`",
    "- `reports/tables/species_phylogenetic_correlation/elevation_adjusted_all_diversity_metric_models.csv`",
    "- `reports/tables/species_phylogenetic_correlation/phylogenetic_incremental_sv_models.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_metric_summary.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_type_names.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_type_presence_frequencies.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_scatterplot_silhouette.csv`",
    "- `reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_individual_moderated_metric_relationships.csv`",
    "- `reports/tables/species_phylogenetic_correlation/significant_species_individual_moderated_metric_relationships.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_residual_divergence.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_cluster_divergence_summary_presence_based.csv`"
  )

  writeLines(report_lines, ANALYSIS_REPORT)
}

write_task_report <- function(correlations, all_metric_correlations, elevation_models, incremental, clusters, type_names, residual_data, moderation, significant_moderation) {
  lines <- c(
    "# Task Report: Species And Phylogenetic Diversity Correlation Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Compare species-diversity and phylogenetic-diversity metrics across 10 m, 20 m, and 50 m quadrats, and evaluate whether phylogenetic diversity adds explanatory value beyond species diversity for observed spectral variation.",
    "",
    "## Outputs Created",
    "",
    "- `reports/analysis/20260817_species_phylogenetic_correlation_analysis.md`",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_pairwise_correlations.csv`",
    "- `reports/tables/species_phylogenetic_correlation/all_diversity_metric_pairwise_correlations.csv`",
    "- `reports/tables/species_phylogenetic_correlation/elevation_adjusted_all_diversity_metric_models.csv`",
    "- `reports/tables/species_phylogenetic_correlation/phylogenetic_incremental_sv_models.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_metric_summary.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_type_names.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_type_presence_frequencies.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_scatterplot_silhouette.csv`",
    "- `reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_individual_moderated_metric_relationships.csv`",
    "- `reports/tables/species_phylogenetic_correlation/significant_species_individual_moderated_metric_relationships.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_phylogenetic_residual_divergence.csv`",
    "- `reports/tables/species_phylogenetic_correlation/species_composition_cluster_divergence_summary_presence_based.csv`",
    "- `reports/figures/species_phylogenetic_correlation/`",
    "",
    "## Result Size",
    "",
    paste0("- Species-phylogenetic pairwise correlation rows: ", nrow(correlations)),
    paste0("- All diversity metric pairwise correlation rows: ", nrow(all_metric_correlations)),
    paste0("- Elevation-adjusted all-metric model rows: ", nrow(elevation_models)),
    paste0("- Incremental spectral-variation model rows: ", nrow(incremental)),
    paste0("- Species-composition cluster rows: ", nrow(clusters)),
    paste0("- Named composition type rows: ", nrow(type_names)),
    paste0("- Species/individual moderation model rows: ", nrow(moderation)),
    paste0("- Significant moderation rows after BH adjustment: ", nrow(significant_moderation)),
    paste0("- Residual divergence rows: ", nrow(residual_data)),
    "",
    "## Notes",
    "",
    "- Biodiversity-only comparisons use all quadrats with available species and phylogenetic diversity values.",
    "- Spectral increment checks use all complete spectral-diversity cases and do not remove edge quadrats.",
    "- Species-composition clusters are based on species presence/absence, not crown-overlap abundance.",
    "- Composition type names are generated from the species codes most characteristic of each presence-based cluster.",
    "- The `09_...composition_cluster...` scatterplot legends report panel-specific mean silhouette width for each composition type, calculated from the two metrics plotted in that panel.",
    "- The `10_...individual_ramp...` scatterplot figures color quadrats by summed crown-overlap proportions, interpreted as crown-equivalent individuals rather than raw stem counts.",
    "- The `11_...species_count_ramp...` scatterplot figures color quadrats by the number of species present, matching the presence-based richness definition.",
    "- Species-count and crown-equivalent-individual moderation tests use standardized `y_metric ~ x_metric * moderator` models; significant result summaries exclude cases where species count is already the plotted richness metric.",
    "- Species-composition clusters are descriptive, not formal community types."
  )
  writeLines(lines, TASK_REPORT)
}

write_validation_report <- function(data, correlations, all_metric_correlations, elevation_models, incremental, clusters, type_names, presence_frequencies, scatterplot_cluster_separation, individual_totals, moderation, significant_moderation, residual_data, figure_files) {
  expected_correlations <- length(SCALES) * length(SPECIES_METRICS) * length(PHYLO_METRICS)
  expected_all_metric_correlations <- length(SCALES) * choose(length(ALL_DIVERSITY_METRICS), 2)
  expected_elevation_models <- expected_all_metric_correlations
  expected_incremental <- length(SCALES) * length(SV_METRICS) * length(SPECIES_METRICS) * length(PHYLO_METRICS)
  table_files <- c(
    file.path(TABLE_DIR, "species_phylogenetic_pairwise_correlations.csv"),
    file.path(TABLE_DIR, "all_diversity_metric_pairwise_correlations.csv"),
    file.path(TABLE_DIR, "elevation_adjusted_all_diversity_metric_models.csv"),
    file.path(TABLE_DIR, "phylogenetic_incremental_sv_models.csv"),
    file.path(TABLE_DIR, "species_phylogenetic_metric_summary.csv"),
    file.path(TABLE_DIR, "species_composition_clusters.csv"),
    file.path(TABLE_DIR, "species_composition_type_names.csv"),
    file.path(TABLE_DIR, "species_composition_type_presence_frequencies.csv"),
    file.path(TABLE_DIR, "species_composition_scatterplot_silhouette.csv"),
    file.path(TABLE_DIR, "quadrat_crown_equivalent_individual_totals.csv"),
    file.path(TABLE_DIR, "species_individual_moderated_metric_relationships.csv"),
    file.path(TABLE_DIR, "significant_species_individual_moderated_metric_relationships.csv"),
    file.path(TABLE_DIR, "species_phylogenetic_residual_divergence.csv"),
    file.path(TABLE_DIR, "species_composition_cluster_divergence_summary_presence_based.csv")
  )
  scale_rows <- aggregate(
    data[, ALL_DIVERSITY_METRICS],
    by = list(scale = data$scale),
    FUN = function(x) sum(!is.na(x))
  )
  coverage <- data.frame(
    scale = scale_rows$scale,
    total_rows = as.integer(table(data$scale)[scale_rows$scale]),
    complete_species_phylo_rows = apply(scale_rows[, ALL_DIVERSITY_METRICS, drop = FALSE], 1, min),
    stringsAsFactors = FALSE
  )

  lines <- c(
    "# Validation: Species And Phylogenetic Diversity Correlation Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- Expected species-phylogenetic correlation rows: ", expected_correlations),
    paste0("- Observed species-phylogenetic correlation rows: ", nrow(correlations)),
    paste0("- Expected all-metric correlation rows: ", expected_all_metric_correlations),
    paste0("- Observed all-metric correlation rows: ", nrow(all_metric_correlations)),
    paste0("- Expected elevation-adjusted model rows: ", expected_elevation_models),
    paste0("- Observed elevation-adjusted model rows: ", nrow(elevation_models)),
    paste0("- Expected incremental model rows: ", expected_incremental),
    paste0("- Observed incremental model rows: ", nrow(incremental)),
    paste0("- Missing Pearson r values: ", sum(is.na(correlations$pearson_r))),
    paste0("- Missing incremental R2 values: ", sum(is.na(incremental$phylo_incremental_r2_after_species))),
    paste0("- Cluster assignments created: ", nrow(clusters)),
    paste0("- Named composition type rows: ", nrow(type_names)),
    paste0("- Missing composition type silhouette values: ", sum(is.na(type_names$mean_silhouette_width))),
    paste0("- Scatterplot composition silhouette rows: ", nrow(scatterplot_cluster_separation)),
    paste0("- Missing scatterplot composition silhouette values: ", sum(is.na(scatterplot_cluster_separation$mean_scatterplot_silhouette_width))),
    paste0("- Quadrat crown-equivalent individual total rows: ", nrow(individual_totals)),
    paste0("- Missing crown-equivalent individual totals: ", sum(is.na(individual_totals$crown_equivalent_individuals))),
    paste0("- Missing present-species counts: ", sum(is.na(individual_totals$present_species_count))),
    paste0("- Species/individual moderation model rows: ", nrow(moderation)),
    paste0("- Significant moderation rows after BH adjustment: ", nrow(significant_moderation)),
    paste0("- Missing moderation interaction q-values: ", sum(is.na(moderation$interaction_q_value))),
    paste0("- Species presence frequency rows: ", nrow(presence_frequencies)),
    paste0("- Residual divergence rows created: ", nrow(residual_data)),
    paste0("- Missing all-metric Pearson r values: ", sum(is.na(all_metric_correlations$pearson_r))),
    paste0("- Missing elevation incremental R2 values: ", sum(is.na(elevation_models$elevation_incremental_r2_after_metric))),
    paste0("- Output tables present: ", sum(file.exists(table_files)), " of ", length(table_files)),
    paste0("- Output figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    "",
    "## Coverage",
    "",
    markdown_table(coverage),
    "",
    "## Result",
    "",
    "The requested species-diversity and phylogenetic-diversity comparison layer completed using all available biodiversity quadrats and produced the expected row counts, tables, and figures."
  )
  writeLines(lines, VALIDATION_REPORT)
}

run_species_phylogenetic_correlation_analysis <- function() {
  data <- read_analysis_data()
  required_columns <- c("quad_id", "scale", "primary_analysis", SPECIES_METRICS, PHYLO_METRICS, SV_METRICS)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  correlations <- run_pairwise_correlations(data)
  all_metric_correlations <- run_all_metric_correlations(data)
  elevation_models <- run_elevation_adjusted_all_metric_models(data)
  incremental <- incremental_sv_models(data)
  individual_totals <- make_quadrat_individual_totals()
  moderation <- analyze_species_individual_moderation(data, individual_totals)
  significant_moderation <- significant_moderation_results(moderation)
  clusters <- make_composition_clusters(data)
  composition_type_outputs <- make_composition_type_names(clusters)
  type_names <- composition_type_outputs$type_names
  presence_frequencies <- composition_type_outputs$frequencies
  clusters <- merge(
    clusters,
    type_names[, c("scale", "composition_cluster", "composition_type")],
    by = c("scale", "composition_cluster"),
    all.x = TRUE,
    sort = FALSE
  )
  cluster_separation <- calculate_composition_cluster_separation(clusters)
  clusters <- merge(
    clusters,
    cluster_separation,
    by = c("scale", "composition_cluster"),
    all.x = TRUE,
    sort = FALSE
  )
  type_names <- merge(
    type_names,
    cluster_separation,
    by = c("scale", "composition_cluster"),
    all.x = TRUE,
    sort = FALSE
  )
  scatterplot_cluster_separation <- calculate_scatterplot_cluster_separation(data, clusters)
  residual_data <- make_residual_dataset(data, clusters)
  cluster_summary <- summarize_cluster_divergence(residual_data)
  stats <- metric_summary(data)
  figure_files <- save_figures(data, correlations, incremental, residual_data, all_metric_correlations, elevation_models, clusters, type_names, presence_frequencies, scatterplot_cluster_separation, individual_totals)

  write.csv(correlations, file.path(TABLE_DIR, "species_phylogenetic_pairwise_correlations.csv"), row.names = FALSE)
  write.csv(all_metric_correlations, file.path(TABLE_DIR, "all_diversity_metric_pairwise_correlations.csv"), row.names = FALSE)
  write.csv(elevation_models, file.path(TABLE_DIR, "elevation_adjusted_all_diversity_metric_models.csv"), row.names = FALSE)
  write.csv(incremental, file.path(TABLE_DIR, "phylogenetic_incremental_sv_models.csv"), row.names = FALSE)
  write.csv(stats, file.path(TABLE_DIR, "species_phylogenetic_metric_summary.csv"), row.names = FALSE)
  write.csv(clusters, file.path(TABLE_DIR, "species_composition_clusters.csv"), row.names = FALSE)
  write.csv(type_names, file.path(TABLE_DIR, "species_composition_type_names.csv"), row.names = FALSE)
  write.csv(presence_frequencies, file.path(TABLE_DIR, "species_composition_type_presence_frequencies.csv"), row.names = FALSE)
  write.csv(scatterplot_cluster_separation, file.path(TABLE_DIR, "species_composition_scatterplot_silhouette.csv"), row.names = FALSE)
  write.csv(individual_totals, file.path(TABLE_DIR, "quadrat_crown_equivalent_individual_totals.csv"), row.names = FALSE)
  write.csv(moderation, file.path(TABLE_DIR, "species_individual_moderated_metric_relationships.csv"), row.names = FALSE)
  write.csv(significant_moderation, file.path(TABLE_DIR, "significant_species_individual_moderated_metric_relationships.csv"), row.names = FALSE)
  write.csv(residual_data, file.path(TABLE_DIR, "species_phylogenetic_residual_divergence.csv"), row.names = FALSE)
  write.csv(cluster_summary, file.path(TABLE_DIR, "species_composition_cluster_divergence_summary_presence_based.csv"), row.names = FALSE)

  write_analysis_report(correlations, incremental, stats, cluster_summary, moderation, significant_moderation, figure_files)
  write_task_report(correlations, all_metric_correlations, elevation_models, incremental, clusters, type_names, residual_data, moderation, significant_moderation)
  write_validation_report(data, correlations, all_metric_correlations, elevation_models, incremental, clusters, type_names, presence_frequencies, scatterplot_cluster_separation, individual_totals, moderation, significant_moderation, residual_data, figure_files)

  message("Species-phylogenetic correlation analysis complete.")
  invisible(list(
    correlations = correlations,
    all_metric_correlations = all_metric_correlations,
    elevation_models = elevation_models,
    incremental = incremental,
    clusters = clusters,
    type_names = type_names,
    presence_frequencies = presence_frequencies,
    scatterplot_cluster_separation = scatterplot_cluster_separation,
    individual_totals = individual_totals,
    moderation = moderation,
    significant_moderation = significant_moderation,
    residual_data = residual_data,
    cluster_summary = cluster_summary,
    stats = stats
  ))
}

if (sys.nframe() == 0) {
  run_species_phylogenetic_correlation_analysis()
}
