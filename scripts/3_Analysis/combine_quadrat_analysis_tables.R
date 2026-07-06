USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c("sf", "dplyr", "readr", "tibble")
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
  library(sf)
  library(dplyr)
  library(readr)
  library(tibble)
})

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
SPECTRAL_DIR <- "Quad_Values/Spectral_diversitySHPs"
DIVERSITY_DIR <- "Quad_Values/Diversity_SHPs"
ENVIRO_DIR <- "Quad_Values/Enviro_SHPs"
TAXA_CSV <- "51sp_taxanomy.csv"
VARIABLE_GUIDE <- "reports/combined_quadrat_variable_guide.md"
VALIDATION_REPORT <- "reports/validation/20260624_combined_quadrat_analysis_tables_validation.md"
TASK_REPORT <- "reports/tasks/20260624_combined_quadrat_analysis_tables.md"

SCALE_CONFIG <- list(
  "10m" = list(
    scale = "10m",
    plant_csv = file.path(DIVERSITY_DIR, "plant_diversity_10m.csv"),
    plant_shp = file.path(DIVERSITY_DIR, "plant_diversity_10m.shp"),
    spectral_csv = file.path(SPECTRAL_DIR, "spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv"),
    enviro_csv = file.path(ENVIRO_DIR, "enviro_variables_10m.csv"),
    output_csv = "quadrat_analysis_10m.csv"
  ),
  "20m" = list(
    scale = "20m",
    plant_csv = file.path(DIVERSITY_DIR, "plant_diversity_20m.csv"),
    plant_shp = file.path(DIVERSITY_DIR, "plant_diversity_20m.shp"),
    spectral_csv = file.path(SPECTRAL_DIR, "spectral_heterogeneity_20m_smooth_masked_5nm_summary.csv"),
    enviro_csv = file.path(ENVIRO_DIR, "enviro_variables_20m.csv"),
    output_csv = "quadrat_analysis_20m.csv"
  ),
  "50m" = list(
    scale = "50m",
    plant_csv = file.path(DIVERSITY_DIR, "plant_diversity_50m.csv"),
    plant_shp = file.path(DIVERSITY_DIR, "plant_diversity_50m.shp"),
    spectral_csv = file.path(SPECTRAL_DIR, "spectral_heterogeneity_50m_smooth_masked_5nm_summary.csv"),
    enviro_csv = file.path(ENVIRO_DIR, "enviro_variables_50m.csv"),
    output_csv = "quadrat_analysis_50m.csv"
  )
)

SPECIES_METRIC_MAP <- c(
  richness = "sp_rich",
  shannon = "sp_shannon",
  simpson = "sp_simpson",
  evenness = "sp_even",
  faith_pd = "phy_faith",
  rao_pd = "phy_rao",
  afaith_pd = "phy_afaith"
)

SPECTRAL_MAP <- c(
  sa_entropy = "spec_sa",
  pca_euclidean_mean = "spec_pca_mean",
  pca_euclidean_median = "spec_pca_med",
  pca_euclidean_sd = "spec_pca_sd",
  rao_q_pca = "spec_rao_q",
  alpha_hull_area = "spec_alpha",
  pca_convex_hull_area = "spec_convex",
  pca_hull_volume_3d = "spec_hull3d_v",
  pca_hull_area_3d = "spec_hull3d_a"
)

ENVIRO_MAP <- c(
  elev_mean = "env_elev",
  tri5_mean = "env_tri5",
  tri11_mean = "env_tri11",
  tri21_mean = "env_tri21"
)

find_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      file.exists(file.path(current_dir, TAXA_CSV)) &&
        dir.exists(file.path(current_dir, SPECTRAL_DIR)) &&
        dir.exists(file.path(current_dir, DIVERSITY_DIR)) &&
        dir.exists(file.path(current_dir, ENVIRO_DIR))
    ) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find the Spectral_Diversity project root.", call. = FALSE)
    }
    current_dir <- parent_dir
  }
}

rename_selected <- function(data, rename_map) {
  available_original <- intersect(names(rename_map), names(data))
  data %>%
    dplyr::select(all_of(available_original)) %>%
    dplyr::rename(!!!setNames(available_original, rename_map[available_original]))
}

read_centers <- function(config) {
  quads <- st_read(config$plant_shp, quiet = TRUE)
  centers <- suppressWarnings(st_centroid(quads))
  coords <- st_coordinates(centers)

  tibble(
    quad_id = as.character(quads$quad_id),
    scale = config$scale,
    center_x = coords[, "X"],
    center_y = coords[, "Y"]
  )
}

read_plant_values <- function(config, taxa_table) {
  plant <- readr::read_csv(config$plant_csv, show_col_types = FALSE) %>%
    mutate(
      quad_id = as.character(.data$quad_id),
      scale = as.character(.data$scale)
    )

  metric_values <- rename_selected(plant, SPECIES_METRIC_MAP)

  bind_cols(
    plant %>% dplyr::select(quad_id, scale),
    metric_values
  )
}

read_spectral_values <- function(config) {
  spectral <- readr::read_csv(config$spectral_csv, show_col_types = FALSE) %>%
    mutate(
      quad_id = as.character(.data$quad_id),
      scale = config$scale
    )

  bind_cols(
    spectral %>% dplyr::select(quad_id, scale),
    rename_selected(spectral, SPECTRAL_MAP)
  )
}

read_enviro_values <- function(config) {
  enviro <- readr::read_csv(config$enviro_csv, show_col_types = FALSE) %>%
    mutate(
      quad_id = as.character(.data$quad_id),
      scale = as.character(.data$scale)
    )

  bind_cols(
    enviro %>% dplyr::select(quad_id, scale),
    rename_selected(enviro, ENVIRO_MAP)
  )
}

combine_scale <- function(config, taxa_table) {
  centers <- read_centers(config)
  plant <- read_plant_values(config, taxa_table)
  spectral <- read_spectral_values(config)
  enviro <- read_enviro_values(config)

  combined <- centers %>%
    left_join(plant, by = c("quad_id", "scale")) %>%
    left_join(spectral, by = c("quad_id", "scale")) %>%
    left_join(enviro, by = c("quad_id", "scale")) %>%
    arrange(quad_id)

  readr::write_csv(combined, config$output_csv)
  combined
}

make_variable_inventory <- function(taxa_table) {
  tibble(
    column = c(
      "quad_id", "scale", "center_x", "center_y",
      unname(SPECTRAL_MAP),
      unname(SPECIES_METRIC_MAP),
      unname(ENVIRO_MAP)
    ),
    group = c(
      rep("Identifier and coordinates", 4),
      rep("Spectral heterogeneity", length(SPECTRAL_MAP)),
      rep("Species and phylogenetic diversity", length(SPECIES_METRIC_MAP)),
      rep("Environmental/topographic", length(ENVIRO_MAP))
    ),
    source_name = c(
      "quad_id", "scale", "centroid X", "centroid Y",
      names(SPECTRAL_MAP),
      names(SPECIES_METRIC_MAP),
      names(ENVIRO_MAP)
    ),
    definition = c(
      "Scale-aware quadrat identifier used for joins across spectral, diversity, and environmental products.",
      "Quadrat grain: 10m, 20m, or 50m.",
      "X coordinate of the quadrat polygon centroid in the plant-diversity shapefile CRS.",
      "Y coordinate of the quadrat polygon centroid in the plant-diversity shapefile CRS.",
      "Primary spectral angle entropy mean from sunlit, shadow-masked smoothed 5 nm spectra; for most quadrats this is the mean of 70 bootstrap iterations using up to 5,000 retained pixels, while the small exact subset uses all retained pixel pairs.",
      "Mean Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Median Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Standard deviation of Euclidean distances of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Rao's Q spectral heterogeneity metric using equal pixel weights and squared Euclidean distance in global PC1-PC3 space.",
      "Alpha-hull area in global PC1-PC2 spectral space.",
      "Convex-hull area in global PC1-PC2 spectral space, used as a supplemental hull summary.",
      "Convex-hull volume in global PC1-PC3 spectral space, used as a supplemental three-axis hull summary.",
      "Convex-hull surface area in global PC1-PC3 spectral space, used as a supplemental three-axis hull summary.",
      "Number of species with positive crown-overlap values in the quadrat.",
      "Shannon species diversity calculated from positive species crown-overlap proportions.",
      "Simpson diversity calculated from positive species crown-overlap proportions as 1 minus the sum of squared proportional abundances.",
      "Species evenness calculated as Shannon diversity divided by log richness when richness is greater than 1.",
      "Faith's phylogenetic diversity calculated as summed branch length for species present in the quadrat.",
      "Rao phylogenetic diversity calculated from species crown-overlap weights and phylogenetic distances.",
      "Abundance-weighted Faith's phylogenetic diversity calculated by weighting phylogenetic branch lengths by species crown-overlap values.",
      "Mean DTM elevation within the quadrat, in meters.",
      "Mean Riley topographic roughness index from a 5x5 DTM moving window.",
      "Mean Riley topographic roughness index from an 11x11 DTM moving window.",
      "Mean Riley topographic roughness index from a 21x21 DTM moving window."
    )
  )
}

write_variable_guide <- function(inventory, center_crs) {
  dir.create(dirname(VARIABLE_GUIDE), recursive = TRUE, showWarnings = FALSE)

  lines <- c(
    "# Combined Quadrat Analysis Tables Variable Guide",
    "",
    "Last updated: 2026-06-24",
    "",
    "This guide documents the shortened columns in the root-level combined CSVs:",
    "",
    "- `quadrat_analysis_10m.csv`",
    "- `quadrat_analysis_20m.csv`",
    "- `quadrat_analysis_50m.csv`",
    "",
    paste0("Center coordinates are polygon centroids from the plant-diversity shapefiles. CRS: `", center_crs, "`."),
    "",
    "Pixel-count, pair-count, bootstrap replicate, method, exclusion, geometry metadata, and per-species composition fields are intentionally excluded from the combined CSVs.",
    "",
    "## Column Definitions",
    "",
    "| Column | Group | Source Name | Definition |",
    "|---|---|---|---|"
  )

  rows <- apply(inventory, 1, function(row) {
    paste0(
      "| `", row[["column"]], "` | ",
      row[["group"]], " | `",
      row[["source_name"]], "` | ",
      row[["definition"]], " |"
    )
  })

  writeLines(c(lines, rows), VARIABLE_GUIDE)
}

validate_outputs <- function(results) {
  validation <- lapply(names(results), function(scale_name) {
    data <- results[[scale_name]]
    tibble(
      scale = scale_name,
      output_csv = SCALE_CONFIG[[scale_name]]$output_csv,
      rows = nrow(data),
      columns = ncol(data),
      duplicate_quad_ids = sum(duplicated(data$quad_id)),
      missing_center_x = sum(is.na(data$center_x)),
      missing_center_y = sum(is.na(data$center_y)),
      missing_spec_sa = sum(is.na(data$spec_sa)),
      missing_spec_pca_mean = sum(is.na(data$spec_pca_mean)),
      missing_env_elev = sum(is.na(data$env_elev))
    )
  }) %>%
    bind_rows()

  forbidden_patterns <- c(
    "pixel", "pixels", "pairs", "boot", "method", "excluded",
    "geometry", "hull_n_points", "alpha_hull_alpha"
  )
  forbidden_columns <- lapply(results, function(data) {
    names(data)[grepl(paste(forbidden_patterns, collapse = "|"), names(data), ignore.case = TRUE)]
  })
  species_code_columns <- lapply(results, function(data) {
    allowed_sp_summary <- c("sp_rich", "sp_shannon", "sp_simpson", "sp_even")
    setdiff(names(data)[grepl("^sp_", names(data))], allowed_sp_summary)
  })

  dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    "# Combined Quadrat Analysis Tables Validation",
    "",
    "Last updated: 2026-06-24",
    "",
    "## Output Summary",
    "",
    "| Scale | Output CSV | Rows | Columns | Duplicate quad IDs | Missing center X | Missing center Y | Missing spectral SA | Missing spectral PCA mean | Missing elevation |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    apply(validation, 1, function(row) {
      paste0(
        "| ", row[["scale"]],
        " | `", row[["output_csv"]],
        "` | ", row[["rows"]],
        " | ", row[["columns"]],
        " | ", row[["duplicate_quad_ids"]],
        " | ", row[["missing_center_x"]],
        " | ", row[["missing_center_y"]],
        " | ", row[["missing_spec_sa"]],
        " | ", row[["missing_spec_pca_mean"]],
        " | ", row[["missing_env_elev"]],
        " |"
      )
    }),
    "",
    "## Metadata Column Check",
    ""
  )

  for (scale_name in names(forbidden_columns)) {
    matches <- forbidden_columns[[scale_name]]
    message <- if (length(matches) == 0) {
      "none detected"
    } else {
      paste(matches, collapse = ", ")
    }
    lines <- c(lines, paste0("- ", scale_name, ": ", message))
  }

  lines <- c(lines, "", "## Species Composition Column Check", "")
  for (scale_name in names(species_code_columns)) {
    matches <- species_code_columns[[scale_name]]
    message <- if (length(matches) == 0) {
      "none detected"
    } else {
      paste(matches, collapse = ", ")
    }
    lines <- c(lines, paste0("- ", scale_name, ": ", message))
  }

  writeLines(lines, VALIDATION_REPORT)
  validation
}

write_task_report <- function(validation) {
  dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)

  lines <- c(
    "# Combined Quadrat Analysis Tables",
    "",
    "Last updated: 2026-06-24",
    "",
    "## Task",
    "",
    "Combined current spectral heterogeneity, species and phylogenetic diversity, environmental/topographic values, and quadrat center coordinates into one root-level CSV for each quadrat scale.",
    "",
    "## Outputs",
    "",
    "- `quadrat_analysis_10m.csv`",
    "- `quadrat_analysis_20m.csv`",
    "- `quadrat_analysis_50m.csv`",
    "- `combined_quadrat_variable_guide.docx`",
    "- `reports/combined_quadrat_variable_guide.md`",
    "- `reports/validation/20260624_combined_quadrat_analysis_tables_validation.md`",
    "",
    "## Notes",
    "",
    "- Pixel-count, pair-count, bootstrap replicate, method, exclusion, and geometry metadata fields were excluded.",
    "- Per-species composition columns were excluded; original composition values remain available in `Quad_Values/Diversity_SHPs/plant_diversity_*m.csv` if needed later.",
    "- Species diversity summaries use the `sp_` prefix.",
    "- Spectral columns use the `spec_` prefix.",
    "- Environmental/topographic columns use the `env_` prefix.",
    "- Spectral gaps remain `NA` where current raster-derived spectral summaries were not available or where PCA-dependent values were excluded upstream.",
    "",
    "## Validation Summary",
    "",
    "| Scale | Rows | Columns | Duplicate quad IDs | Missing center X | Missing center Y | Missing spectral SA | Missing spectral PCA mean | Missing elevation |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    apply(validation, 1, function(row) {
      paste0(
        "| ", row[["scale"]],
        " | ", row[["rows"]],
        " | ", row[["columns"]],
        " | ", row[["duplicate_quad_ids"]],
        " | ", row[["missing_center_x"]],
        " | ", row[["missing_center_y"]],
        " | ", row[["missing_spec_sa"]],
        " | ", row[["missing_spec_pca_mean"]],
        " | ", row[["missing_env_elev"]],
        " |"
      )
    })
  )

  writeLines(lines, TASK_REPORT)
}

run_combined_quadrat_analysis_tables <- function() {
  project_root <- find_project_root(PROJECT_DIR)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  taxa_table <- read.csv(TAXA_CSV, stringsAsFactors = FALSE)
  results <- lapply(SCALE_CONFIG, combine_scale, taxa_table = taxa_table)
  inventory <- make_variable_inventory(taxa_table)
  center_crs <- st_crs(st_read(SCALE_CONFIG[["10m"]]$plant_shp, quiet = TRUE))$input

  write_variable_guide(inventory, center_crs)
  validation <- validate_outputs(results)
  write_task_report(validation)

  print(validation)
  invisible(list(results = results, validation = validation))
}

if (sys.nframe() == 0) {
  run_combined_quadrat_analysis_tables()
}
