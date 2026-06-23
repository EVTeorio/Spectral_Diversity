USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

required_packages <- c(
  "sf", "dplyr", "tidyr", "readr", "ape", "picante", "V.PhyloMaker2"
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
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ape)
  library(picante)
  library(V.PhyloMaker2)
})

CRS_PROJ <- 26916
TREE_CSV <- "PR_tree_DL.csv"
TAXA_CSV <- "51sp_taxanomy.csv"
OUTPUT_DIR <- "Indices_SHPs/Diversity_SHPs"

SCALE_CONFIG <- list(
  "10m" = list(
    scale = "10m",
    quadrat_path = "Quad_Scale_SHPs/PR_10m.shp",
    id_col = "sub_id",
    output_stem = "plant_diversity_10m"
  ),
  "20m" = list(
    scale = "20m",
    quadrat_path = "Quad_Scale_SHPs/PR_20m.shp",
    id_col = "Name",
    output_stem = "plant_diversity_20m"
  ),
  "50m" = list(
    scale = "50m",
    quadrat_path = "Quad_Scale_SHPs/PR_50m.shp",
    id_col = "Name",
    output_stem = "plant_diversity_50m"
  )
)

find_project_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      file.exists(file.path(current_dir, TREE_CSV)) &&
        file.exists(file.path(current_dir, TAXA_CSV)) &&
        dir.exists(file.path(current_dir, "Quad_Scale_SHPs"))
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

read_quadrats <- function(config) {
  quads <- st_read(config$quadrat_path, quiet = TRUE) %>%
    dplyr::select(-matches("^Dscrptn")) %>%
    st_transform(CRS_PROJ) %>%
    mutate(
      quad_id = as.character(.data[[config$id_col]]),
      scale = config$scale
    )

  quads %>%
    relocate(quad_id, scale)
}

read_tree_crowns <- function(species_cols) {
  tree_df <- readr::read_csv(TREE_CSV, show_col_types = FALSE) %>%
    filter(
      .data[["DBH.2024"]] >= 200,
      .data[["crown.position"]] %in% c(3, 4, 5),
      .data[["cluster_status"]] %in% c("A", "R"),
      !is.na(.data[["UTMX_CURRENT"]]),
      !is.na(.data[["UTMY_CURRENT"]]),
      !is.na(.data[["cw_m_2025"]]),
      .data[["cw_m_2025"]] > 0,
      .data[["sp"]] %in% species_cols
    ) %>%
    mutate(
      crown_id = row_number(),
      radius_m = .data[["cw_m_2025"]] / 2
    )

  points_sf <- st_as_sf(
    tree_df,
    coords = c("UTMX_CURRENT", "UTMY_CURRENT"),
    crs = CRS_PROJ,
    remove = FALSE
  )

  st_buffer(points_sf, dist = points_sf$radius_m, endCapStyle = "ROUND") %>%
    dplyr::select(
      crown_id, sp, quadrat, tag, DBH.2024, crown.position,
      UTMX_CURRENT, UTMY_CURRENT, radius_m
    )
}

make_species_matrix <- function(quads, crowns, species_cols) {
  intersections <- suppressWarnings(st_intersection(crowns, quads)) %>%
    mutate(intersect_area = st_area(.))

  crown_areas <- crowns %>%
    mutate(total_area = st_area(.)) %>%
    st_drop_geometry() %>%
    dplyr::select(crown_id, total_area)

  quad_species <- intersections %>%
    left_join(crown_areas, by = "crown_id") %>%
    mutate(prop_area = as.numeric(intersect_area / total_area)) %>%
    st_drop_geometry() %>%
    group_by(quad_id, sp) %>%
    summarise(prop_sum = sum(prop_area, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(
      quad_id = quads$quad_id,
      sp = species_cols,
      fill = list(prop_sum = 0)
    )

  quad_species %>%
    pivot_wider(
      names_from = sp,
      values_from = prop_sum,
      values_fill = 0
    )
}

calculate_species_diversity <- function(abundance_matrix) {
  metrics <- t(apply(abundance_matrix, 1, function(abundances) {
    abundances <- as.numeric(abundances)
    abundances <- abundances[abundances > 0]

    if (length(abundances) == 0) {
      return(c(richness = 0, shannon = 0, simpson = 0, evenness = 0))
    }

    p <- abundances / sum(abundances)
    richness <- length(abundances)
    shannon <- -sum(p * log(p))
    simpson <- 1 - sum(p^2)
    evenness <- ifelse(richness > 1, shannon / log(richness), 0)

    c(
      richness = richness,
      shannon = shannon,
      simpson = simpson,
      evenness = evenness
    )
  }))

  as.data.frame(metrics)
}

build_phylogeny <- function(taxa_table) {
  phylo_input <- taxa_table %>%
    mutate(species = paste(genus, species, sep = " ")) %>%
    dplyr::select(species, genus, family) %>%
    filter(species != "Carya sp") %>%
    distinct(species, .keep_all = TRUE)

  phylo.maker(
    sp.list = phylo_input,
    tree = GBOTB.extended.TPL,
    nodes = nodes.info.1.TPL,
    scenarios = "S3"
  )$scenario.3
}

make_tip_abundances <- function(abundances, species_cols, sp_to_tip, all_tips) {
  names(abundances) <- species_cols
  present_species <- names(abundances)[abundances > 0]
  tip_labels <- sp_to_tip[present_species]
  valid <- !is.na(tip_labels) & tip_labels %in% all_tips

  if (!any(valid)) {
    return(numeric(0))
  }

  tip_abundances <- abundances[present_species][valid]
  names(tip_abundances) <- tip_labels[valid]
  tapply(tip_abundances, names(tip_abundances), sum)
}

calculate_faith_pd <- function(tip_abundances, phylo_tree) {
  if (length(tip_abundances) == 0) {
    return(0)
  }
  if (length(tip_abundances) == 1) {
    return(1)
  }

  sub_tree <- ape::keep.tip(phylo_tree, names(tip_abundances))
  sum(sub_tree$edge.length)
}

calculate_rao_pd <- function(tip_abundances, phylo_tree) {
  if (length(tip_abundances) <= 1) {
    return(0)
  }

  sub_tree <- ape::keep.tip(phylo_tree, names(tip_abundances))
  dist_mat <- cophenetic(sub_tree)

  p <- as.numeric(tip_abundances[sub_tree$tip.label])
  p <- p / sum(p)

  sum(dist_mat * outer(p, p))
}

calculate_abundance_faith_pd <- function(tip_abundances, phylo_tree) {
  if (length(tip_abundances) == 0) {
    return(0)
  }
  if (length(tip_abundances) == 1) {
    return(as.numeric(tip_abundances))
  }

  sub_tree <- ape::keep.tip(phylo_tree, names(tip_abundances))
  edge_lengths <- sub_tree$edge.length
  weights <- rep(1, length(edge_lengths))

  for (tip_label in sub_tree$tip.label) {
    tip_index <- match(tip_label, sub_tree$tip.label)
    edge_index <- which(sub_tree$edge[, 2] == tip_index)
    weights[edge_index] <- tip_abundances[tip_label]
  }

  sum(edge_lengths * weights)
}

calculate_phylogenetic_diversity <- function(
  abundance_matrix,
  species_cols,
  phylo_tree,
  sp_to_tip
) {
  all_tips <- phylo_tree$tip.label

  metrics <- t(apply(abundance_matrix, 1, function(abundances) {
    tip_abundances <- make_tip_abundances(
      abundances = as.numeric(abundances),
      species_cols = species_cols,
      sp_to_tip = sp_to_tip,
      all_tips = all_tips
    )

    c(
      faith_pd = calculate_faith_pd(tip_abundances, phylo_tree),
      rao_pd = calculate_rao_pd(tip_abundances, phylo_tree),
      afaith_pd = calculate_abundance_faith_pd(tip_abundances, phylo_tree)
    )
  }))

  as.data.frame(metrics)
}

write_scale_outputs <- function(diversity_sf, config) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(OUTPUT_DIR, paste0(config$output_stem, ".csv"))
  shp_path <- file.path(OUTPUT_DIR, paste0(config$output_stem, ".shp"))

  readr::write_csv(st_drop_geometry(diversity_sf), csv_path)
  st_write(diversity_sf, shp_path, delete_layer = TRUE, quiet = TRUE)

  list(csv = csv_path, shapefile = shp_path)
}

process_scale <- function(config, crowns, species_cols, phylo_tree, sp_to_tip) {
  message("Processing plant diversity for ", config$scale, " quadrats")

  quads <- read_quadrats(config)
  species_matrix <- make_species_matrix(quads, crowns, species_cols)

  abundance_matrix <- species_matrix %>%
    dplyr::select(all_of(species_cols)) %>%
    as.data.frame()

  species_metrics <- calculate_species_diversity(abundance_matrix)
  phylo_metrics <- calculate_phylogenetic_diversity(
    abundance_matrix = abundance_matrix,
    species_cols = species_cols,
    phylo_tree = phylo_tree,
    sp_to_tip = sp_to_tip
  )

  tabular_values <- species_matrix %>%
    left_join(
      bind_cols(
        quad_id = species_matrix$quad_id,
        species_metrics,
        phylo_metrics
      ),
      by = "quad_id"
    )

  diversity_sf <- quads %>%
    left_join(tabular_values, by = "quad_id") %>%
    mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))

  outputs <- write_scale_outputs(diversity_sf, config)

  message(
    "Wrote ", config$scale, " outputs: ",
    outputs$csv, " and ", outputs$shapefile
  )

  diversity_sf
}

run_plant_diversity_workflow <- function(scales = names(SCALE_CONFIG)) {
  project_root <- find_project_root()
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  taxa_table <- read.csv(TAXA_CSV, stringsAsFactors = FALSE) %>%
    filter(sp_code != "COOB2")

  species_cols <- taxa_table$sp_code
  crowns <- read_tree_crowns(species_cols)
  phylo_tree <- build_phylogeny(taxa_table)
  sp_to_tip <- setNames(
    paste(taxa_table$genus, taxa_table$species, sep = "_"),
    taxa_table$sp_code
  )

  results <- lapply(scales, function(scale_name) {
    process_scale(
      config = SCALE_CONFIG[[scale_name]],
      crowns = crowns,
      species_cols = species_cols,
      phylo_tree = phylo_tree,
      sp_to_tip = sp_to_tip
    )
  })
  names(results) <- scales

  invisible(results)
}

if (sys.nframe() == 0) {
  run_plant_diversity_workflow()
}
