USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

if (!requireNamespace("officer", quietly = TRUE)) {
  stop("Missing required R package: officer", call. = FALSE)
}

suppressPackageStartupMessages(library(officer))

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
OUTPUT_DOCX <- "combined_quadrat_variable_guide.docx"

add_definition_list <- function(doc, rows) {
  for (i in seq_len(nrow(rows))) {
    doc <- body_add_par(
      doc,
      value = paste0(rows$column[i], " - ", rows$meaning[i]),
      style = "heading 2"
    )
    doc <- body_add_par(
      doc,
      value = paste0("How calculated: ", rows$calculation[i]),
      style = "Normal"
    )
  }
  doc
}

add_bullets <- function(doc, bullets) {
  for (bullet in bullets) {
    doc <- body_add_par(doc, value = paste0("- ", bullet), style = "Normal")
  }
  doc
}

make_rows <- function(column, meaning, calculation) {
  data.frame(
    column = column,
    meaning = meaning,
    calculation = calculation,
    stringsAsFactors = FALSE
  )
}

write_combined_variable_guide_docx <- function() {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(PROJECT_DIR)

  doc <- read_docx()
  doc <- body_add_par(doc, "Combined Quadrat Analysis Tables", style = "heading 1")
  doc <- body_add_par(
    doc,
    "User-friendly variable guide for quadrat_analysis_10m.csv, quadrat_analysis_20m.csv, and quadrat_analysis_50m.csv",
    style = "Normal"
  )
  doc <- body_add_par(doc, "Last updated: 2026-06-24", style = "Normal")

  doc <- body_add_par(doc, "What Is Included", style = "heading 1")
  doc <- body_add_par(
    doc,
    "Each CSV contains one row per quadrat and 32 value columns. The tables retain identifiers, quadrat center coordinates, species diversity summaries, phylogenetic diversity summaries, spectral heterogeneity metrics, standardized_PCA spectral heterogeneity metrics, and environmental/topographic metrics.",
    style = "Normal"
  )
  doc <- add_bullets(doc, c(
    "10 m table: 2,000 quadrats",
    "20 m table: 500 quadrats",
    "50 m table: 80 quadrats",
    "Center coordinates use NAD83 / UTM zone 16N and come from plant-diversity shapefile polygon centroids."
  ))

  doc <- body_add_par(doc, "What Is Excluded", style = "heading 1")
  doc <- body_add_par(
    doc,
    "The combined CSVs intentionally omit fields that are useful for processing diagnostics but not requested for the value-only analysis tables.",
    style = "Normal"
  )
  doc <- add_bullets(doc, c(
    "Per-species composition columns using species codes",
    "Pixel counts, pair counts, and bootstrap replicate columns",
    "Processing method fields, manual exclusion flags, and geometry columns",
    "Original source identifier duplicates such as Name and sub_id"
  ))
  doc <- body_add_par(
    doc,
    "If per-species composition values are needed later, use the original plant-diversity files in Quad_Values/Diversity_SHPs/plant_diversity_*m.csv.",
    style = "Normal"
  )

  doc <- body_add_par(doc, "Identifier and Coordinate Columns", style = "heading 1")
  doc <- add_definition_list(doc, make_rows(
    column = c("quad_id", "scale", "center_x", "center_y"),
    meaning = c(
      "Scale-aware quadrat identifier.",
      "Quadrat grain: 10m, 20m, or 50m.",
      "X coordinate of the quadrat center.",
      "Y coordinate of the quadrat center."
    ),
    calculation = c(
      "Carried from the aligned source tables and used for all joins.",
      "Assigned from the scale-specific source table.",
      "Calculated as the polygon centroid X coordinate from the plant-diversity shapefile.",
      "Calculated as the polygon centroid Y coordinate from the plant-diversity shapefile."
    )
  ))

  doc <- body_add_par(doc, "Species and Phylogenetic Diversity Columns", style = "heading 1")
  doc <- add_definition_list(doc, make_rows(
    column = c("sp_rich", "sp_shannon", "sp_simpson", "sp_even", "phy_faith", "phy_rao", "phy_afaith"),
    meaning = c(
      "Species richness.",
      "Shannon species diversity.",
      "Simpson species diversity.",
      "Species evenness.",
      "Faith's phylogenetic diversity.",
      "Rao phylogenetic diversity.",
      "Abundance-weighted Faith's phylogenetic diversity."
    ),
    calculation = c(
      "Number of species with positive crown-overlap values in the quadrat.",
      "Calculated from positive species crown-overlap proportions.",
      "Calculated as 1 minus the sum of squared proportional species abundances.",
      "Calculated as Shannon diversity divided by log richness when richness is greater than 1.",
      "Calculated as summed branch length for species present in the quadrat.",
      "Calculated from species crown-overlap weights and phylogenetic distances.",
      "Calculated by weighting phylogenetic branch lengths by species crown-overlap values."
    )
  ))

  doc <- body_add_par(doc, "Spectral Heterogeneity Columns", style = "heading 1")
  doc <- add_definition_list(doc, make_rows(
    column = c(
      "spec_sa", "spec_pca_mean", "spec_pca_med", "spec_pca_sd",
      "spec_rao_q", "spec_alpha", "spec_convex",
      "spec_hull3d_v", "spec_hull3d_a",
      "spec_spca_mean", "spec_spca_med", "spec_spca_sd",
      "spec_spca_rao", "spec_spca_alpha", "spec_spca_convex",
      "spec_spca_hull3d_v", "spec_spca_hull3d_a"
    ),
    meaning = c(
      "Spectral angle entropy mean.",
      "Mean PCA-space spectral distance.",
      "Median PCA-space spectral distance.",
      "Standard deviation of PCA-space spectral distances.",
      "Rao's Q in spectral PCA space.",
      "Alpha-hull area in PC1-PC2 space.",
      "Convex-hull area in PC1-PC2 space.",
      "Convex-hull volume in PC1-PC3 space.",
      "Convex-hull surface area in PC1-PC3 space.",
      "Mean standardized_PCA-space spectral distance.",
      "Median standardized_PCA-space spectral distance.",
      "Standard deviation of standardized_PCA-space spectral distances.",
      "Rao's Q in standardized_PCA spectral space.",
      "Alpha-hull area in standardized_PCA PC1-PC2 space.",
      "Convex-hull area in standardized_PCA PC1-PC2 space.",
      "Convex-hull volume in standardized_PCA PC1-PC3 space.",
      "Convex-hull surface area in standardized_PCA PC1-PC3 space."
    ),
    calculation = c(
      "Primary entropy mean from sunlit, shadow-masked smoothed 5 nm spectra; for most quadrats this is the mean of 70 bootstrap iterations using up to 5,000 retained pixels, while the small exact subset uses all retained pixel pairs.",
      "Mean Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Median Euclidean distance of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Standard deviation of Euclidean distances of retained pixels from the quadrat centroid in global PC1-PC3 spectral space.",
      "Calculated with equal pixel weights and squared Euclidean distance in global PC1-PC3 space.",
      "Calculated as alpha-hull area in global PC1-PC2 spectral space.",
      "Calculated as convex-hull area in global PC1-PC2 spectral space.",
      "Calculated as convex-hull volume in global PC1-PC3 spectral space.",
      "Calculated as convex-hull surface area in global PC1-PC3 spectral space.",
      "Mean Euclidean distance of retained vector-normalized spectra from the quadrat centroid in standardized_PCA PC1-PC3 spectral space.",
      "Median Euclidean distance of retained vector-normalized spectra from the quadrat centroid in standardized_PCA PC1-PC3 spectral space.",
      "Standard deviation of Euclidean distances of retained vector-normalized spectra from the quadrat centroid in standardized_PCA PC1-PC3 spectral space.",
      "Calculated with equal pixel weights and squared Euclidean distance in standardized_PCA PC1-PC3 space.",
      "Calculated as alpha-hull area in standardized_PCA PC1-PC2 spectral space.",
      "Calculated as convex-hull area in standardized_PCA PC1-PC2 spectral space.",
      "Calculated as convex-hull volume in standardized_PCA PC1-PC3 spectral space.",
      "Calculated as convex-hull surface area in standardized_PCA PC1-PC3 spectral space."
    )
  ))

  doc <- body_add_par(doc, "Environmental and Topographic Columns", style = "heading 1")
  doc <- add_definition_list(doc, make_rows(
    column = c("env_elev", "env_tri5", "env_tri11", "env_tri21"),
    meaning = c(
      "Mean elevation.",
      "Mean fine-window topographic roughness.",
      "Mean medium-window topographic roughness.",
      "Mean broad-window topographic roughness."
    ),
    calculation = c(
      "Mean DTM elevation within each quadrat, in meters.",
      "Mean Riley topographic roughness index from a 5x5 DTM moving window.",
      "Mean Riley topographic roughness index from an 11x11 DTM moving window.",
      "Mean Riley topographic roughness index from a 21x21 DTM moving window."
    )
  ))

  doc <- body_add_par(doc, "Validation Snapshot", style = "heading 1")
  doc <- add_bullets(doc, c(
    "10 m: 2,000 rows, 32 columns, 0 duplicate quad IDs, 97 missing spec_sa values, 256 missing spec_pca_mean values, and 0 missing elevation values.",
    "20 m: 500 rows, 32 columns, 0 duplicate quad IDs, 15 missing spec_sa values, 64 missing spec_pca_mean values, and 0 missing elevation values.",
    "50 m: 80 rows, 32 columns, 0 duplicate quad IDs, 0 missing spec_sa values, 6 missing spec_pca_mean values, and 0 missing elevation values.",
    "No per-species composition columns or metadata-style processing columns were detected in the regenerated CSVs."
  ))

  doc <- body_add_par(doc, "Source Files", style = "heading 1")
  doc <- add_bullets(doc, c(
    "Spectral heterogeneity: Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_*_smooth_masked_5nm_summary.csv",
    "Species and phylogenetic diversity: Quad_Values/Diversity_SHPs/plant_diversity_*m.csv",
    "Environmental/topographic values: Quad_Values/Enviro_SHPs/enviro_variables_*m.csv",
    "Reproducible join script: scripts/3_Analysis/combine_quadrat_analysis_tables.R"
  ))

  print(doc, target = OUTPUT_DOCX)
  normalizePath(OUTPUT_DOCX, winslash = "/", mustWork = TRUE)
}

if (sys.nframe() == 0) {
  cat(write_combined_variable_guide_docx(), "\n")
}
