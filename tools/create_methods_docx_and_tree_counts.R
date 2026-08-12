options(stringsAsFactors = FALSE)

USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root, "Documents", "Paper")
figure_dir <- file.path(root, "reports", "figures", "methods_tree_counts")
table_dir <- file.path(root, "reports", "tables", "methods_tree_counts")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

today_stamp <- "20260811"
output_docx <- file.path(paper_dir, paste0(today_stamp, "_methods_draft_spectral_phylogenetic_diversity.docx"))

required_packages <- c("sf", "dplyr", "tidyr", "readr", "ape", "picante", "V.PhyloMaker2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
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

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

ps_quote <- function(path) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "\\", mustWork = FALSE), fixed = TRUE), "'")
}

expand_archive <- function(zip_path, destination) {
  cmd <- paste("Expand-Archive", "-LiteralPath", ps_quote(zip_path), "-DestinationPath", ps_quote(destination), "-Force")
  status <- system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(destination, "word", "document.xml"))) {
    stop("Archive extraction failed for ", zip_path, "\n", paste(status, collapse = "\n"))
  }
}

compress_archive <- function(source_dir, out_zip) {
  if (file.exists(out_zip)) unlink(out_zip)
  cmd <- paste("Compress-Archive", "-Path", ps_quote(file.path(source_dir, "*")), "-DestinationPath", ps_quote(out_zip), "-Force")
  status <- system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE)
  if (!file.exists(out_zip)) stop("Archive compression failed.\n", paste(status, collapse = "\n"))
}

extract_docx_text <- function(docx_path) {
  if (!file.exists(docx_path)) return(character())
  tmp <- tempfile("docx_extract_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  zip_copy <- tempfile(fileext = ".zip")
  copied <- tryCatch(file.copy(docx_path, zip_copy, overwrite = TRUE), error = function(e) FALSE)
  if (!isTRUE(copied)) return(character())
  ok <- tryCatch({ expand_archive(zip_copy, tmp); TRUE }, error = function(e) FALSE)
  if (!ok) return(character())
  xml <- paste(readLines(file.path(tmp, "word", "document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
  paras <- unlist(strsplit(xml, "</w:p>", fixed = TRUE), use.names = FALSE)
  out <- vapply(paras, function(p) {
    texts <- regmatches(p, gregexpr("<w:t[^>]*>.*?</w:t>", p, perl = TRUE))[[1]]
    if (length(texts) == 0 || identical(texts, character(0))) return("")
    texts <- gsub("^<w:t[^>]*>", "", texts, perl = TRUE)
    texts <- gsub("</w:t>$", "", texts, perl = TRUE)
    texts <- gsub("&amp;", "&", texts, fixed = TRUE)
    texts <- gsub("&lt;", "<", texts, fixed = TRUE)
    texts <- gsub("&gt;", ">", texts, fixed = TRUE)
    texts <- gsub("&quot;", "\"", texts, fixed = TRUE)
    trimws(paste(texts, collapse = ""))
  }, character(1))
  out[nzchar(out)]
}

style_sources <- list(
  methods_sample = extract_docx_text(file.path(root, "Documents", "Methods-SampleT.docx")),
  spectral_variation_paper = extract_docx_text(file.path(root, "Documents", "Spectral Variation Paper.docx")),
  svh_paper_read_only = extract_docx_text(file.path(root, "Documents", "Paper", "SVH_Paper.docx"))
)

md_files <- c(
  file.path(root, "reports", "project_state.md"),
  file.path(root, "reports", "directory_map.md"),
  file.path(root, "reports", "analysis", "20260710_sv_diversity_pairwise_correlations.md"),
  file.path(root, "reports", "validation", "20260622_spectral_heterogeneity_all_metrics_validation.md"),
  file.path(root, "reports", "validation", "20260624_combined_quadrat_analysis_tables_validation.md"),
  file.path(root, "reports", "validation", "20260624_multiscale_spectral_biodiversity_analysis_validation.md"),
  file.path(root, "Documents", "Paper", "20260803_writing_style_notes.md"),
  file.path(root, "Documents", "Paper", "20260803_remote_sensing_mdpi_author_requirements.md")
)
md_context <- unlist(lapply(md_files[file.exists(md_files)], readLines, warn = FALSE), use.names = FALSE)

CRS_PROJ <- 26916
taxa <- read.csv(file.path(root, "51sp_taxanomy.csv"), stringsAsFactors = FALSE) %>%
  filter(sp_code != "COOB2") %>%
  mutate(
    Genus = genus,
    Species = species
  )
species_cols <- taxa$sp_code

tree_df <- readr::read_csv(file.path(root, "PR_tree_DL.csv"), show_col_types = FALSE) %>%
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
  ) %>%
  select(-any_of(c("Genus", "Species", "Family", "family"))) %>%
  left_join(taxa %>% select(sp_code, Genus, Species, family), by = c("sp" = "sp_code"))

tree_points <- st_as_sf(
  tree_df,
  coords = c("UTMX_CURRENT", "UTMY_CURRENT"),
  crs = CRS_PROJ,
  remove = FALSE
)
tree_crowns <- st_buffer(tree_points, dist = tree_points$radius_m, endCapStyle = "ROUND") %>%
  select(crown_id, sp, Genus, Species, family, quadrat, tag, DBH.2024, crown.position)

read_quads <- function(scale_name) {
  shp_path <- file.path(root, "Quad_Scale_SHPs", paste0("PR_", scale_name, ".shp"))
  id_col <- if (scale_name == "10m") "sub_id" else "Name"
  st_read(shp_path, quiet = TRUE) %>%
    st_transform(CRS_PROJ) %>%
    mutate(
      quad_id = as.character(.data[[id_col]]),
      scale = scale_name
    ) %>%
    select(quad_id, scale)
}

analysis_data <- read.csv(
  file.path(root, "reports", "tables", "multiscale_spectral_biodiversity", "sv_diversity_analysis_dataset.csv"),
  stringsAsFactors = FALSE
)

included_quads <- analysis_data %>%
  filter(primary_analysis, !is.na(spec_spca_mean), !is.na(phy_afaith)) %>%
  select(scale, quad_id)

count_species <- function(crowns) {
  crowns %>%
    st_drop_geometry() %>%
    distinct(crown_id, .keep_all = TRUE) %>%
    count(Genus, Species, name = "Number of individuals") %>%
    arrange(Genus, Species)
}

all_plot_counts <- count_species(tree_crowns)

count_included_scale <- function(scale_name) {
  quads <- read_quads(scale_name) %>%
    inner_join(included_quads %>% filter(scale == scale_name), by = c("scale", "quad_id"))
  intersections <- suppressWarnings(st_intersection(tree_crowns, quads))
  count_species(intersections)
}

included_counts <- list(
  "10m" = count_included_scale("10m"),
  "20m" = count_included_scale("20m"),
  "50m" = count_included_scale("50m")
)

write_count_csv <- function(data, filename) {
  path <- file.path(table_dir, filename)
  readr::write_csv(data, path)
  path
}

count_paths <- c(
  all_plot = write_count_csv(all_plot_counts, paste0(today_stamp, "_all_plot_tree_counts_by_species.csv")),
  included_10m = write_count_csv(included_counts[["10m"]], paste0(today_stamp, "_included_10m_tree_counts_by_species.csv")),
  included_20m = write_count_csv(included_counts[["20m"]], paste0(today_stamp, "_included_20m_tree_counts_by_species.csv")),
  included_50m = write_count_csv(included_counts[["50m"]], paste0(today_stamp, "_included_50m_tree_counts_by_species.csv"))
)

combined_counts <- all_plot_counts %>%
  rename(all_plot_individuals = `Number of individuals`) %>%
  full_join(included_counts[["10m"]] %>% rename(included_10m_individuals = `Number of individuals`), by = c("Genus", "Species")) %>%
  full_join(included_counts[["20m"]] %>% rename(included_20m_individuals = `Number of individuals`), by = c("Genus", "Species")) %>%
  full_join(included_counts[["50m"]] %>% rename(included_50m_individuals = `Number of individuals`), by = c("Genus", "Species")) %>%
  mutate(across(ends_with("individuals"), ~ tidyr::replace_na(.x, 0L))) %>%
  arrange(Genus, Species)
combined_count_path <- write_count_csv(combined_counts, paste0(today_stamp, "_tree_counts_summary_by_species_and_analysis_set.csv"))

save_count_plot <- function(data, filename, title) {
  png_path <- file.path(figure_dir, filename)
  data <- data %>%
    mutate(label = paste(Genus, Species)) %>%
    arrange(`Number of individuals`)
  height_px <- max(1400, 90 * nrow(data) + 300)
  grDevices::png(png_path, width = 2400, height = height_px, res = 220)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mar = c(5, 12, 4, 2))
  barplot(
    data$`Number of individuals`,
    names.arg = data$label,
    horiz = TRUE,
    las = 1,
    col = "#4C78A8",
    border = NA,
    main = title,
    xlab = "Number of individuals",
    cex.names = 0.72
  )
  png_path
}

figure_paths <- c(
  all_plot = save_count_plot(all_plot_counts, paste0(today_stamp, "_all_plot_tree_counts_by_species.png"), "All plot trees used for phylogenetic diversity"),
  included_10m = save_count_plot(included_counts[["10m"]], paste0(today_stamp, "_included_10m_tree_counts_by_species.png"), "Trees intersecting included 10 m analysis quadrats"),
  included_20m = save_count_plot(included_counts[["20m"]], paste0(today_stamp, "_included_20m_tree_counts_by_species.png"), "Trees intersecting included 20 m analysis quadrats"),
  included_50m = save_count_plot(included_counts[["50m"]], paste0(today_stamp, "_included_50m_tree_counts_by_species.png"), "Trees intersecting included 50 m analysis quadrats")
)

build_phylogeny <- function(taxa_table) {
  phylo_input <- taxa_table %>%
    mutate(species_name = paste(genus, species, sep = " ")) %>%
    select(species = species_name, genus, family) %>%
    filter(species != "Carya sp") %>%
    distinct(species, .keep_all = TRUE)

  phylo.maker(
    sp.list = phylo_input,
    tree = GBOTB.extended.TPL,
    nodes = nodes.info.1.TPL,
    scenarios = "S3"
  )$scenario.3
}

phylo_tree <- build_phylogeny(taxa)
analysis_summary <- analysis_data %>%
  group_by(scale) %>%
  summarise(
    total_quadrats = n(),
    primary_quadrats = sum(primary_analysis),
    focal_complete_quadrats = sum(primary_analysis & !is.na(spec_spca_mean) & !is.na(phy_afaith)),
    complete_sa_quadrats = sum(primary_analysis & !is.na(spec_sa)),
    .groups = "drop"
  )

tree_total <- nrow(tree_df)
species_total <- nrow(all_plot_counts)
included_tree_totals <- sapply(included_counts, function(x) sum(x$`Number of individuals`))

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

styles_xml <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:pPr><w:jc w:val="both"/><w:spacing w:after="160" w:line="320" w:lineRule="auto"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:pPr><w:spacing w:after="160"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:b/><w:color w:val="0B2545"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:pPr><w:spacing w:after="200"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:color w:val="555555"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="360" w:after="200"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="120"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="26"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Note"><w:name w:val="Note"/><w:pPr><w:spacing w:before="80" w:after="160"/><w:ind w:left="360" w:right="360"/><w:shd w:fill="F4F6F9"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:i/><w:color w:val="1F3A5F"/><w:sz w:val="22"/></w:rPr></w:style>
</w:styles>'

make_para <- function(text, style = "Normal", bold = FALSE, italic = FALSE, color = NULL, size = NULL) {
  props <- c()
  if (bold) props <- c(props, "<w:b/>")
  if (italic) props <- c(props, "<w:i/>")
  if (!is.null(color)) props <- c(props, paste0("<w:color w:val=\"", color, "\"/>"))
  if (!is.null(size)) props <- c(props, paste0("<w:sz w:val=\"", as.integer(size * 2), "\"/>"))
  rpr <- if (length(props)) paste0("<w:rPr>", paste(props, collapse = ""), "</w:rPr>") else ""
  paste0("<w:p><w:pPr><w:pStyle w:val=\"", style, "\"/></w:pPr><w:r>", rpr, "<w:t xml:space=\"preserve\">", xml_escape(text), "</w:t></w:r></w:p>")
}

method_paragraphs <- c(
  make_para("Methods Draft: Spectral Heterogeneity and Phylogenetic Diversity Workflow", "Title"),
  make_para("Created 2026-08-11; draft Methods section for manuscript development.", "Subtitle"),
  make_para("Context note: this draft was generated from the active R workflow scripts, current Markdown project reports, `Methods-SampleT.docx`, and `Spectral Variation Paper.docx`. Species-count tables and charts were generated as external files under `reports/figures/methods_tree_counts/` and `reports/tables/methods_tree_counts/`; they are intentionally not embedded in this document.", "Note"),
  make_para("2. Materials and Methods", "Heading1"),
  make_para("2.1 Study Site and Plot Data", "Heading2"),
  make_para(paste0("This study was conducted in the Paint Rock Forest Dynamics Plot using georeferenced tree-census data, quadrat boundary layers, and drone-based hyperspectral imagery. Site and background language should be refined against the existing Spectral Variation Paper text during final editing. The analytical workflow used trees from the 2024 tree table that met the crown and phylogenetic-diversity inclusion criteria: DBH in 2024 greater than or equal to 200 mm, crown position classes 3-5, cluster status A or R, valid current UTM coordinates, valid 2025 crown-width estimates, and species codes present in the 51-species taxonomy table after excluding `COOB2`. Under these criteria, ", fmt_int(tree_total), " tree crowns representing ", fmt_int(species_total), " species were available for plot-level phylogenetic-diversity calculations.")),
  make_para("Tree crowns were represented as circular crown buffers centered on the current UTM coordinates. Crown radius was calculated as one half of the 2025 crown-width estimate. These crown buffers were intersected with 10 m, 20 m, and 50 m quadrat polygons to estimate species contributions within each quadrat. Species abundances were therefore based on summed crown-overlap proportions rather than stem counts alone, preserving the canopy-weighted structure of the original plant-diversity workflow."),
  make_para("External tree-count outputs: all-plot and included-analysis species-count charts were written to `reports/figures/methods_tree_counts/`; matching CSV tables were written to `reports/tables/methods_tree_counts/`. These files provide Genus, Species, and Number of individuals summaries and are not embedded in this Methods draft.", "Note"),
  make_para("2.2 Quadrat Scales and Analysis Inclusion", "Heading2"),
  make_para("Analyses were conducted at three spatial grains: 10 m, 20 m, and 50 m. Quadrat identifiers were harmonized across plant-diversity, spectral-heterogeneity, environmental, and analysis-ready tables. The 10 m workflow used `sub_id` identifiers, the 20 m workflow used `Name` identifiers, and the 50 m workflow used `Name` identifiers such as `sub50_1`. Combined analysis tables were generated for each scale and retained quadrat centroid coordinates, species and phylogenetic diversity metrics, spectral heterogeneity metrics, environmental values, and topographic values."),
  make_para(paste0("Documented edge quadrats were excluded from the primary 10 m and 20 m analyses. No separate 50 m edge rule was documented, so all 50 m quadrats were initially retained before metric-specific completeness checks. For the focal standardized PCA mean-distance and abundance-weighted Faith's PD results, complete primary-analysis sample sizes were 1,440 quadrats at 10 m, 360 quadrats at 20 m, and 74 quadrats at 50 m. The corresponding tree-count summaries include ", fmt_int(included_tree_totals[['10m']]), " trees intersecting included 10 m quadrats, ", fmt_int(included_tree_totals[['20m']]), " trees intersecting included 20 m quadrats, and ", fmt_int(included_tree_totals[['50m']]), " trees intersecting included 50 m quadrats.")),
  make_para("2.3 Species and Phylogenetic Diversity Metrics", "Heading2"),
  make_para("Species matrices were calculated separately for each quadrat scale by intersecting buffered crowns with quadrat polygons, dividing each crown-intersection area by the full buffered crown area, and summing proportional crown overlap by species and quadrat. The resulting species-by-quadrat abundance matrix was used to calculate species richness, Shannon diversity, Simpson diversity, and evenness. Shannon diversity was calculated as the negative sum of p log p across species with positive crown-overlap values; Simpson diversity was calculated as one minus the sum of squared proportional abundances; and evenness was calculated as Shannon diversity divided by log richness when richness exceeded one."),
  make_para("Phylogenetic diversity was calculated from the 51-species taxonomy table using `V.PhyloMaker2` with the GBOTB.extended.TPL backbone and scenario S3. The workflow generated a plot-level phylogeny after excluding the unresolved `Carya sp` entry from tree construction. Species codes were mapped to phylogenetic tip labels using genus and species names. For each quadrat, species crown-overlap values were mapped to phylogenetic tips before calculating Faith's PD, phylogenetic Rao's Q, and abundance-weighted Faith's PD. Faith's PD was calculated as summed branch length for species present in a quadrat. Phylogenetic Rao's Q was calculated from the cophenetic phylogenetic distance matrix and normalized species weights. Abundance-weighted Faith's PD was calculated by weighting terminal branch contributions by the crown-overlap abundance of each tip before summing branch lengths."),
  make_para("2.4 Hyperspectral Quadrat Processing", "Heading2"),
  make_para("The spectral workflow used the confirmed partitioned quadrat spectra from `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m`. Spectra were smoothed using a Savitzky-Golay workflow and then resampled to 5 nm spacing. Current analysis inputs are the smoothed 5 nm quadrat spectra in `Quad_Spectra/10m_smooth_5nm`, `Quad_Spectra/20m_smooth_5nm`, and `Quad_Spectra/50m_smooth_5nm`. A 563 nm shadow mask was applied using threshold 0.0305476, retaining pixels above the threshold under the current workflow direction."),
  make_para("Spectral angle entropy was calculated from retained sunlit pixels within each quadrat. Exact all-pixel spectral angle entropy was used only for sufficiently small quadrats because pairwise angle counts scale quadratically with pixel number. Larger quadrats used a bootstrap-mean fallback, with 70 bootstrap replicates, up to 5,000 sampled pixels per replicate, and 10,000 sampled pixel pairs per large bootstrap replicate. Bootstrap standard deviation, coefficient of variation, confidence intervals, and method fields were retained as quality-control metadata."),
  make_para("2.5 PCA-Based Spectral Heterogeneity Metrics", "Heading2"),
  make_para("PCA-based spectral heterogeneity metrics were calculated from global PCA coordinates. The current PCA basis was rebuilt from the 10 m spectral footprint only, avoiding nested multiscale resampling of the same footprint. Up to 450 illuminated pixels were sampled from each current 10 m smoothed 5 nm raster after shadow masking. The raw global PCA used 854,767 sampled pixels and 121 spectral bands. A second standardized PCA basis used the same sampling design after vector-normalizing each retained spectrum. The standardized PCA basis explained 45.46% of variation on PC1, 22.54% on PC2, and 6.80% on PC3, for a cumulative 74.80% across PC1-PC3."),
  make_para("The focal PCA metric in the Results is standardized PCA mean Euclidean distance. For each quadrat, retained spectra were projected into the standardized PCA basis, a quadrat centroid was calculated in PC1-PC3 space, each retained pixel's Euclidean distance from that centroid was calculated, and the mean of those distances was used as the quadrat-level spectral heterogeneity value. This metric asks how far the average retained pixel lies from the quadrat's spectral centroid in standardized PCA space."),
  make_para("Spectral Rao's Q was also calculated in PCA space. The active script calculates the quadrat centroid, computes each retained pixel's squared Euclidean distance from that centroid, and stores Rao's Q as `2 * mean(squared_radius)`. With equal pixel weights, this centroid-based expression is mathematically equivalent to the mean pairwise squared Euclidean distance among retained pixels. Thus, spectral Rao's Q is conceptually an equal-weight pairwise squared Euclidean diversity metric, but computationally it is implemented through the equivalent centroid formula."),
  make_para("2.6 Combined Quadrat Tables and Statistical Analysis", "Heading2"),
  make_para("Scale-specific plant-diversity, spectral-heterogeneity, environmental, and quadrat-center data were joined into `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`. Per-species composition columns, pixel-count fields, pair-count fields, bootstrap replicate fields, method fields, manual-exclusion fields, and geometry columns were excluded from these compact combined CSVs. Metadata and sampling provenance were retained in companion spectral-heterogeneity and multiscale-analysis output tables."),
  make_para("The focal first-layer statistical analysis was intentionally direct. For each scale, standardized PCA mean distance and spectral angle entropy were independently paired with species and phylogenetic diversity measures using simple linear models. The model form was `SV_measure ~ diversity_measure`, fit separately within each scale. Reported statistics include Pearson correlation, R2, F statistic, F-test p-value, slope, intercept, and Spearman rank correlation. In one-predictor models, R2 is the squared Pearson correlation. The current focal manuscript emphasis is on abundance-weighted Faith's PD and Shannon diversity, with the strongest current support appearing for phylogenetic diversity metrics rather than species-diversity metrics."),
  make_para("2.7 Reproducibility and Software", "Heading2"),
  make_para("All processing and analysis workflows were implemented in R. The active project scripts use R 4.2.3 from `C:/Program Files/R/R-4.2.3/bin/Rscript.exe` and project-accessible packages including `sf`, `terra`, `dplyr`, `tidyr`, `readr`, `ape`, `picante`, `V.PhyloMaker2`, `alphahull`, and `geometry`. Primary workflows write intermediate and final outputs to `Quad_Values/`, `reports/analysis/`, `reports/tables/`, and `reports/figures/`. Manuscript-preparation outputs are stored under `Documents/Paper/` with dated filenames.")
)

document_xml <- paste0(
  '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
  '<w:document xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" mc:Ignorable=""><w:body>',
  paste(method_paragraphs, collapse = ""),
  '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/></w:sectPr>',
  '</w:body></w:document>'
)

tmp_docx <- tempfile("docx_build_")
dir.create(file.path(tmp_docx, "_rels"), recursive = TRUE)
dir.create(file.path(tmp_docx, "word", "_rels"), recursive = TRUE)
dir.create(file.path(tmp_docx, "docProps"), recursive = TRUE)

writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/><Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/></Types>', file.path(tmp_docx, "[Content_Types].xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/><Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/></Relationships>', file.path(tmp_docx, "_rels", ".rels"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>', file.path(tmp_docx, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
writeLines(document_xml, file.path(tmp_docx, "word", "document.xml"), useBytes = TRUE)
writeLines(styles_xml, file.path(tmp_docx, "word", "styles.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Methods Draft Spectral Phylogenetic Diversity</dc:title><dc:creator>Codex via R</dc:creator><cp:lastModifiedBy>Codex via R</cp:lastModifiedBy><dcterms:created xsi:type="dcterms:W3CDTF">2026-08-11T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-11T00:00:00Z</dcterms:modified></cp:coreProperties>', file.path(tmp_docx, "docProps", "core.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>R base OOXML builder</Application></Properties>', file.path(tmp_docx, "docProps", "app.xml"), useBytes = TRUE)

zipfile <- tempfile(fileext = ".zip")
compress_archive(tmp_docx, zipfile)
file.copy(zipfile, output_docx, overwrite = TRUE)

cat("Created:", output_docx, "\n")
cat("Style/context paragraphs read: Methods-SampleT =", length(style_sources$methods_sample),
    "; Spectral Variation Paper =", length(style_sources$spectral_variation_paper),
    "; SVH_Paper read-only =", length(style_sources$svh_paper_read_only), "\n")
cat("Markdown/context lines read:", length(md_context), "\n")
cat("Tree total:", tree_total, "Species total:", species_total, "\n")
print(analysis_summary)
cat("CSV outputs:\n", paste(c(count_paths, combined_count_path), collapse = "\n"), "\n")
cat("Figure outputs:\n", paste(figure_paths, collapse = "\n"), "\n")
