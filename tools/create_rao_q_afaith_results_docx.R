options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root, "Documents", "Paper")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

today_stamp <- "20260805"
output_docx <- file.path(paper_dir, paste0(today_stamp, "_results_spectral_rao_q_afaith_centroid_clarified.docx"))

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
  spectral_variation_paper = extract_docx_text(file.path(root, "Documents", "Spectral Variation Paper.docx")),
  svh_outline = extract_docx_text(file.path(root, "Documents", "SVH Outline.docx")),
  results_sample = extract_docx_text(file.path(root, "Documents", "RESULTS-SampleT.docx")),
  svh_results_read_only = extract_docx_text(file.path(root, "Documents", "Paper", "SVH_Results.docx"))
)

md_files <- c(
  file.path(root, "reports", "project_state.md"),
  file.path(root, "reports", "directory_map.md"),
  file.path(root, "reports", "analysis", "20260710_sv_diversity_pairwise_correlations.md"),
  file.path(root, "reports", "analysis", "20260725_edge_bootstrap_sensitivity_sv_diversity.md"),
  file.path(root, "Documents", "Paper", "20260730_paper_workspace_index.md"),
  file.path(root, "Documents", "Paper", "20260803_writing_style_notes.md"),
  file.path(root, "Documents", "Paper", "20260803_remote_sensing_mdpi_author_requirements.md")
)
md_context <- unlist(lapply(md_files[file.exists(md_files)], readLines, warn = FALSE), use.names = FALSE)

analysis_dataset <- read.csv(file.path(root, "reports", "tables", "multiscale_spectral_biodiversity", "analysis_dataset_with_flags.csv"))
analysis_dataset <- subset(analysis_dataset, primary_analysis == TRUE)

fit_one_scale <- function(scale_name) {
  d <- subset(analysis_dataset, scale == scale_name & !is.na(spec_rao_q) & !is.na(phy_afaith))
  m <- lm(spec_rao_q ~ phy_afaith, data = d)
  sm <- summary(m)
  f <- sm$fstatistic
  spearman <- suppressWarnings(cor.test(d$spec_rao_q, d$phy_afaith, method = "spearman", exact = FALSE))
  data.frame(
    scale = scale_name,
    n = nrow(d),
    pearson_r = unname(cor(d$spec_rao_q, d$phy_afaith)),
    r_squared = unname(sm$r.squared),
    f_statistic = unname(f[["value"]]),
    f_df1 = unname(f[["numdf"]]),
    f_df2 = unname(f[["dendf"]]),
    f_p_value = pf(f[["value"]], f[["numdf"]], f[["dendf"]], lower.tail = FALSE),
    slope = unname(coef(m)[["phy_afaith"]]),
    intercept = unname(coef(m)[["(Intercept)"]]),
    spearman_r = unname(spearman$estimate),
    spearman_p_value = spearman$p.value
  )
}

focal <- do.call(rbind, lapply(c("10m", "20m", "50m"), fit_one_scale))

fmt_p <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x))
}

row_text <- function(scale) {
  row <- focal[focal$scale == scale, ]
  paste0(
    scale, ": n = ", row$n,
    ", Pearson r = ", sprintf("%.3f", row$pearson_r),
    ", R2 = ", sprintf("%.3f", row$r_squared),
    ", F(", row$f_df1, ", ", row$f_df2, ") = ", sprintf("%.2f", row$f_statistic),
    ", p = ", fmt_p(row$f_p_value),
    ", slope = ", sprintf("%.6f", row$slope),
    ", Spearman rho = ", sprintf("%.3f", row$spearman_r),
    "."
  )
}

styles_xml <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:pPr><w:jc w:val="both"/><w:spacing w:after="160" w:line="320" w:lineRule="auto"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:pPr><w:spacing w:after="160"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:b/><w:color w:val="0B2545"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:pPr><w:spacing w:after="200"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:color w:val="555555"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="360" w:after="200"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="120"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="26"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Placeholder"><w:name w:val="Placeholder"/><w:pPr><w:spacing w:before="80" w:after="160"/><w:ind w:left="360" w:right="360"/><w:shd w:fill="F4F6F9"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:i/><w:color w:val="1F3A5F"/><w:sz w:val="22"/></w:rPr></w:style>
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

body <- c(
  make_para("Results Draft: Spectral Rao's Q and Abundance-Weighted Faith's PD", "Title"),
  make_para("Created 2026-08-05; focused Results-section prose for the Remote Sensing manuscript.", "Subtitle"),
  make_para("Style/context note: this draft was generated after reading current project Markdown, the multiscale analysis outputs, the existing writing-style notes, and accessible text from the provided Word references. `SVH_Results.docx` was read only for context and was not edited. The prose is intentionally limited to direct spectral diversity and abundance-weighted phylogenetic diversity relationships, without environmental-factor interpretation.", "Normal", italic = TRUE),
  make_para("3. Results", "Heading1"),
  make_para("3.1 Scale-dependent relationship between spectral Rao's Q and abundance-weighted Faith's PD", "Heading2"),
  make_para("Spectral Rao's Q showed a weak to moderate positive relationship with abundance-weighted Faith's phylogenetic diversity, with the clearest association emerging at the 50 m grain. This metric summarizes spectral dispersion in global PC1-PC3 space using equal pixel weights and squared Euclidean distance, so the result is framed as a relationship between within-quadrat spectral dissimilarity and abundance-weighted phylogenetic breadth."),
  make_para("Calculation note: in the current R workflow, spectral Rao's Q is not calculated by explicitly constructing every pixel-pair distance. The script first calculates the quadrat centroid in PC1-PC3 space, then computes squared distances from each retained pixel to that centroid. Rao's Q is stored as `2 * mean(squared_radius)`, where `squared_radius` is the squared Euclidean distance from each pixel to the quadrat centroid. For equal pixel weights, this centroid-based expression is mathematically equivalent to the mean pairwise squared Euclidean distance among pixels. Thus, the metric is conceptually pairwise squared-distance Rao's Q, but computationally it is implemented through the equal-weight centroid shortcut."),
  make_para("[Figure placeholder: scale-specific scatterplots of spectral Rao's Q versus abundance-weighted Faith's PD at 10 m, 20 m, and 50 m. Include regression lines, confidence bands, point counts, and consistent axis labels. In text, cite this figure immediately after the opening scale-pattern sentence.]", "Placeholder"),
  make_para("At 10 m, spectral Rao's Q was not meaningfully associated with abundance-weighted Faith's PD. The primary analysis included 1,440 quadrats and produced a near-zero Pearson correlation (r = 0.003, R2 < 0.001; F(1, 1438) = 0.01, p = 0.907). The rank correlation was similarly weak, indicating that fine-grain spectral Rao's Q did not track abundance-weighted phylogenetic diversity in a consistent monotonic way at this scale."),
  make_para("At 20 m, the relationship became positive but remained modest. Across 360 primary quadrats, spectral Rao's Q increased with abundance-weighted Faith's PD (r = 0.112, R2 = 0.013; F(1, 358) = 4.57, p = 0.033). The effect was weaker than the standardized PCA mean-distance result at the same scale, but it suggests some recovery of the spectral-phylogenetic signal as quadrats integrate over a broader neighborhood."),
  make_para("At 50 m, spectral Rao's Q showed the strongest relationship with abundance-weighted Faith's PD. Across 74 complete quadrats, the correlation was positive and statistically detectable (r = 0.345, R2 = 0.119; F(1, 72) = 9.74, p = 0.003). This scale therefore produced the clearest evidence that greater abundance-weighted phylogenetic diversity corresponded to greater within-quadrat spectral dissimilarity as measured by spectral Rao's Q."),
  make_para("[Table placeholder: primary pairwise results for spectral Rao's Q versus abundance-weighted Faith's PD by scale. Columns should include scale, n, Pearson r, R2, F statistic, p-value, slope, intercept, Spearman rho, and Spearman p-value. In text, cite this table after reporting the three scale-specific models.]", "Placeholder"),
  make_para("3.2 Interpretation of the scale pattern within the Results section", "Heading2"),
  make_para("The scale pattern for spectral Rao's Q differs from the stronger standardized PCA mean-distance pattern. Spectral Rao's Q was essentially absent at 10 m, modest at 20 m, and most apparent at 50 m. This makes spectral Rao's Q useful as a complementary result: it supports the broader scale-dependence of the spectral-phylogenetic relationship, but it does not appear to be the strongest spectral metric at finer grains."),
  make_para("[Figure placeholder: comparison panel showing Pearson r or R2 for spectral Rao's Q across the three scales, optionally paired with the standardized PCA mean-distance values for context. Keep the main visual focused on spectral Rao's Q if this remains a separate Results subsection.]", "Placeholder"),
  make_para("In the final manuscript, this subsection can be positioned after the standardized PCA mean-distance results. The recommended structure is to present spectral Rao's Q as a secondary spectral diversity measure that reinforces the strongest scale-level conclusion at 50 m while showing that not all spectral heterogeneity metrics recover the relationship equally at 10 m and 20 m. When describing the metric, use the wording 'equal-weight pairwise squared Euclidean Rao's Q, computed via the equivalent centroid-based formula' to avoid implying that the script directly calculated an explicit pairwise distance matrix."),
  make_para("No environmental covariates are included in this Results subsection. The intent is to present the direct bivariate relationship between spectral Rao's Q and abundance-weighted Faith's PD across spatial scales. Environmental and spatial sensitivity results should remain separate unless the manuscript later adds a robustness subsection or supplement."),
  make_para("Draft result values checked from `reports/tables/multiscale_spectral_biodiversity/analysis_dataset_with_flags.csv`, using `primary_analysis == TRUE` rows only:", "Normal", italic = TRUE),
  make_para(row_text("10m"), "Normal", italic = TRUE),
  make_para(row_text("20m"), "Normal", italic = TRUE),
  make_para(row_text("50m"), "Normal", italic = TRUE)
)

document_xml <- paste0(
  '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
  '<w:document xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" mc:Ignorable=""><w:body>',
  paste(body, collapse = ""),
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
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Results Draft Spectral Rao Q and Abundance-Weighted Faith PD</dc:title><dc:creator>Codex via R</dc:creator><cp:lastModifiedBy>Codex via R</cp:lastModifiedBy><dcterms:created xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:modified></cp:coreProperties>', file.path(tmp_docx, "docProps", "core.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>R base OOXML builder</Application></Properties>', file.path(tmp_docx, "docProps", "app.xml"), useBytes = TRUE)

zipfile <- tempfile(fileext = ".zip")
compress_archive(tmp_docx, zipfile)
file.copy(zipfile, output_docx, overwrite = TRUE)

cat("Created:", output_docx, "\n")
cat("Style source paragraphs read:",
    "Spectral Variation Paper =", length(style_sources$spectral_variation_paper),
    "; SVH Outline =", length(style_sources$svh_outline),
    "; RESULTS-SampleT =", length(style_sources$results_sample),
    "; SVH_Results read-only =", length(style_sources$svh_results_read_only), "\n")
cat("Markdown/context lines read:", length(md_context), "\n")
print(focal)
