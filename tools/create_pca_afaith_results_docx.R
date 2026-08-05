options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root, "Documents", "Paper")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

today_stamp <- "20260805"
output_docx <- file.path(paper_dir, paste0(today_stamp, "_results_pca_mean_distance_afaith.docx"))

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
  results_sample = extract_docx_text(file.path(root, "Documents", "RESULTS-SampleT.docx"))
)

md_context <- c(
  readLines(file.path(root, "reports", "project_state.md"), warn = FALSE),
  readLines(file.path(root, "reports", "analysis", "20260710_sv_diversity_pairwise_correlations.md"), warn = FALSE),
  readLines(file.path(root, "Documents", "Paper", "20260803_writing_style_notes.md"), warn = FALSE),
  readLines(file.path(root, "Documents", "Paper", "20260803_remote_sensing_mdpi_author_requirements.md"), warn = FALSE)
)

results <- read.csv(file.path(root, "reports", "tables", "multiscale_spectral_biodiversity", "sv_diversity_pairwise_correlations.csv"))
focal <- subset(results, sv_measure == "spec_spca_mean" & diversity_measure == "phy_afaith")
focal <- focal[match(c("10m", "20m", "50m"), focal$scale), ]

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
  make_para("Results Draft: Standardized PCA Mean Distance and Abundance-Weighted Faith's PD", "Title"),
  make_para("Created 2026-08-05; focused Results-section prose for the Remote Sensing manuscript.", "Subtitle"),
  make_para("Style/context note: this draft was generated after reading the current project-state Markdown, the July 10 pairwise correlation report, the paper style notes, Remote Sensing author requirements, and accessible text from the provided Word references. The prose is intentionally results-focused and avoids discussion of mechanisms until the Discussion section.", "Normal", italic = TRUE),
  make_para("3. Results", "Heading1"),
  make_para("3.1 Scale-dependent relationship between standardized PCA mean distance and abundance-weighted Faith's PD", "Heading2"),
  make_para("Standardized PCA mean Euclidean distance showed a positive relationship with abundance-weighted Faith's phylogenetic diversity at each quadrat scale. The association was weakest at the 10 m grain, strengthened at 20 m, and was strongest at the 50 m grain. This pattern indicates that the spectral heterogeneity signal associated with abundance-weighted phylogenetic diversity became clearer as quadrats integrated broader canopy and community structure."),
  make_para("[Figure 3 placeholder: scale-specific scatterplots of standardized PCA mean Euclidean distance versus abundance-weighted Faith's PD at 10 m, 20 m, and 50 m. Include regression lines, confidence bands, point counts, and consistent axis labels across panels. In text, cite as Fig. 3 after the first sentence describing the scale trend.]", "Placeholder"),
  make_para("At 10 m, standardized PCA mean distance was positively but weakly correlated with abundance-weighted Faith's PD. The primary non-edge analysis included 1,440 quadrats and produced a Pearson correlation of r = 0.106, explaining approximately 1.1% of the variance in standardized PCA mean distance (R2 = 0.011). The relationship was statistically detectable given the large number of quadrats, F(1, 1438) = 16.36, p = 5.51e-05, but the effect size was small."),
  make_para("At 20 m, the relationship was stronger. Across 360 primary quadrats, standardized PCA mean distance increased with abundance-weighted Faith's PD (r = 0.291, R2 = 0.085, F(1, 358) = 33.17, p = 1.82e-08). This intermediate scale therefore captured a more pronounced spectral-phylogenetic association than the 10 m grain, while retaining substantially more observations than the 50 m grain."),
  make_para("At 50 m, the association was strongest among the three scales. Across 74 complete quadrats, standardized PCA mean distance was positively correlated with abundance-weighted Faith's PD (r = 0.419, R2 = 0.175, F(1, 72) = 15.31, p = 2.05e-04). Although the 50 m analysis had fewer quadrats, the larger effect size indicates that broader spatial aggregation increased the apparent coupling between spectral heterogeneity and abundance-weighted phylogenetic diversity."),
  make_para("[Table 2 placeholder: primary pairwise results for standardized PCA mean distance versus abundance-weighted Faith's PD by scale. Columns should include scale, n, Pearson r, R2, F statistic, p-value, slope, intercept, and Spearman r. In text, cite as Table 2 at the end of the paragraph summarizing all three scale-specific models.]", "Placeholder"),
  make_para("The slope estimates were positive at every scale, consistent with the correlation results. The fitted slope was highest at 10 m (0.001053), similar but slightly lower at 20 m (0.000948), and lower at 50 m (0.000245). Because the scale-specific predictor and response distributions differ across quadrat grains, the standardized correlation and R2 values are the clearest basis for comparing relationship strength among scales."),
  make_para("3.2 Pattern of increasing explanatory strength across scale", "Heading2"),
  make_para("The increase in correlation strength from 10 m to 50 m suggests that the relationship between standardized PCA mean distance and abundance-weighted Faith's PD is scale dependent. The 10 m relationship was detectable but modest, the 20 m relationship was intermediate, and the 50 m relationship accounted for the largest fraction of spectral heterogeneity variation. This scale pattern should be presented as the central result for the standardized PCA mean-distance analysis."),
  make_para("[Figure 4 placeholder: scale-comparison summary for the standardized PCA mean distance versus abundance-weighted Faith's PD relationship. Recommended design: bar or point-range plot showing Pearson r and/or R2 by scale, with sample sizes annotated. In text, cite as Fig. 4 when describing the monotonic increase in association strength across scales.]", "Placeholder"),
  make_para("A concise results statement for the manuscript can therefore be built around the following sequence: first, all scale-specific relationships were positive; second, the magnitude of the relationship increased with quadrat grain; third, the 50 m grain produced the largest effect size but had the smallest number of complete observations. This framing keeps the interpretation anchored to the observed correlations while leaving mechanistic explanation for the Discussion."),
  make_para("3.3 Results text to connect with planned figures and tables", "Heading2"),
  make_para("When the figure set is assembled, the Results section should first direct readers to the three-panel scatterplot, then report the numerical model summaries in the table. The figure should carry the visual scale comparison, while the table should carry exact r, R2, F, p-value, and slope values. This division will keep the Results section readable and avoid overloading the prose with repeated statistics."),
  make_para("[Table 3 placeholder: optional compact comparison table for focal biodiversity contrasts. If later included, add Shannon species diversity beside abundance-weighted Faith's PD to show that the abundance-weighted phylogenetic relationship is the clearer signal. Keep this table supplemental unless it becomes central to the Results narrative.]", "Placeholder"),
  make_para("No environmental covariates are included in this Results subsection. The current purpose is to present the direct correlation between standardized PCA mean distance and abundance-weighted Faith's PD across spatial scales. Environmental or spatial sensitivity models should be introduced later only if the manuscript adds a separate robustness subsection or supplement."),
  make_para("Draft result values checked from `reports/tables/multiscale_spectral_biodiversity/sv_diversity_pairwise_correlations.csv`:", "Normal", italic = TRUE),
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
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Results Draft PCA Mean Distance and Abundance-Weighted Faith PD</dc:title><dc:creator>Codex via R</dc:creator><cp:lastModifiedBy>Codex via R</cp:lastModifiedBy><dcterms:created xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:modified></cp:coreProperties>', file.path(tmp_docx, "docProps", "core.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>R base OOXML builder</Application></Properties>', file.path(tmp_docx, "docProps", "app.xml"), useBytes = TRUE)

zipfile <- tempfile(fileext = ".zip")
compress_archive(tmp_docx, zipfile)
file.copy(zipfile, output_docx, overwrite = TRUE)

cat("Created:", output_docx, "\n")
cat("Style source paragraphs read:",
    "Spectral Variation Paper =", length(style_sources$spectral_variation_paper),
    "; SVH Outline =", length(style_sources$svh_outline),
    "; RESULTS-SampleT =", length(style_sources$results_sample), "\n")
cat("Markdown/context lines read:", length(md_context), "\n")
