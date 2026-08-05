options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root, "Documents", "Paper")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

today_stamp <- "20260803"
style_note_path <- file.path(paper_dir, paste0(today_stamp, "_writing_style_notes.md"))
outline_docx_path <- file.path(paper_dir, paste0(today_stamp, "_detailed_manuscript_outline.docx"))
md_context_path <- file.path(paper_dir, paste0(today_stamp, "_markdown_context_digest.md"))

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
  cmd <- paste(
    "Expand-Archive",
    "-LiteralPath", ps_quote(zip_path),
    "-DestinationPath", ps_quote(destination),
    "-Force"
  )
  status <- system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(destination, "word", "document.xml"))) {
    stop("PowerShell archive extraction failed for ", zip_path, "\n", paste(status, collapse = "\n"))
  }
}

compress_archive <- function(source_dir, out_zip) {
  if (file.exists(out_zip)) unlink(out_zip)
  cmd <- paste(
    "Compress-Archive",
    "-Path", ps_quote(file.path(source_dir, "*")),
    "-DestinationPath", ps_quote(out_zip),
    "-Force"
  )
  status <- system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE)
  if (!file.exists(out_zip)) {
    stop("PowerShell archive compression failed.\n", paste(status, collapse = "\n"))
  }
}

read_text <- function(path) {
  if (!file.exists(path)) return(character())
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

extract_docx_text <- function(docx_path) {
  tmp <- tempfile("docx_extract_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  zip_copy <- tempfile(fileext = ".zip")
  copied <- tryCatch(file.copy(docx_path, zip_copy, overwrite = TRUE), error = function(e) FALSE)
  if (!isTRUE(copied)) {
    warning("Could not copy locked or unavailable DOCX: ", docx_path)
    return(character())
  }
  ok <- tryCatch({
    expand_archive(zip_copy, tmp)
    TRUE
  }, error = function(e) {
    warning(conditionMessage(e))
    FALSE
  })
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

md_files <- list.files(root, pattern = "\\.md$", recursive = TRUE, full.names = TRUE)
md_files <- md_files[!grepl("/\\.git/", normalizePath(md_files, winslash = "/", mustWork = FALSE), fixed = TRUE)]
md_rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", root), "/?"), "", normalizePath(md_files, winslash = "/", mustWork = FALSE))
md_headings <- lapply(seq_along(md_files), function(i) {
  lines <- read_text(md_files[i])
  headings <- grep("^#{1,4} ", lines, value = TRUE)
  headings <- head(headings, 12)
  c(paste0("- `", md_rel[i], "`"), if (length(headings)) paste0("  - ", headings) else "  - No Markdown headings detected.")
})

writeLines(c(
  "# Markdown Context Digest",
  "",
  "Created: 2026-08-03",
  "",
  "Workflow note: this digest was produced by `tools/create_paper_outline_docx.R` using R 4.2.3. It records that the paper-outline workflow read the repository Markdown corpus before generating the style notes and Word outline.",
  "",
  paste0("Markdown files read: ", length(md_files)),
  "",
  "## File And Heading Inventory",
  "",
  unlist(md_headings, use.names = FALSE)
), md_context_path, useBytes = TRUE)

previous_paper <- extract_docx_text(file.path(root, "Documents", "Spectral Variation Paper.docx"))
previous_outline <- extract_docx_text(file.path(root, "Documents", "SVH Outline.docx"))

sv_report <- read_text(file.path(root, "reports", "analysis", "20260710_sv_diversity_pairwise_correlations.md"))
project_state <- read_text(file.path(root, "reports", "project_state.md"))
variable_guide <- read_text(file.path(root, "reports", "combined_quadrat_variable_guide.md"))

style_note <- c(
  "# Writing Style Notes From Prior Drafts",
  "",
  "Created: 2026-08-03",
  "",
  "Workflow note: the prior Word documents and repository Markdown context were read with `tools/create_paper_outline_docx.R` using R 4.2.3. Continue using R for Word-output generation in this paper workflow unless a later instruction changes that.",
  "",
  "Source documents used for style context:",
  "",
  paste0("- `Documents/Spectral Variation Paper.docx`: ", if (length(previous_paper)) paste0(length(previous_paper), " paragraphs read.") else "not readable during this pass because the file was locked by another process."),
  paste0("- `Documents/SVH Outline.docx`: ", length(previous_outline), " paragraphs read."),
  "- Repository Markdown context summarized in `Documents/Paper/20260803_markdown_context_digest.md`",
  "",
  "## Observed Writing Style",
  "",
  "- The prior writing is exploratory but hypothesis-driven. It tends to frame the work around a clear research question, then build toward methods and interpretation through staged explanation.",
  "- The prose favors direct scientific framing over ornament: define the ecological problem, introduce spectral variation as the possible proxy, then specify the scale and diversity metrics.",
  "- Section-level organization matters. The existing outline style uses manuscript sections and nested conceptual prompts rather than a continuous essay.",
  "- Methods should be concrete and reproducible. The strongest prior material names data sources, quadrat scales, filtering decisions, spectral processing steps, and expected outputs.",
  "- The voice is cautious about interpretation. Claims should be framed as relationships, correlations, or evidence for/against the spectral variation hypothesis, not as broad proof.",
  "- The writing benefits from explicit bridges between remote sensing and ecology: spectral heterogeneity should be tied back to canopy, species, phylogenetic, and scale concepts.",
  "",
  "## Style To Preserve In The New Outline",
  "",
  "- Use working headings that can later become manuscript headings.",
  "- Keep bullets detailed enough to draft from, but avoid writing final manuscript paragraphs too early.",
  "- State methodological updates plainly and in chronological workflow order.",
  "- Keep technical terms visible: spectral angle entropy, standardized PCA mean distance, spectral Rao's Q, abundance-weighted Faith's PD, Shannon species diversity, 10 m, 20 m, and 50 m quadrats.",
  "- Keep environmental factors out of the main outline for this stage because the current focus is the correlation between species/phylogenetic diversity and spectral diversity across scales.",
  "",
  "## Key Methodological Updates To Carry Forward",
  "",
  "- Species and phylogenetic diversity are now generated at 10 m, 20 m, and 50 m scales using the current all-scale plant-diversity workflow and quadrat IDs aligned with spectral outputs.",
  "- Species diversity emphasis for this outline should be Shannon diversity only.",
  "- Phylogenetic diversity emphasis for this outline should be abundance-weighted Faith's PD only.",
  "- Spectral heterogeneity emphasis should include spectral angle entropy, standardized PCA mean distance, and standardized PCA Rao's Q as the spectral Rao-style measure.",
  "- Spectral angle entropy is calculated from sunlit, shadow-masked, smoothed 5 nm spectra, with most quadrats using bootstrap means rather than exhaustive all-pair calculations.",
  "- Standardized PCA metrics use vector-normalized spectra before projection, reducing broad brightness dominance compared with regular PCA.",
  "- The manuscript question should be narrowed to whether spectral-diversity relationships differ between abundance-weighted Faith's PD and Shannon species diversity, and how these relationships change across 10 m, 20 m, and 50 m scales.",
  "",
  "## Result Emphasis For This Outline",
  "",
  "- Abundance-weighted Faith's PD shows positive relationships with the primary spectral heterogeneity metrics across scales, with the standardized PCA mean distance relationship strongest at 20 m and 50 m.",
  "- Shannon diversity relationships are weak, near zero, or negative in the primary non-edge analysis, making it useful as the species-diversity contrast.",
  "- Scale should be treated as a central result, not a footnote: the phylogenetic relationship becomes more interpretable at broader grains, while species Shannon does not show a comparable pattern.",
  "- Results should be framed as correlation and scale-dependence, not environmental explanation."
)
writeLines(style_note, style_note_path, useBytes = TRUE)

result_rows <- data.frame(
  scale = c("10 m", "10 m", "10 m", "10 m", "20 m", "20 m", "20 m", "20 m", "50 m", "50 m", "50 m", "50 m"),
  spectral_metric = c(
    "Spectral angle entropy", "Spectral angle entropy", "Standardized PCA mean distance", "Standardized PCA mean distance",
    "Spectral angle entropy", "Spectral angle entropy", "Standardized PCA mean distance", "Standardized PCA mean distance",
    "Spectral angle entropy", "Spectral angle entropy", "Standardized PCA mean distance", "Standardized PCA mean distance"
  ),
  diversity_metric = rep(c("Abundance-weighted Faith's PD", "Shannon species diversity"), 6),
  n = c(1587, 1587, 1440, 1440, 405, 405, 360, 360, 80, 80, 74, 74),
  r = c(0.102, -0.051, 0.106, -0.049, 0.166, -0.065, 0.291, -0.075, 0.228, 0.047, 0.419, 0.009),
  r2 = c(0.010, 0.003, 0.011, 0.002, 0.027, 0.004, 0.085, 0.006, 0.052, 0.002, 0.175, 0.000),
  p = c("4.47e-05", "0.043", "5.51e-05", "0.062", "8.12e-04", "0.195", "1.82e-08", "0.155", "0.042", "0.678", "2.05e-04", "0.942")
)

outline <- list(
  list(level = 1, title = "Working Title", body = c(
    "Evaluating scale-dependent relationships between drone hyperspectral spectral heterogeneity, abundance-weighted phylogenetic diversity, and Shannon species diversity in the Paint Rock Forest Dynamics Plot.",
    "Purpose of this document: provide a detailed outline to build upon; it is not a full manuscript draft."
  )),
  list(level = 1, title = "Core Manuscript Question", body = c(
    "Do spectral heterogeneity metrics derived from drone hyperspectral imagery correlate more strongly with abundance-weighted Faith's phylogenetic diversity than with Shannon species diversity?",
    "How does that relationship change across 10 m, 20 m, and 50 m quadrat scales?"
  )),
  list(level = 1, title = "Working Argument", body = c(
    "The paper should argue that spectral heterogeneity shows a clearer relationship with abundance-weighted phylogenetic diversity than with Shannon species diversity, and that this contrast changes with spatial grain.",
    "The strongest current evidence is the positive relationship between standardized PCA mean distance and abundance-weighted Faith's PD, especially at 20 m and 50 m. Spectral angle entropy supports a weaker but directionally similar pattern.",
    "Shannon diversity should serve as the species-diversity contrast because its relationships with the same spectral metrics are weak, near zero, or negative in the primary analysis."
  )),
  list(level = 1, title = "Introduction Outline", body = c(
    "Opening frame: forest biodiversity is difficult to monitor continuously, while hyperspectral imagery offers spatially extensive spectral information that may reflect canopy and compositional variation.",
    "Introduce the spectral variation hypothesis as the conceptual bridge between spectral heterogeneity and biodiversity.",
    "Narrow the biodiversity concept: species diversity and phylogenetic diversity can behave differently because phylogenetic diversity incorporates evolutionary dissimilarity and abundance-weighted lineage structure.",
    "Identify the key contrast for this manuscript: abundance-weighted Faith's PD versus Shannon species diversity.",
    "Introduce the scale problem: a 10 m quadrat may capture fine canopy texture and individual crown structure, whereas 20 m and 50 m quadrats integrate broader community and canopy composition.",
    "End the introduction with two focused questions: first, whether spectral heterogeneity correlates more strongly with abundance-weighted Faith's PD than Shannon diversity; second, whether the strength or direction of those relationships changes with scale."
  )),
  list(level = 1, title = "Methods Outline", body = character()),
  list(level = 2, title = "Study System And Spatial Design", body = c(
    "Describe the Paint Rock Forest Dynamics Plot and its mapped forest-census context.",
    "Define the quadrat framework: 10 m, 20 m, and 50 m analysis grains from the quadrat shapefiles.",
    "State that quadrat IDs were harmonized across spectral, species-diversity, and phylogenetic-diversity products before analysis."
  )),
  list(level = 2, title = "Species And Phylogenetic Diversity Workflow", body = c(
    "Use the current all-scale plant-diversity workflow as the methods basis.",
    "Tree census data were filtered to the canopy-focused analysis population using the current workflow logic, then summarized within quadrats.",
    "Species composition was represented by summed crown-overlap proportions by species and quadrat.",
    "Shannon species diversity was calculated from positive species crown-overlap proportions and retained as the species-diversity comparison metric.",
    "Abundance-weighted Faith's PD was calculated by weighting phylogenetic branch lengths by species crown-overlap values and retained as the focal phylogenetic-diversity metric.",
    "Mention other diversity products only as workflow context if needed; do not let richness, Simpson, evenness, Faith's PD, or phylogenetic Rao's Q become main manuscript outcomes in this outline."
  )),
  list(level = 2, title = "Spectral Data Processing Workflow", body = c(
    "Use confirmed quadrat spectra in `Quad_Spectra/10m`, `Quad_Spectra/20m`, and `Quad_Spectra/50m` as the basis.",
    "State that spectra were smoothed with the current Savitzky-Golay workflow and resampled to 5 nm spacing, producing `Quad_Spectra/*_smooth_5nm` inputs.",
    "Apply the current sunlit/shadow mask using the 563 nm reflectance threshold of 0.0305476, retaining pixels greater than the threshold.",
    "Describe manual atmospheric/cloud exclusions for PCA-dependent metrics only as a data-quality step."
  )),
  list(level = 2, title = "Spectral Heterogeneity Metrics", body = c(
    "Spectral angle entropy (`spec_sa`): calculate angular heterogeneity from sunlit, shadow-masked spectra. Most quadrats use the mean of 70 bootstrap iterations with up to 5,000 retained pixels per iteration; only small quadrats below the pair-count threshold use exact all-pixel pairwise angles.",
    "Standardized PCA mean distance (`spec_spca_mean`): vector-normalize retained spectra, project into the standardized PCA basis, and summarize the mean Euclidean distance of retained pixels from the quadrat centroid in PC1-PC3 space.",
    "Standardized PCA Rao's Q (`spec_spca_rao`): use equal pixel weights and squared Euclidean distances in standardized PCA space as the spectral Rao-style heterogeneity measure.",
    "Keep raw PCA, alpha-hull, convex hull, and environmental variables out of the main outline unless they become supplemental methods later."
  )),
  list(level = 2, title = "Correlation Analysis", body = c(
    "Analyze each scale separately: 10 m, 20 m, and 50 m.",
    "Fit direct pairwise relationships between each focal spectral heterogeneity metric and each focal biodiversity metric.",
    "Focal biodiversity metrics: abundance-weighted Faith's PD and Shannon species diversity.",
    "Focal spectral metrics: spectral angle entropy, standardized PCA mean distance, and standardized PCA Rao's Q.",
    "Report Pearson r, R2, F-test p-value, sample size, slope, intercept, and Spearman rank correlation where available.",
    "Exclude documented 10 m and 20 m edge quadrats from primary inference; retain all 50 m quadrats because no separate 50 m edge rule is documented.",
    "Do not include environmental factors in the main analysis for this manuscript stage."
  )),
  list(level = 1, title = "Results Outline", body = character()),
  list(level = 2, title = "Analysis Coverage By Scale", body = c(
    "Report total and primary-analysis quadrat counts by scale.",
    "From the current primary dataset: 10 m has 1,656 primary non-edge quadrats, 20 m has 414 primary non-edge quadrats, and 50 m has 80 quadrats.",
    "Report metric-specific complete cases: SA entropy has 1,587, 405, and 80 complete primary values at 10 m, 20 m, and 50 m; standardized PCA mean distance has 1,440, 360, and 74 complete primary values."
  )),
  list(level = 2, title = "Abundance-Weighted Faith's PD Relationships", body = c(
    "Lead with standardized PCA mean distance because it gives the clearest focal phylogenetic signal.",
    "Standardized PCA mean distance versus abundance-weighted Faith's PD: r = 0.106 at 10 m, r = 0.291 at 20 m, and r = 0.419 at 50 m.",
    "Spectral angle entropy versus abundance-weighted Faith's PD: r = 0.102 at 10 m, r = 0.166 at 20 m, and r = 0.228 at 50 m.",
    "Interpretation to develop later: both focal spectral metrics show positive phylogenetic relationships, with the standardized PCA metric strengthening more clearly across scale."
  )),
  list(level = 2, title = "Shannon Species Diversity Relationships", body = c(
    "Use Shannon diversity as the contrast case.",
    "Standardized PCA mean distance versus Shannon diversity: r = -0.049 at 10 m, r = -0.075 at 20 m, and r = 0.009 at 50 m.",
    "Spectral angle entropy versus Shannon diversity: r = -0.051 at 10 m, r = -0.065 at 20 m, and r = 0.047 at 50 m.",
    "Interpretation to develop later: Shannon diversity does not show the same positive scale-dependent relationship as abundance-weighted Faith's PD."
  )),
  list(level = 2, title = "Spectral Rao's Q As A Supporting Spectral Heterogeneity Measure", body = c(
    "Add standardized PCA Rao's Q as a supporting spectral heterogeneity metric once its focal pairwise results are pulled from the current tables.",
    "Use this subsection to ask whether a Rao-style spectral dispersion measure tells the same story as standardized PCA mean distance.",
    "Keep this as supporting evidence unless the final table shows it materially changes the interpretation."
  )),
  list(level = 2, title = "Scale Dependence", body = c(
    "Organize the scale result around the contrast between abundance-weighted Faith's PD and Shannon diversity.",
    "The phylogenetic relationship strengthens from 10 m toward broader grains, especially for standardized PCA mean distance.",
    "The Shannon relationship remains weak and does not show a comparable positive scale pattern.",
    "This section should become the manuscript's main results narrative rather than a long list of pairwise tests."
  )),
  list(level = 1, title = "Discussion Outline", body = c(
    "Start with the central contrast: spectral heterogeneity aligns more clearly with abundance-weighted phylogenetic diversity than with Shannon species diversity.",
    "Discuss why abundance-weighted phylogenetic diversity may better match canopy spectral variation: it combines abundance structure with evolutionary dissimilarity, which may better capture trait and canopy differences that affect spectra.",
    "Discuss why Shannon diversity may be weak: species diversity can be high without strong spectral separation if species are spectrally similar, occur below the dominant canopy, or vary in abundance in ways not visible to the sensor.",
    "Discuss scale: coarser grains may integrate more canopy community structure, while fine grains may be noisier because they are influenced by individual crowns, shadows, and local spectral variability.",
    "Keep limitations focused on correlation, sample size by scale, edge/missing quadrats, spectral masking, and spectral metric uncertainty. Do not broaden into environmental explanations in the main outline.",
    "End with a practical implication: drone hyperspectral spectral heterogeneity may be more useful for detecting phylogenetic structure than conventional species-diversity summaries, but the relationship is scale-dependent and modest."
  )),
  list(level = 1, title = "Figures And Tables To Build", body = c(
    "Figure 1: study system and multiscale quadrat design.",
    "Figure 2: methods workflow from tree census and hyperspectral imagery to diversity and spectral heterogeneity metrics.",
    "Figure 3: scale-by-scale scatterplots for abundance-weighted Faith's PD against spectral angle entropy and standardized PCA mean distance.",
    "Figure 4: scale-by-scale scatterplots for Shannon diversity against the same spectral metrics.",
    "Figure 5 or supplement: standardized PCA Rao's Q relationships if they clarify the spectral Rao-style result.",
    "Table 1: focal variables and their calculation notes.",
    "Table 2: primary pairwise results by scale, spectral metric, and diversity metric.",
    "Supplementary table: coverage and missingness by scale."
  )),
  list(level = 1, title = "Immediate Gaps Before Drafting", body = c(
    "Pull or calculate the focused standardized PCA Rao's Q pairwise correlations against abundance-weighted Faith's PD and Shannon diversity.",
    "Decide whether the main text should use only non-edge quadrats or show all quadrats with edge points flagged in figures.",
    "Decide whether 50 m should be framed as the strongest result or as a coarse-scale endpoint with lower sample size.",
    "Choose target journal before converting this outline into full prose."
  ))
)

make_para <- function(text, style = "Normal", bold = FALSE, italic = FALSE, size = NULL, color = NULL) {
  props <- c()
  if (bold) props <- c(props, "<w:b/>")
  if (italic) props <- c(props, "<w:i/>")
  if (!is.null(color)) props <- c(props, paste0("<w:color w:val=\"", color, "\"/>"))
  if (!is.null(size)) props <- c(props, paste0("<w:sz w:val=\"", as.integer(size * 2), "\"/>"))
  rpr <- if (length(props)) paste0("<w:rPr>", paste(props, collapse = ""), "</w:rPr>") else ""
  paste0("<w:p><w:pPr><w:pStyle w:val=\"", style, "\"/></w:pPr><w:r>", rpr, "<w:t xml:space=\"preserve\">", xml_escape(text), "</w:t></w:r></w:p>")
}

make_bullet <- function(text) {
  paste0("<w:p><w:pPr><w:pStyle w:val=\"ListParagraph\"/><w:numPr><w:ilvl w:val=\"0\"/><w:numId w:val=\"1\"/></w:numPr></w:pPr><w:r><w:t xml:space=\"preserve\">", xml_escape(text), "</w:t></w:r></w:p>")
}

make_table <- function(df) {
  widths <- c(900, 2400, 2500, 700, 700, 700, 900)
  rows <- list()
  header <- names(df)
  cell <- function(value, width, header = FALSE) {
    fill <- if (header) "<w:shd w:fill=\"F2F4F7\"/>" else ""
    bold <- if (header) "<w:rPr><w:b/></w:rPr>" else ""
    paste0("<w:tc><w:tcPr><w:tcW w:w=\"", width, "\" w:type=\"dxa\"/>", fill, "</w:tcPr><w:p><w:r>", bold, "<w:t xml:space=\"preserve\">", xml_escape(as.character(value)), "</w:t></w:r></w:p></w:tc>")
  }
  rows[[1]] <- paste0("<w:tr>", paste(mapply(cell, header, widths, MoreArgs = list(header = TRUE)), collapse = ""), "</w:tr>")
  for (i in seq_len(nrow(df))) {
    vals <- unname(unlist(df[i, ], use.names = FALSE))
    rows[[length(rows) + 1]] <- paste0("<w:tr>", paste(mapply(cell, vals, widths, MoreArgs = list(header = FALSE)), collapse = ""), "</w:tr>")
  }
  paste0(
    "<w:tbl><w:tblPr><w:tblW w:w=\"9360\" w:type=\"dxa\"/><w:tblBorders>",
    "<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>",
    "<w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>",
    "<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>",
    "<w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"BFBFBF\"/>",
    "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/>",
    "<w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"D9D9D9\"/>",
    "</w:tblBorders><w:tblCellMar><w:top w:w=\"80\" w:type=\"dxa\"/><w:bottom w:w=\"80\" w:type=\"dxa\"/><w:start w:w=\"120\" w:type=\"dxa\"/><w:end w:w=\"120\" w:type=\"dxa\"/></w:tblCellMar></w:tblPr>",
    "<w:tblGrid>", paste0("<w:gridCol w:w=\"", widths, "\"/>", collapse = ""), "</w:tblGrid>",
    paste(rows, collapse = ""),
    "</w:tbl>"
  )
}

styles_xml <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:pPr><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:pPr><w:spacing w:after="160"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:b/><w:color w:val="0B2545"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:pPr><w:spacing w:after="200"/></w:pPr><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri"/><w:color w:val="555555"/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="320" w:after="160"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="120"/></w:pPr><w:rPr><w:b/><w:color w:val="2E74B5"/><w:sz w:val="26"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="ListParagraph"><w:name w:val="List Paragraph"/><w:basedOn w:val="Normal"/><w:pPr><w:ind w:left="720" w:hanging="360"/><w:spacing w:after="80" w:line="280" w:lineRule="auto"/></w:pPr></w:style>
</w:styles>'

numbering_xml <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:numbering xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:abstractNum w:abstractNumId="0"><w:lvl w:ilvl="0"><w:start w:val="1"/><w:numFmt w:val="bullet"/><w:lvlText w:val="•"/><w:lvlJc w:val="left"/><w:pPr><w:tabs><w:tab w:val="num" w:pos="720"/></w:tabs><w:ind w:left="720" w:hanging="360"/></w:pPr></w:lvl></w:abstractNum>
<w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num>
</w:numbering>'

content <- c(
  make_para("Detailed Manuscript Outline", "Title"),
  make_para("Spectral heterogeneity, abundance-weighted phylogenetic diversity, and Shannon diversity across scale", "Subtitle"),
  make_para("Created: 2026-08-03. Generated with R 4.2.3 from prior Word drafts, repository Markdown context, and current analysis reports.", "Normal", italic = TRUE)
)

for (item in outline) {
  style <- if (item$level == 1) "Heading1" else "Heading2"
  content <- c(content, make_para(item$title, style))
  for (line in item$body) content <- c(content, make_bullet(line))
  if (identical(item$title, "Results Outline")) {
    content <- c(content, make_para("Focused primary pairwise results currently available for the two requested biodiversity metrics:", "Normal", italic = TRUE), make_table(result_rows))
  }
}

document_xml <- paste0(
  '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
  '<w:document xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas" xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math" xmlns:v="urn:schemas-microsoft-com:vml" xmlns:wp14="http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing" xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" xmlns:w10="urn:schemas-microsoft-com:office:word" xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml" xmlns:wpg="http://schemas.microsoft.com/office/word/2010/wordprocessingGroup" xmlns:wpi="http://schemas.microsoft.com/office/word/2010/wordprocessingInk" xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml" xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape" mc:Ignorable="w14 wp14"><w:body>',
  paste(content, collapse = ""),
  '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/></w:sectPr>',
  '</w:body></w:document>'
)

tmp_docx <- tempfile("docx_build_")
dir.create(file.path(tmp_docx, "_rels"), recursive = TRUE)
dir.create(file.path(tmp_docx, "word", "_rels"), recursive = TRUE)
dir.create(file.path(tmp_docx, "docProps"), recursive = TRUE)

writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/><Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/></Types>', file.path(tmp_docx, "[Content_Types].xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/><Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/></Relationships>', file.path(tmp_docx, "_rels", ".rels"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering" Target="numbering.xml"/></Relationships>', file.path(tmp_docx, "word", "_rels", "document.xml.rels"), useBytes = TRUE)
writeLines(document_xml, file.path(tmp_docx, "word", "document.xml"), useBytes = TRUE)
writeLines(styles_xml, file.path(tmp_docx, "word", "styles.xml"), useBytes = TRUE)
writeLines(numbering_xml, file.path(tmp_docx, "word", "numbering.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Detailed Manuscript Outline</dc:title><dc:creator>Codex via R</dc:creator><cp:lastModifiedBy>Codex via R</cp:lastModifiedBy><dcterms:created xsi:type="dcterms:W3CDTF">2026-08-03T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-03T00:00:00Z</dcterms:modified></cp:coreProperties>', file.path(tmp_docx, "docProps", "core.xml"), useBytes = TRUE)
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>R base OOXML builder</Application></Properties>', file.path(tmp_docx, "docProps", "app.xml"), useBytes = TRUE)

zipfile <- tempfile(fileext = ".zip")
compress_archive(tmp_docx, zipfile)
file.copy(zipfile, outline_docx_path, overwrite = TRUE)

cat("Created style note:", style_note_path, "\n")
cat("Created Markdown context digest:", md_context_path, "\n")
cat("Created outline DOCX:", outline_docx_path, "\n")
cat("Prior Word paragraphs read:", length(previous_paper), "from Spectral Variation Paper;", length(previous_outline), "from SVH Outline\n")
cat("Markdown files read:", length(md_files), "\n")
