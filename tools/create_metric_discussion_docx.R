options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root, "Documents", "Paper")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

today_stamp <- "20260805"
output_docx <- file.path(paper_dir, paste0(today_stamp, "_discussion_pca_distance_vs_rao_q_scale_centroid_clarified.docx"))

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
  file.path(root, "reports", "analysis", "20260707_standardized_pca_metric_sample_size_effects.md"),
  file.path(root, "Documents", "Paper", "20260730_paper_workspace_index.md"),
  file.path(root, "Documents", "Paper", "20260803_writing_style_notes.md"),
  file.path(root, "Documents", "Paper", "20260803_remote_sensing_mdpi_author_requirements.md")
)
md_context <- unlist(lapply(md_files[file.exists(md_files)], readLines, warn = FALSE), use.names = FALSE)

pairs <- read.csv(file.path(root, "reports", "tables", "multiscale_spectral_biodiversity", "sv_diversity_pairwise_correlations.csv"))
pca_afaith <- subset(pairs, sv_measure == "spec_spca_mean" & diversity_measure == "phy_afaith")
pca_afaith <- pca_afaith[match(c("10m", "20m", "50m"), pca_afaith$scale), ]

analysis_dataset <- read.csv(file.path(root, "reports", "tables", "multiscale_spectral_biodiversity", "analysis_dataset_with_flags.csv"))
analysis_dataset <- subset(analysis_dataset, primary_analysis == TRUE)

fit_rao <- function(scale_name) {
  d <- subset(analysis_dataset, scale == scale_name & !is.na(spec_rao_q) & !is.na(phy_afaith))
  m <- lm(spec_rao_q ~ phy_afaith, data = d)
  data.frame(
    scale = scale_name,
    n = nrow(d),
    pearson_r = unname(cor(d$spec_rao_q, d$phy_afaith)),
    r_squared = unname(summary(m)$r.squared),
    slope = unname(coef(m)[["phy_afaith"]])
  )
}
rao_afaith <- do.call(rbind, lapply(c("10m", "20m", "50m"), fit_rao))

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

pca_summary <- paste0(
  "For abundance-weighted Faith's PD, standardized PCA mean distance increased from r = ",
  sprintf("%.3f", pca_afaith$pearson_r[1]), " at 10 m to r = ",
  sprintf("%.3f", pca_afaith$pearson_r[2]), " at 20 m and r = ",
  sprintf("%.3f", pca_afaith$pearson_r[3]), " at 50 m."
)

rao_summary <- paste0(
  "By comparison, spectral Rao's Q increased from r = ",
  sprintf("%.3f", rao_afaith$pearson_r[1]), " at 10 m to r = ",
  sprintf("%.3f", rao_afaith$pearson_r[2]), " at 20 m and r = ",
  sprintf("%.3f", rao_afaith$pearson_r[3]), " at 50 m."
)

`%+%` <- paste0

body <- c(
  make_para("Discussion Draft: Why PCA Mean Distance May Outperform Spectral Rao's Q", "Title"),
  make_para("Created 2026-08-05; manuscript discussion prose for the Remote Sensing paper workspace.", "Subtitle"),
  make_para("Context note: this draft expands the user's conceptual notes using the current project reports and result tables. `SVH_Results.docx` was treated as read-only context and was not edited.", "Note"),
  make_para("4. Discussion", "Heading1"),
  make_para("4.1 Why the centroid-based metric may align more strongly with abundance-weighted phylogenetic diversity", "Heading2"),
  make_para("The stronger performance of standardized PCA mean Euclidean distance suggests that the spectral signal associated with abundance-weighted Faith's phylogenetic diversity is expressed more as broad dispersion around a quadrat-level spectral center than as the full distribution of pairwise spectral differences among pixels. The mean-distance metric asks how atypical, on average, pixels are relative to the quadrat's own spectral centroid. In this sense, it summarizes the average departure of canopy pixels from the local spectral condition of the quadrat. Because the spectra were vector-normalized before the standardized PCA basis was fit, this metric is also less directly dominated by overall brightness than a raw reflectance-space distance would be."),
  make_para("Spectral Rao's Q is conceptually closer to a formal diversity measure because it represents average pairwise dissimilarity among spectral units. In the current workflow, however, the script does not explicitly construct all pixel-pair distances. It calculates the quadrat centroid in global PC1-PC3 space, computes each retained pixel's squared Euclidean distance from that centroid, and stores spectral Rao's Q as `2 * mean(squared_radius)`. With equal pixel weights, this centroid-based expression is mathematically equivalent to the mean pairwise squared Euclidean distance among pixels. Thus, the metric should be described as equal-weight pairwise squared Euclidean Rao's Q computed through the equivalent centroid formula."),
  make_para("This implementation detail is important for interpretation. The mean Euclidean distance metric uses the average unsquared distance from each pixel to the quadrat centroid, whereas spectral Rao's Q uses twice the average squared distance from that same centroid. Both metrics therefore use the quadrat centroid computationally, but they summarize dispersion differently: mean distance preserves the original Euclidean distance scale, while Rao's Q emphasizes variance-like squared dispersion that is equivalent to pairwise squared dissimilarity under equal weights."),
  make_para("Because Rao's Q is based on squared distances, rare but extreme pixels can contribute disproportionately to the final value even though the script uses the centroid shortcut. Those extreme pixels may represent meaningful canopy differences, but they may also reflect residual shadow, sunlit crown edges, gaps, mixed pixels, registration mismatch, or atmospheric artifacts."),
  make_para("This distinction matters because abundance-weighted Faith's PD is itself a community-level summary that gives more influence to dominant taxa than to rare taxa. If the biodiversity signal is most strongly carried by the dominant canopy assemblage, a centroid-based measure may be better matched to the ecological weighting of the field metric. Mean Euclidean distance can capture the dominant spread of spectral conditions without allowing a relatively small number of unusual pixels to dominate the quadrat score. Spectral Rao's Q may remain more faithful to the concept of spectral diversity in the abstract, but it may be less tightly matched to the specific ecological and sampling structure of abundance-weighted Faith's PD."),
  make_para("4.2 Why a conceptually richer spectral diversity metric may be noisier", "Heading2"),
  make_para("The contrast between the two metrics is therefore not a failure of Rao's Q as a diversity concept. Rather, it suggests that conceptual richness can come with greater sensitivity to sources of variation that are not the target ecological signal. Although the current script computes Rao's Q through the equal-weight centroid identity, the value still corresponds to pairwise squared Euclidean dissimilarity. That can be valuable when the research question is about total internal spectral differentiation, but it can be problematic when the target response is a field-derived community index that is weighted by species abundance and measured at the quadrat scale."),
  make_para("In a forest canopy, pixels are not independent biological individuals. They are a mixture of leaf angles, crown positions, sunlit and shaded surfaces, sub-canopy gaps, branch and background exposure, and spatial resampling artifacts. A pairwise squared-distance metric can amplify these non-taxonomic contrasts. Mean Euclidean distance from the quadrat centroid is not immune to these effects, but it compresses them into a simpler summary of overall dispersion. This may explain why standardized PCA mean distance produced stronger abundance-weighted Faith correlations than spectral Rao's Q, particularly at the 20 m and 50 m scales."),
  make_para("4.3 Evolutionary and physiological constraints on the spectral-phylogenetic relationship", "Heading2"),
  make_para("The positive but relatively shallow relationships also suggest that spectral diversity should not be expected to increase linearly with phylogenetic diversity without limit. Tree lineages differ in leaf chemistry, structure, water content, pigment concentrations, phenology, and canopy architecture, all of which can influence reflectance. However, leaves are also constrained by shared physiological requirements. Photosynthetic tissues must absorb and manage solar radiation, maintain water balance, and operate within biochemical and structural limits. These shared constraints may place an upper bound on how different leaf and canopy spectra can become, even as species become more phylogenetically distinct."),
  make_para("This provides one possible interpretation for the relatively flat but positive slopes. Quadrats with greater abundance-weighted phylogenetic breadth may contain more spectrally differentiated canopy assemblages, but the spectral space occupied by trees is not unlimited. As taxa branch evolutionarily, some traits may diverge while others converge because all species must solve similar problems of light capture, photoprotection, water regulation, and canopy deployment. The result could be a detectable increase in spectral heterogeneity with phylogenetic diversity, but one that is weaker than a one-to-one mapping between phylogenetic distance and spectral distance would imply."),
  make_para("This interpretation should be framed cautiously. The current data support a positive spectral-phylogenetic association, but they do not by themselves identify which leaf or canopy traits are responsible. A useful discussion point is that phylogenetic diversity may be linked to spectral heterogeneity through a combination of conserved lineage differences, convergent photosynthetic constraints, and canopy-structural expression of dominant species. The observed relationships are consistent with this framework, but trait measurements or species-level spectral signatures would be needed to separate these mechanisms directly."),
  make_para("4.4 Why correlations increase with quadrat scale", "Heading2"),
  make_para("The increase in correlation strength with quadrat size indicates that the spectral-phylogenetic relationship becomes clearer as the spatial grain better matches the scale at which canopy composition and abundance-weighted phylogenetic structure are expressed. At 10 m, spectral heterogeneity can be strongly affected by individual crowns, crown margins, within-crown variation, small gaps, local shadow geometry, and positional mismatch between field quadrats and remotely sensed pixels. These fine-grain sources of variation can obscure the community-level relationship between spectral heterogeneity and phylogenetic diversity."),
  make_para("As quadrats expand to 20 m and 50 m, the spectral metric integrates across more canopy individuals and a larger portion of the local species assemblage. This aggregation can average out idiosyncratic crown-level noise while retaining differences associated with dominant species composition. Larger quadrats are also more likely to capture the mixture of lineages represented in abundance-weighted Faith's PD, which may explain why the strongest standardized PCA mean-distance relationship occurred at 50 m. " %+% pca_summary),
  make_para("The spectral Rao's Q results point in the same broad direction, although more weakly. " %+% rao_summary %+% " This reinforces the idea that scale is central to the spectral diversity hypothesis in this system. The metric does not merely need to be theoretically appropriate; it must also be calculated at a spatial grain that corresponds to the ecological organization of the canopy and the biodiversity measure being tested."),
  make_para("4.5 Implications for manuscript framing", "Heading2"),
  make_para("A strong way to frame the discussion is that standardized PCA mean Euclidean distance may be the more useful operational metric for detecting abundance-weighted phylogenetic structure, while spectral Rao's Q remains an important conceptual comparison. Mean distance appears to recover the broad canopy-level dispersion associated with dominant phylogenetic composition, whereas Rao's Q summarizes a variance-like squared dispersion that is equivalent to equal-weight pairwise squared dissimilarity. This does not mean Rao's Q is wrong; it means that the best metric depends on whether the goal is formal spectral diversity, sensitivity to all within-quadrat dissimilarities, or robust detection of a biodiversity signal."),
  make_para("The scale pattern should be treated as one of the central findings. The relationship between spectral heterogeneity and abundance-weighted Faith's PD is not fixed across grain size; it strengthens as the quadrat becomes large enough to represent canopy assemblage structure more reliably. This supports a multiscale interpretation of the spectral diversity hypothesis and argues against treating a single quadrat grain as universally optimal. In the final paper, the discussion can emphasize that metric choice and spatial scale interact: a stable centroid-based metric at broader grain may be better aligned with abundance-weighted phylogenetic diversity than a pairwise diversity metric at fine grain."),
  make_para("Potential caveat for final editing: environmental and spatial sensitivity results should be discussed separately if included. The direct metric comparison here focuses on the bivariate spectral-phylogenetic relationship, while the broader manuscript can later ask whether topography, spatial autocorrelation, or edge effects modify that relationship.", "Note"),
  make_para("Current result anchors used in this draft:", "Heading2"),
  make_para(pca_summary, "Note"),
  make_para(rao_summary, "Note")
)

body <- c(
  make_para("Discussion Draft: Why PCA Mean Distance May Outperform Spectral Rao's Q", "Title"),
  make_para("Created 2026-08-05; manuscript discussion prose for the Remote Sensing paper workspace.", "Subtitle"),
  make_para("Context note: this draft expands the user's conceptual notes using the current project reports and result tables. `SVH_Results.docx` was treated as read-only context and was not edited.", "Note"),
  make_para("4. Discussion", "Heading1"),
  make_para("4.1 Why the centroid-based metric may align more strongly with abundance-weighted phylogenetic diversity", "Heading2"),
  make_para("The stronger performance of standardized PCA mean Euclidean distance suggests that the spectral signal associated with abundance-weighted Faith's phylogenetic diversity is expressed more as broad dispersion around a quadrat-level spectral center than as the full distribution of pairwise spectral differences among pixels. The mean-distance metric asks how atypical, on average, pixels are relative to the quadrat's own spectral centroid. In this sense, it summarizes the average departure of canopy pixels from the local spectral condition of the quadrat. Because the spectra were vector-normalized before the standardized PCA basis was fit, this metric is also less directly dominated by overall brightness than a raw reflectance-space distance would be."),
  make_para("Spectral Rao's Q is conceptually closer to a formal diversity measure because it represents average pairwise dissimilarity among spectral units. In the current workflow, however, the script does not explicitly construct all pixel-pair distances. It calculates the quadrat centroid in global PC1-PC3 space, computes each retained pixel's squared Euclidean distance from that centroid, and stores spectral Rao's Q as `2 * mean(squared_radius)`. With equal pixel weights, this centroid-based expression is mathematically equivalent to the mean pairwise squared Euclidean distance among pixels. Thus, the metric should be described as equal-weight pairwise squared Euclidean Rao's Q computed through the equivalent centroid formula."),
  make_para("This implementation detail is important for interpretation. The mean Euclidean distance metric uses the average unsquared distance from each pixel to the quadrat centroid, whereas spectral Rao's Q uses twice the average squared distance from that same centroid. Both metrics therefore use the quadrat centroid computationally, but they summarize dispersion differently: mean distance preserves the original Euclidean distance scale, while Rao's Q emphasizes variance-like squared dispersion that is equivalent to pairwise squared dissimilarity under equal weights."),
  make_para("Because Rao's Q is based on squared distances, rare but extreme pixels can contribute disproportionately to the final value even though the script uses the centroid shortcut. Those extreme pixels may represent meaningful canopy differences, but they may also reflect residual shadow, sunlit crown edges, gaps, mixed pixels, registration mismatch, or atmospheric artifacts."),
  make_para("This distinction matters because abundance-weighted Faith's PD is itself a community-level summary that gives more influence to dominant taxa than to rare taxa. If the biodiversity signal is most strongly carried by the dominant canopy assemblage, a centroid-based measure may be better matched to the ecological weighting of the field metric. Mean Euclidean distance can capture the dominant spread of spectral conditions without allowing a relatively small number of unusual pixels to dominate the quadrat score. Spectral Rao's Q may remain more faithful to the concept of spectral diversity in the abstract, but it may be less tightly matched to the specific ecological and sampling structure of abundance-weighted Faith's PD."),
  make_para("4.2 Why a conceptually richer spectral diversity metric may be noisier", "Heading2"),
  make_para("The contrast between the two metrics is therefore not a failure of Rao's Q as a diversity concept. Rather, it suggests that conceptual richness can come with greater sensitivity to sources of variation that are not the target ecological signal. Although the current script computes Rao's Q through the equal-weight centroid identity, the value still corresponds to pairwise squared Euclidean dissimilarity. That can be valuable when the research question is about total internal spectral differentiation, but it can be problematic when the target response is a field-derived community index that is weighted by species abundance and measured at the quadrat scale."),
  make_para("In a forest canopy, pixels are not independent biological individuals. They are a mixture of leaf angles, crown positions, sunlit and shaded surfaces, sub-canopy gaps, branch and background exposure, and spatial resampling artifacts. A pairwise squared-distance metric can amplify these non-taxonomic contrasts. Mean Euclidean distance from the quadrat centroid is not immune to these effects, but it compresses them into a simpler summary of overall dispersion. This may explain why standardized PCA mean distance produced stronger abundance-weighted Faith correlations than spectral Rao's Q, particularly at the 20 m and 50 m scales."),
  make_para("4.3 Evolutionary and physiological constraints on the spectral-phylogenetic relationship", "Heading2"),
  make_para("The positive but relatively shallow relationships also suggest that spectral diversity should not be expected to increase linearly with phylogenetic diversity without limit. Tree lineages differ in leaf chemistry, structure, water content, pigment concentrations, phenology, and canopy architecture, all of which can influence reflectance. However, leaves are also constrained by shared physiological requirements. Photosynthetic tissues must absorb and manage solar radiation, maintain water balance, and operate within biochemical and structural limits. These shared constraints may place an upper bound on how different leaf and canopy spectra can become, even as species become more phylogenetically distinct."),
  make_para("This provides one possible interpretation for the relatively flat but positive slopes. Quadrats with greater abundance-weighted phylogenetic breadth may contain more spectrally differentiated canopy assemblages, but the spectral space occupied by trees is not unlimited. As taxa branch evolutionarily, some traits may diverge while others converge because all species must solve similar problems of light capture, photoprotection, water regulation, and canopy deployment. The result could be a detectable increase in spectral heterogeneity with phylogenetic diversity, but one that is weaker than a one-to-one mapping between phylogenetic distance and spectral distance would imply."),
  make_para("This interpretation should be framed cautiously. The current data support a positive spectral-phylogenetic association, but they do not by themselves identify which leaf or canopy traits are responsible. A useful discussion point is that phylogenetic diversity may be linked to spectral heterogeneity through a combination of conserved lineage differences, convergent photosynthetic constraints, and canopy-structural expression of dominant species. The observed relationships are consistent with this framework, but trait measurements or species-level spectral signatures would be needed to separate these mechanisms directly."),
  make_para("4.4 Why correlations increase with quadrat scale", "Heading2"),
  make_para("The increase in correlation strength with quadrat size indicates that the spectral-phylogenetic relationship becomes clearer as the spatial grain better matches the scale at which canopy composition and abundance-weighted phylogenetic structure are expressed. At 10 m, spectral heterogeneity can be strongly affected by individual crowns, crown margins, within-crown variation, small gaps, local shadow geometry, and positional mismatch between field quadrats and remotely sensed pixels. These fine-grain sources of variation can obscure the community-level relationship between spectral heterogeneity and phylogenetic diversity."),
  make_para(paste0("As quadrats expand to 20 m and 50 m, the spectral metric integrates across more canopy individuals and a larger portion of the local species assemblage. This aggregation can average out idiosyncratic crown-level noise while retaining differences associated with dominant species composition. Larger quadrats are also more likely to capture the mixture of lineages represented in abundance-weighted Faith's PD, which may explain why the strongest standardized PCA mean-distance relationship occurred at 50 m. ", pca_summary)),
  make_para(paste0("The spectral Rao's Q results point in the same broad direction, although more weakly. ", rao_summary, " This reinforces the idea that scale is central to the spectral diversity hypothesis in this system. The metric does not merely need to be theoretically appropriate; it must also be calculated at a spatial grain that corresponds to the ecological organization of the canopy and the biodiversity measure being tested.")),
  make_para("4.5 Implications for manuscript framing", "Heading2"),
  make_para("A strong way to frame the discussion is that standardized PCA mean Euclidean distance may be the more useful operational metric for detecting abundance-weighted phylogenetic structure, while spectral Rao's Q remains an important conceptual comparison. Mean distance appears to recover the broad canopy-level dispersion associated with dominant phylogenetic composition, whereas Rao's Q summarizes a variance-like squared dispersion that is equivalent to equal-weight pairwise squared dissimilarity. This does not mean Rao's Q is wrong; it means that the best metric depends on whether the goal is formal spectral diversity, sensitivity to all within-quadrat dissimilarities, or robust detection of a biodiversity signal."),
  make_para("The scale pattern should be treated as one of the central findings. The relationship between spectral heterogeneity and abundance-weighted Faith's PD is not fixed across grain size; it strengthens as the quadrat becomes large enough to represent canopy assemblage structure more reliably. This supports a multiscale interpretation of the spectral diversity hypothesis and argues against treating a single quadrat grain as universally optimal. In the final paper, the discussion can emphasize that metric choice and spatial scale interact: a stable centroid-based metric at broader grain may be better aligned with abundance-weighted phylogenetic diversity than a pairwise diversity metric at fine grain."),
  make_para("Potential caveat for final editing: environmental and spatial sensitivity results should be discussed separately if included. The direct metric comparison here focuses on the bivariate spectral-phylogenetic relationship, while the broader manuscript can later ask whether topography, spatial autocorrelation, or edge effects modify that relationship.", "Note"),
  make_para("Current result anchors used in this draft", "Heading2"),
  make_para(pca_summary, "Note"),
  make_para(rao_summary, "Note")
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
writeLines('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>Discussion Draft PCA Distance Versus Spectral Rao Q</dc:title><dc:creator>Codex via R</dc:creator><cp:lastModifiedBy>Codex via R</cp:lastModifiedBy><dcterms:created xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-05T00:00:00Z</dcterms:modified></cp:coreProperties>', file.path(tmp_docx, "docProps", "core.xml"), useBytes = TRUE)
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
print(pca_afaith[, c("scale", "n", "pearson_r", "r_squared", "slope")])
print(rao_afaith)
