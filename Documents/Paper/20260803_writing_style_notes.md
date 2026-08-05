# Writing Style Notes From Prior Drafts

Created: 2026-08-03

Workflow note: the prior Word documents and repository Markdown context were read with `tools/create_paper_outline_docx.R` using R 4.2.3. Continue using R for Word-output generation in this paper workflow unless a later instruction changes that.

Source documents used for style context:

- `Documents/Spectral Variation Paper.docx`: not readable during this pass because the file was locked by another process.
- `Documents/SVH Outline.docx`: 66 paragraphs read.
- Repository Markdown context summarized in `Documents/Paper/20260803_markdown_context_digest.md`

## Observed Writing Style

- The prior writing is exploratory but hypothesis-driven. It tends to frame the work around a clear research question, then build toward methods and interpretation through staged explanation.
- The prose favors direct scientific framing over ornament: define the ecological problem, introduce spectral variation as the possible proxy, then specify the scale and diversity metrics.
- Section-level organization matters. The existing outline style uses manuscript sections and nested conceptual prompts rather than a continuous essay.
- Methods should be concrete and reproducible. The strongest prior material names data sources, quadrat scales, filtering decisions, spectral processing steps, and expected outputs.
- The voice is cautious about interpretation. Claims should be framed as relationships, correlations, or evidence for/against the spectral variation hypothesis, not as broad proof.
- The writing benefits from explicit bridges between remote sensing and ecology: spectral heterogeneity should be tied back to canopy, species, phylogenetic, and scale concepts.

## Style To Preserve In The New Outline

- Use working headings that can later become manuscript headings.
- Keep bullets detailed enough to draft from, but avoid writing final manuscript paragraphs too early.
- State methodological updates plainly and in chronological workflow order.
- Keep technical terms visible: spectral angle entropy, standardized PCA mean distance, spectral Rao's Q, abundance-weighted Faith's PD, Shannon species diversity, 10 m, 20 m, and 50 m quadrats.
- Keep environmental factors out of the main outline for this stage because the current focus is the correlation between species/phylogenetic diversity and spectral diversity across scales.

## Key Methodological Updates To Carry Forward

- Species and phylogenetic diversity are now generated at 10 m, 20 m, and 50 m scales using the current all-scale plant-diversity workflow and quadrat IDs aligned with spectral outputs.
- Species diversity emphasis for this outline should be Shannon diversity only.
- Phylogenetic diversity emphasis for this outline should be abundance-weighted Faith's PD only.
- Spectral heterogeneity emphasis should include spectral angle entropy, standardized PCA mean distance, and standardized PCA Rao's Q as the spectral Rao-style measure.
- Spectral angle entropy is calculated from sunlit, shadow-masked, smoothed 5 nm spectra, with most quadrats using bootstrap means rather than exhaustive all-pair calculations.
- Standardized PCA metrics use vector-normalized spectra before projection, reducing broad brightness dominance compared with regular PCA.
- The manuscript question should be narrowed to whether spectral-diversity relationships differ between abundance-weighted Faith's PD and Shannon species diversity, and how these relationships change across 10 m, 20 m, and 50 m scales.

## Result Emphasis For This Outline

- Abundance-weighted Faith's PD shows positive relationships with the primary spectral heterogeneity metrics across scales, with the standardized PCA mean distance relationship strongest at 20 m and 50 m.
- Shannon diversity relationships are weak, near zero, or negative in the primary non-edge analysis, making it useful as the species-diversity contrast.
- Scale should be treated as a central result, not a footnote: the phylogenetic relationship becomes more interpretable at broader grains, while species Shannon does not show a comparable pattern.
- Results should be framed as correlation and scale-dependence, not environmental explanation.
