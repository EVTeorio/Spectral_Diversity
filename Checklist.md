# Checklist

Created: 2026-08-19

Purpose: identify what is missing from `Documents/Paper/SVH_Paper.docx` relative to the revised research direction: **Assessing the Spectral Variation Hypothesis Across Spatial Scales Using PCA-based Canopy Spectral Heterogeneity**.

[A] = I've added an answer and or context to be considered.
[X] = Task is complete. Please update as tasks are completed.
[O] = Tasks is disgarded and will not be completed or stored for a later date.


## Analysis Needed

- [X] Run a transformation analysis for spectral-biodiversity correlations.
  - Test log, `log1p`, square-root, Box-Cox, Yeo-Johnson, and simple power transforms where appropriate.
  - Apply to PCA-based spectral metrics and biodiversity metrics by scale.
  - Report whether transformations improve Pearson r, R2, residual structure, and interpretability.
- [O] Evaluate the optimal quadrat scale for spectral-phylogenetic correlation.
  - Compare 10 m, 20 m, and 50 m, 100m using correlation strength, complete-case sample size, uncertainty, and ecological interpretability.
- [X] Add spatial autocorrelation diagnostics.
  - Minimum: Moran's I or variograms for key response variables and model residuals.
  - Priority models: `spec_spca_alpha`, `spec_spca_mean`, and `spec_spca_rao` against `phy_rao`, `phy_afaith`, and `sp_shannon`.
- [X] Summarize all-metric spectral-biodiversity relationships under the revised PCA-centered framing.
  - Include vector-normalized PCA metrics as the main set.
  - Include raw PCA metrics as brightness-sensitive comparison.
  - Include spectral angle entropy only as a supporting non-PCA comparison.
  - Use existing spectral heterogeneity relationship outputs as the base.
- [X] Compare spectral heterogeneity metrics against each other by scale.
  - Explain when metrics are redundant, when they diverge, and what those differences imply.
  - Use existing spectral heterogeneity relationship outputs as the base.
- [X] Compare biodiversity metrics against each other by scale.
  - Show why Faith's PD tracks richness, why `phy_rao` and `phy_afaith` diverge, and why species diversity often fails to track spectral heterogeneity.
- [X] Add an analysis layer for what drives metric values.
  - Spectral drivers: elevation, retained-pixel brightness, blue/green/red/NIR brightness, pixel count, edge status, and composition type.
  - Biodiversity drivers: species count, crown-equivalent individuals, evenness, dominance, and species-composition clusters.
  - Based on current outputs, only species/ indivdual counts and use of evolutionary history drives these metrics from this list
  - Based on current outputs only BGR/NIR and overall reflectance matter.
- [X] Decide whether `spec_spca_alpha` should become the lead spectral result.
  - Current all-metric results show it has the strongest relationships at 50 m and strong 20 m performance.
  - We should. We will treat all pca methods as equals
- [X] Decide whether all alpha hull including 3d metrics should be included
  - We should for now and will be able to take out later if need be.
- [X] Decide whether `spec_spca_rao` remains a main metric or becomes supporting.
  - It is theoretically attractive but appears more sensitive to squared distances and possibly to extreme pixels.
- [X] Clarify whether final PCA metrics were calculated from all retained pixels or from a 4,000-pixel sample.
  - The manuscript currently leans heavily on the 4,000-pixel sample-size section, but current documentation says main PCA metric values used all retained pixels for included quadrats.
  - Final PCA mean-distance and Rao-style metric fields are marked `all_pixels` for included quadrats. Alpha-hull and 3D hull fields use `all_pixels` for many 10 m quadrats but use sampled-pixel/fallback methods for larger 20 m and 50 m quadrats. The fixed 4,000-pixel outputs remain sample-size sensitivity results, not the final main PCA metric values.

## Missing Manuscript Sections

- [ ] Add a full title page or front matter using the revised working title.
- [ ] Add an abstract that foregrounds scale, PCA-based spectral heterogeneity, and metric comparison.
- [ ] Add keywords.
- [ ] Add the Introduction.
  - Introduce the Spectral Variation Hypothesis.
  - Explain why scale matters.
  - Explain why PCA-based spectral heterogeneity is useful.
  - Explain why biodiversity metrics can disagree.
  - State the revised research questions.
- [ ] Add a clear objectives/hypotheses section.
  - Include metric correlation, metric mechanism, scale dependence, transformation testing, and optimal quadrat scale.
- [ ] Add a table or subsection defining all spectral and biodiversity metrics.
- [ ] Add a Results section on spectral metric relationships.
  - Raw PCA versus vector-normalized PCA.
  - Mean distance versus alpha hull versus Rao's Q.
  - Metric redundancy and divergence across scale.
- [ ] Add a Results section on biodiversity metric relationships.
  - Richness versus Faith's PD.
  - Shannon/Simpson/evenness versus phylogenetic metrics.
  - Why abundance weighting changes the phylogenetic metrics.
- [ ] Add a Results section on all spectral-biodiversity metric correlations by scale.
  - This should replace the current narrow focus on only standardized PCA mean distance and spectral Rao's Q versus abundance-weighted Faith's PD.
- [ ] Add a Results section on transformations if the analysis is completed.
- [ ] Add a Results section on spatial autocorrelation or spatial sensitivity.
- [ ] Add a Discussion section specifically about scale.
  - Why 10 m is weak.
  - Why 20 m may be a practical compromise.
  - Why 50 m currently has the strongest correlations but needs sample-size caution.
- [ ] Add a Discussion section on metric construction.
  - Linear distance versus squared distance.
  - Area/hull metrics versus centroid-distance metrics.
  - Raw reflectance PCA versus vector-normalized PCA.
- [ ] Add a Discussion section on why phylogenetic diversity can correlate when Shannon diversity does not.
- [ ] Add limitations.
  - Spatial non-independence.
  - Single forest plot.
  - 50 m sample size.
  - Illumination and topographic effects.
  - Crown-buffer approximation.
  - Phylogenetic tree construction assumptions.
- [ ] Add conclusions.
- [ ] Add data availability, code availability, author contributions, funding, acknowledgments, conflicts of interest, and references for MDPI-style completion.

## Methods That Need Expansion Or Correction

- [ ] Expand crown-radius calculation.
  - Current text says "Crown radius was calculated by..." and is incomplete.
  - This was also flagged by an existing Word comment.
- [ ] Document the crown-overlap abundance procedure more fully.
  - State whether overlap proportions are summed by species.
  - Clarify whether these values are proportionalized before diversity metrics.
  - Explain why canopy-weighted abundance matches the remote sensing signal better than stem counts.
- [ ] Document abundance-weighted Faith's PD more precisely.
  - Explain which branches are weighted and how species abundance enters the calculation.
  - This was also flagged in Evan's notes.
- [ ] Clarify tree inclusion filters.
  - Current manuscript says DBH >= 200 mm and crown position class 3-5.
  - Current project documentation elsewhere says CPI = 4 or 5 OR DBH > 200 mm.
  - Verify the actual script logic and make the manuscript match it exactly.
- [ ] Clarify species count.
  - Current manuscript says 4,007 tree crowns representing 50 species out of 79.
  - Earlier objective text expected approximately 4,519 upper-canopy individuals and 51 species.
  - Use the current validated all-plot methods counts if they are final.
- [ ] Clarify coordinate reference system.
  - Current manuscript says WGS 1984 UTM Zone 16.
  - Current combined-table documentation says centroid coordinates are NAD83 / UTM zone 16N.
  - Verify imagery, quadrat, tree, and centroid CRS language separately.
- [ ] Clarify PCA basis construction.
  - State that the current PCA basis samples current 10 m smoothed 5 nm rasters after shadow masking.
  - State that raw PCA and vector-normalized PCA are separate bases.
- [ ] Clarify raw PCA brightness dominance.
  - Raw PCA PC1 is strongly correlated with mean reflectance.
  - Standardized PCA reduces but does not eliminate brightness relationships.
- [ ] Clarify PCA metric computation.
  - Current documentation indicates final PCA Euclidean metrics use all retained pixels for included quadrats.
  - Do not imply that the final main PCA values are fixed-4,000-pixel bootstrap estimates unless scripts confirm that.
- [ ] Add `spec_spca_alpha` methods if it remains a lead result.
  - The current manuscript does not adequately describe alpha-hull area.
- [ ] Add raw PCA metric methods as supporting comparison.
  - Mention raw mean distance, raw Rao's Q, and raw alpha hull without letting them dominate.
- [A] Decide whether sample-size results belong in main Methods/Results or supplement.
  - They belong in results
  - The revised narrative excludes sample-size runs from the concise metric list and centers scale and metric behavior instead.

## Results That Need Updating

- [ ] Replace the current single-biodiversity-metric focus with all-metric comparison.
  - Current Results emphasize abundance-weighted Faith's PD.
  - Revised direction needs `phy_rao`, `phy_afaith`, and `sp_shannon` as the main interpretive contrast, with other metrics summarized.
- [ ] Add the August 19 all-metric result.
  - `spec_spca_alpha` versus `phy_afaith` at 50 m: r = 0.497, R2 = 0.247.
  - `spec_spca_alpha` versus `phy_rao` at 50 m: r = 0.487, R2 = 0.238.
  - `spec_spca_mean` versus `phy_rao` at 50 m: r = 0.444, R2 = 0.197.
- [ ] Add 20 m strength for `spec_spca_alpha`.
  - `spec_spca_alpha` versus `phy_afaith`: r = 0.384, R2 = 0.147.
  - `spec_spca_alpha` versus `phy_rao`: r = 0.336, R2 = 0.113.
- [ ] Add species-diversity contrast results.
  - Shannon diversity and related species metrics are weak, near zero, or inconsistent compared with phylogenetic metrics.
- [ ] Add biodiversity metric concordance results.
  - Faith's PD tracks richness strongly.
  - `phy_rao` and `phy_afaith` diverge more from richness and Shannon.
- [ ] Add spectral metric correlation results.
  - Standardized PCA mean distance, alpha hull, and Rao's Q are strongly correlated, but their calculation differences matter.
- [ ] Add scale comparison across all selected metrics, not just two metric pairings.
- [X] Add transformation results once completed.
- [X] Add spatial autocorrelation results once completed.

## Possible Inaccuracies Or Mistakes In `SVH_Paper.docx`

- [ ] Species name typo: "Caryacordiformis" should likely be `Carya cordiformis`.
- [ ] Formatting issue: "50 (Table)species" needs spacing and table reference cleanup.
- [ ] Incomplete sentence: "Crown radius was calculated by..." must be finished.
- [ ] Filter logic may be inaccurate.
  - Manuscript: DBH >= 200 mm and crown position 3-5.
  - Project docs: CPI = 4 or 5 OR DBH > 200 mm.
  - Verify against the active plant-diversity script.
- [ ] Species/count language may be inconsistent.
  - Manuscript: 4,007 crowns, 50 species, 79 available.
  - Earlier project objective: approximately 4,519 individuals, 51 species.
  - Use the latest validated tree-count outputs and explain why species were removed.
- [ ] The manuscript currently says "upper and intermediate canopy" with crown position class 3-5.
  - Verify whether class 3 was actually included in final diversity outputs.
- [ ] The current Results imply the main PCA analysis uses 4,000 pixels.
  - Current project documentation suggests the final PCA metrics use all retained pixels, while 4,000-pixel values belong to the sample-size sensitivity branch.
- [ ] "PCAbased" needs a hyphen or spacing cleanup.
- [ ] "50m" should be standardized to "50 m" throughout.
- [ ] "spectral-angle entropy" and "spectral angle entropy" should be made consistent.
- [ ] "Spectral ROA's Q" in comments appears to be a typo for Rao's Q.
- [ ] "cloudshadow" should be "cloud-shadow" or "cloud shadow" consistently.
- [ ] "World Geodetic System 1984 Universal Transverse Mercator coordinate system, Zone 16" should be checked against the actual CRS used for imagery and quadrat centroids.
- [ ] "This preprocessing removed approximately 90% of the illumination signal" is too strong unless supported by a specific calculation.
  - Better: "substantially reduced brightness dominance" and cite standardized PCA mean-reflectance diagnostics.
- [ ] The manuscript frames spectral Rao's Q as one of two main results, but the revised direction may make it supporting behind `spec_spca_alpha` and `spec_spca_mean`.
- [ ] The manuscript does not yet address existing Word comment 1: "Check accuracy."
- [ ] The manuscript does not yet address existing Word comment 2: "Expound in crown calculations."
- [ ] The manuscript does not yet resolve existing comments 3 and 4 about narrowing spectral and biodiversity metrics.
- [ ] Existing Word comment 5 asks about wording around "reasonable" or "feasible"; revise computational feasibility language in the spectral Rao's Q section.

## Figure And Table Needs

- [ ] Figure: study site, quadrat scales, and workflow.
- [ ] Figure: spectral metric relationships by scale.
- [ ] Figure: biodiversity metric relationships by scale.
- [ ] Figure: all spectral-biodiversity correlations by scale, centered on vector-normalized PCA metrics.
- [ ] Figure: scale comparison of strongest correlations.
- [ ] Figure or supplement: raw PCA versus vector-normalized PCA brightness/illumination comparison.
- [ ] Figure or supplement: spatial autocorrelation diagnostics.
- [ ] Figure or supplement: transformation comparison.
- [ ] Table: concise metric definitions and calculation formulas.
- [ ] Table: quadrat coverage and complete-case counts.
- [ ] Table: top metric correlations by scale.
- [ ] Table: biodiversity metric concordance.
- [ ] Table: final main-versus-supplement figure inventory.

## Priority Order

- [ ] 1. Verify tree filtering, species counts, CRS language, and final PCA pixel-use details.
- [A] 2. Decide final main spectral metrics: likely `spec_spca_alpha`, `spec_spca_mean`, and supporting `spec_spca_rao`.
- I need to run correlation to the other alpha hulls metrics first
- [X] 3. Run transformation analysis.
- [X] 4. Run spatial autocorrelation diagnostics.
- [ ] 5. Rewrite Results around scale and metric behavior.
- [ ] 6. Rewrite Discussion around why metrics correlate, why they differ, and how scale changes the SVH signal.
- [ ] 7. Add missing Introduction, Abstract, Conclusions, and journal front/back matter.
