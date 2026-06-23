# ExecPlan: Plant Diversity Values at All Scales

Date: 2026-06-22

## Objective

Update the plant-diversity workflow so species composition, species diversity, and phylogenetic diversity values are generated consistently for 10 m, 20 m, and 50 m quadrats.

## Requested Task

- Align quadrat identifiers with the spectral diversity outputs.
- Preserve the existing crown-buffer and crown-overlap logic used to build the species matrix.
- Create one CSV and one shapefile per scale under `Indices_SHPs/Diversity_SHPs`.

## Files To Review

- `scripts/2_Indices Creation/Plant_diversity/sp_weighted_matrix.R`
- `scripts/2_Indices Creation/Plant_diversity/species_diversity.R`
- `scripts/2_Indices Creation/Plant_diversity/phylogenetic_diversity.R`
- `Quad_Scale_SHPs/PR_10m.shp`
- `Quad_Scale_SHPs/PR_20m.shp`
- `Quad_Scale_SHPs/PR_50m.shp`
- `PR_tree_DL.csv`
- `51sp_taxanomy.csv`

## Proposed Changes

- Add a reusable all-scale workflow script with shared functions for:
  - reading scale-specific quadrats;
  - assigning canonical `quad_id` values;
  - intersecting buffered tree crowns with quadrats;
  - building species-by-quadrat matrices;
  - calculating richness, Shannon, Simpson, evenness, Faith PD, Rao PD, and abundance-weighted Faith PD;
  - writing per-scale CSVs and shapefiles.
- Update the existing plant-diversity scripts to use the shared workflow instead of separate hard-coded 20 m implementations.
- Keep output names stable and scale explicit.

## Expected Outputs

- `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_20m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_50m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.shp`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_20m.shp`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_50m.shp`

## Validation Plan

- Confirm scripts parse under R 4.2.3.
- Run the all-scale workflow if required R packages are visible.
- Confirm outputs have expected row counts: 2,000 for 10 m, 500 for 20 m, and 80 for 50 m.
- Confirm `quad_id` values align to spectral output IDs.
- Confirm no missing expected diversity fields in the CSVs.

## Risks

- `sf`, `dplyr`, `readr`, `tidyr`, `ape`, `picante`, and `V.PhyloMaker2` may require access to the user's R library outside the workspace sandbox.
- Full crown-polygon intersections may take time, especially for the 10 m scale.
- Shapefile field-name limits can truncate long species or metric names, so the CSVs should be treated as the full-fidelity tabular outputs.
