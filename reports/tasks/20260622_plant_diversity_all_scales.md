# Task Report: Plant Diversity Values at All Scales

Date: 2026-06-22

## Objective

Generate species composition, species diversity, and phylogenetic diversity outputs for 10 m, 20 m, and 50 m quadrats with quadrat identifiers aligned to the current spectral diversity outputs.

## Scripts Updated

- `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`
- `scripts/2_Indices Creation/Plant_diversity/sp_weighted_matrix.R`
- `scripts/2_Indices Creation/Plant_diversity/species_diversity.R`
- `scripts/2_Indices Creation/Plant_diversity/phylogenetic_diversity.R`

## Key Changes

- Added a shared all-scale plant-diversity workflow.
- Preserved the existing crown-buffer and quadrat-intersection logic:
  - tree census input: `PR_tree_DL.csv`;
  - retained trees: `DBH.2024 >= 200`, `crown.position` in `3, 4, 5`, and `cluster_status` in `A, R`;
  - crown radius: `cw_m_2025 / 2`;
  - species matrix values: summed proportion of each tree crown intersecting each quadrat.
- Added scale-aware quadrat ID handling:
  - 10 m uses `sub_id` values such as `0_a`;
  - 20 m uses `Name` values such as `0`;
  - 50 m uses `Name` values such as `sub50_1`.
- Added per-quadrat species metrics:
  - `richness`
  - `shannon`
  - `simpson`
  - `evenness`
- Added per-quadrat phylogenetic metrics:
  - `faith_pd`
  - `rao_pd`
  - `afaith_pd`
- Replaced the previous 20 m-only scripts with wrappers that call the shared all-scale workflow.

## Outputs Created

- `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_20m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_50m.csv`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_10m.shp`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_20m.shp`
- `Indices_SHPs/Diversity_SHPs/plant_diversity_50m.shp`

## Notes

- The CSVs should be treated as the full-fidelity tabular outputs.
- The shapefiles contain the same analytical fields where possible, but future additions should remember the DBF field-name length limits.
- Existing older 20 m diversity shapefiles were preserved after an initial shapefile overwrite attempt displaced them into a temporary backup folder.
