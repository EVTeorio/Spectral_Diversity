# Validation Report: Plant Diversity Values at All Scales

Date: 2026-06-22

## Validation Scope

Validate the updated plant-diversity scripts and the 10 m, 20 m, and 50 m CSV/shapefile outputs.

## Checks Performed

| Check | Result | Notes |
|---|---|---|
| Script parse | Pass | All four plant-diversity scripts parse under R 4.2.3. |
| Package-backed execution | Pass | `plant_diversity_all_scales.R` completed using the R 4.2 user library. |
| 10 m output rows | Pass | CSV and DBF each contain 2,000 rows. |
| 20 m output rows | Pass | CSV and DBF each contain 500 rows. |
| 50 m output rows | Pass | CSV and DBF each contain 80 rows. |
| Required metric fields | Pass | `richness`, `shannon`, `simpson`, `evenness`, `faith_pd`, `rao_pd`, and `afaith_pd` have no missing values in all three CSVs. |
| Spectral ID overlap | Pass | Every current spectral `quad_id` exists in the matching plant-diversity CSV. |
| Legacy 20 m files preserved | Pass | Existing `species_diversity_20m.*` and `Phylogenetic_Diversity_20m.*` files were restored after the first write attempt displaced them. |

## Row Counts

| Scale | CSV Rows | Shapefile DBF Rows | Expected Rows |
|---|---:|---:|---:|
| 10 m | 2,000 | 2,000 | 2,000 |
| 20 m | 500 | 500 | 500 |
| 50 m | 80 | 80 | 80 |

## ID Alignment With Current Spectral Summaries

| Scale | Spectral IDs Missing From Plant Output | Plant IDs Missing From Spectral Output | Interpretation |
|---|---:|---:|---|
| 10 m | 0 | 91 | Matches known spectral-raster missing quadrats. |
| 20 m | 0 | 15 | Matches known spectral-raster missing quadrats. |
| 50 m | 0 | 0 | Complete overlap. |

## Notes

- The new plant outputs deliberately include every quadrat polygon at each scale. Downstream joins to spectral entropy should therefore produce missing spectral fields for the known 10 m and 20 m quadrats without spectral raster summaries.
- The workflow uses the same quadrat ID convention as current spectral outputs: 10 m `sub_id`, 20 m `Name`, and 50 m `Name`.
