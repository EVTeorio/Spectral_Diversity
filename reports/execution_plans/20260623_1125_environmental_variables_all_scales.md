# ExecPlan: Environmental Variables at All Quadrat Scales

## Objective

Create reproducible environmental covariate outputs for 10 m, 20 m, and 50 m quadrats using the local DTM raster.

## Requested Task

Calculate mean elevation and mean Riley topographic roughness index for each quadrat. Riley TRI should be calculated from the DTM using 5x5, 11x11, and 21x21 moving windows. Outputs should be CSVs and shapefiles written to `Quad_Values/Enviro_SHPs/`.

## Files To Review

- `RESEARCH_OBJECTIVES.md`
- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `reports/directory_map.md`
- `reports/data_dictionary.md`
- `scripts/2_Indices Creation/Enviro_Variables/elevation.R`
- `scripts/2_Indices Creation/Plant_diversity/plant_diversity_all_scales.R`

## Proposed Changes

- Add a new all-scale R workflow under `scripts/2_Indices Creation/Enviro_Variables/`.
- Preserve the existing single-scale `elevation.R` as historical reference.
- Add focused tests for the reusable Riley TRI helper logic.
- Update project documentation after outputs are generated.

## Expected Outputs

- `Quad_Values/Enviro_SHPs/enviro_variables_10m.csv`
- `Quad_Values/Enviro_SHPs/enviro_variables_20m.csv`
- `Quad_Values/Enviro_SHPs/enviro_variables_50m.csv`
- Matching `.shp`, `.dbf`, `.prj`, `.shx`, and related shapefile sidecars for each scale.

## Validation Plan

- Confirm all required R packages are available.
- Run the workflow with `Rscript`.
- Confirm output row counts match quadrat shapefiles: 2,000 for 10 m, 500 for 20 m, and 80 for 50 m.
- Confirm required columns are present and non-missing for valid plot quadrats.
- Run the focused test file for Riley TRI helper behavior.

## Risks

- Shapefile field names are limited to 10 characters, so output metric columns will use compact names: `elev_mean`, `tri5_mean`, `tri11_mean`, and `tri21_mean`.
- The DTM CRS is a compound NAD83(CSRS) UTM 16N CRS while quadrat shapefiles are WGS84; the workflow will transform quadrats to the DTM CRS before extraction.
- Mean extraction uses raster cells intersecting each quadrat polygon; at 1 m DTM resolution this should be appropriate for 10 m, 20 m, and 50 m quadrats.
