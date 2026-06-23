# Builds the scale-aware species matrix and diversity outputs.
source(file.path(
  "scripts",
  "2_Indices Creation",
  "Plant_diversity",
  "plant_diversity_all_scales.R"
))

run_plant_diversity_workflow()
