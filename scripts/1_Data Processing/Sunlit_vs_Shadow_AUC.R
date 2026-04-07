
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(pROC)
library(ggplot2)
library(beepr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/")

illumination <- st_read("Shadow_vs_sunlit_SHP/Shadow_vs_sunlit_Rev.shp")

tile_dir <- "Quad_Spectra/20m_smoothed_5nm"
allfiles <- list.files(tile_dir, full.names = TRUE)
tile_files <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]

# illumination object must already be loaded (sf polygon with Type field)

# ----------------------------- Function ----------------------------------
extract_tile_pixels <- function(tile_path) {
  
  tile_name <- tools::file_path_sans_ext(basename(tile_path))
  cat("Processing:", tile_name, "\n")
  
  tile <- rast(tile_path)
  
  # ------- FIX: Standardize band names -----For use only with VegIndex below
  original_names <- names(tile)
  
  # Extract only the VI_# portion from each band name
  clean_names <- sub(".*(VI_\\d+)$", "\\1", original_names)
  
  names(tile) <- clean_names
  #---------------------------- only with VegIndex above _-------#
  
  # Reproject polygons
  illum_proj <- st_transform(illumination, crs(tile))
  
  # Convert tile extent to sf polygon
  tile_extent_sf <- st_as_sf(as.polygons(ext(tile)))
  st_crs(tile_extent_sf) <- st_crs(illum_proj)
  
  # Intersect
  illum_tile <- st_intersection(illum_proj, tile_extent_sf)
  
  if (nrow(illum_tile) == 0) {
    return(NULL)
  }
  
  illum_vect <- vect(illum_tile)
  
  vals <- terra::extract(tile, illum_vect)
  
  vals$class <- tolower(illum_tile$Type)[vals$ID]
  vals$tile  <- tile_name
  
  vals %>%
    select(-ID)
}


# ----------------------------- Extraction --------------------------------
training_data <- map_dfr(tile_files, extract_tile_pixels)
beep()

# Remove duplicate pixels if tiles overlap
training_data <- distinct(training_data)

# ----------------------------- Long Format --------------------------------
long_vals <- training_data %>%
  pivot_longer(cols = -class: -tile,
               names_to = "band",
               values_to = "value") %>%
  drop_na()

# ----------------------------- ROC AUC -----------------------------------
roc_scores <- long_vals %>%
  group_by(band) %>%
  summarise(
    auc = tryCatch({
      roc_obj <- roc(response = class,
                     predictor = value,
                     levels = c("shadow", "sunlit"),
                     quiet = TRUE)
      as.numeric(auc(roc_obj))
    }, error = function(e) NA_real_)
  ) %>%
  arrange(desc(auc))

cat("\nTop Performing Bands:\n")
print(head(roc_scores, 10))

# ------------------------ Optimal Threshold -------------------------------
best_band <- roc_scores$band[1]

best_vals <- long_vals %>%
  filter(band == best_band)

roc_obj <- roc(response = best_vals$class,
               predictor = best_vals$value,
               levels = c("shadow", "sunlit"))

opt_threshold <- coords(roc_obj,
                        x = "best",
                        ret = c("threshold", "sensitivity", "specificity"))
best_band
roc_obj$direction
print(opt_threshold)
##################################################################################

plot_data <- roc_scores %>%
  mutate(
    wavelength = as.numeric(gsub(" nm", "", band))
  ) %>%
  filter(!is.na(wavelength))

#for VIs
plot_data <- roc_scores %>%
  mutate(
    wavelength = band)


# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

ggplot(plot_data, aes(x = wavelength, y = auc)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "black", size = 2) +
  geom_vline(xintercept = as.numeric(gsub(" nm", "", best_band)),
             linetype = "dashed",
             color = "red") +
  labs(
    title = "Wavelength vs ROC AUC",
    x = "Wavelength (nm)",
    y = "AUC"
  ) +
  theme_minimal(base_size = 8)

################################################################################

# Count pixels per tile per class
tile_summary <- training_data %>%
  count(tile, class) %>%
  tidyr::pivot_wider(
    names_from = class,
    values_from = n,
    values_fill = 0
  )

# Add totals per tile
tile_summary <- tile_summary %>%
  mutate(total = rowSums(across(where(is.numeric))))

# Add grand total row
grand_total <- tile_summary %>%
  summarise(
    tile = "TOTAL",
    across(where(is.numeric), sum)
  )

tile_summary <- bind_rows(tile_summary, grand_total)

print(tile_summary)

###########                Visual Test                ###################################
#Go to "Spectral_Diversity/scripts/visuals/shadow_mask.R"




