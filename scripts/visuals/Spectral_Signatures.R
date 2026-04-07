
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Parameters (same as your masking script) ---
n_sample  <- 2000
best_band <- "603 nm"
threshold <- 0.01898618
direction <- "<"

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/")

rast1 <- rast("Spectral_Diversity/Quad_Spectra/20m_resampled_5nm/503")
rast2 <- rast("Spectral_Diversity/Quad_Spectra/20m_resampled_5nm/504")
rast3 <- rast("Spectral_Diversity/Quad_Spectra/20m_resampled_5nm/603")
rast4 <- rast("Spectral_Diversity/Quad_Spectra/20m_resampled_5nm/604")

raster_list <- list(
  "Quad 0503" = rast1,
  "Quad 0504" = rast2,
  "Quad 0603" = rast3,
  "Quad 0604" = rast4
)

# --- Extract masked spectra from each raster ---
extract_spectra <- function(r, label) {
  
  # Apply sunlit mask
  band_mask   <- r[[best_band]]
  sunlit_mask <- if (direction == ">") {
    band_mask < threshold
  } else {
    band_mask > threshold
  }
  r_masked <- mask(r, sunlit_mask, maskvalues = 0)
  
  # Extract pixel matrix
  X <- values(r_masked, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  X <- X[rowSums(X) > 0, , drop = FALSE]
  
  if (nrow(X) == 0) return(NULL)
  
  # Subsample if needed
  if (nrow(X) > n_sample) X <- X[sample(nrow(X), n_sample), , drop = FALSE]
  
  # Parse wavelengths from band names
  band_names  <- names(r)
  wavelengths <- as.numeric(gsub("[^0-9.]", "", band_names))
  
  # Build long format dataframe
  df <- as.data.frame(X)
  names(df) <- wavelengths
  df$pixel_id <- seq_len(nrow(df))
  
  df_long <- df %>%
    pivot_longer(
      cols      = -pixel_id,
      names_to  = "Wavelength",
      values_to = "Reflectance"
    ) %>%
    mutate(
      Wavelength = as.numeric(Wavelength),
      Raster     = label
    )
  
  df_long
}

# --- Run extraction ---
all_spectra <- bind_rows(
  mapply(extract_spectra, raster_list, names(raster_list), SIMPLIFY = FALSE)
)

# --- Summarise per raster per wavelength ---
summary_df <- all_spectra %>%
  filter(!is.na(Wavelength), Wavelength <= 850) %>%
  group_by(Raster, Wavelength) %>%
  summarise(
    mean_reflectance = mean(Reflectance, na.rm = TRUE),
    p25    = quantile(Reflectance, 0.25,  na.rm = TRUE),
    p75    = quantile(Reflectance, 0.75,  na.rm = TRUE),
    p12.5  = quantile(Reflectance, 0.125, na.rm = TRUE),
    p87.5  = quantile(Reflectance, 0.875, na.rm = TRUE),
    .groups = "drop"
  )

# --- Visible spectrum bands ---
visible_bands <- tibble(
  xmin = c(400, 495, 570),
  xmax = c(495, 570, 700),
  fill = c("blue", "green", "red")
)

# --- Plot: one panel per raster ---
raster_ids <- unique(summary_df$Raster)

y_min <- 0
y_max <- .3

for (ras in raster_ids) {
  
  plot_df  <- summary_df %>% filter(Raster == ras)
  n_pixels <- all_spectra %>%
    filter(Raster == ras) %>%
    distinct(pixel_id) %>%
    nrow()
  
  p <- ggplot() +
    # Visible spectrum background
    geom_rect(
      data         = visible_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
      alpha        = 0.35,
      inherit.aes  = FALSE
    ) +
    scale_fill_manual(
      values = c("blue" = "#ADD8E6", "green" = "#90EE90", "red" = "#FFC0CB"),
      guide  = "none"
    ) +
    
    # 75th percentile shaded band
    geom_ribbon(
      data  = plot_df,
      aes(x = Wavelength, ymin = p12.5, ymax = p87.5),
      fill  = "lightgray",
      alpha = 0.5
    ) +
    # 50th percentile shaded band
    geom_ribbon(
      data  = plot_df,
      aes(x = Wavelength, ymin = p25, ymax = p75),
      fill  = "darkgray",
      alpha = 0.6
    ) +
    # Mean line
    geom_line(
      data  = plot_df,
      aes(x = Wavelength, y = mean_reflectance),
      color = "black",
      size  = 1.1
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      title    = paste0("Spectral Signature: ", ras,
                        " (n = ", n_pixels, "pixels)"),
      subtitle = "Dark gray = IQR (25–75%), Light gray = 12.5–87.5%",
      x        = "Wavelength (nm)",
      y        = "Reflectance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  print(p)
  
  # #Optionally save
  # ggsave(filename = paste0("Spectral_Diversity/Graphs", ras, "_spectrum.png"),
  #        plot = p, width = 10, height = 6)
}

