library(dplyr)

het_df <- read.csv("Indices_SHPs/20m_SA_smooth_masked_7_11.csv")
csv <- read.csv("Indices_SHPs/20m_spectral_sp.csv")

# Ensure consistent column name (case-sensitive!)
het_df_clean <- het_df %>%
  rename(SA_entropy_smooth_masked = spectral_entropy)
  

# Join and relocate column
csv_updated <- csv %>%
  left_join(
    het_df %>%
      mutate(Name = as.integer(Name)) %>%
      rename(SA_entropy_smooth_masked_711 = spectral_entropy),
    by = "Name"
  ) %>%
  relocate(SA_entropy_smooth_masked_711, .after = Name) %>%
  select(-X)

write.csv(csv_updated, "Indices_SHPs/20m_spectral_sp.csv")
