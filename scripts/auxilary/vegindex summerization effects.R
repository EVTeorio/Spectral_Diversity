setwd("C:/lecospec")
source("C:/lecospec/Functions/lecospectR.R")

library(spectrolab)
library(beepr)
library(dplyr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/2025_Biomass_LiDAR/SOFOR")

summary <- read.csv("big_data/PaintRock_Spectra_Summary.CSV")

full <- read.csv("big_data/PaintRock_Spectra.CSV")
beep()


# Calculate Vegindices
trees_image_spectra <- full
trees_image_spectra_df<- speclib_to_df(trees_image_spectra)

#Calculate vegetation indices for the pixels
trees_image_spectra_VIs<-get_vegetation_indices(trees_image_spectra_df, NULL) 
beep()

VIs_full <-cbind(as.data.frame(trees_image_spectra)[,3:5],trees_image_spectra_VIs)
###################################################################################

VIs_summary <- VIs_full %>%
  group_by(ImageID, TreeID, SpeciesID) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(VIs_summary, "big_data/summaryof_VIs.CSV")
write.csv(VIs_summarized, "big_data/presummarized_VIs.CSV")
############# comparing data sets #########
df_a <- VIs_summary
df_b <- VIs_summarized
df_a <- df_a %>% dplyr::select(-SpeciesID)

# Numeric comparison tolerance
tolerance <- 1e-8

# Composite key columns
key_cols <- c("ImageID", "TreeID")

# ---- SORT BY KEY ----
df_a <- df_a %>% arrange(across(all_of(key_cols)))
df_b <- df_b %>% arrange(across(all_of(key_cols)))


# ---- VALIDATE KEYS ----
if (!identical(df_a[key_cols], df_b[key_cols])) {
  stop("ERROR: ImageID + TreeID keys do not match between data frames.")
}


# ---- COMPARE DATA ----
comparison_result <- all.equal(df_a, df_b, tolerance = tolerance)

if (isTRUE(comparison_result)) {
  cat("SUCCESS: Data frames are equal within tolerance.\n")
} else {
  cat("WARNING: Data frames differ.\n\n")
  print(comparison_result)
}


# ---- LOCATE DIFFERENCES ----
# ---- Column-wise difference summary ----
col_diff_summary <- data.frame(
  Column    = character(),
  Difference = logical(),
  NA_Diff   = logical(),
  Inf_Diff  = logical(),
  stringsAsFactors = FALSE
)

for (col in names(df_a)) {
  
  # Skip key columns
  if (col %in% key_cols) next
  
  a_vals <- df_a[[col]]
  b_vals <- df_b[[col]]
  
  # Any difference, ignoring NAs for numeric comparison
  diff_exists <- any(a_vals != b_vals, na.rm = TRUE)
  
  # Difference caused by NA mismatches
  na_diff <- any(is.na(a_vals) != is.na(b_vals))
  
  # Difference caused by Inf mismatches
  inf_diff <- any(is.infinite(a_vals) != is.infinite(b_vals))
  
  # Add row to summary table
  col_diff_summary <- rbind(
    col_diff_summary,
    data.frame(
      Column = col,
      Difference = diff_exists,
      NA_Diff = na_diff,
      Inf_Diff = inf_diff,
      stringsAsFactors = FALSE
    )
  )
}
col_diff_summary


#################################

# ---- Identify columns with any NA or Inf ----
cols_with_na_or_inf <- sapply(names(df_a), function(col) {
  if (col %in% key_cols) return(TRUE)  # exclude key columns
  any(is.na(df_a[[col]]) | is.na(df_b[[col]]) |
        is.infinite(df_a[[col]]) | is.infinite(df_b[[col]]))
})

# Columns to compare (exclude NAs/Inf)
cols_to_compare <- names(df_a)[!cols_with_na_or_inf]

cat("Columns excluded due to NA or Inf:\n")
print(names(df_a)[cols_with_na_or_inf])

# ---- Column-wise numeric difference summary (clean columns only) ----
col_numeric_diff <- data.frame(
  Column = character(),
  Num_Diff = integer(),
  Prop_Diff = numeric(),
  Mean_Abs_Diff = numeric(),
  Max_Abs_Diff = numeric(),
  stringsAsFactors = FALSE
)

for (col in cols_to_compare) {
  
  a_vals <- df_a[[col]]
  b_vals <- df_b[[col]]
  
  # safe numeric comparison since we excluded NA/Inf columns
  diffs <- a_vals - b_vals
  
  col_numeric_diff <- rbind(
    col_numeric_diff,
    data.frame(
      Column = col,
      Num_Diff = sum(diffs != 0),
      Prop_Diff = mean(diffs != 0),
      Mean_Abs_Diff = mean(abs(diffs)),
      Max_Abs_Diff = max(abs(diffs)),
      stringsAsFactors = FALSE
    )
  )
}

# Sort by proportion of difference
col_numeric_diff <- col_numeric_diff[order(-col_numeric_diff$Prop_Diff), ]

# View summary
col_numeric_diff

################################################################3

# ---- Compute relative difference ----
col_numeric_combined_stats$Relative_Diff <- NA

for (i in 1:nrow(col_numeric_combined_stats)) {
  
  min_val <- col_numeric_combined_stats$Combined_Min[i]
  max_val <- col_numeric_combined_stats$Combined_Max[i]
  max_diff <- col_numeric_combined_stats$Max_Abs_Diff[i]
  
  # Only compute if range is non-zero
  if (!is.na(min_val) & !is.na(max_val) & (max_val - min_val) != 0) {
    col_numeric_combined_stats$Relative_Diff[i] <- max_diff / (max_val - min_val)
  } else {
    col_numeric_combined_stats$Relative_Diff[i] <- NA
  }
}

# ---- Reorder columns: put Relative_Diff first ----
col_numeric_combined_stats <- col_numeric_combined_stats[, c(
  "Relative_Diff", setdiff(names(col_numeric_combined_stats), "Relative_Diff")
)]

write.csv(col_numeric_combined_stats,"big_data/VI_Summary_Diffrences.CSV")
# View the updated table
col_numeric_combined_stats
