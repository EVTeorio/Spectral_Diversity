setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens")


# ---- Define directories ----
dir1 <- "Spectral_Diversity/Quad_Spectra/20m_VegIndex"
dir2 <- "Spectral_Diversity/Quad_Spectra/20m_RGB_old"

# ---- Function to list filtered files ----
get_filtered_files <- function(directory) {
  allfiles <- list.files(directory, full.names = TRUE)
  
  # Remove unwanted extensions
  filtered <- allfiles[!grepl("\\.hdr$|\\.aux$|\\.xml$|\\.enp$|\\.sta$", allfiles)]
  
  # Return only file names (not full paths)
  basename(filtered)
}

# ---- Get file names from both directories ----
files_dir1 <- get_filtered_files(dir1)
files_dir2 <- get_filtered_files(dir2)
files_dir2 <- sub("\\.tif$", "", files_dir2, ignore.case = TRUE)
# ---- Compare file names ----
missing_in_dir2 <- setdiff(files_dir1, files_dir2)
missing_in_dir1 <- setdiff(files_dir2, files_dir1)

# ---- Results ----
cat("Files in dir1 but missing in dir2:\n")
print(missing_in_dir2)

cat("\nFiles in dir2 but missing in dir1:\n")
print(missing_in_dir1)

cat("\nSummary:\n")
cat("Total files in dir1:", length(files_dir1), "\n")
cat("Total files in dir2:", length(files_dir2), "\n")
cat("Missing in dir2:", length(missing_in_dir2), "\n")
cat("Missing in dir1:", length(missing_in_dir1), "\n")
