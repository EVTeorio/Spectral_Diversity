
saveRDS(global_pca, file = "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity/global_pca.rds")


library(terra)
library(snow)

library(terra)
library(snow)

# -----------------------------
# Parallel global PCA build
# -----------------------------

n_cores  <- max(1, parallel::detectCores() - 2)
n_sample <- 5000
seed     <- 123

cl <- makeCluster(n_cores, type = "SOCK")

# Load terra and set deterministic seeds per worker
clusterEvalQ(cl, library(terra))
clusterSetRNGStream(cl, seed)

# Export variables needed by workers
clusterExport(cl, list = "n_sample")

X_list <- parLapply(cl, ras_files, function(f) {
  
  r <- rast(f)
  
  X <- values(r, mat = TRUE)
  X <- X[rowSums(is.na(X)) == 0, ]
  X <- X[rowSums(X) > 0, ]
  
  if (nrow(X) > n_sample) {
    X <- X[sample(nrow(X), n_sample), , drop = FALSE]
  }
  
  X
})

stopCluster(cl)

X_global <- do.call(rbind, X_list)

X_global_scaled <- scale(X_global)

global_center <- attr(X_global_scaled, "scaled:center")
global_scale  <- attr(X_global_scaled, "scaled:scale")

global_pca <- prcomp(X_global_scaled, center = FALSE, scale. = FALSE)
