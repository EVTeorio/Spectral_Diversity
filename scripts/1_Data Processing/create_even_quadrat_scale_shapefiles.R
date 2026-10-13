library(sf)

OUT_DIR <- "Quad_Scale_SHPs"
BASE_SHP <- file.path(OUT_DIR, "PR_50m.shp")
TARGET_CRS <- 4326
WORK_CRS <- 26916

delete_shapefile <- function(path_no_ext) {
  sidecars <- paste0(path_no_ext, c(".shp", ".shx", ".dbf", ".prj", ".cpg"))
  existing <- sidecars[file.exists(sidecars)]
  if (length(existing) > 0) {
    unlink(existing)
  }
}

unit_vector <- function(x) {
  x / sqrt(sum(x^2))
}

make_oriented_grid <- function(footprint, cell_size, id_label) {
  footprint_utm <- st_transform(st_union(st_geometry(footprint)), WORK_CRS)
  obb <- st_minimum_rotated_rectangle(footprint_utm)
  coords <- st_coordinates(obb)[, 1:2, drop = FALSE]
  coords <- coords[!duplicated(coords), , drop = FALSE]

  if (nrow(coords) != 4) {
    stop("Expected four unique coordinates from the oriented bounding box.", call. = FALSE)
  }

  edges <- rbind(coords[c(2, 3, 4, 1), , drop = FALSE] - coords)
  edge_lengths <- sqrt(rowSums(edges^2))
  start_idx <- which.max(edge_lengths)

  origin <- coords[start_idx, ]
  x_vec <- unit_vector(edges[start_idx, ])
  y_vec <- unit_vector(edges[((start_idx) %% 4) + 1, ])
  width <- edge_lengths[start_idx]
  height <- edge_lengths[((start_idx) %% 4) + 1]

  n_col <- as.integer(round(width / cell_size))
  n_row <- as.integer(round(height / cell_size))

  if (!isTRUE(all.equal(width, n_col * cell_size, tolerance = 1))) {
    stop("Footprint width does not divide cleanly by requested cell size.", call. = FALSE)
  }
  if (!isTRUE(all.equal(height, n_row * cell_size, tolerance = 1))) {
    stop("Footprint height does not divide cleanly by requested cell size.", call. = FALSE)
  }

  make_cell <- function(col, row) {
    p0 <- origin + x_vec * ((col - 1) * cell_size) + y_vec * ((row - 1) * cell_size)
    p1 <- p0 + x_vec * cell_size
    p2 <- p1 + y_vec * cell_size
    p3 <- p0 + y_vec * cell_size
    st_polygon(list(rbind(p0, p1, p2, p3, p0)))
  }

  geometries <- vector("list", n_col * n_row)
  names <- character(n_col * n_row)
  index <- 1
  for (row in seq_len(n_row)) {
    for (col in seq_len(n_col)) {
      geometries[[index]] <- make_cell(col, row)
      names[index] <- paste0("sub", id_label, "_", index)
      index <- index + 1
    }
  }

  grid <- st_sf(
    Name = names,
    Dscrptn = NA_character_,
    geometry = st_sfc(geometries, crs = WORK_CRS)
  )

  st_transform(grid, TARGET_CRS)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
base <- st_read(BASE_SHP, quiet = TRUE)

targets <- data.frame(
  scale_label = c("25m", "100m"),
  id_label = c("25", "100"),
  cell_size = c(25, 100),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(targets))) {
  scale_label <- targets$scale_label[i]
  id_label <- targets$id_label[i]
  cell_size <- targets$cell_size[i]
  out_base <- file.path(OUT_DIR, paste0("PR_", scale_label))
  out_shp <- paste0(out_base, ".shp")

  grid <- make_oriented_grid(base, cell_size, id_label)
  delete_shapefile(out_base)
  st_write(grid, out_shp, driver = "ESRI Shapefile", quiet = TRUE)

  message("Wrote ", out_shp, " with ", nrow(grid), " polygons.")
}
