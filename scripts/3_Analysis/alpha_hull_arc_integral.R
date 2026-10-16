# Exact area of the same circular-arc boundary produced by alphahull::ahull.
# Each concave boundary arc is traversed clockwise about its circle center so
# that the hull interior remains on the left. Green's theorem gives
# 1/2 * integral(x dy - y dx), summed over all arcs. Closed components and holes
# are handled by their orientation; no component ordering or polygon sampling
# is necessary. This avoids areaahulleval's endpoint-reordering assumptions.
alpha_hull_arc_integral <- function(hull) {
  stopifnot(inherits(hull, "ahull"))
  a <- hull$arcs
  if (is.null(a) || nrow(a) == 0) return(0)
  a <- a[a[, 3] > 0, , drop = FALSE]
  if (nrow(a) == 0) return(0)
  stopifnot(all(is.finite(a[, 1:6])), all(a[, 6] >= 0), all(a[, 6] <= pi + 1e-10),
            all(abs(rowSums(a[, 4:5, drop = FALSE]^2) - 1) < 1e-8))
  center <- a[, 1:2, drop = FALSE]
  radius <- a[, 3]
  vx <- a[, 4]; vy <- a[, 5]; theta <- a[, 6]
  start <- center + cbind(radius * (vx * cos(theta) - vy * sin(theta)),
                          radius * (vx * sin(theta) + vy * cos(theta)))
  end <- center + cbind(radius * (vx * cos(theta) + vy * sin(theta)),
                        radius * (-vx * sin(theta) + vy * cos(theta)))
  closure_error <- sqrt(sum(colSums(end - start)^2))
  stopifnot(closure_error < 1e-8 * max(1, sum(2 * radius * theta)))
  # A translation of the origin leaves a closed line integral unchanged and
  # reduces cancellation when coordinates are far from zero.
  center <- sweep(center, 2, colMeans(hull$xahull), "-")
  area <- sum((center[, 1] * (end[, 2] - start[, 2]) -
                 center[, 2] * (end[, 1] - start[, 1]) - 2 * radius^2 * theta) / 2)
  stopifnot(is.finite(area), area >= 0)
  area
}

# Equivalent endpoint/chord form of the same line integral. Exact vertex-ID
# balance establishes closure without reconstructing endpoints from rounded
# center/bisector parameters. Each clockwise arc contributes its oriented chord
# area minus r^2 * (2 * theta - sin(2 * theta)) / 2.
alpha_hull_endpoint_integral <- function(hull) {
  stopifnot(inherits(hull, "ahull"))
  a <- hull$arcs
  if (is.null(a) || nrow(a) == 0) return(0)
  a <- a[a[, 3] > 0, , drop = FALSE]
  if (nrow(a) == 0) return(0)
  stopifnot(all(is.finite(a)), all(a[, 6] >= 0), all(a[, 6] <= pi + 1e-10))
  i <- as.integer(a[, 7]); j <- as.integer(a[, 8])
  stopifnot(all(i >= 1), all(j >= 1), all(i <= nrow(hull$xahull)), all(j <= nrow(hull$xahull)))
  p <- hull$xahull[i, , drop = FALSE]
  q <- hull$xahull[j, , drop = FALSE]
  center <- a[, 1:2, drop = FALSE]
  tolerance <- 1e-7 * pmax(1, a[, 3])
  stopifnot(all(abs(sqrt(rowSums((p - center)^2)) - a[, 3]) <= tolerance),
            all(abs(sqrt(rowSums((q - center)^2)) - a[, 3]) <= tolerance))
  cross <- a[, 4] * (p[, 2] - center[, 2]) - a[, 5] * (p[, 1] - center[, 1])
  flip <- cross < 0
  start_id <- i; end_id <- j
  start_id[flip] <- j[flip]; end_id[flip] <- i[flip]
  balance <- tabulate(start_id, nrow(hull$xahull)) - tabulate(end_id, nrow(hull$xahull))
  if (any(balance != 0)) stop(paste("Unbalanced boundary vertices:", sum(balance != 0)))
  origin <- colMeans(hull$xahull)
  start <- sweep(hull$xahull[start_id, , drop = FALSE], 2, origin, "-")
  end <- sweep(hull$xahull[end_id, , drop = FALSE], 2, origin, "-")
  area <- sum((start[, 1] * end[, 2] - start[, 2] * end[, 1] -
                 a[, 3]^2 * (2 * a[, 6] - sin(2 * a[, 6]))) / 2)
  stopifnot(is.finite(area), area >= 0)
  area
}
