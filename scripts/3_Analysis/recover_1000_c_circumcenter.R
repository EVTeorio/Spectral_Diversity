# Test a local numerical correction; installed packages and point set stay intact.
.libPaths(c('C:/Users/PaintRock/AppData/Local/R/win-library/4.2', .libPaths()))
source('scripts/3_Analysis/alpha_hull_arc_integral.R')
out <- 'reports/tables/20261015_alpha_half_small_scales/1000_c_numeric_geometry'
saved <- readRDS(file.path(out, 'triangulation.rds'))
p <- saved$points
ns <- asNamespace('alphahull')
standard_circum <- get('circum', ns)
records <- list()
for (method in c('translated_determinant', 'translated_linear_solve')) {
  private <- new.env(parent = ns)
  private$tri.mesh <- function(x, y = NULL, ...) saved$triangulation
  replacements <- 0L
  private$circum <- function(x, y) {
    ans <- standard_circum(x, y)
    if (all(is.finite(c(ans$x, ans$y)))) return(ans)
    replacements <<- replacements + 1L
    xy <- cbind(x, y)
    u <- xy[2, ] - xy[1, ]; v <- xy[3, ] - xy[1, ]
    determinant <- u[1] * v[2] - u[2] * v[1]
    stopifnot(is.finite(determinant), determinant != 0)
    offset <- if (method == 'translated_determinant') {
      c(sum(u^2) * v[2] - sum(v^2) * u[2],
        u[1] * sum(v^2) - v[1] * sum(u^2)) / (2 * determinant)
    } else solve(2 * rbind(u, v), c(sum(u^2), sum(v^2)), tol = 1e-14)
    center <- xy[1, ] + offset
    radii <- sqrt(rowSums(sweep(xy, 2, center)^2))
    stopifnot(all(is.finite(center)), diff(range(radii)) < 1e-12 * max(radii))
    ans$x <- center[1]; ans$y <- center[2]; ans$radius <- mean(radii)
    ans
  }
  local_delvor <- get('delvor', ns)
  environment(local_delvor) <- private
  record <- data.frame(method = method, replacements = NA_integer_,
    mesh_nonfinite = NA_integer_, package_area = NA_real_, arc_area = NA_real_,
    endpoint_area = NA_real_, warnings = '', error = '', arc_error = '', endpoint_error = '')
  withCallingHandlers(tryCatch({
    dd <- local_delvor(p[, 1], p[, 2])
    stopifnot(identical(unname(dd$x), unname(p)), identical(dd$tri.obj, saved$triangulation))
    record$replacements <- replacements
    record$mesh_nonfinite <- sum(!is.finite(dd$mesh))
    stopifnot(replacements == 1, record$mesh_nonfinite == 0)
    saveRDS(dd, file.path(out, paste0('corrected_delvor_', method, '.rds')))
    h <- alphahull::ahull(dd, alpha = .5)
    saveRDS(h, file.path(out, paste0('corrected_hull_', method, '.rds')))
    record$package_area <- tryCatch(alphahull::areaahull(h, timeout = 120),
      error = function(e) {record$error <<- conditionMessage(e); NA_real_})
    record$arc_area <- tryCatch(alpha_hull_arc_integral(h),
      error = function(e) {record$arc_error <<- conditionMessage(e); NA_real_})
    record$endpoint_area <- tryCatch(alpha_hull_endpoint_integral(h),
      error = function(e) {record$endpoint_error <<- conditionMessage(e); NA_real_})
  }, error = function(e) record$error <<- conditionMessage(e)), warning = function(w) {
    record$warnings <<- paste(record$warnings, conditionMessage(w), sep = ' | ')
    invokeRestart('muffleWarning')
  })
  records[[method]] <- record
  write.csv(do.call(rbind, records), file.path(out, 'circumcenter_recovery_trials.csv'), row.names = FALSE)
  print(record, digits = 16); flush.console()
}
