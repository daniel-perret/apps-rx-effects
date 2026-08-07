#### Multi-scale terrain roughness and rugosity metrics
# 
# Computes standard deviation of elevation (rugosity) and Terrain Ruggedness Index (TRI)
# at multiple scales across the AOI. 
#
# Data: Elevation band ("elev") from combined topo_stack
# Output: Multi-scale roughness raster stack


library(terra)

# ===== Load Data =====

# Load topo_stack and extract elevation and slope bands
topo_stack <- rast("data/topography/topo_stack.tif")
dem <- topo_stack[["elev"]]
slope <- topo_stack[["slp.dg"]]

cat("DEM (elevation) loaded from topo_stack band 'elev'\n")
cat("  CRS:", crs(dem), "\n")
cat("  Resolution:", res(dem), "\n")

# ===== Configuration =====

# Window sizes (must be odd integers, in cells)
# Scales chosen for landscape heterogeneity relevant to fire/forest responses
scales <- c(3, 9, 35)

# ===== Helper Function: Focal Standard Deviation =====

focal_sd <- function(x, w) {
  # x: SpatRaster
  # w: odd integer window size (in cells)
  
  if (w %% 2 == 0) stop("Window size must be odd")
  
  m <- matrix(1, nrow = w, ncol = w)
  
  focal(x, w = m, fun = function(v, ...) {
    # v is a numeric vector of cell values (in focal order: top-left to bottom-right)
    if (all(is.na(v))) return(NA_real_)
    stats::sd(v, na.rm = TRUE)
  }, na.policy = "omit", fillvalue = NA_real_)
}

# ===== Compute Multi-Scale Roughness (Elevation-Based) =====

cat("\nComputing multi-scale elevation roughness...\n")

rough_list <- list()

for (w in scales) {
  cat("  Window size:", w, "cells (~", round(w * res(dem)[1], 0), "m)\n")

  # Elevation roughness: standard deviation of elevation in focal window
  # Captures vertical variation across the window
  r_dem <- focal_sd(dem, w)
  names(r_dem) <- paste0("rough_dem_w", w)

  # Slope roughness: standard deviation of slope angle in focal window
  # Captures variation in slope gradient within window
  r_slope <- focal_sd(slope, w)
  names(r_slope) <- paste0("rough_slope_w", w)

  rough_list[[paste0("w", w)]] <- c(r_dem, r_slope)
}

rough_stack <- do.call(c, rough_list)

# ===== Compute Multi-Scale Terrain Ruggedness Index (TRI) =====
#
# TRI = mean absolute difference between center cell and neighboring cells
# Independent of window size interpretation (just captures elevation variation)
# Lower computational load than SD-based approaches

cat("\nComputing multi-scale TRI...\n")

tri_list <- list()

for (w in scales) {
  cat("  Window size:", w, "cells (~", round(w * res(dem)[1], 0), "m)\n")
  
  m <- matrix(1, nrow = w, ncol = w)
  
  tri <- focal(dem, w = m, fun = function(v, ...) {
    # v is flattened focal window (length = w * w)
    # For square window, center is at position ceiling(w^2 / 2)
    
    k <- ceiling(length(v) / 2)  # center index
    c0 <- v[k]
    
    if (is.na(c0)) return(NA_real_)
    
    # Exclude center; compute mean absolute deviation
    nb <- v[-k]
    if (all(is.na(nb))) return(NA_real_)
    
    mean(abs(nb - c0), na.rm = TRUE)
  }, na.policy = "omit", fillvalue = NA_real_)
  
  names(tri) <- paste0("tri_w", w)
  tri_list[[paste0("w", w)]] <- tri
}

tri_stack <- do.call(c, tri_list)

# ===== Combine All Outputs =====

roughness_stack <- c(rough_stack, tri_stack)

cat("\nRoughness stack summary:\n")
print(roughness_stack)

# ===== Write Output =====

output_path <- "data/topography/roughness_multiscale.tif"
writeRaster(roughness_stack, output_path, overwrite = TRUE)

cat("\nOutput written to:", output_path, "\n")
