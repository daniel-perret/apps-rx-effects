### gathers, processes, summarizes topographic and terrain-related variables

# AOI shapefile

aoi <- vect("data/aoi.shp")

# CHILI

chili <- rast("data/topography/CHILI_clipped.tif")

# TPI

tpi <- rast("data/topography/TPI_clip_all.tif")

# Slope angle

slp.dg <- rast("../../SHARED_DATA/LANDFIRE/LF2020_SlpD_CONUS/LF2020_SlpD_CONUS/Tif/LF2020_SlpD_CONUS.tif")
slp.dg <- crop(slp.dg, 
               aoi %>% 
                 terra::project(crs(slp.dg)),
               mask = T)
slp.dg <- terra::project(slp.dg, crs(tpi))
slp.dg <- crop(slp.dg, tpi)
slp.dg <- resample(slp.dg, chili, method = "bilinear")

writeRaster(slp.dg, "data/topography/slope_degrees.tif", overwrite = T)

# DEM

dem <- rast("../../SHARED_DATA/LANDFIRE/LF2020_Elev_CONUS/LF2020_Elev_CONUS/Tif/LF2020_Elev_CONUS.tif")
dem <- crop(dem, 
               aoi %>% 
                 terra::project(crs(dem)),
               mask = T)
dem <- terra::project(dem, crs(tpi))
dem <- crop(dem, tpi)
dem <- resample(dem, chili, method = "bilinear")

writeRaster(dem, "data/topography/slope_degrees.tif", overwrite = T)

# Stack 'em
topo_stack <- c(chili, tpi, slp.dg, dem)
names(topo_stack) <- c("chili","tpi","slp.dg", "elev")

writeRaster(topo_stack, "data/topography/topo_stack.tif", overwrite = TRUE)
