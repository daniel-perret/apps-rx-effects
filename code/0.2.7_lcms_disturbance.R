#### Mosaic LCMS disturbance tiles (GEE export)
#
# data/LCMS contains tiled GeoTIFF exports from Google Earth Engine
# (LCMS_Change_2000_2025_Stack-<row>-<col>.tif). Tiles are adjacent,
# non-overlapping, and share a common grid, so they are combined with
# terra::merge() rather than exactextractr-style averaging.
#
# Output: data/LCMS/LCMS_Change_2000_2025_mosaic.tif

library(terra)

tile_files <- list.files("data/LCMS/lcms_v2", pattern = "\\.tif$", full.names = TRUE)

tiles <- lapply(tile_files, rast)


lcms_mosaic <- sprc(tiles) %>%
  mosaic()

aoi <- vect("data/aoi.shp") %>% 
  terra::project(., crs(lcms_mosaic))

lcms_mosaic_clip <- crop(lcms_mosaic, aoi, mask = T)

writeRaster(lcms_mosaic_clip, "data/LCMS/lcms_v2/LCMS_Change_2000_2025_mosaic.tif", overwrite = TRUE)

#### Annual binary insect/disease disturbance masks
#
# lcms_mosaic_clip is the LCMS "Change" product (16 classes per year; see
# LCMS codebook). We reclassify each annual layer into a binary mask of
# insect/disease-related disturbance, keeping one layer per year. The 5-year
# "has this point been disturbed in the past 5 years" rollup is NOT done
# here -- it happens downstream (in the modeling script) after point-level
# extraction, once burn year is known, mirroring how ndvi_winter_prior_year
# and canopy_cover_prior_year are resolved in 0.4_exploratory_model.R.
#
# Insect/disease-related Change classes:
#   10 = Defoliation
#   11 = Southern Pine Beetle
#   12 = Insect and Disease or Drought Stress
#
# Output: data/LCMS/lcms_v2/LCMS_InsectDisease_binary_2000_2025_mosaic.tif
# (one band per year, named insdis_YYYY; 1 = insect/disease disturbance,
# 0 = no insect/disease disturbance that year, NA preserved from source,
# including the Non-Processing Area Mask class 16)

rcl <- matrix(c(1,0,
                2,0,
                3,0,
                4,0,
                5,0,
                6,0,
                7,0,
                8,0,
                9,0,
                10,1,
                11,1,
                12,1,
                13,0,
                14,0,
                15,0,
                16,0), ncol = 2, byrow = T)

lcms_insdis_binary <- terra::classify(lcms_mosaic_clip, rcl)

names(lcms_insdis_binary) <- gsub("^Change_", "insdis_", names(lcms_mosaic_clip))

terrain_stack <- rast("data/topography/topo_stack.tif")

lcms_insdis_binary <- terra::project(lcms_insdis_binary, terrain_stack, method = "near")

writeRaster(
  lcms_insdis_binary,
  "data/LCMS/lcms_v2/LCMS_InsectDisease_binary_2000_2025_mosaic.tif",
  overwrite = TRUE
)
