#### Burn unit geometry: perimeter and area
#
# Computes the perimeter and area of each burn unit from the outcome
# polygons (1burn_cnpy.shp). That shapefile has multiple sub-polygons per
# unit (one per canopy severity class), so units are dissolved to a single
# polygon per burn_name before measuring.
#
# Data: data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel/1burn_cnpy.shp
# Output: data/unit_geometry.csv (burn_name, perimeter_m, area_m2, area_ha)

library(tidyverse)
library(sf)

# ---- 1. load + dissolve ----

rx.shp <- sf::read_sf("data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel/1burn_cnpy.shp")

rx_units <- rx.shp %>%
  group_by(burn_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ---- 2. project to an equal-area CRS for accurate length/area ----
# Albers Equal Area (CONUS) keeps area/perimeter distortion low for VA/WV extent.

rx_units_proj <- st_transform(rx_units, crs = 5070)

# ---- 3. compute perimeter and area ----

unit_geometry <- rx_units_proj %>%
  st_drop_geometry() %>%
  mutate(
    perimeter_m = as.numeric(st_length(st_boundary(rx_units_proj))),
    area_m2 = as.numeric(st_area(rx_units_proj)),
    area_ha = area_m2 / 10000,
    perimeter_mi = perimeter_m / 1609.344,
    area_acres = area_m2 / 4046.8564224
  ) %>%
  select(burn_name, perimeter_m, perimeter_mi, area_m2, area_ha, area_acres)

# ---- 4. write output ----

write.csv(unit_geometry, "data/unit_geometry.csv", row.names = FALSE)
