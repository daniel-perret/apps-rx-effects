#### Burn polygon prep: harmonize 1burn-4burn canopy-outcome shapefiles into a
#### single multi-entry object, plus unit boundaries from Rxb_tomodel

library(tidyverse)
library(sf)

shp_dir <- "data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel"
out_dir <- "data/burn_polygons"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# common CRS (matches the original 1burn_cnpy.shp)
target_crs <- 26917 # NAD83 / UTM 17N

# ---- 1. name crosswalk (EDIT HERE) ----
# Maps polygon names (names = shapefile) to Fire_name_fixed in the burn history
# (values). Applied to all entries. Remaining mismatches are written to the
# audit table below for manual review.
name_fixes <- c(
  # spelling / case
  "Short mtn"    = "Short Mtn",
  "Gauley ridge" = "Gauley Ridge",
  # same unit under a different name; year agrees with hist
  "Cole Mtn"     = "Cole Mtn A1",
  # tomodel notes: High Butte 4 Wl was formerly Wells Branch
  "Wells Branch" = "High Butte 4 Wl"
)

# ---- 2. read + harmonize each entry's canopy shapefile ----
# Schemas differ across files (name, year, acres, zone columns), so rename
# explicitly to one schema.

read_entry <- function(entry) {
  read_sf(file.path(shp_dir, paste0(entry, "burn_cnpy.shp"))) %>%
    rename(any_of(c(
      burn_name  = "Fire_name",
      Burn_Yea_1 = "burn_yea_1",
      gis_acres  = "gisAcres",
      gis_acres  = "acres",
      zone       = "ZONE",
      zone       = "Zone"
    ))) %>%
    st_transform(target_crs) %>%
    st_make_valid() %>%
    transmute(
      burn_name  = burn_name %>% str_squish() %>% recode(!!!name_fixes),
      burn_entry = as.integer(entry),
      Burn_Yea_1 = as.integer(Burn_Yea_1),
      Severity_1 = as.integer(Severity_1),
      gis_acres  = gis_acres,
      zone       = zone,
      geometry   = geometry
    )
}

cnpy_list <- map(1:4, read_entry)

rx_all <- bind_rows(cnpy_list) %>%
  filter(!is.na(burn_name),
         Severity_1 %in% 0:2)

# ---- 3. check against burn history ----
# Severity_1 codes outside 0:2 (3, 9, 91, 92, 999) are not canopy classes used
# in the model; they are kept here but flagged.

unit.hist <- read.csv("data/Rxb_unit_history_092126.csv", stringsAsFactors = FALSE) %>%
  mutate(keep. = trimws(keep.),
         burn_entry = trimws(burn_entry),
         Fire_name_fixed = str_squish(Fire_name_fixed)) %>%
  filter(keep. == "yes", burn_entry %in% as.character(1:4)) %>%
  mutate(burn_entry = as.integer(burn_entry),
         hist_year = as.integer(format(as.Date(date_fixed, format = "%m/%d/%Y"), "%Y"))) %>%
  group_by(burn_name = Fire_name_fixed, burn_entry) %>%
  summarise(hist_year = min(hist_year, na.rm = TRUE), .groups = "drop") %>%
  mutate(hist_year = if_else(is.infinite(hist_year), NA_integer_, hist_year))

rx_all <- rx_all %>%
  left_join(unit.hist, by = c("burn_name", "burn_entry")) %>%
  mutate(
    severity_valid = Severity_1 %in% 0:2,
    in_hist        = !is.na(hist_year),
    year_check     = case_when(
      !in_hist ~ "no_hist_match",
      Burn_Yea_1 == hist_year ~ "ok",
      .default = "year_mismatch"
    )
  )

# audit tables (unit x entry level)
name_audit <- rx_all %>%
  st_drop_geometry() %>%
  group_by(burn_name, burn_entry) %>%
  summarise(shp_year = paste(sort(unique(Burn_Yea_1)), collapse = "/"),
            hist_year = first(hist_year),
            year_check = first(year_check),
            n_polys = n(),
            acres = sum(gis_acres, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(year_check, burn_name, burn_entry)

name_audit %>% count(year_check)

# canopy-class polygons present in hist (keep == yes) but missing from shapefiles
missing_polys <- unit.hist %>%
  anti_join(st_drop_geometry(rx_all), by = c("burn_name", "burn_entry"))

write.csv(name_audit, file.path(out_dir, "name_audit_cnpy_vs_hist.csv"), row.names = FALSE)
write.csv(missing_polys, file.path(out_dir, "hist_entries_without_polygons.csv"), row.names = FALSE)

# ---- 4. unit boundaries (Rxb_tomodel) ----
# Fire_name is the original name; burn_name has 5 NAs and 15 duplicates, so
# fall back to Fire_name where burn_name is missing.
rx_unit_bounds <- read_sf(file.path(shp_dir, "Rxb_tomodel.shp")) %>%
  st_transform(target_crs) %>%
  st_make_valid() %>%
  mutate(burn_name = coalesce(str_squish(burn_name), str_squish(Fire_name)),
         in_cnpy = burn_name %in% rx_all$burn_name)

rx_unit_bounds %>%
  st_drop_geometry() %>%
  filter(!in_cnpy) %>%
  select(Fire_name, burn_name, burncount)

# ---- 5. checks and write ----
rx_all %>%
  st_drop_geometry() %>%
  count(burn_entry, Severity_1) %>%
  pivot_wider(names_from = burn_entry, values_from = n, values_fill = 0)

st_write(rx_all, file.path(out_dir, "rx_canopy_all_entries.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
st_write(rx_unit_bounds, file.path(out_dir, "rx_unit_boundaries.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
