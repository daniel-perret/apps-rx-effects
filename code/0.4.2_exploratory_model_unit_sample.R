#### Exploratory random forest, UNIT-SAMPLED training design
#
# Variant of 0.4_exploratory_model.R with a different training sample:
#
#   0.4  (old): n points per severity class x burn entry stratum, i.e. classes
#               are forced to equal size and units contribute in proportion to
#               the area of each class they hold.
#   0.4.2 (new): n points per burn EVENT (burn_name x burn_entry), drawn at
#               random over the whole event footprint. Every event gets the same
#               weight, and class mix within an event is NOT balanced -- it
#               follows polygon area -- so between-unit variation in class
#               composition is preserved and unit-level predictors can carry
#               signal.
#
# Everything downstream (feature construction, predictor set, rf settings) is
# kept identical to 0.4 so results can be compared. New objects carry a
# `_unit` suffix so the old pipeline's objects (rf_data, rf_unit, ...) are not
# overwritten:
#   unit_sample, model_data_unit, rf_data_unit, predictors_unit, rf_unit
# Shared inputs (rx_pool, rasters, final.weather, unit_tpi_sd) are re-created
# here and are needed by 0.4.3_blocked_cv_unit_validation.R.
#
# Note: the class mix of this sample is imbalanced (natural area shares), so
# OOB accuracy here is not comparable to 0.4's balanced OOB accuracy. Judge
# models with balanced accuracy / kappa and the unit-level checks in 0.4.3.

library(tidyverse)
library(sf)
library(terra)
library(randomForest)

# ---- 1. load data ----

# fire weather summary (used here for burn date and for the event filter)
final.weather <- read.csv("data/fire_weather_summaries/fireweather_summary_RxAllEntries_110226.csv",
                          header = TRUE, stringsAsFactors = FALSE) %>%
  filter(date_confirm == "y")

# Rx outcome polygons (canopy severity), all burn entries
# (built by 0.1_burn_polygon_prep.R)
rx.shp <- sf::read_sf("data/burn_polygons/rx_canopy_all_entries.gpkg")

# strict year filter (EDIT HERE): keep only polygons whose shapefile burn year
# equals the burn history year, so weather is joined to the right burn.
require_year_match <- TRUE
if (require_year_match) {
  rx.shp <- rx.shp %>% filter(year_check == "ok")
}

terrain_stack <- rast("data/topography/topo_stack.tif")
eric_stack    <- rast("data/winter_NDVI/NDVI_bins_mosaic.tif")
cc_stack      <- rast("data/RAP/TreeCover_mosaic.tif")
insdis_stack  <- rast("data/LCMS/lcms_v2/LCMS_InsectDisease_binary_2000_2025_mosaic.tif")
ecozones      <- rast("data/ecozones/merged_ecozones.tif") %>%
  terra::project(., terrain_stack, method = "near")

# unit geometry summary (area_ha, perimeter_m, ...)
unit.geo <- read.csv("data/unit_geometry.csv")

# ---- 1b. unit-level TPI SD (per burn unit) ----

tpi_layer <- terrain_stack[["tpi"]]

rx_units <- rx.shp %>%
  group_by(burn_name) %>%
  summarise(do_union = TRUE, .groups = "drop")

unit_tpi_sd <- exactextractr::exact_extract(tpi_layer, rx_units, fun = "stdev") %>%
  tibble(tpi_sd = ., burn_name = rx_units$burn_name) %>%
  select(burn_name, tpi_sd)

# ---- 2. equal-n point sample per burn event ----

set.seed(3917)

# points per burn event = burn_name x burn_entry (EDIT HERE)
n_per_event <- 100

# polygons with a weather summary for the same unit AND entry, in the terrain
# raster CRS so point extraction below is CRS-consistent
rx_pool <- rx.shp %>%
  filter(Severity_1 %in% 0:2,
         paste(burn_name, burn_entry) %in%
           paste(final.weather$fire_name, final.weather$burn_entry)) %>%
  st_transform(terra::crs(terrain_stack))

# Sample over each event's whole footprint (all severity polygons together), so
# class shares within the event follow polygon area. Attributes are joined only
# against that event's polygons (entries overlap spatially, so a join against
# the full layer would duplicate points). Points landing in the same terrain
# cell are collapsed so small events don't repeat identical pixels.
unit_sample <- rx_pool %>%
  group_split(burn_name, burn_entry) %>%
  map(function(g) {
    pts <- st_sample(g, size = n_per_event, type = "random") %>%
      st_as_sf() %>%
      st_join(g %>% select(burn_name, Burn_Yea_1, Severity_1, burn_entry),
              join = st_within, left = TRUE)
    pts$cell <- terra::cellFromXY(terrain_stack, st_coordinates(pts))
    pts %>% distinct(cell, .keep_all = TRUE)
  }) %>%
  bind_rows()

unit_sample <- unit_sample %>%
  bind_cols(
    terra::extract(terrain_stack, vect(unit_sample), ID = FALSE),
    terra::extract(eric_stack, vect(unit_sample), ID = FALSE),
    terra::extract(cc_stack, vect(unit_sample), ID = FALSE),
    terra::extract(insdis_stack, vect(unit_sample), ID = FALSE),
    terra::extract(ecozones, vect(unit_sample), ID = FALSE)
  )

# ---- 3. join burn date and build annual predictors ----

model_data_unit <- unit_sample %>%
  sf::st_drop_geometry() %>%
  # all fire weather / fuels variables, joined by unit AND entry
  left_join(final.weather, by = c("burn_name" = "fire_name", "burn_entry")) %>%
  mutate(
    Severity_1 = case_when(
      Severity_1 == 0 ~ "Closed",
      Severity_1 == 1 ~ "Open",
      Severity_1 == 2 ~ "Early Seral"
    ) %>%
      factor(levels = c("Closed", "Open", "Early Seral")),
    burn_year = as.integer(format(as.Date(date), "%Y")),
    burn_doy = as.integer(format(as.Date(date), "%j")),
    ndvi_col = paste0("y_", burn_year - 1),
    cc_col = paste0(burn_year - 1, "_TRE_", burn_year - 1),
    # factor (not character) so randomForest treats it as categorical
    ecozone = factor(as.character(merged_ecozones))
  )

# NDVI value from the year before the burn (categorical bins)
ndvi_year_cols <- grep("^y_[0-9]{4}$", names(model_data_unit), value = TRUE)
ndvi_matrix <- as.matrix(model_data_unit[ndvi_year_cols])
model_data_unit$ndvi_winter_prior_year <- ndvi_matrix[
  cbind(seq_len(nrow(model_data_unit)), match(model_data_unit$ndvi_col, ndvi_year_cols))
] %>%
  as.character() %>%
  factor()

# canopy cover value from the year before the burn
cc_year_cols <- grep("^[0-9]{4}_TRE_[0-9]{4}$", names(model_data_unit), value = TRUE)
cc_matrix <- as.matrix(model_data_unit[cc_year_cols])
model_data_unit$canopy_cover_prior_year <- cc_matrix[
  cbind(seq_len(nrow(model_data_unit)), match(model_data_unit$cc_col, cc_year_cols))
]

# insect/disease flag: 1 if any disturbance in the 5 years before the burn
# (uses available years if the window is truncated by the data start)
insdis_year_cols <- grep("^insdis_[0-9]{4}$", names(model_data_unit), value = TRUE)
insdis_matrix <- as.matrix(model_data_unit[insdis_year_cols])
insdis_col_years <- as.integer(gsub("^insdis_", "", insdis_year_cols))

model_data_unit$insect_disease_5yr <- vapply(seq_len(nrow(model_data_unit)), function(i) {
  if (is.na(model_data_unit$burn_year[i])) {
    return(NA_integer_)
  }
  window_years <- (model_data_unit$burn_year[i] - 5):(model_data_unit$burn_year[i] - 1)
  cols_in_window <- which(insdis_col_years %in% window_years)
  if (length(cols_in_window) == 0) {
    return(NA_integer_)
  }
  as.integer(any(insdis_matrix[i, cols_in_window] == 1, na.rm = TRUE))
}, integer(1))

model_data_unit <- model_data_unit %>%
  select(-all_of(ndvi_year_cols), -ndvi_col, -burn_year,
         -all_of(cc_year_cols), -cc_col, -all_of(insdis_year_cols),
         -merged_ecozones) %>%
  left_join(unit.geo, by = "burn_name")

# ---- 4. random forest ----

# predictors used in the model (edit here only). All candidates are available
# in model_data_unit; comment lines out to drop them. Unit/event-level variables
# (tpi_sd, burn_doy, burn_entry, area_ha, perimeter_m and all fire weather /
# fuels variables) are constant within a burn event; 0.4.3 picks them up
# automatically from this list.
predictors_unit <- c(
  # pixel-level
  "chili", "tpi", "slp.dg", "elev",
  "ecozone", "burn_doy", "area_ha", "X10HrFM",
  "ndvi_winter_prior_year", "canopy_cover_prior_year",
  # unit/event-level: geometry and timing
  "tpi_sd", "burn_entry", 
  # unit/event-level: fire weather and fuels
  "days_since_precip", 
  "daily_mean_temp_aft", "daily_meanRH_aft",
  "KBDI")

#predictors_unit <- predictors

rf_data_unit <- model_data_unit %>%
  left_join(unit_tpi_sd, by = "burn_name") %>%
  select(Severity_1, burn_name, date, all_of(predictors_unit)) %>%
  na.omit() %>%
  as.data.frame()

# rows per event after NA removal -- should be close to n_per_event
rf_data_unit %>%
  count(burn_name, burn_entry) %>%
  summarise(events = n(), min_n = min(n), median_n = median(n), max_n = max(n))

### MODEL -----

samp.size <- c("Closed" = 397,
               "Open" = 245,
               "Early Seral" = 258)

set.seed(2764)
rf_unit <- randomForest::randomForest(reformulate(predictors_unit, response = "Severity_1"),
                                      data = rf_data_unit,
                                      ntree = 1000,
                                      mtry = 6,
                                      sampsize = samp.size,
                                      type = "classification")

rf_unit

randomForest::varImpPlot(rf_unit)

# OOB confusion matrix + per-class error
confusion_table <- rf_unit$confusion %>%
  as.data.frame() %>%
  rownames_to_column("Predicted") %>%
  rename(Error_Rate = class.error) %>%
  mutate(Error_Rate = round(Error_Rate, 4))

confusion_table

oob_accuracy <- 1 - rf_unit$err.rate[nrow(rf_unit$err.rate), "OOB"]
oob_accuracy

# OOB error vs. number of trees

err_df <- as.data.frame(rf_unit$err.rate) %>% 
  mutate(trees = row_number()) %>% 
  pivot_longer(-trees, names_to = "class", values_to = "error")

ggplot(err_df, aes(trees, error, color = class)) +
  geom_line() +
  labs(title = "OOB error rate vs. number of trees",
       x = "Number of trees", y = "OOB error rate", color = NULL)

# variable importance
importance_df <- importance(rf_unit) %>% 
  as.data.frame() %>% 
  rownames_to_column("variable") %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  mutate(
    type = case_when(
      str_detect(variable, "FM") ~ "fuels",
      str_detect(variable, "prior_year") ~ "fuels",
      str_detect(variable, "trend") ~ "weather",
      str_detect(variable, "MAP") ~ "weather",
      str_detect(variable, "daily") ~ "weather",
      str_detect(variable, "precip") ~ "weather",
      variable %in% c("ERC", "KBDI") ~ "weather",
      variable %in% c("chili", "elev", "slp.dg", "tpi", "tpi_sd") ~ "topo",
      .default = "other"
    )
  )

ggplot(importance_df, aes(y = MeanDecreaseGini, 
                          x = fct_reorder(variable, -MeanDecreaseGini))) +
  geom_col(aes(fill = type)) +
  labs(title = "Variable importance",
       y = "Mean decrease in Gini impurity", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# partial dependence for top predictors, one panel set per class
top_vars <- head(importance_df$variable, 8)

for (cl in levels(rf_data$Severity_1)) {
  par(mfrow = c(4,2))
  for (v in top_vars) {
    do.call(partialPlot, list(x = rf_unit, pred.data = rf_data, x.var = v,
                              which.class = cl,
                              main = paste0(v, " (class: ", cl, ")")))
  }
}
par(mfrow = c(1, 1))

# ---- 6. partial dependence with ICE (pdp) ----

library(pdp)

# ICE curves per top variable and class
pdp_plots <- map_df(top_vars, function(var) {
  map_df(levels(rf_data$Severity_1), function(cl) {
    partial(
      rf_unit,
      pred.var = var,
      pred.data = rf_data,
      train = rf_data,
      which.class = cl,
      grid.resolution = 50,
      ice = TRUE,
      prob = TRUE
    ) %>%
      as_tibble() %>%
      mutate(variable = var, class = cl)
  })
})

# one plot per variable, faceted by class
library(gridExtra)
pdp_plot_list <- map(top_vars, function(var) {
  data_for_var <- pdp_plots %>% 
    filter(variable == var) %>%
    mutate(
      x_val = .data[[var]],
      class = factor(class, levels = c("Closed", "Open", "Early Seral"))
    )
  
  # average ICE curves to get the PDP
  pdp_data <<- data_for_var %>%
    group_by(across(all_of(c(var, "class")))) %>%
    summarise(yhat = mean(yhat), .groups = "drop") %>%
    mutate(x_val = .data[[var]])
  
  ggplot(data_for_var, aes(x = x_val, y = yhat)) +
    # geom_line(aes(group = interaction(class, yhat.id)), 
    #           alpha = 0.1, color = "gray50") +
    geom_line(data = pdp_data, aes(color = class), linewidth = 1) +
    facet_wrap(~ class, ncol = 1) +
    labs(title = var, y = "Predicted probability", x = NULL) +
    theme(legend.position = "none") 
})

do.call(grid.arrange, c(pdp_plot_list, ncol = 4))

