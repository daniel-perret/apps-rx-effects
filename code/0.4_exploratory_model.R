#### Exploratory random forest: canopy severity (Severity_1) from terrain,
#### fire-weather, and vegetation predictors

#roughness_stack <- rast("data/topography/")

# combined raster predictor stack
#terrain_stack <- c(topo_stack, roughness_stack)

# fire weather summary
final.weather <- read.csv("data/fire_weather_summaries/fireweather_summary_RxAllEntries_110226.csv", header = T, stringsAsFactors = F) %>% 
  filter(date_confirm == "y")

# ---- 1. load data ----

# Rx outcome polygons (canopy severity), all burn entries
# (built by 0.1_burn_polygon_prep.R)
rx.shp <- sf::read_sf("data/burn_polygons/rx_canopy_all_entries.gpkg")

# strict year filter (EDIT HERE): keep only polygons whose shapefile burn year
# equals the burn history year, so weather is joined to the right burn.
# Set FALSE to keep mismatches.
require_year_match <- TRUE
if (require_year_match) {
  rx.shp <- rx.shp %>% filter(year_check == "ok")
}

# terrain/topographic predictors
terrain_stack <- rast("data/topography/topo_stack.tif")

# winter NDVI predictors
eric_stack <- rast("data/winter_NDVI/NDVI_bins_mosaic.tif")

# canopy cover
cc_stack <- rast("data/RAP/TreeCover_mosaic.tif")

# insect/disease disturbance (binary, one layer per year)
insdis_stack <- rast("data/LCMS/lcms_v2/LCMS_InsectDisease_binary_2000_2025_mosaic.tif")

# ecozones
ecozones <- rast('data/ecozones/merged_ecozones.tif') %>% 
  terra::project(., terrain_stack, method = "near")

# multi-scale roughness predictors

# unit geometry summary
unit.geo <- read.csv("data/unit_geometry.csv")

# ---- 1b. unit-level TPI SD (per burn unit) ----
tpi_layer <- terrain_stack[["tpi"]]

# dissolve to one polygon per unit
rx_units <- rx.shp %>%
  group_by(burn_name) %>%
  summarise(do_union = TRUE, .groups = "drop")

unit_tpi_sd <- exactextractr::exact_extract(tpi_layer, rx_units, fun = "stdev") %>%
  tibble(tpi_sd = ., burn_name = rx_units$burn_name) %>%
  select(burn_name, tpi_sd)

# ---- 2. balanced point sample per severity class ----

set.seed(7984)

# points per severity class x burn entry stratum (EDIT HERE)
n_per_stratum <- 500

# polygons with a weather summary for the same unit AND entry
rx_pool <- rx.shp %>%
  filter(Severity_1 %in% 0:2,
         paste(burn_name, burn_entry) %in%
           paste(final.weather$fire_name, final.weather$burn_entry))

# sample within each class x entry stratum and join attributes only against
# that stratum's polygons (entries overlap spatially, so a join against the
# full layer would duplicate points)
balanced_sample <- rx_pool %>%
  group_split(Severity_1, burn_entry) %>%
  map(function(g) {
    st_sample(g, size = n_per_stratum, type = "random") %>%
      st_as_sf() %>%
      st_join(g %>% select(burn_name, Burn_Yea_1, Severity_1, burn_entry),
              join = st_within, left = TRUE)
  }) %>%
  bind_rows()

# extract all raster layers (annual layers are filtered after the join) and
# coordinates at each point
balanced_sample <- balanced_sample %>% 
  bind_cols(
    terra::extract(terrain_stack, vect(balanced_sample), ID = F),
    terra::extract(eric_stack, vect(balanced_sample), ID = F),
    terra::extract(cc_stack, vect(balanced_sample), ID = F),
    terra::extract(insdis_stack, vect(balanced_sample), ID = F),
    terra::extract(ecozones, vect(balanced_sample), ID = F),
    st_coordinates(balanced_sample) %>% as_tibble()
  )

# ---- 3. join fire weather vars by burn name and filter annual vars ----

model_data <- balanced_sample %>% 
  sf::st_drop_geometry() %>% 
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
    ndvi_prior_year = burn_year - 1,
    ndvi_col = paste0("y_", ndvi_prior_year),
    cc_prior_year = burn_year - 1,
    cc_col = paste0(cc_prior_year, "_TRE_", cc_prior_year),
    # factor (not character) so randomForest treats it as categorical
    ecozone = factor(as.character(merged_ecozones))
  )

# NDVI value from the year before the burn
ndvi_year_cols <- grep("^y_[0-9]{4}$", names(model_data), value = TRUE)
ndvi_matrix <- as.matrix(model_data[ndvi_year_cols])
model_data$ndvi_winter_prior_year <- ndvi_matrix[
  cbind(seq_len(nrow(model_data)), match(model_data$ndvi_col, ndvi_year_cols))
]

# NDVI values are categorical bins
model_data$ndvi_winter_prior_year <- factor(as.character(model_data$ndvi_winter_prior_year))

# NDVI values are categorical bins
model_data$ndvi_winter_prior_year <- factor(as.character(model_data$ndvi_winter_prior_year))

# canopy cover value from the year before the burn
cc_year_cols <- grep("^[0-9]{4}_TRE_[0-9]{4}$", names(model_data), value = TRUE)
cc_matrix <- as.matrix(model_data[cc_year_cols])
model_data$canopy_cover_prior_year <- cc_matrix[
  cbind(seq_len(nrow(model_data)), match(model_data$cc_col, cc_year_cols))
]

# insect/disease flag: 1 if any disturbance in the 5 years before the burn
# (uses available years if the window is truncated by the data start)
insdis_year_cols <- grep("^insdis_[0-9]{4}$", names(model_data), value = TRUE)
insdis_matrix <- as.matrix(model_data[insdis_year_cols])
insdis_col_years <- as.integer(gsub("^insdis_", "", insdis_year_cols))

model_data$insect_disease_5yr <- vapply(seq_len(nrow(model_data)), function(i) {
  if (is.na(model_data$burn_year[i])) {
    return(NA_integer_)
  }
  window_years <- (model_data$burn_year[i] - 5):(model_data$burn_year[i] - 1)
  cols_in_window <- which(insdis_col_years %in% window_years)
  if (length(cols_in_window) == 0) {
    return(NA_integer_)
  }
  as.integer(any(insdis_matrix[i, cols_in_window] == 1, na.rm = TRUE))
}, integer(1))

model_data <- model_data %>% 
  select(-all_of(ndvi_year_cols), -ndvi_col, -burn_year, -ndvi_prior_year,
         -all_of(cc_year_cols), -cc_col, -cc_prior_year, -all_of(insdis_year_cols),
         - merged_ecozones) %>% 
  left_join(unit.geo,by = "burn_name")

# inspect join match rate / missingness before modeling
model_data %>% filter(if_any(everything(), is.na)) %>% view()


# ---- 4. random forest ----

library(randomForest)

# predictors used in the model (edit here only)
predictors <- c(
  "chili", "tpi", "tpi_sd", "slp.dg", "elev",
  "burn_doy", "burn_entry",
  "area_ha", "perimeter_m",
  "ecozone", "insect_disease_5yr",
  "daily_max_temp", "daily_minRH",
  "daily_mean_temp", "daily_meanRH",
  "days_since_precip", "ts_MAP",
  "daily_mean_temp_aft", "daily_minRH_aft",
  "days_since_precip_025",
  "X1HrFM", "X10HrFM", "X100HrFM", "X1000HrFM",
  "HerbFM", "WoodyFM",
  "KBDI", "ERC",
  "rh_trend", "temp_trend", "KBDI_trend", "ERC_trend",
  "ndvi_winter_prior_year",
  "canopy_cover_prior_year"
)

# predictors used in the model (edit here only)
predictors <- c(
  "chili", "tpi", "tpi_sd", "slp.dg", "elev",
  "burn_doy", "burn_entry",
  "area_ha", "perimeter_m",
  "ecozone", "insect_disease_5yr",
  "daily_max_temp", "daily_minRH",
  "days_since_precip_025",
  "X10HrFM","KBDI",
  "rh_trend", "temp_trend", "KBDI_trend",
  "ndvi_winter_prior_year",
  "canopy_cover_prior_year"
)

predictors <- c(
  "chili", "tpi", "tpi_sd", "slp.dg", "elev",
  "burn_entry",
  "ecozone", "insect_disease_5yr",
  "ndvi_winter_prior_year",
  "canopy_cover_prior_year"
)


# model data
rf_data <- model_data %>% 
  left_join(unit_tpi_sd, by = "burn_name") %>% 
  select(Severity_1, burn_name, date, all_of(predictors)) %>% 
  na.omit() %>% 
  as.data.frame()

# model
set.seed(117)
rf1 <- randomForest::randomForest(reformulate(predictors, response = "Severity_1"),
                                  data = rf_data,
                                  ntree = 1000,
                                  mtry = 6,
                                  type = "classification")

rf1
summary(rf1)

randomForest::varImpPlot(rf1)

# ---- 5. model fit diagnostics ----

# OOB confusion matrix + per-class error
confusion_table <- rf1$confusion %>%
  as.data.frame() %>%
  rownames_to_column("Predicted") %>%
  rename(Error_Rate = class.error) %>%
  mutate(Error_Rate = round(Error_Rate, 4))

confusion_table

oob_accuracy <- 1 - rf1$err.rate[nrow(rf1$err.rate), "OOB"]
oob_accuracy

# OOB error vs. number of trees

err_df <- as.data.frame(rf1$err.rate) %>% 
  mutate(trees = row_number()) %>% 
  pivot_longer(-trees, names_to = "class", values_to = "error")

ggplot(err_df, aes(trees, error, color = class)) +
  geom_line() +
  labs(title = "OOB error rate vs. number of trees",
       x = "Number of trees", y = "OOB error rate", color = NULL)

# variable importance
importance_df <- importance(rf1) %>% 
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
    do.call(partialPlot, list(x = rf1, pred.data = rf_data, x.var = v,
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
      rf1,
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

