#### Exploratory random forest: predicting canopy severity outcomes from
#### terrain and fire-weather predictors
####
#### Combines:
####   - terrain/roughness rasters (0.2.3_topo_terrain.R, 0.2.3_terrain_roughness.R)
####   - fire weather summary + trends (0.2.2_raws_timeseries.R)
####   - 1burn canopy severity outcome polygons
#### into a spatially balanced point sample, then fits a random forest
#### classifier for Severity_1.

library(tidyverse)
library(sf)
library(terra)

# ---- 1. load data ----

# Rx outcome polygons (canopy severity)
rx.shp <- sf::read_sf("data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel/1burn_cnpy.shp")

# terrain/topographic predictors
terrain_stack <- rast("data/topography/topo_stack.tif")

# winter NDVI predictors
eric_stack <- rast("data/winter_NDVI/NDVI_bins_mosaic.tif")

# canopy cover
cc_stack <- rast("data/RAP/TreeCover_mosaic.tif")

# multi-scale roughness predictors
#roughness_stack <- rast("data/topography/")

# combined raster predictor stack
#terrain_stack <- c(topo_stack, roughness_stack)

# fire weather summary (per-fire, day-of-burn conditions + week-prior trends)
final.weather <- read.csv("data/fire_weather_summaries/fireweather_summary_Rx1b_082126.csv", header = T, stringsAsFactors = F) %>% 
  filter(date_confirm == "y")

# ---- 1b. unit-level TPI SD (per burn unit) ----
# SD of TPI across each burn unit (one polygon per unit, summarised in terra)
tpi_layer <- terrain_stack[["tpi"]]

rx_units <- rx.shp %>%
  group_by(burn_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

unit_tpi_sd <- terra::extract(tpi_layer, vect(rx_units), 
                              fun = sd, na.rm = TRUE) %>%
  rename(tpi_sd = 2) %>%
  mutate(burn_name = rx_units$burn_name) %>%
  select(burn_name, tpi_sd)


# ---- 2. spatially balanced sample of outcome polygons ----
# (mirrors 0.3.1_exploratory_topo.R)

set.seed(7984)

n_per_class <- 1000

balanced_sample <- rx.shp %>%
  filter(Severity_1 %in% 0:2,
         burn_name %in% final.weather$fire_name) %>%
  nest(data = -Severity_1) %>%
  mutate(
    points = map(data, ~st_sample(.x, size = n_per_class, type = "random") %>%
                   st_as_sf())
  ) %>%
  select(-data) %>%
  unnest(points) %>% 
  st_as_sf() %>% 
  # recover burn identity for each sampled point
  st_join(rx.shp %>% 
            select(burn_name, Burn_Yea_1), 
          join = st_within, left = TRUE)

# extract terrain/roughness values, all winter-NDVI (eric_stack) layers, all
# canopy-cover (cc_stack) layers, and point coordinates at each sampled
# point. eric_stack and cc_stack each have one layer per year; we don't yet
# know which layer is relevant (that depends on burn year, joined in next
# step), so pull them all for now and narrow down after the join.
balanced_sample <- balanced_sample %>% 
  bind_cols(
    terra::extract(terrain_stack, vect(balanced_sample), ID = F),
    terra::extract(eric_stack, vect(balanced_sample), ID = F),
    terra::extract(cc_stack, vect(balanced_sample), ID = F),
    st_coordinates(balanced_sample) %>% as_tibble()
  )

# ---- 3. join fire weather vars by burn name ----

model_data <- balanced_sample %>% 
  sf::st_drop_geometry() %>% 
  left_join(final.weather, by = c("burn_name" = "fire_name")) %>% 
  mutate(
    Severity_1 = case_when(
      Severity_1 == 0 ~ "Closed",
      Severity_1 == 1 ~ "Open",
      Severity_1 == 2 ~ "Early Seral"
    ) %>% 
      factor(levels = c("Closed", "Open", "Early Seral")),
    # burn year comes from the fire-weather table's burn date, not the
    # shapefile's Burn_Yea_1 field
    burn_year = as.integer(format(as.Date(date), "%Y")),
    burn_doy = as.integer(format(as.Date(date), "%j")),
    ndvi_prior_year = burn_year - 1,
    ndvi_col = paste0("y_", ndvi_prior_year),
    cc_prior_year = burn_year - 1,
    cc_col = paste0(cc_prior_year, "_TRE_", cc_prior_year)
  )

# pull the single eric_stack layer value that corresponds to each row's
# year-prior-to-fire, then drop the now-unneeded per-year NDVI columns
ndvi_year_cols <- grep("^y_[0-9]{4}$", names(model_data), value = TRUE)
ndvi_matrix <- as.matrix(model_data[ndvi_year_cols])
model_data$ndvi_winter_prior_year <- ndvi_matrix[
  cbind(seq_len(nrow(model_data)), match(model_data$ndvi_col, ndvi_year_cols))
]

# pull the single cc_stack (canopy cover) layer value that corresponds to
# each row's year-prior-to-fire, then drop the now-unneeded per-year columns
cc_year_cols <- grep("^[0-9]{4}_TRE_[0-9]{4}$", names(model_data), value = TRUE)
cc_matrix <- as.matrix(model_data[cc_year_cols])
model_data$canopy_cover_prior_year <- cc_matrix[
  cbind(seq_len(nrow(model_data)), match(model_data$cc_col, cc_year_cols))
]

model_data <- model_data %>% 
  select(-all_of(ndvi_year_cols), -ndvi_col, -burn_year, -ndvi_prior_year,
         -all_of(cc_year_cols), -cc_col, -cc_prior_year)

# TODO: inspect join match rate / missingness before modeling
# model_data %>% summarise(across(everything(), ~sum(is.na(.))))

# ---- 4. random forest: predict Severity_1 from terrain, fuels, weather vars ----

library(randomForest)

# keep the na.omit()'d frame used for fitting so diagnostics (partial
# dependence, etc.) reference the exact same rows/factor levels as rf1.
# as.data.frame() matters here: partialPlot() does pred.data[, xname] to
# pull a single column, which on a tibble returns a one-column tibble
# (breaking is.finite()/min()/max() inside partialPlot) instead of a vector.
rf_data <- model_data %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  left_join(unit.hist %>% 
              filter(burn_entry==1,
                     date_confirm=="y") %>% 
              select(burn_name = Fire_name_fixed,
                     ACRES)) %>%
  left_join(unit_tpi_sd, by = "burn_name")

set.seed(117)
rf1 <- randomForest::randomForest(Severity_1 ~
                                    chili + tpi + tpi_sd + slp.dg + elev +
                                    burn_doy +
                                    ACRES +
                                    daily_max_temp + daily_minRH +
                                    daily_mean_temp + daily_meanRH +
                                    days_since_precip + ts_MAP +
                                    daily_mean_temp_aft + daily_minRH_aft +
                                    days_since_precip_025 +
                                    X1HrFM + X10HrFM + X100HrFM + X1000HrFM +
                                    HerbFM + WoodyFM +
                                    KBDI + ERC +
                                    rh_trend + temp_trend + KBDI_trend + ERC_trend +
                                    ndvi_winter_prior_year +
                                    canopy_cover_prior_year,
                                  data = rf_data,
                                  ntree = 1000,
                                  mtry = 6,
                                  type = "classification")
#                                   

rf1 <- randomForest(Severity_1 ~ .,
                    data = rf_data,
                    ntree = 1000,
                    mtry = 6, 
                    type = "classification") 


rf1
summary(rf1)

randomForest::varImpPlot(rf1)

# ---- 5. model fit diagnostics ----
# rf1 has no held-out test set (fit on all of rf_data), so these diagnostics
# lean on the out-of-bag (OOB) estimates randomForest computes internally
# (each tree's OOB rows are ones it wasn't trained on), which substitute for
# cross-validated performance here.

# OOB confusion matrix + per-class error (last column of rf1$confusion)
confusion_table <- rf1$confusion %>%
  as.data.frame() %>%
  rownames_to_column("Predicted") %>%
  rename(Error_Rate = class.error) %>%
  mutate(Error_Rate = round(Error_Rate, 4))

confusion_table

oob_accuracy <- 1 - rf1$err.rate[nrow(rf1$err.rate), "OOB"]
oob_accuracy

# OOB error rate vs. number of trees, to check ntree = 1000 is sufficient
# for the error estimate to stabilize

err_df <- as.data.frame(rf1$err.rate) %>% 
  mutate(trees = row_number()) %>% 
  pivot_longer(-trees, names_to = "class", values_to = "error")

ggplot(err_df, aes(trees, error, color = class)) +
  geom_line() +
  labs(title = "OOB error rate vs. number of trees",
       x = "Number of trees", y = "OOB error rate", color = NULL)

# variable importance, tidied (mean decrease in Gini impurity)
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

# partial dependence for the top predictors, to check the direction/shape
# of each predictor's marginal effect on the predicted class probability.
# One 2x2 panel of base plots per Severity_1 class.
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

# ---- 5. partial dependencies with ICE lines (pdp) ----

library(pdp)

# Generate PDP plots with ICE lines for each top variable and each class
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

# Create individual plots per variable, faceted by class
library(gridExtra)
pdp_plot_list <- map(top_vars, function(var) {
  data_for_var <- pdp_plots %>% 
    filter(variable == var) %>%
    mutate(
      x_val = .data[[var]],
      class = factor(class, levels = c("Closed", "Open", "Early Seral"))
    )
  
  # Compute PDP by averaging ICE curves across all observations at each grid point
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

