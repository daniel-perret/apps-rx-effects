#### Spatial predictions for each burn unit
####
#### Applies rf1 (from 0.4_exploratory_model.R) to generate pixel-level
#### canopy severity predictions across all confirmed burn units.
####
#### Predictions are restricted to burns represented in rf_data (the exact
#### training set), so unit-level constants — fire weather, burn_doy, ACRES —
#### are pulled directly from rf_data rather than re-joined from source files.
####
#### Pixel-varying predictors  : terrain (chili, tpi, slp.dg, elev),
####                             winter NDVI (year-before-burn),
####                             canopy cover (year-before-burn)
#### Unit-level constants       : all other model predictors (from rf_data)
####
#### Outputs
####   output/prelim_model/unit_predictions/<burn>_pred.tif   — 4-layer raster per unit
####     (prob_Closed, prob_Open, prob_EarlySeral, pred_class)
####   output/prelim_model/unit_predictions/unit_summary.csv  — polygon-level summary

library(tidyverse)
library(sf)
library(terra)
library(randomForest)

# ---- 0. Prerequisites --------------------------------------------------------
# rf1, rf_data, terrain_stack, eric_stack, cc_stack must be present.
# Run 0.4_exploratory_model.R if needed, or load a saved model:
#   rf1    <- readRDS("outputs/rf1.rds")
#   rf_data <- readRDS("outputs/rf_data.rds")

stopifnot(
  "rf1 not found — run 0.4_exploratory_model.R first"    = exists("rf1"),
  "rf_data not found — run 0.4_exploratory_model.R first" = exists("rf_data"),
  "terrain_stack not found" = exists("terrain_stack"),
  "eric_stack not found"    = exists("eric_stack"),
  "cc_stack not found"      = exists("cc_stack")
)

# Predictors the model was actually trained on
model_preds <- rownames(rf1$importance)

# Raster-level (pixel-varying) predictors; everything else is unit-level
raster_preds <- c("chili", "tpi", "slp.dg", "elev",
                  "ndvi_winter_prior_year", "canopy_cover_prior_year")
unit_preds   <- setdiff(model_preds, raster_preds)

# ---- 1. Burns in training data only ------------------------------------------

# One row per burn — take the first sampled point (all unit-level vars are
# constant within a burn in rf_data)
unit_constants <- rf_data %>%
  group_by(burn_name) %>%
  slice(1) %>%
  ungroup() %>%
  select(burn_name, all_of(unit_preds))

confirmed_burns <- unit_constants$burn_name
cat("Burns in training data:", length(confirmed_burns), "\n")

# ---- 2. Burn unit boundaries (dissolve multi-polygon burns) ------------------

rx.shp <- sf::read_sf(
  "data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel/1burn_cnpy.shp"
)

burn_bounds <- rx.shp %>%
  filter(burn_name %in% confirmed_burns) %>%
  group_by(burn_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_transform(crs(terrain_stack))

# Burn year for selecting year-specific raster layers — derive from rf_data
burn_years <- rf_data %>%
  group_by(burn_name) %>%
  summarise(burn_year = as.integer(format(as.Date(date[1]), "%Y")), .groups = "drop")

burn_bounds <- burn_bounds %>%
  left_join(burn_years, by = "burn_name")

# ---- 3. Layer name lookups for year-specific rasters -------------------------

ndvi_layer_names <- grep("^y_[0-9]{4}$",             names(eric_stack), value = TRUE)
cc_layer_names   <- grep("^[0-9]{4}_TRE_[0-9]{4}$", names(cc_stack),   value = TRUE)

# ---- 4. Predict for each burn unit ------------------------------------------

pred_list    <- list()
unit_summary <- list()

dir.create("output/prelim_model/unit_predictions", showWarnings = FALSE, recursive = TRUE)

for (burn in confirmed_burns) {

  unit_poly  <- burn_bounds %>% filter(burn_name == burn)
  unit_const <- unit_constants %>% filter(burn_name == burn)
  burn_year  <- unit_poly$burn_year

  ndvi_col <- paste0("y_", burn_year - 1)
  cc_col   <- paste0(burn_year - 1, "_TRE_", burn_year - 1)

  # -- 4a. Crop rasters to burn boundary ---------------------------------------

  topo_crop <- crop(terrain_stack, vect(unit_poly)) %>% mask(vect(unit_poly))

  if (ndvi_col %in% ndvi_layer_names) {
    ndvi_crop <- crop(eric_stack[[ndvi_col]], vect(unit_poly)) %>%
      mask(vect(unit_poly)) %>%
      resample(topo_crop, method = "bilinear")
  } else {
    ndvi_crop <- topo_crop[[1]]; values(ndvi_crop) <- NA_real_
    warning("NDVI layer ", ndvi_col, " not found for burn: ", burn)
  }
  names(ndvi_crop) <- "ndvi_winter_prior_year"

  if (cc_col %in% cc_layer_names) {
    cc_crop <- crop(cc_stack[[cc_col]], vect(unit_poly)) %>%
      mask(vect(unit_poly)) %>%
      resample(topo_crop, method = "bilinear")
  } else {
    cc_crop <- topo_crop[[1]]; values(cc_crop) <- NA_real_
    warning("CC layer ", cc_col, " not found for burn: ", burn)
  }
  names(cc_crop) <- "canopy_cover_prior_year"

  # -- 4b. Build pixel-level data frame ----------------------------------------

  topo_df <- as.data.frame(topo_crop, xy = TRUE, na.rm = FALSE)
  ndvi_df <- as.data.frame(ndvi_crop, xy = FALSE, na.rm = FALSE)
  cc_df   <- as.data.frame(cc_crop,   xy = FALSE, na.rm = FALSE)

  pix_df     <- bind_cols(topo_df, ndvi_df, cc_df)
  valid_mask <- !is.na(pix_df$chili)
  n_total    <- nrow(pix_df)
  n_valid    <- sum(valid_mask)

  if (n_valid == 0) {
    message("  No valid pixels for: ", burn, " — skipping")
    next
  }

  # Add unit-level constants (same value for every pixel in this burn)
  pix_valid <- pix_df[valid_mask, ] %>%
    bind_cols(unit_const[rep(1, n_valid), unit_preds])

  # -- 4c. Predict -------------------------------------------------------------

  pred_input <- pix_valid[, model_preds] %>% as.data.frame()

  pred_class <- predict(rf1, newdata = pred_input, type = "response")
  pred_probs <- predict(rf1, newdata = pred_input, type = "prob") %>%
    as.data.frame()

  class_num  <- as.integer(pred_class)

  # -- 4d. Map predictions back into rasters -----------------------------------

  prob_closed_vec <- rep(NA_real_, n_total)
  prob_open_vec   <- rep(NA_real_, n_total)
  prob_es_vec     <- rep(NA_real_, n_total)
  pred_class_vec  <- rep(NA_real_, n_total)

  prob_closed_vec[valid_mask] <- pred_probs$Closed
  prob_open_vec[valid_mask]   <- pred_probs$Open
  prob_es_vec[valid_mask]     <- pred_probs[["Early Seral"]]
  pred_class_vec[valid_mask]  <- class_num

  out_ras <- c(topo_crop[[1]], topo_crop[[1]], topo_crop[[1]], topo_crop[[1]])
  names(out_ras) <- c("prob_Closed", "prob_Open", "prob_EarlySeral", "pred_class")
  values(out_ras[[1]]) <- prob_closed_vec
  values(out_ras[[2]]) <- prob_open_vec
  values(out_ras[[3]]) <- prob_es_vec
  values(out_ras[[4]]) <- pred_class_vec

  levels(out_ras[[4]]) <- data.frame(id = 1:3,
                                     label = c("Closed", "Open", "Early Seral"))

  # -- 4e. Save raster ---------------------------------------------------------

  out_path <- file.path(
    "output/prelim_model/unit_predictions",
    paste0(gsub("[^A-Za-z0-9_-]", "_", burn), "_pred.tif")
  )
  writeRaster(out_ras, out_path, overwrite = TRUE)
  pred_list[[burn]] <- out_ras

  # -- 4f. Unit-level summary --------------------------------------------------

  dominant <- names(which.max(colMeans(pred_probs)))

  unit_summary[[burn]] <- tibble(
    burn_name            = burn,
    n_pixels             = n_valid,
    mean_prob_Closed     = mean(pred_probs$Closed),
    mean_prob_Open       = mean(pred_probs$Open),
    mean_prob_EarlySeral = mean(pred_probs[["Early Seral"]]),
    dominant_class       = dominant
  )

  message("Predicted: ", burn,
          "  [", n_valid, " px | dominant: ", dominant, "]")
}

# ---- 5. Compile and save unit summary ----------------------------------------

unit_summary_df <- bind_rows(unit_summary)

write.csv(unit_summary_df,
          "output/prelim_model/unit_predictions/unit_summary.csv",
          row.names = FALSE)

print(unit_summary_df %>%
        select(burn_name, n_pixels, dominant_class,
               mean_prob_Closed, mean_prob_Open, mean_prob_EarlySeral))

# ---- 6. Map: dominant predicted class per burn unit --------------------------

rx.pred <- burn_bounds %>%
  left_join(unit_summary_df, by = "burn_name") %>%
  mutate(dominant_class = factor(dominant_class,
                                 levels = c("Closed", "Open", "Early Seral")))

ggplot(rx.pred) +
  geom_sf(aes(fill = dominant_class), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    name   = "Dominant predicted class",
    values = c("Closed"      = "forestgreen",
               "Open"        = "goldenrod3",
               "Early Seral" = "firebrick2"),
    na.value = "grey80"
  ) +
  labs(title = "Predicted canopy severity — dominant class per burn unit",
       subtitle = "Conditioned on actual burn-day weather; restricted to model training burns")

# ---- 7. Mosaics: one per predicted layer ------------------------------------
# Each burn unit raster is cropped to its own extent; mosaic tiles them into
# study-area-wide surfaces. Units burned in different years carry different
# weather conditions, so each mosaic reflects actual per-unit conditions.
# For probability layers, overlapping pixels (if any) are averaged.
# For pred_class, the modal value is used.

layer_specs <- list(
  list(idx = 1, name = "prob_Closed",     fun = "mean"),
  list(idx = 2, name = "prob_Open",       fun = "mean"),
  list(idx = 3, name = "prob_EarlySeral", fun = "mean"),
  list(idx = 4, name = "pred_class",      fun = "first")
)

mosaic_list <- lapply(layer_specs, function(spec) {
  layer_rasters <- sprc(lapply(pred_list, \(r) r[[spec$idx]]))
  mos <- mosaic(layer_rasters, fun = spec$fun)
  names(mos) <- spec$name
  out_path <- file.path("output/prelim_model/unit_predictions",
                        paste0("mosaic_", spec$name, ".tif"))
  writeRaster(mos, out_path, overwrite = TRUE)
  message("Saved mosaic: ", spec$name)
  mos
})
names(mosaic_list) <- sapply(layer_specs, `[[`, "name")

# ----- 8. Simple comparisons at unit-level -----------------------------

comp.df <- rx.shp %>%
  sf::st_drop_geometry() %>% 
  filter(burn_name %in% confirmed_burns) %>%
  mutate(Severity_class = case_when(Severity_1 == 0 ~ "Closed",
                                Severity_1 == 1 ~ "Open",
                                Severity_1 == 2 ~ "Early",
                                TRUE ~ "Other")) %>%
  filter(!Severity_class == "Other") %>% 
  group_by(burn_name, Severity_class) %>% 
  summarise(area = sum(gisAcres)) %>% 
  group_by(burn_name) %>% 
  mutate(observed.prop = area/sum(area)) %>%
  ungroup() %>%
  select(-area) %>%
  pivot_wider(names_from = Severity_class, 
              values_from = observed.prop, 
              names_prefix = "obs_") %>% 
  left_join(unit_summary_df,
            by = "burn_name")
  
ggplot(comp.df,
       aes(x = obs_Early,
           y = mean_prob_EarlySeral)) +
  geom_point(aes(col = "Early"),
             pch = 19,
             size = 3) +
  geom_point(aes(x = obs_Closed,
                 y = mean_prob_Closed,
                 col = "Closed"),
             pch = 19,
             size = 3) +
  geom_point(aes(x = obs_Open,
                 y = mean_prob_Open,
                 col = "Open"),
             size = 3, pch = 19) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed", y = "Predicted")
































