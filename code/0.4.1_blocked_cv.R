#### Burn-unit-blocked k-fold cross-validation for rf1
#
# OOB error from randomForest() resamples rows, not burn units. Since many
# sample points share a burn unit (and therefore identical unit-level
# weather/fuels/geometry predictors, plus spatial autocorrelation), a
# point's OOB trees are often still trained on other points from the same
# unit -- OOB is a reasonable estimate of accuracy at new pixels *within*
# already-sampled units, but optimistic for accuracy at new, unseen burn
# units. This script re-evaluates rf1 with cross-validation blocked on 
# burn_name, so no unit straddles a train/test split.
#
# rf_data now spans burn entries 1-4. Entries of the same unit share (nearly)
# the same footprint, so blocking is on burn_name only: ALL entries of a unit
# go to the same fold. Blocking on burn_name x entry would leak the unit's
# terrain/spatial signal across train and test.
#
# Depends on objects from:
#   0.4_exploratory_model.R
#
# Outputs:
#   unit_folds        -- burn_name -> fold assignment (stratified by each
#                         unit's dominant Severity_1 class, pooled over entries)
#   blocked_cv         -- list(unit_folds, predictions, fold_metrics,
#                         confusion, models)
#     predictions       one row per rf_data row, with its held-out fold and
#                        predicted class
#     fold_metrics       per-fold accuracy, kappa, n units/rows, sparse_fold
#                        flag (fold's test units don't cover all 3 classes)
#     confusion           pooled confusion matrix across all folds (each row
#                        tested exactly once, in its one held-out fold)
#     models              the fold-f randomForest fit, for predicting rasters
#                        of fold f's held-out (test) units
#
#   output/prelim_v2/blocked_cv_predictions/<burn>_entry<n>_blockedCV_pred.tif
#                       -- genuine out-of-unit prediction raster for each
#                          burn entry of a sample of held-out test units
#                          (section 9), i.e. each unit's surface comes from the
#                          fold model that never trained on that unit. Compare
#                          against 0.5_prelim_prediction.R's <burn>_pred.tif
#                          (rf1, trained on all units) to see how much the
#                          surface degrades once unit-level leakage is removed.

library(tidyverse)
library(randomForest)
library(sf)
library(terra)

if (!exists("rf_data")) {
  stop("rf_data not found -- run code/0.4_exploratory_model.R first.")
}

# Same predictor set as rf1 in 0.4_exploratory_model.R
cv_formula <- reformulate(predictors, response = "Severity_1")

# ---- helper: Cohen's kappa from a confusion matrix (observed x predicted) ----

cohens_kappa <- function(cm) {
  n <- sum(cm)
  po <- sum(diag(cm)) / n
  pe <- sum((rowSums(cm) / n) * (colSums(cm) / n))
  (po - pe) / (1 - pe)
}

# ---- helper: burn-unit-blocked k-fold CV ----
# Whole burn units (not rows) are assigned to folds, so every point from a
# given unit -- across all of its burn entries -- is always on the same side
# of the train/test split. Units are stratified by their *dominant*
# Severity_1 class before assignment so folds aren't skewed toward one class;
# an unstratified random split could easily leave a fold with no test units of
# a given class.

run_blocked_cv <- function(dat, formula, group_var = "burn_name",
                            response_var = "Severity_1", k = 5,
                            ntree = 1000, mtry = 6, seed) {
  set.seed(seed)

  unit_dominant_class <- dat %>%
    count(.data[[group_var]], .data[[response_var]]) %>%
    group_by(across(all_of(group_var))) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup()

  unit_folds <- unit_dominant_class %>%
    group_by(across(all_of(response_var))) %>%
    mutate(fold = sample(rep_len(seq_len(k), n()))) %>%
    ungroup() %>%
    select(all_of(group_var), fold)

  dat <- dat %>% left_join(unit_folds, by = group_var)
  n_classes <- n_distinct(dat[[response_var]])

  fold_results <- purrr::map(seq_len(k), function(f) {
    train <- dat %>% filter(fold != f) %>% select(-fold)
    test  <- dat %>% filter(fold == f) %>% select(-fold)

    fit  <- randomForest(formula, data = train, ntree = ntree, mtry = mtry)
    pred <- predict(fit, newdata = test)
    obs  <- test[[response_var]]

    list(
      predictions = tibble(
        !!group_var := test[[group_var]],
        fold = f,
        observed = obs,
        predicted = pred
      ),
      metrics = tibble(
        fold = f,
        n_train_units = n_distinct(train[[group_var]]),
        n_test_units  = n_distinct(test[[group_var]]),
        n_train       = nrow(train),
        n_test        = nrow(test),
        sparse_fold   = n_distinct(obs) < n_classes,
        accuracy      = mean(pred == obs),
        kappa         = cohens_kappa(table(obs, pred))
      ),
      fit = fit
    )
  })

  predictions  <- purrr::map_dfr(fold_results, "predictions")
  fold_metrics <- purrr::map_dfr(fold_results, "metrics")
  confusion    <- table(observed = predictions$observed,
                         predicted = predictions$predicted)
  # Named by fold number so blocked_cv$models[[f]] gives fold f's fit -- the
  # model trained on every unit *except* fold f's test units.
  models <- purrr::map(fold_results, "fit") %>% setNames(seq_len(k))

  list(unit_folds = unit_folds, predictions = predictions,
       fold_metrics = fold_metrics, confusion = confusion, models = models)
}

# ---- run ----

blocked_cv <- run_blocked_cv(rf_data, cv_formula, k = 5, seed = 6203)

cat("Burn-unit-blocked 5-fold CV:\n")
print(blocked_cv$fold_metrics)
cat("  mean accuracy =", round(mean(blocked_cv$fold_metrics$accuracy), 3),
    " mean kappa =", round(mean(blocked_cv$fold_metrics$kappa), 3),
    " (", sum(blocked_cv$fold_metrics$sparse_fold), "/", nrow(blocked_cv$fold_metrics),
    "folds flagged sparse -- test units don't cover all", n_distinct(rf_data$Severity_1),
    "classes)\n")

cat("\nPooled confusion matrix (each row tested once, in its held-out fold):\n")
print(blocked_cv$confusion)

# per-class error rate, pooled across folds -- compare against rf1$confusion
# (OOB, in 0.4_exploratory_model.R) to see how much OOB understates error
# for predicting genuinely new burn units
pooled_class_error <- blocked_cv$confusion %>%
  as.data.frame() %>%
  group_by(observed) %>%
  summarise(
    n = sum(Freq),
    correct = sum(Freq[predicted == observed]),
    error_rate = round(1 - correct / n, 4),
    .groups = "drop"
  )

pooled_class_error

# Note: unit_dominant_class/unit_folds assign a whole unit by its majority
# class; units with no majority (e.g. an even split across classes) get
# whichever class `slice_max()` returns first -- fine for balancing fold
# composition, but not a claim about the unit's "true" class.

# ---- 9. Raster predictions for a sample of held-out test units ---------------
# Mirrors the prediction pipeline in 0.5_prelim_prediction.R, but for each
# rendered unit uses blocked_cv$models[[f]] -- the fold model that did NOT
# see that unit during training -- instead of rf1 (trained on all units).
# These are genuine out-of-unit prediction surfaces; compare visually against
# 0.5's <burn>_pred.tif for the same unit to see what blocked CV's accuracy
# drop actually looks like on the ground.
#
# One raster is made per burn entry of each rendered unit, using that entry's
# own weather, burn_doy, burn_entry, and burn year (prior-year NDVI / canopy
# cover / insect-disease window).
#
# Predictor set matches cv_formula, which includes insect_disease_5yr --
# resolved per pixel from insdis_stack, same 5-year pre-burn window logic as
# 0.5_prelim_prediction.R.
#
# Requires the same raster prerequisites as 0.5_prelim_prediction.R:
#   terrain_stack, eric_stack, cc_stack, ecozones, insdis_stack

stopifnot(
  "terrain_stack not found -- run 0.4_exploratory_model.R first" = exists("terrain_stack"),
  "eric_stack not found -- run 0.4_exploratory_model.R first"    = exists("eric_stack"),
  "cc_stack not found -- run 0.4_exploratory_model.R first"      = exists("cc_stack"),
  "ecozones not found -- run 0.4_exploratory_model.R first"      = exists("ecozones"),
  "insdis_stack not found -- run 0.4_exploratory_model.R first"  = exists("insdis_stack")
)

cv_model_preds  <- all.vars(cv_formula)[-1]
raster_preds_cv <- c("chili", "tpi", "slp.dg", "elev",
                      "ndvi_winter_prior_year", "canopy_cover_prior_year",
                      "ecozone", "insect_disease_5yr")
# includes burn_entry, which comes from the per-entry unit constants
unit_preds_cv   <- setdiff(cv_model_preds, raster_preds_cv)

insdis_layer_names_cv <- grep("^insdis_[0-9]{4}$", names(insdis_stack), value = TRUE)
insdis_layer_years_cv <- as.integer(gsub("^insdis_", "", insdis_layer_names_cv))

# Recode a pixel-level categorical column to the factor levels a given fold's
# model actually saw at fit time, warning about (and NA-ing out) any code
# that fold's training data never saw.
recode_to_training_levels <- function(x, levels, label, burn) {
  recoded <- factor(as.character(x), levels = levels)
  n_new <- sum(is.na(recoded) & !is.na(x))
  if (n_new > 0) {
    warning(n_new, " pixel(s) in ", burn, " have a ", label,
            " code not seen in training -- set to NA")
  }
  recoded
}

# same polygons as 0.4 (combined entries, strict year match)
if (!exists("rx.shp")) {
  rx.shp <- sf::read_sf("data/burn_polygons/rx_canopy_all_entries.gpkg") %>%
    filter(year_check == "ok")
}

# unit-level constants, one row per burn unit x entry (weather and burn_doy
# differ by entry)
unit_constants_cv <- rf_data %>%
  group_by(burn_name, burn_entry) %>%
  slice(1) %>%
  ungroup() %>%
  select(burn_name, burn_entry, all_of(unit_preds_cv))

# one dissolved polygon per unit x entry
burn_bounds_cv <- rx.shp %>%
  filter(paste(burn_name, burn_entry) %in%
           paste(unit_constants_cv$burn_name, unit_constants_cv$burn_entry)) %>%
  group_by(burn_name, burn_entry) %>%
  # sf's summarise unions geometries itself (gpkg geometry column is "geom")
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_transform(crs(terrain_stack)) %>%
  left_join(
    rf_data %>%
      group_by(burn_name, burn_entry) %>%
      summarise(burn_year = as.integer(format(as.Date(date[1]), "%Y")), .groups = "drop"),
    by = c("burn_name", "burn_entry")
  )

predict_blocked_raster <- function(burn, entry, fit, model_preds, out_dir) {

  unit_poly  <- burn_bounds_cv %>% filter(burn_name == burn, burn_entry == entry)
  unit_const <- unit_constants_cv %>% filter(burn_name == burn, burn_entry == entry)
  burn_year  <- unit_poly$burn_year

  ndvi_col <- paste0("y_", burn_year - 1)
  cc_col   <- paste0(burn_year - 1, "_TRE_", burn_year - 1)

  topo_crop <- crop(terrain_stack, vect(unit_poly)) %>% mask(vect(unit_poly))

  if (ndvi_col %in% names(eric_stack)) {
    ndvi_crop <- crop(eric_stack[[ndvi_col]], vect(unit_poly)) %>%
      mask(vect(unit_poly)) %>%
      resample(topo_crop, method = "near")
  } else {
    ndvi_crop <- topo_crop[[1]]; values(ndvi_crop) <- NA_real_
    warning("NDVI layer ", ndvi_col, " not found for burn: ", burn, " entry ", entry)
  }
  names(ndvi_crop) <- "ndvi_winter_prior_year"

  if (cc_col %in% names(cc_stack)) {
    cc_crop <- crop(cc_stack[[cc_col]], vect(unit_poly)) %>%
      mask(vect(unit_poly)) %>%
      resample(topo_crop, method = "bilinear")
  } else {
    cc_crop <- topo_crop[[1]]; values(cc_crop) <- NA_real_
    warning("CC layer ", cc_col, " not found for burn: ", burn, " entry ", entry)
  }
  names(cc_crop) <- "canopy_cover_prior_year"

  ecozone_crop <- crop(ecozones, vect(unit_poly)) %>%
    mask(vect(unit_poly)) %>%
    resample(topo_crop, method = "near")
  names(ecozone_crop) <- "ecozone"

  # Insect/disease disturbance, 5 yrs prior to burn. Binary 0/1 layers, so
  # "any year in the window shows disturbance" reduces to a cellwise max
  # across whichever window years are available in insdis_stack.
  insdis_window_years <- (burn_year - 5):(burn_year - 1)
  insdis_cols_in_window <- which(insdis_layer_years_cv %in% insdis_window_years)

  if (length(insdis_cols_in_window) > 0) {
    insdis_window_stack <- crop(insdis_stack[[insdis_layer_names_cv[insdis_cols_in_window]]],
                                 vect(unit_poly)) %>%
      mask(vect(unit_poly))
    insdis_crop <- if (nlyr(insdis_window_stack) > 1) {
      app(insdis_window_stack, fun = function(x) {
        if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
      })
    } else {
      insdis_window_stack
    }
    insdis_crop <- resample(insdis_crop, topo_crop, method = "near")
  } else {
    insdis_crop <- topo_crop[[1]]; values(insdis_crop) <- NA_real_
    warning("No insect/disease layers in [", burn_year - 5, ", ", burn_year - 1,
            "] for burn: ", burn, " entry ", entry)
  }
  names(insdis_crop) <- "insect_disease_5yr"

  topo_df    <- as.data.frame(topo_crop,    xy = TRUE,  na.rm = FALSE)
  ndvi_df    <- as.data.frame(ndvi_crop,    xy = FALSE, na.rm = FALSE)
  cc_df      <- as.data.frame(cc_crop,      xy = FALSE, na.rm = FALSE)
  ecozone_df <- as.data.frame(ecozone_crop, xy = FALSE, na.rm = FALSE)
  insdis_df  <- as.data.frame(insdis_crop,  xy = FALSE, na.rm = FALSE)

  pix_df     <- bind_cols(topo_df, ndvi_df, cc_df, ecozone_df, insdis_df)
  valid_mask <- !is.na(pix_df$chili)
  n_total    <- nrow(pix_df)
  n_valid    <- sum(valid_mask)

  if (n_valid == 0) {
    message("  No valid pixels for: ", burn, " entry ", entry, " -- skipping")
    return(invisible(NULL))
  }

  # burn_entry is a model predictor and stays in; only burn_name is dropped
  pix_valid <- pix_df[valid_mask, ] %>%
    bind_cols(unit_const[rep(1, n_valid), setdiff(names(unit_const), "burn_name")])

  ecozone_levels <- fit$forest$xlevels$ecozone
  ndvi_levels    <- fit$forest$xlevels$ndvi_winter_prior_year

  pix_valid$ecozone <- recode_to_training_levels(
    pix_valid$ecozone, ecozone_levels, "ecozone", burn)
  pix_valid$ndvi_winter_prior_year <- recode_to_training_levels(
    pix_valid$ndvi_winter_prior_year, ndvi_levels, "ndvi_winter_prior_year", burn)

  pred_input <- pix_valid[, model_preds] %>% as.data.frame()

  pred_class <- predict(fit, newdata = pred_input, type = "response")
  pred_probs <- predict(fit, newdata = pred_input, type = "prob") %>% as.data.frame()
  class_num  <- as.integer(pred_class)

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

  out_path <- file.path(
    out_dir,
    paste0(gsub("[^A-Za-z0-9_-]", "_", burn), "_entry", entry, "_blockedCV_pred.tif")
  )
  writeRaster(out_ras, out_path, overwrite = TRUE)
  message("Saved blocked-CV prediction raster: ", burn, " entry ", entry)

  out_ras
}

# Which held-out units to render. Default: the first 2 test units per fold
# (10 units total), each rendered for every burn entry present in rf_data --
# enough to inspect visually without rebuilding rasters for every unit. Set
# to a larger number, or Inf, to render all test units.
n_units_per_fold <- 2

test_units_to_render <- blocked_cv$unit_folds %>%
  group_by(fold) %>%
  slice_head(n = n_units_per_fold) %>%
  ungroup()

units_to_render <- rf_data %>%
  distinct(burn_name, burn_entry) %>%
  inner_join(test_units_to_render, by = "burn_name") %>%
  arrange(fold, burn_name, burn_entry)

out_dir_cv <- "output/prelim_v3/blocked_cv_predictions"
dir.create(out_dir_cv, showWarnings = FALSE, recursive = TRUE)

blocked_rasters <- purrr::pmap(
  list(units_to_render$burn_name, units_to_render$burn_entry, units_to_render$fold),
  function(burn, entry, f) {
    predict_blocked_raster(
      burn        = burn,
      entry       = entry,
      fit         = blocked_cv$models[[as.character(f)]],
      model_preds = cv_model_preds,
      out_dir     = out_dir_cv
    )
  }
) %>% setNames(paste0(units_to_render$burn_name, "_entry", units_to_render$burn_entry))
