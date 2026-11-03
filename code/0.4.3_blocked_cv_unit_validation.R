#### Burn-unit-blocked CV with whole-unit validation (pixel-level + class
#### proportions per burn event)
#
# Companion to 0.4.2_exploratory_model_unit_sample.R; replaces the validation in
# 0.4.1_blocked_cv.R. Differences from 0.4.1:
#
#   1. Folds are blocked on burn_name (all entries of a unit in the same fold,
#      as in 0.4.1), but units are stratified by their dominant class computed
#      from polygon AREA, not from the noisy sample counts.
#   2. Held-out units are not scored on the sampled training-style points.
#      Instead the fold model predicts EVERY pixel of each held-out burn event
#      (burn_name x burn_entry), and the predictions are compared with the
#      observed severity polygons rasterized to the same 30 m grid.
#   3. Per event we report pixel accuracy / kappa AND the predicted vs observed
#      share of the event in each severity class (hard-classified pixel shares
#      and mean predicted probabilities), which is what operational planning
#      needs. Each event is also compared with a naive baseline that assigns
#      every event the training-fold average class mix.
#
# Pixel predictors are taken at pixel centres from the native layers
# (nearest cell, as in training-point extraction), and any pixel with a missing
# predictor is excluded from scoring and counted (n_excluded), so predictions
# can never misalign with pixels.
#
# Depends on objects from: 0.4.2_exploratory_model_unit_sample.R
#   rf_data_unit, predictors_unit, rx_pool, unit_tpi_sd, terrain_stack,
#   eric_stack, cc_stack, ecozones, insdis_stack
#
# Outputs (object `unit_validation`, plus CSVs in output/unit_validation/):
#   unit_folds          burn_name -> fold, dominant class (by area)
#   event_metrics       one row per held-out event: pixel accuracy, kappa,
#                       total-variation distance (TVD) of class shares for the
#                       model vs the baseline, dominant-class match
#   event_props         one row per event x class: observed share, predicted
#                       share (hard), predicted share (soft), baseline share
#   class_proportion_summary   bias / MAE / RMSE / r of class shares across
#                       events, share of events within `prop_tol`
#   pooled_confusion, class_accuracy, pooled_summary, event_summary
#   fold_models

library(tidyverse)
library(randomForest)
library(sf)
library(terra)

needed <- c("rf_data_unit", "predictors_unit", "rx_pool", "unit_tpi_sd",
            "terrain_stack", "eric_stack", "cc_stack", "ecozones", "insdis_stack")
missing_objs <- needed[!vapply(needed, exists, logical(1))]
if (length(missing_objs) > 0) {
  stop("Missing: ", paste(missing_objs, collapse = ", "),
       " -- run code/0.4.2_exploratory_model_unit_sample.R first.")
}

# ---- settings (EDIT HERE) ----

k        <- 5
seed_cv  <- 5481
ntree    <- 1000
mtry     <- 6
# a predicted class share counts as "close" if within this many proportion
# points of the observed share (0.10 = 10 percentage points)
prop_tol <- 0.10

# Held-out prediction rasters for visual comparison: the first
# `n_render_per_fold` events of each fold (sorted by name) are written as
# GeoTIFFs (use Inf for all events). Each event is predicted by the fold model
# that never saw its unit. Set to 0 to skip writing rasters.
n_render_per_fold <- 10
raster_dir <- "output/unit_validation/rasters"
if (n_render_per_fold > 0) dir.create(raster_dir, showWarnings = FALSE, recursive = TRUE)

class_levels <- levels(rf_data_unit$Severity_1) # Closed, Open, Early Seral = 0, 1, 2
cv_formula   <- reformulate(predictors_unit, response = "Severity_1")

cohens_kappa <- function(cm) {
  n <- sum(cm)
  po <- sum(diag(cm)) / n
  pe <- sum((rowSums(cm) / n) * (colSums(cm) / n))
  if (isTRUE(all.equal(pe, 1))) return(NA_real_)
  (po - pe) / (1 - pe)
}

# ---- 1. folds: blocked on burn_name, stratified by area-dominant class ----

pool_area <- tibble(
  burn_name  = rx_pool$burn_name,
  Severity_1 = rx_pool$Severity_1,
  area_m2    = as.numeric(st_area(rx_pool))
)

set.seed(seed_cv)

unit_folds <- pool_area %>%
  filter(burn_name %in% rf_data_unit$burn_name) %>%
  group_by(burn_name, Severity_1) %>%
  summarise(area_m2 = sum(area_m2), .groups = "drop") %>%
  group_by(burn_name) %>%
  slice_max(area_m2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(dominant_class = factor(Severity_1, levels = 0:2, labels = class_levels)) %>%
  group_by(dominant_class) %>%
  mutate(fold = sample(rep_len(seq_len(k), n()))) %>%
  ungroup() %>%
  select(burn_name, dominant_class, fold)

train_data <- rf_data_unit %>%
  left_join(unit_folds %>% select(burn_name, fold), by = "burn_name")

# units per fold and dominant class
unit_folds %>% count(fold, dominant_class) %>% pivot_wider(names_from = dominant_class,
                                                            values_from = n, values_fill = 0)

# ---- 2. fold models (trained on the equal-n-per-event sample) ----

fold_models <- map(seq_len(k), function(f) {
  set.seed(seed_cv + f)
  randomForest(cv_formula,
               data = train_data %>% filter(fold != f) %>% select(-fold),
               ntree = ntree, 
               mtry = mtry,
               sampsize = samp.size, # pulled from 0.4.2
               type = "classification")
}) %>% setNames(seq_len(k))

# baseline class mix per fold = class shares of that fold's training sample
# (equal points per event, so this is an event-averaged mix)
baseline_props <- map(seq_len(k), function(f) {
  tab <- table(train_data$Severity_1[train_data$fold != f])
  as.numeric(prop.table(tab))[match(class_levels, names(tab))]
}) %>% setNames(seq_len(k))

# ---- 3. pixel-level prediction + scoring for each held-out event ----

events <- rf_data_unit %>%
  group_by(burn_name, burn_entry) %>%
  summarise(burn_year = as.integer(format(as.Date(dplyr::first(date)), "%Y")),
            .groups = "drop") %>%
  left_join(unit_folds %>% select(burn_name, fold), by = "burn_name")

# Predictors split by how they are obtained for each pixel:
#   pixel_preds -- taken from rasters at each pixel centre
#   unit_preds  -- constant within a burn event (tpi_sd, burn_entry, burn_doy,
#                  geometry, fire weather, fuels, ...); read from the event's
#                  rows in rf_data_unit, so any predictor added to
#                  predictors_unit in 0.4.2 is handled automatically
pixel_preds <- c("chili", "tpi", "slp.dg", "elev", "ecozone", "insect_disease_5yr",
                 "ndvi_winter_prior_year", "canopy_cover_prior_year")
unit_preds  <- setdiff(predictors_unit, pixel_preds)

unit_constants <- rf_data_unit %>%
  group_by(burn_name, burn_entry) %>%
  slice(1) %>%
  ungroup() %>%
  select(burn_name, burn_entry, all_of(setdiff(unit_preds, "burn_entry")))

insdis_layer_names <- grep("^insdis_[0-9]{4}$", names(insdis_stack), value = TRUE)
insdis_layer_years <- as.integer(gsub("^insdis_", "", insdis_layer_names))

# extract raster values at pixel-centre coordinates, dropping terra's ID column
extract_xy <- function(r, xy) {
  out <- terra::extract(r, xy)
  out[, setdiff(names(out), "ID"), drop = FALSE]
}

extract_layer_or_na <- function(stack, layer, xy) {
  if (layer %in% names(stack)) {
    extract_xy(stack[[layer]], xy)[[1]]
  } else {
    rep(NA_real_, nrow(xy))
  }
}

validate_event <- function(burn, entry, burn_year, fold, write_raster = FALSE) {
  unit_const <- unit_constants %>%
    filter(burn_name == burn, burn_entry == entry) %>%
    select(-burn_name)

  fit <- fold_models[[as.character(fold)]]

  ev   <- rx_pool %>% filter(burn_name == burn, burn_entry == entry)
  ev_v <- terra::vect(ev)

  topo  <- crop(terrain_stack, ev_v) %>% mask(ev_v)
  obs_r <- terra::rasterize(ev_v, topo[[1]], field = "Severity_1")

  pix <- as.data.frame(topo, xy = TRUE, cells = TRUE, na.rm = FALSE)
  pix$obs <- factor(terra::values(obs_r)[, 1], levels = 0:2, labels = class_levels)
  # scored population: pixels with an observed severity class and terrain data
  pix <- pix[!is.na(pix$obs) & !is.na(pix$chili), ]
  if (nrow(pix) == 0) return(NULL)

  xy <- as.matrix(pix[, c("x", "y")])

  # insect/disease: same rule as training (any disturbance in the 5 prior
  # years; all-NA counts as 0; no layers in window -> NA)
  window_cols <- insdis_layer_names[insdis_layer_years %in% (burn_year - 5):(burn_year - 1)]
  insect <- if (length(window_cols) == 0) {
    rep(NA_integer_, nrow(xy))
  } else {
    as.integer(rowSums(as.matrix(extract_xy(insdis_stack[[window_cols]], xy)) == 1,
                       na.rm = TRUE) > 0)
  }

  # categorical predictors recoded to the levels the fold model saw
  ecozone <- factor(as.character(extract_xy(ecozones[[1]], xy)[[1]]),
                    levels = fit$forest$xlevels$ecozone)
  ndvi <- factor(as.character(extract_layer_or_na(eric_stack, paste0("y_", burn_year - 1), xy)),
                 levels = fit$forest$xlevels$ndvi_winter_prior_year)

  cc_layer <- paste0(burn_year - 1, "_TRE_", burn_year - 1)

  pred_df <- data.frame(
    chili = pix$chili, tpi = pix$tpi, slp.dg = pix$slp.dg, elev = pix$elev,
    # event-level constants are appended below
    ecozone = ecozone,
    insect_disease_5yr = insect,
    ndvi_winter_prior_year = ndvi,
    canopy_cover_prior_year = extract_layer_or_na(cc_stack, cc_layer, xy)
  ) %>%
    bind_cols(unit_const[rep(1, nrow(pix)), ]) %>%
    .[, predictors_unit]

  ok <- complete.cases(pred_df)
  if (!any(ok)) return(NULL)

  probs <- predict(fit, newdata = pred_df[ok, ], type = "prob")[, class_levels, drop = FALSE]
  pred_class <- factor(class_levels[max.col(probs, ties.method = "first")],
                       levels = class_levels)
  obs <- pix$obs[ok]

  cm <- table(observed = obs, predicted = pred_class)

  # held-out prediction raster on the terrain grid (cells not scored stay NA)
  if (write_raster) {
    cell_ids <- pix$cell[ok]
    layer_from <- function(vals) {
      r <- terra::rast(topo[[1]])
      terra::values(r) <- NA_real_
      r[cell_ids] <- vals
      r
    }
    class_tbl <- data.frame(id = seq_along(class_levels), label = class_levels)
    obs_layer  <- layer_from(as.integer(obs))
    pred_layer <- layer_from(as.integer(pred_class))
    levels(obs_layer)  <- class_tbl
    levels(pred_layer) <- class_tbl

    out_ras <- c(obs_layer, pred_layer,
                 layer_from(as.integer(obs == pred_class)),
                 layer_from(probs[, 1]), layer_from(probs[, 2]), layer_from(probs[, 3]))
    names(out_ras) <- c("observed", "predicted", "correct",
                        paste0("prob_", gsub(" ", "", class_levels)))

    out_path <- file.path(
      raster_dir,
      paste0(gsub("[^A-Za-z0-9_-]", "_", burn), "_entry", entry, "_fold", fold, "_unitCV_pred.tif")
    )
    writeRaster(out_ras, out_path, overwrite = TRUE)
    message("  Saved held-out prediction raster: ", basename(out_path))
  }

  obs_prop  <- as.numeric(prop.table(table(obs)))
  pred_prop <- as.numeric(prop.table(table(pred_class)))
  soft_prop <- as.numeric(colMeans(probs))
  base_prop <- baseline_props[[as.character(fold)]]

  list(
    metrics = tibble(
      burn_name  = burn, burn_entry = entry, fold = fold,
      n_pix      = length(obs),
      n_excluded = sum(!ok),
      accuracy   = sum(diag(cm)) / sum(cm),
      kappa      = cohens_kappa(cm),
      tvd_model     = 0.5 * sum(abs(pred_prop - obs_prop)),
      tvd_soft      = 0.5 * sum(abs(soft_prop - obs_prop)),
      tvd_baseline  = 0.5 * sum(abs(base_prop - obs_prop)),
      obs_dominant  = class_levels[which.max(obs_prop)],
      pred_dominant = class_levels[which.max(pred_prop)]
    ),
    props = tibble(
      burn_name = burn, burn_entry = entry, fold = fold, n_pix = length(obs),
      class = factor(class_levels, levels = class_levels),
      obs_prop = obs_prop, pred_prop = pred_prop,
      soft_prop = soft_prop, baseline_prop = base_prop
    ),
    cm = cm
  )
}

# which events get a prediction raster
render_keys <- events %>%
  arrange(fold, burn_name, burn_entry) %>%
  group_by(fold) %>%
  slice_head(n = n_render_per_fold) %>%
  ungroup() %>%
  mutate(key = paste(burn_name, burn_entry)) %>%
  pull(key)
events$render <- paste(events$burn_name, events$burn_entry) %in% render_keys

results <- pmap(
  list(events$burn_name, events$burn_entry, events$burn_year, events$fold, events$render),
  function(burn, entry, burn_year, fold, render) {
    message("Validating ", burn, " entry ", entry, " (fold ", fold, ")")
    tryCatch(
      validate_event(burn, entry, burn_year, fold, write_raster = render),
      error = function(e) {
        warning("FAILED ", burn, " entry ", entry, ": ", conditionMessage(e))
        NULL
      }
    )
  }
) %>% setNames(paste0(events$burn_name, "_entry", events$burn_entry))

failed_events <- names(results)[vapply(results, is.null, logical(1))]
failed_events
results <- compact(results)

# ---- 4. summaries ----

event_metrics <- map_dfr(results, "metrics") %>%
  mutate(
    dominant_match = obs_dominant == pred_dominant,
    # > 0 means the model's class shares are closer to truth than the baseline
    tvd_skill = 1 - tvd_model / tvd_baseline
  )

event_props <- map_dfr(results, "props") %>%
  mutate(err_hard = pred_prop - obs_prop,
         err_soft = soft_prop - obs_prop,
         err_base = baseline_prop - obs_prop)

pooled_confusion <- Reduce(`+`, map(results, "cm"))

# pixel-level per-class accuracy, pooled over all held-out pixels
class_accuracy <- tibble(
  class = class_levels,
  n_obs_pix = as.numeric(rowSums(pooled_confusion)),
  pooled_obs_prop  = n_obs_pix / sum(n_obs_pix),
  pooled_pred_prop = as.numeric(colSums(pooled_confusion)) / sum(pooled_confusion),
  recall    = diag(pooled_confusion) / rowSums(pooled_confusion),
  precision = diag(pooled_confusion) / colSums(pooled_confusion)
)

pooled_summary <- tibble(
  n_events = nrow(event_metrics),
  n_pixels = sum(pooled_confusion),
  pixel_accuracy    = sum(diag(pooled_confusion)) / sum(pooled_confusion),
  balanced_accuracy = mean(class_accuracy$recall, na.rm = TRUE),
  kappa = cohens_kappa(pooled_confusion)
)

# event-level (each event weighted equally)
event_summary <- event_metrics %>%
  summarise(
    n_events = n(),
    mean_accuracy   = mean(accuracy),
    median_accuracy = median(accuracy),
    mean_kappa      = mean(kappa, na.rm = TRUE),
    mean_tvd_model    = mean(tvd_model),
    mean_tvd_soft     = mean(tvd_soft),
    mean_tvd_baseline = mean(tvd_baseline),
    share_events_beating_baseline = mean(tvd_model < tvd_baseline),
    share_dominant_class_correct  = mean(dominant_match)
  )

# class-share agreement across events (the operational-planning metric)
class_proportion_summary <- event_props %>%
  group_by(class) %>%
  summarise(
    n_events      = n(),
    mean_obs      = mean(obs_prop),
    mean_pred     = mean(pred_prop),
    bias_hard     = mean(err_hard),
    mae_hard      = mean(abs(err_hard)),
    rmse_hard     = sqrt(mean(err_hard^2)),
    r_hard        = suppressWarnings(cor(obs_prop, pred_prop)),
    within_tol_hard = mean(abs(err_hard) <= prop_tol),
    bias_soft     = mean(err_soft),
    mae_soft      = mean(abs(err_soft)),
    r_soft        = suppressWarnings(cor(obs_prop, soft_prop)),
    within_tol_soft = mean(abs(err_soft) <= prop_tol),
    mae_baseline  = mean(abs(err_base)),
    .groups = "drop"
  )

fold_summary <- event_metrics %>%
  group_by(fold) %>%
  summarise(n_events = n(), n_pix = sum(n_pix),
            mean_accuracy = mean(accuracy), mean_tvd_model = mean(tvd_model),
            mean_tvd_baseline = mean(tvd_baseline), .groups = "drop")

cat("\nPooled held-out pixels:\n"); print(pooled_summary)
cat("\nEvent-level (events weighted equally):\n"); print(t(event_summary))
cat("\nFold summary:\n"); print(fold_summary)
cat("\nPooled confusion (rows observed, cols predicted):\n"); print(pooled_confusion)
cat("\nPer-class pixel accuracy:\n"); print(class_accuracy)
cat("\nClass-share agreement across events:\n"); print(t(class_proportion_summary))

event_metrics %>% arrange(desc(tvd_model))

# predicted vs observed class share per event, 1:1 line dashed
ggplot(event_props, aes(obs_prop, pred_prop)) +
  geom_abline(linetype = "dashed") +
  geom_point(alpha = 0.6) +
  facet_wrap(~ class) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Observed share of event", y = "Predicted share of event (hard classes)")

# ---- 5. save ----

unit_validation <- list(
  unit_folds = unit_folds, event_metrics = event_metrics, event_props = event_props,
  class_proportion_summary = class_proportion_summary, pooled_confusion = pooled_confusion,
  class_accuracy = class_accuracy, pooled_summary = pooled_summary,
  event_summary = event_summary, fold_summary = fold_summary,
  fold_models = fold_models, failed_events = failed_events
)

out_dir_val <- "output/unit_validation"
dir.create(out_dir_val, showWarnings = FALSE, recursive = TRUE)
write_csv(event_metrics, file.path(out_dir_val, "event_metrics.csv"))
write_csv(event_props, file.path(out_dir_val, "event_class_proportions.csv"))
write_csv(class_proportion_summary, file.path(out_dir_val, "class_proportion_summary.csv"))

# ---- 6. observed vs predicted maps, one JPG per rendered event ----
# Reads the GeoTIFFs written in section 3 (so only events selected by
# n_render_per_fold get a map; set it to Inf above to map every event).

jpg_dir <- file.path(raster_dir, "jpg")
dir.create(jpg_dir, showWarnings = FALSE, recursive = TRUE)

class_cols <- c("Closed" = "#2c7a3f", "Open" = "#f2c14e", "Early Seral" = "#c8553d")

# class layers read back as labelled factors or as raw 1:3 ids, depending on
# whether terra restored the category table
to_class <- function(x) {
  if (is.factor(x)) factor(as.character(x), levels = class_levels)
  else factor(x, levels = seq_along(class_levels), labels = class_levels)
}

for (i in seq_len(nrow(event_metrics))) {
  m    <- event_metrics[i, ]
  stem <- paste0(gsub("[^A-Za-z0-9_-]", "_", m$burn_name),
                 "_entry", m$burn_entry, "_fold", m$fold, "_unitCV_pred")
  tif  <- file.path(raster_dir, paste0(stem, ".tif"))
  if (!file.exists(tif)) next

  map_df <- terra::rast(tif)[[c("observed", "predicted")]] %>%
    as.data.frame(xy = TRUE) %>%
    mutate(across(c(observed, predicted), to_class)) %>%
    pivot_longer(c(observed, predicted), names_to = "layer", values_to = "class") %>%
    mutate(layer = factor(layer, levels = c("observed", "predicted")))

  p <- ggplot(map_df, aes(x, y, fill = class)) +
    geom_raster() +
    facet_wrap(~ layer) +
    scale_fill_manual(values = class_cols, drop = FALSE, na.value = "grey90") +
    coord_quickmap() +
    labs(title = paste0(m$burn_name, " (entry ", m$burn_entry, ", fold ", m$fold, ")"),
         subtitle = sprintf("pixel accuracy %.2f, kappa %.2f", m$accuracy, m$kappa),
         x = NULL, y = NULL, fill = NULL)

  ggsave(file.path(jpg_dir, paste0(stem, ".jpg")), p, width = 10, height = 5, dpi = 150)
}
