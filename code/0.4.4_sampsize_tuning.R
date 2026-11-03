#### Tune partial class balancing (sampsize) for out-of-unit class-share fidelity
#
# Full balancing (equal sampsize per class) fixes the collapse to "Closed" but
# over-predicts Open and Early Seral. This script scans partial balancing
# between the natural class mix and equal classes, using the same
# burn-unit-blocked folds as 0.4.3, and scores each setting on how well the
# held-out events' predicted class shares match their observed shares.
#
# Balancing parameter alpha in [0, 1] (per fold, from that fold's training rows):
#   q_c       = (1 - alpha) * natural_share_c + alpha / n_classes
#   sampsize_c = round(total * q_c)
# alpha = 0 -> natural mix, alpha = 1 -> equal classes. `total` is the number
# of rows drawn per tree (sum of sampsize). A reference row with randomForest's
# default sampling (no sampsize) is also included.
#
# Selection criterion: mean over held-out events of the total-variation
# distance (TVD) between hard-classified predicted class shares and observed
# class shares. Per-class bias, soft (probability) TVD, balanced accuracy and
# kappa are reported alongside. "best" minimises mean TVD; "parsimonious" is
# the least-balanced setting within 1 SE of the best (a one-standard-error
# rule; SE is across events, which are not fully independent).
#
# Speed: instead of predicting every pixel for every setting, a random sample
# of up to `n_pix_cap` scored pixels per event is cached once and reused for
# all settings. Class shares from that sample are unbiased estimates of the
# event's full-pixel shares.
#
# Caveats: tuning and scoring use the same folds, so the best setting's
# numbers are somewhat optimistic; re-score the chosen setting with a new
# fold seed in 0.4.3 before reporting. Use a fresh fold seed there to confirm.
#
# Depends on objects from 0.4.2 and 0.4.3 (run both first):
#   rf_data_unit, predictors_unit, rx_pool, terrain_stack, eric_stack,
#   cc_stack, ecozones, insdis_stack, unit_folds, events, unit_constants,
#   insdis_layer_names, insdis_layer_years, extract_xy, extract_layer_or_na,
#   cohens_kappa, class_levels, k
#
# Outputs: tuning_results (one row per setting), tune_events (event-level TVD
# per setting), chosen_sampsize; CSV in output/unit_validation/.

library(tidyverse)
library(randomForest)
library(sf)
library(terra)

needed <- c("rf_data_unit", "predictors_unit", "rx_pool", "terrain_stack", "eric_stack",
            "cc_stack", "ecozones", "insdis_stack", "unit_folds", "events",
            "unit_constants", "insdis_layer_names", "insdis_layer_years",
            "extract_xy", "extract_layer_or_na", "cohens_kappa", "class_levels", "k")
missing_objs <- needed[!vapply(needed, exists, logical(1))]
if (length(missing_objs) > 0) {
  stop("Missing: ", paste(missing_objs, collapse = ", "),
       " -- run 0.4.2 and 0.4.3 first.")
}

# ---- settings (EDIT HERE) ----

alpha_grid  <- c(0, 0.25, 0.5, 0.75, 1)  # 0 = natural mix, 1 = equal classes
total_grid  <- c(600, 900)               # rows drawn per tree (sum of sampsize)
ntree_tune  <- 500
mtry_tune   <- 6
n_pix_cap   <- 5000                      # scored pixels cached per event
seed_tune   <- 8153

cv_formula <- reformulate(predictors_unit, response = "Severity_1")
n_classes  <- length(class_levels)
bias_names <- paste0("bias_", gsub(" ", "", class_levels))

train_data <- rf_data_unit %>%
  left_join(unit_folds %>% select(burn_name, fold), by = "burn_name")

fold_class_counts <- map(seq_len(k), function(f) {
  table(factor(train_data$Severity_1[train_data$fold != f], levels = class_levels))
})
fold_natural_props <- map(fold_class_counts, ~ as.numeric(prop.table(.x)))

make_sampsize <- function(natural_props, alpha, total) {
  q <- (1 - alpha) * natural_props + alpha / n_classes
  setNames(pmax(round(total * q), 1), class_levels)
}

# ---- 1. cache a pixel sample for every held-out event ----

set.seed(seed_tune)

build_event_pixels <- function(burn, entry, burn_year, fold) {
  unit_const <- unit_constants %>%
    filter(burn_name == burn, burn_entry == entry) %>%
    select(-burn_name)

  ev_v <- rx_pool %>% filter(burn_name == burn, burn_entry == entry) %>% terra::vect()

  topo  <- crop(terrain_stack, ev_v) %>% mask(ev_v)
  obs_r <- terra::rasterize(ev_v, topo[[1]], field = "Severity_1")

  pix <- as.data.frame(topo, xy = TRUE, na.rm = FALSE)
  pix$obs <- factor(terra::values(obs_r)[, 1], levels = 0:2, labels = class_levels)
  pix <- pix[!is.na(pix$obs) & !is.na(pix$chili), ]
  if (nrow(pix) == 0) return(NULL)

  xy <- as.matrix(pix[, c("x", "y")])

  # same predictor rules as 0.4.3 / training
  window_cols <- insdis_layer_names[insdis_layer_years %in% (burn_year - 5):(burn_year - 1)]
  insect <- if (length(window_cols) == 0) {
    rep(NA_integer_, nrow(xy))
  } else {
    as.integer(rowSums(as.matrix(extract_xy(insdis_stack[[window_cols]], xy)) == 1,
                       na.rm = TRUE) > 0)
  }

  df <- data.frame(
    chili = pix$chili, tpi = pix$tpi, slp.dg = pix$slp.dg, elev = pix$elev,
    ecozone = factor(as.character(extract_xy(ecozones[[1]], xy)[[1]]),
                     levels = levels(rf_data_unit$ecozone)),
    insect_disease_5yr = insect,
    ndvi_winter_prior_year = factor(
      as.character(extract_layer_or_na(eric_stack, paste0("y_", burn_year - 1), xy)),
      levels = levels(rf_data_unit$ndvi_winter_prior_year)),
    canopy_cover_prior_year = extract_layer_or_na(
      cc_stack, paste0(burn_year - 1, "_TRE_", burn_year - 1), xy)
  ) %>%
    bind_cols(unit_const[rep(1, nrow(pix)), ])

  df$obs <- pix$obs
  df <- df[complete.cases(df[, predictors_unit]), ]
  if (nrow(df) == 0) return(NULL)
  if (nrow(df) > n_pix_cap) df <- df[sample.int(nrow(df), n_pix_cap), ]

  df$burn_name <- burn
  df$fold <- fold
  df
}

pix_cache <- pmap(
  list(events$burn_name, events$burn_entry, events$burn_year, events$fold),
  function(burn, entry, burn_year, fold) {
    message("Caching pixels: ", burn, " entry ", entry)
    tryCatch(build_event_pixels(burn, entry, burn_year, fold),
             error = function(e) {
               warning("FAILED ", burn, " entry ", entry, ": ", conditionMessage(e))
               NULL
             })
  }
) %>%
  compact() %>%
  bind_rows()

# ---- 2. fit fold models and score each setting ----

evaluate_config <- function(alpha, total) {
  use_default <- is.na(alpha)

  preds <- map_dfr(seq_len(k), function(f) {
    train_f <- train_data %>% filter(fold != f) %>% select(-fold)

    set.seed(seed_tune + f)
    fit <- if (use_default) {
      randomForest(cv_formula, data = train_f, ntree = ntree_tune, mtry = mtry_tune)
    } else {
      randomForest(cv_formula, data = train_f, ntree = ntree_tune, mtry = mtry_tune,
                   sampsize = make_sampsize(fold_natural_props[[f]], alpha, total))
    }

    test  <- pix_cache %>% filter(fold == f)
    probs <- predict(fit, newdata = test, type = "prob")[, class_levels, drop = FALSE]
    colnames(probs) <- paste0("p", seq_len(n_classes))

    tibble(burn_name = test$burn_name, burn_entry = test$burn_entry, obs = test$obs,
           pred = factor(class_levels[max.col(probs, ties.method = "first")],
                         levels = class_levels)) %>%
      bind_cols(as_tibble(probs))
  })

  # per-event class shares
  key    <- paste(preds$burn_name, preds$burn_entry, sep = "|")
  n_ev   <- rowsum(rep(1, nrow(preds)), key)[, 1]
  share  <- function(m) rowsum(m * 1, key) / n_ev
  obs_p  <- share(vapply(class_levels, function(cl) preds$obs == cl, logical(nrow(preds))))
  pred_p <- share(vapply(class_levels, function(cl) preds$pred == cl, logical(nrow(preds))))
  soft_p <- share(as.matrix(preds[, paste0("p", seq_len(n_classes))]))

  tvd_hard <- 0.5 * rowSums(abs(pred_p - obs_p))
  tvd_soft <- 0.5 * rowSums(abs(soft_p - obs_p))
  bias <- colMeans(pred_p - obs_p)

  cm <- table(observed = preds$obs, predicted = preds$pred)

  max_draw_ratio <- if (use_default) NA_real_ else {
    max(map_dbl(seq_len(k), function(f) {
      max(make_sampsize(fold_natural_props[[f]], alpha, total) / as.numeric(fold_class_counts[[f]]))
    }))
  }

  list(
    summary = tibble(
      alpha = alpha, total = total,
      mean_tvd_hard = mean(tvd_hard),
      se_tvd_hard   = sd(tvd_hard) / sqrt(length(tvd_hard)),
      mean_tvd_soft = mean(tvd_soft),
      mean_abs_bias = mean(abs(bias)),
      balanced_accuracy = mean(diag(cm) / rowSums(cm), na.rm = TRUE),
      kappa = cohens_kappa(cm),
      pixel_accuracy = sum(diag(cm)) / sum(cm),
      max_draw_ratio = max_draw_ratio
    ) %>% bind_cols(as_tibble(as.list(setNames(bias, bias_names)))),
    events = tibble(event = names(tvd_hard), tvd_hard = as.numeric(tvd_hard))
  )
}

configs <- bind_rows(
  tibble(alpha = NA_real_, total = NA_real_),
  expand_grid(alpha = alpha_grid, total = total_grid)
)

tune_out <- pmap(configs, function(alpha, total) {
  message("Setting: alpha = ", alpha, ", total = ", total)
  evaluate_config(alpha, total)
})

tuning_results <- map_dfr(tune_out, "summary")
tune_events <- map(tune_out, "events") %>%
  setNames(paste0("alpha", configs$alpha, "_total", configs$total))

# ---- 3. choose a setting ----

tuned <- tuning_results %>% filter(!is.na(alpha))

best <- tuned %>% slice_min(mean_tvd_hard, n = 1, with_ties = FALSE)

parsimonious <- tuned %>%
  filter(mean_tvd_hard <= best$mean_tvd_hard + best$se_tvd_hard) %>%
  arrange(alpha, total) %>%
  slice(1)

cat("\nAll settings (default sampling row has alpha = NA):\n")
tuning_results %>% arrange(mean_tvd_hard) %>% print(width = Inf)
cat("\nBest (min mean TVD):\n"); print(best, width = Inf)
cat("\nParsimonious (least balanced within 1 SE of best):\n"); print(parsimonious, width = Inf)

# sampsize to use for the final all-units model in 0.4.2 for a chosen setting
chosen <- parsimonious
chosen_sampsize <- make_sampsize(as.numeric(prop.table(table(rf_data_unit$Severity_1))),
                                 chosen$alpha, chosen$total)
chosen_sampsize

# ---- 4. plots ----

# class-share error vs balancing
ggplot(tuned, aes(alpha, mean_tvd_hard, color = factor(total), group = total)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = tuning_results$mean_tvd_hard[is.na(tuning_results$alpha)],
             linetype = "dashed") +
  labs(x = "Balancing alpha (0 = natural mix, 1 = equal classes)",
       y = "Mean TVD of class shares (hard)", color = "Rows per tree",
       caption = "Dashed line: randomForest default sampling")

# per-class bias vs balancing: where does over-prediction cross zero?
tuned %>%
  pivot_longer(all_of(bias_names), names_to = "class", values_to = "bias") %>%
  mutate(class = sub("^bias_", "", class)) %>%
  ggplot(aes(alpha, bias, color = class)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ total, labeller = label_both) +
  labs(x = "Balancing alpha (0 = natural mix, 1 = equal classes)",
       y = "Mean predicted - observed class share")

out_dir_val <- "output/unit_validation"
dir.create(out_dir_val, showWarnings = FALSE, recursive = TRUE)
write_csv(tuning_results, file.path(out_dir_val, "sampsize_tuning.csv"))
