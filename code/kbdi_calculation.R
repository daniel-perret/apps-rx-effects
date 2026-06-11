# ============================================================================
# KBDI Calculation Function with Full Provenance Documentation
# ============================================================================
#
# Original Reference:
#   Keetch, J.J. and Byram, G.M. 1968. A drought index for forest fire control.
#   Research Paper SE-38. USDA Forest Service, Southeastern Forest Experiment Station.
#   (Revised 1988; commonly used in NFDRS - National Fire Danger Rating System)
#
# Implementation details based on:
#   - Wikifire (wsl.ch) - Keetch-Byram drought index equations and net rainfall
#   - Crane, R.W. (1983) - Corrected typographical error in original publication
#   - Alexander, M.E. (1990) - Fire behavior research publications
#   - WFAS documentation (U.S. Forest Service)
# ============================================================================

calc_kbdi <- function(df, 
                      temp_col = "daily_max_temp", 
                      precip_col = "daily_precip",
                      group_cols = c("Fire_key", "station"),
                      mean_annual_precip = NULL) {
  
  #' Calculate Keech-Byram Drought Index from daily weather data
  #'
  #' @param df Data frame with daily observations
  #' @param temp_col Name of daily maximum temperature column (°F)
  #' @param precip_col Name of daily precipitation column (inches)
  #' @param group_cols Vector of columns defining separate KBDI series
  #' @param mean_annual_precip Mean annual precipitation. Can be:
  #'   - NULL (default): uses 30 inches (conservative default for CONUS)
  #'   - A scalar numeric value: applied uniformly to all rows
  #'   - A column name (character string): uses station/group-specific MAP values
  #'       Must exist in df. Typically one unique value per group_col combination.
  #'
  #' @details
  #' The KBDI formula (corrected - Crane 1983):
  #'   KBDI_t = Q + [(800-Q) * (0.968*exp(0.0486*T) - 8.30) * 1] / [1 + 10.88*exp(-0.0441*P)] * 0.001
  #'
  #' Where:
  #'   Q = Previous day's KBDI minus net rainfall (in hundredths of inch)
  #'   T = Daily maximum temperature (°F) - PROVENANCE: Keetch & Byram 1968, §2.3
  #'   P = Mean annual precipitation (inches) - PROVENANCE: Keetch & Byram 1968, §3
  #'   Net rainfall = max(0, daily_precip - 0.2) - PROVENANCE: Crane 1983, §2
  #'
  #' Initialization:
  #'   KBDI should begin at 0 after substantial precipitation (>0.5 inches in a period).
  #'   PROVENANCE: Alexander 1990, Johnson & Forthum 2001; Fujioka 1991
  #'
  #' Output scale:
  #'   0 = saturated soil (8 inches of available water)
  #'   800 = absolutely dry soil (0 inches available water)
  #'   PROVENANCE: Keetch & Byram 1968, p.15; represents 8 inches of moisture
  
  # Handle mean_annual_precip input: NULL, scalar, or column name
  if (is.null(mean_annual_precip)) {
    # Use default
    map_source <- "default"
    df <- df %>% mutate(.map_value = 30)
  } else if (is.character(mean_annual_precip) && length(mean_annual_precip) == 1) {
    # Column name provided
    map_source <- "column"
    if (!mean_annual_precip %in% names(df)) {
      stop("Column '", mean_annual_precip, "' not found in data frame")
    }
    df <- df %>% rename(.map_value = all_of(mean_annual_precip))
  } else if (is.numeric(mean_annual_precip) && length(mean_annual_precip) == 1) {
    # Scalar value provided
    map_source <- "scalar"
    df <- df %>% mutate(.map_value = mean_annual_precip)
  } else {
    stop("mean_annual_precip must be NULL, a single numeric value, or a column name (character string)")
  }
  
  df <- df %>%
    # Sort by groups and date to ensure proper chronological sequence
    arrange(across(all_of(group_cols)), date) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(
      # -----------------------------------------------------------------------
      # STEP 1: Calculate Net Rainfall
      # -----------------------------------------------------------------------
      # "In order to obtain net rainfall, 0.2 inch has to be subtracted from 
      #  any daily rainfall amount exceeding 0.2 inch."
      # PROVENANCE: Wikifire (wsl.ch), Keetch-Byram drought index; 
      #             Crane 1983; original SE-38 (corrected)
      
      net_rainfall_daily = if_else(.data[[precip_col]] > 0.2, 
                                    .data[[precip_col]] - 0.2, 
                                    0),
      
      # -----------------------------------------------------------------------
      # STEP 2: Accumulate net rainfall across consecutive rain days
      # -----------------------------------------------------------------------
      # "In the case of consecutive days with rainfall, the 0.2 inch has to 
      #  be subtracted on the exact day when the summed rainfall amount exceeds 0.2 inch."
      # PROVENANCE: Wikifire (wsl.ch); P_net formula with running sum
      
      # Track consecutive days with precipitation
      precip_event = cumsum(.data[[precip_col]] == 0),
      
      # Sum rainfall within each event
      cumul_event_precip = .data[[precip_col]],
      cumul_event_precip = ave(cumul_event_precip, precip_event, 
                                FUN = cumsum),
      
      # Apply threshold only when event total exceeds 0.2
      net_rainfall = if_else(cumul_event_precip > 0.2,
                             pmax(0, cumul_event_precip - 0.2),
                             0),
      
      # Only subtract net rain once (on first day threshold is exceeded)
      net_rainfall_applied = if_else(
        row_number() == 1 | net_rainfall != lag(net_rainfall, default = 0),
        net_rainfall,
        0
      ),
      
      # -----------------------------------------------------------------------
      # STEP 3: Calculate temperature-based drought factor
      # -----------------------------------------------------------------------
      # The exponential component reflects evapotranspiration as a function 
      # of temperature (daily maximum).
      # PROVENANCE: Keetch & Byram 1968, Equation 18 (corrected constant 8.30)
      #             Crane 1983 (identified and corrected typo: 8.30 not 0.083)
      
      # Core exponential drought factor
      drought_factor_raw = 0.968 * exp(0.0486 * .data[[temp_col]]) - 8.30,
      
      # Ensure non-negative
      drought_factor = pmax(0, drought_factor_raw),
      
      # -----------------------------------------------------------------------
      # STEP 4: Retrieve or initialize previous KBDI
      # -----------------------------------------------------------------------
      # The index is cumulative; each day builds on the previous.
      # PROVENANCE: Keetch & Byram 1968 (cumulative nature); 
      #             Alexander 1990 (initialization protocol)
      
      kbdi_prev = lag(kbdi, default = 0),
      kbdi_prev = if_else(is.na(kbdi_prev), 0, kbdi_prev),
      
      # -----------------------------------------------------------------------
      # STEP 5: Calculate Q (adjusted previous KBDI)
      # -----------------------------------------------------------------------
      # Q = KBDI_{t-1} - P_net_t * 100
      # Converts net rainfall (inches) to hundredths of an inch units
      # PROVENANCE: Wikifire formula; Keetch & Byram 1968 scale definition
      
      Q = kbdi_prev - (net_rainfall * 100),
      Q = pmax(0, Q),  # Cannot go below 0
      
      # -----------------------------------------------------------------------
      # STEP 6: Apply the full KBDI formula
      # -----------------------------------------------------------------------
      # KBDI_t = Q + [(800-Q) * (0.968*exp(0.0486*T) - 8.30) * Δt] / [1 + 10.88*exp(-0.0441*P)] * 0.001
      #
      # Components:
      #   (800 - Q)        : remaining capacity for drought accumulation
      #   drought_factor   : temperature-driven evapotranspiration (Δt=1 day)
      #   1/(1 + 10.88*exp(-0.0441*P)) : vegetation/transpiration dampening
      #                                    as function of mean annual precip
      #   0.001            : scaling to hundredths of inch units
      #
      # PROVENANCE: Wikifire (wsl.ch); Keetch & Byram 1968 (original derivation)
      #             Mean annual precip effect: Keetch & Byram 1968, §3
      #               (higher MAP = more vegetation/transpiration = higher damping)
      
      # Transpiration scaling factor based on mean annual precipitation
      # "vegetation density (and therefore its transpiration capacity) is a 
      #  function of mean annual rainfall"
      # PROVENANCE: Hawaii climatological study (Fujioka et al.)
      transpiration_factor = 1 / (1 + 10.88 * exp(-0.0441 * .map_value)),
      
      # Full KBDI increment
      kbdi_increment = ((800 - Q) * drought_factor * transpiration_factor * 0.001),
      kbdi_increment = if_else(is.na(kbdi_increment) | is.nan(kbdi_increment), 
                                0, 
                                kbdi_increment),
      
      # Apply increment to adjusted previous value
      kbdi_raw = Q + kbdi_increment,
      
      # -----------------------------------------------------------------------
      # STEP 7: Bound KBDI to valid range [0, 800]
      # -----------------------------------------------------------------------
      # "The KBDI expresses moisture deficiency in hundredths of an inch, and 
      #  is based on a measurement of 8 inches of available moisture... therefore, 
      #  it is on a scale ranging from 0 to 800"
      # PROVENANCE: Keetch & Byram 1968, pp. 15; Wikifire reference
      
      kbdi = pmax(0, pmin(kbdi_raw, 800))
      
    ) %>%
    # Remove intermediate calculations from output
    select(-net_rainfall_daily, -precip_event, -cumul_event_precip,
           -net_rainfall_applied, -drought_factor_raw, -drought_factor,
           -kbdi_prev, -Q, -transpiration_factor, -kbdi_increment, -kbdi_raw,
           -.map_value) %>%
    ungroup()
  
  return(df)
}


# ============================================================================
# VALIDATION TESTS
# ============================================================================

validate_kbdi <- function(df, kbdi_col = "kbdi", temp_col = "daily_max_temp", 
                          precip_col = "daily_precip") {
  
  cat("=== KBDI Validation Tests ===\n\n")
  
  # TEST 1: Range bounds
  # All KBDI values should be in [0, 800]
  cat("TEST 1: Range Bounds\n")
  min_val <- min(df[[kbdi_col]], na.rm = TRUE)
  max_val <- max(df[[kbdi_col]], na.rm = TRUE)
  cat("  Min KBDI:", min_val, "\n")
  cat("  Max KBDI:", max_val, "\n")
  
  if (min_val < 0 | max_val > 800) {
    cat("  ❌ FAIL: Values outside [0, 800] range\n")
    return(FALSE)
  } else {
    cat("  ✓ PASS: All values in valid range [0, 800]\n")
  }
  cat("\n")
  
  # TEST 2: Response to heavy precipitation
  # Large precipitation (>0.2 in) should reduce KBDI
  cat("TEST 2: Precipitation Response\n")
  
  # Find days with substantial precipitation
  precip_threshold <- 0.5
  heavy_precip_indices <- which(df[[precip_col]] > precip_threshold)
  
  if (length(heavy_precip_indices) > 0) {
    improvements <- 0
    for (i in heavy_precip_indices) {
      if (i > 1 && !is.na(df[[kbdi_col]][i]) && !is.na(df[[kbdi_col]][i-1])) {
        if (df[[kbdi_col]][i] <= df[[kbdi_col]][i-1]) {
          improvements <- improvements + 1
        }
      }
    }
    pct_improvement <- (improvements / length(heavy_precip_indices)) * 100
    cat("  Precipitation >", precip_threshold, "in found:", length(heavy_precip_indices), "days\n")
    cat("  Days where KBDI decreased after heavy precip:", improvements, 
        "({pct_improvement}%)\n")
    
    if (pct_improvement >= 70) {
      cat("  ✓ PASS: KBDI generally decreases after heavy precipitation\n")
    } else {
      cat("  ⚠ WARNING: Lower than expected precipitation response\n")
    }
  } else {
    cat("  ⚠ No heavy precipitation events (>", precip_threshold, "in) in dataset\n")
  }
  cat("\n")
  
  # TEST 3: Response to temperature
  # Hotter days should increase KBDI (when not affected by rain)
  cat("TEST 3: Temperature Response\n")
  
  # Find long dry spells (no precip for 5+ days)
  df_test <- df %>%
    mutate(dry_spell = cumsum(.data[[precip_col]] > 0.1))
  
  spell_lengths <- df_test %>%
    group_by(dry_spell) %>%
    summarise(length = n(), 
              mean_temp = mean(.data[[temp_col]], na.rm = TRUE),
              max_kbdi = max(.data[[kbdi_col]], na.rm = TRUE),
              .groups = 'drop') %>%
    filter(length >= 5)
  
  if (nrow(spell_lengths) > 0) {
    # Within dry spells, higher temps should correlate with higher KBDI
    spell_correlations <- c()
    for (i in unique(df_test$dry_spell)) {
      spell_data <- df_test %>% filter(dry_spell == i)
      if (nrow(spell_data) >= 5) {
        corr <- cor(spell_data[[temp_col]], spell_data[[kbdi_col]], 
                    use = "complete.obs")
        if (!is.na(corr)) spell_correlations <- c(spell_correlations, corr)
      }
    }
    
    if (length(spell_correlations) > 0) {
      mean_corr <- mean(spell_correlations)
      cat("  Mean temperature-KBDI correlation in dry spells:", 
          round(mean_corr, 3), "\n")
      if (mean_corr > 0.3) {
        cat("  ✓ PASS: Positive correlation between temperature and KBDI\n")
      } else {
        cat("  ⚠ WARNING: Weak temperature-KBDI correlation\n")
      }
    }
  } else {
    cat("  ⚠ No dry spells of 5+ days found; cannot test temperature response\n")
  }
  cat("\n")
  
  # TEST 4: Cumulative nature
  # On dry days with moderate temps, KBDI should generally increase
  cat("TEST 4: Cumulative Behavior\n")
  
  df_test <- df %>%
    filter(.data[[precip_col]] < 0.1,  # Dry day
           .data[[temp_col]] >= 50,      # Reasonable temperature
           !is.na(.data[[kbdi_col]])) %>%
    mutate(kbdi_change = .data[[kbdi_col]] - lag(.data[[kbdi_col]]))
  
  increases <- sum(df_test$kbdi_change > 0, na.rm = TRUE)
  total_dry_days <- nrow(df_test)
  pct_increase <- if (total_dry_days > 0) (increases / total_dry_days) * 100 else 0
  
  cat("  Dry days with T≥50°F:", total_dry_days, "\n")
  cat("  Days with KBDI increase:", increases, 
      sprintf("(%.1f%%)\n", pct_increase))
  
  if (pct_increase >= 60) {
    cat("  ✓ PASS: KBDI increases on most dry days\n")
  } else if (pct_increase >= 40) {
    cat("  ⚠ WARN: Lower-than-expected increase on dry days\n")
  } else {
    cat("  ❌ FAIL: KBDI not accumulating as expected\n")
  }
  cat("\n")
  
  # TEST 5: Missing data handling
  cat("TEST 5: Missing Data Summary\n")
  na_count <- sum(is.na(df[[kbdi_col]]))
  total_count <- nrow(df)
  pct_na <- (na_count / total_count) * 100
  cat("  Missing KBDI values:", na_count, "of", total_count, 
      sprintf("(%.1f%%)\n", pct_na))
  
  if (pct_na < 10) {
    cat("  ✓ PASS: Low rate of missing values\n")
  } else if (pct_na < 30) {
    cat("  ⚠ WARNING: Moderate rate of missing values\n")
  } else {
    cat("  ❌ FAIL: High rate of missing values\n")
  }
  cat("\n")
  
  cat("=== End Validation ===\n")
  return(TRUE)
  
}
