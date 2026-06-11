# ============================================================================
# KBDI Vector Function - Returns a vector for use within mutate()
# ============================================================================
#
# This function calculates KBDI as a vector, allowing it to be called
# within dplyr::mutate(). It assumes the data is already grouped and ordered
# correctly by the calling code.
#
# Original Reference:
#   Keetch, J.J. and Byram, G.M. 1968. A drought index for forest fire control.
#   Research Paper SE-38. USDA Forest Service, Southeastern Forest Experiment Station.
#   Corrected: Crane, R.W. 1983; Wikifire.ch

calc_kbdi_vector <- function(temp_vec,
                             precip_vec,
                             map_vec = 30) {
  
  #' Calculate KBDI as a vector
  #'
  #' @param temp_vec Vector of daily maximum temperatures (°F)
  #' @param precip_vec Vector of daily precipitation (inches)
  #' @param map_vec Vector or scalar of mean annual precipitation (inches).
  #'   Can be a single value or a vector matching the length of temp_vec.
  #'
  #' @return A numeric vector of KBDI values (0-800 scale)
  #'
  #' @details
  #' Assumes input vectors are in chronological order and represent a
  #' continuous time series. If called within dplyr::mutate() with group_by(),
  #' ensure the data is arranged by date within each group.
  #'
  #' KBDI formula (corrected - Crane 1983):
  #'   KBDI_t = Q + [(800-Q) * (0.968*exp(0.0486*T) - 8.30) * 1] / [1 + 10.88*exp(-0.0441*P)] * 0.001
  #'
  #' Where:
  #'   Q = KBDI_{t-1} - (net_rainfall * 100)
  #'   T = daily maximum temperature (°F)
  #'   P = mean annual precipitation (inches)
  #'   
  #' Net rainfall calculation:
  #'   - New rain event begins when precip > 0 after dry day(s)
  #'   - Within event, accumulate daily precipitation
  #'   - When cumulative event precip > 0.2 in, net_rainfall = cumul - 0.2
  #'   - Otherwise net_rainfall = 0
  
  n <- length(temp_vec)
  
  # Validate inputs
  if (length(precip_vec) != n) {
    stop("temp_vec and precip_vec must have the same length")
  }
  
  # Handle map_vec
  if (length(map_vec) == 1) {
    map_vec <- rep(map_vec, n)
  } else if (length(map_vec) != n) {
    stop("map_vec must be length 1 or match the length of temp_vec")
  }
  
  # Handle NAs in input vectors
  # Missing precip: treat as 0 (no recorded rainfall)
  # Missing temp: drought increment will be 0 for that day (KBDI carried forward)
  n_na_precip <- sum(is.na(precip_vec))
  n_na_temp   <- sum(is.na(temp_vec))
  if (n_na_precip > 0) {
    warning(n_na_precip, " NA(s) in precip_vec replaced with 0")
    precip_vec[is.na(precip_vec)] <- 0
  }
  if (n_na_temp > 0) {
    warning(n_na_temp, " NA(s) in temp_vec — drought increment will be 0 for those days")
  }
  
  # Initialize output vector
  kbdi <- numeric(n)
  
  # =========================================================================
  # Pre-calculate net rainfall for each day
  # =========================================================================
  # Identify rain events: new event starts when precip > 0 after dry day(s).
  # Within a rain event, the first 0.2" does not reduce drought (threshold).
  # Once the threshold is crossed, all additional rain in the event reduces
  # drought incrementally. This prevents double-counting in multi-day events.
  is_rain_day <- precip_vec > 0
  
  net_rainfall <- numeric(n)
  cumul_in_event <- 0  # running cumulative precip for the current rain event
  
  for (i in seq_along(precip_vec)) {
    if (is_rain_day[i]) {
      # Accumulate within rain event; reset cumulative when a new event begins
      if (i == 1 || !is_rain_day[i - 1]) {
        cumul_in_event <- precip_vec[i]
      } else {
        cumul_in_event <- cumul_in_event + precip_vec[i]
      }
      
      # Net for entire event up to today, and up to yesterday
      net_cumul_today <- max(0, cumul_in_event - 0.2)
      net_cumul_prev  <- max(0, (cumul_in_event - precip_vec[i]) - 0.2)
      
      # Today's incremental net rainfall (avoids double-counting multi-day events)
      net_rainfall[i] <- net_cumul_today - net_cumul_prev
      
    } else {
      cumul_in_event <- 0  # reset on dry days
    }
  }
  
  # =========================================================================
  # Main KBDI calculation loop
  # =========================================================================
  for (i in seq_along(temp_vec)) {
    
    # Get previous KBDI (0 for first day of each group)
    kbdi_prev <- if (i == 1) 0 else kbdi[i - 1]
    
    # -----------------------------------------------------------------------
    # Step 1: Adjust previous KBDI for net rainfall from today
    # -----------------------------------------------------------------------
    # Q = KBDI_{t-1} - net_rainfall_{t} * 100
    # (converts inches to hundredths of inch)
    Q <- kbdi_prev - (net_rainfall[i] * 100)
    Q <- max(0, Q)
    
    # -----------------------------------------------------------------------
    # Step 2: Calculate temperature-based drought factor
    # -----------------------------------------------------------------------
    # PROVENANCE: Keetch & Byram 1968, Equation 18 (corrected by Crane 1983)
    # The exponential component models evapotranspiration
    drought_factor_raw <- 0.968 * exp(0.0486 * temp_vec[i]) - 8.30
    drought_factor <- max(0, drought_factor_raw)
    
    # -----------------------------------------------------------------------
    # Step 3: Calculate transpiration dampening factor
    # -----------------------------------------------------------------------
    # Higher mean annual precip = more vegetation = higher damping
    # PROVENANCE: Keetch & Byram 1968, §3
    transpiration_factor <- 1 / (1 + 10.88 * exp(-0.0441 * map_vec[i]))
    
    # -----------------------------------------------------------------------
    # Step 4: Calculate daily KBDI increment
    # -----------------------------------------------------------------------
    # KBDI_t = Q + [(800-Q) * drought_factor * transpiration_factor] * 0.001
    kbdi_increment <- ((800 - Q) * drought_factor * transpiration_factor * 0.001)
    
    # Handle NaN/Inf (shouldn't happen, but be safe)
    if (!is.finite(kbdi_increment)) {
      kbdi_increment <- 0
    }
    
    # -----------------------------------------------------------------------
    # Step 5: Calculate KBDI for today and bound to [0, 800]
    # -----------------------------------------------------------------------
    kbdi[i] <- max(0, min(Q + kbdi_increment, 800))
  }
  
  return(kbdi)
}
