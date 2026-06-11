# ============================================================================
# Apply KBDI calculation to raws.burn.ts and save results
# ============================================================================

source("code/kbdi_calculation.R")

# Apply KBDI calculation with your data
# Using default MAP=30 inches (typical for much of the region)
# If you have station-specific data, pass mean_annual_precip as a named vector
raws.burn.ts <- calc_kbdi(raws.burn.ts, mean_annual_precip = 30)

# Run validation tests
cat("\n=== Running Validation Tests ===\n")
validate_kbdi(raws.burn.ts, 
              kbdi_col = "kbdi", 
              temp_col = "daily_max_temp",
              precip_col = "daily_precip")

# Summary statistics
cat("\n=== KBDI Summary Statistics ===\n")
summary_stats <- raws.burn.ts %>%
  summarise(
    n_records = n(),
    n_fires = n_distinct(Fire_key),
    date_range = paste(min(date), "to", max(date)),
    kbdi_min = min(kbdi, na.rm = TRUE),
    kbdi_q25 = quantile(kbdi, 0.25, na.rm = TRUE),
    kbdi_median = median(kbdi, na.rm = TRUE),
    kbdi_mean = mean(kbdi, na.rm = TRUE),
    kbdi_q75 = quantile(kbdi, 0.75, na.rm = TRUE),
    kbdi_max = max(kbdi, na.rm = TRUE)
  )
print(summary_stats)

# By-fire summary
cat("\n=== KBDI by Fire ===\n")
fire_summary <- raws.burn.ts %>%
  group_by(Fire_key, station) %>%
  summarise(
    n_days = n(),
    date_start = min(date),
    date_end = max(date),
    mean_kbdi = mean(kbdi, na.rm = TRUE),
    max_kbdi = max(kbdi, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(max_kbdi))

print(fire_summary, n = 10)

# Optional: Save to CSV if needed
# write.csv(raws.burn.ts, "data/raws_burn_ts_with_kbdi.csv", row.names = FALSE)

cat("\nKBDI calculation complete!\n")
