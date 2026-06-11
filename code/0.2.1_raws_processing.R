#### Standardize WRCC RAWS xlsx files into a single combined CSV
#### Reads 10 station files from data/RAWS/, handles two column schemas,
#### and writes data/RAWS/RAWS_all_stations.csv

source("code/00_setup.R")

# ---- helper function ----

read_raws_file <- function(path) {
  # Extract station name from row 1
  station_name <- readxl::read_excel(path, n_max = 1, col_names = FALSE)[[1, 1]] %>%
    sub("\\s+Virginia.*", "", x = .) %>%
    trimws()

  df <- readxl::read_excel(path, skip = 3)

  # Two schemas exist across stations:
  #   Schema A (CraigValley, FortValley, LimeKiln, Marlinton):
  #     cols 9-13: fuel_moisture, temp_air_max, temp_air_min, rh_max, rh_min
  #   Schema B (Flatwoods, Glenped, HQ, SawmillRidge, StonyFork, UpperTract):
  #     cols 9-13: voltage (dropped), fuel_moisture, gust_dir (dropped),
  #                gust_speed (dropped), solar_rad (dropped)
  schema_b <- grepl("volt", colnames(df)[9], ignore.case = TRUE)

  if (schema_b) {
    df <- df %>%
      select(date = 1, time = 2, precip = 3,
             wind_speed = 4, wind_dir = 5,
             temp_air_avg = 6, temp_fuel = 7, rh = 8,
             fuel_moisture = 10) %>%
      mutate(temp_air_max = NA_real_, temp_air_min = NA_real_,
             rh_max = NA_real_, rh_min = NA_real_)
  } else {
    df <- df %>%
      select(date = 1, time = 2, precip = 3,
             wind_speed = 4, wind_dir = 5,
             temp_air_avg = 6, temp_fuel = 7, rh = 8,
             fuel_moisture = 9, temp_air_max = 10, temp_air_min = 11,
             rh_max = 12, rh_min = 13)
  }

  df %>%
    mutate(
      station     = station_name,
      date        = as.Date(date),
      time        = format(time, "%H:%M"),
      .before     = 1
    ) %>%
    mutate(
      station      = as.character(station),
      time         = as.character(time),
      across(c(precip, wind_speed, wind_dir, temp_air_avg, temp_fuel,
               rh, fuel_moisture, temp_air_max, temp_air_min,
               rh_max, rh_min),
             as.numeric)
    )
}

# ---- read and combine all stations ----

files <- list.files("data/RAWS", pattern = "^[^~].*\\.xlsx$", full.names = TRUE)

all_raws_clean <- purrr::map(files, function(f) {
  message("Reading: ", basename(f))
  read_raws_file(f)
}) %>%
  bind_rows()

# ---- write output ----

readr::write_csv(all_raws_clean, "data/RAWS/RAWS_all_stations.csv")
message("Done. Wrote ", nrow(all_raws_clean), " rows across ",
        n_distinct(all_raws_clean$station), " stations.")

