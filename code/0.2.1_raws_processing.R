#### Standardize WRCC RAWS xlsx files into a single combined CSV
#### Reads 11 station files from data/RAWS/, handles three column schemas,
#### and writes data/RAWS/RAWS_all_stations.csv

source("code/00_setup.R")

# ---- helper function ----

read_raws_file <- function(path) {
  # Extract station name from row 1
  head_rows <- readxl::read_excel(path, n_max = 4, col_names = FALSE,
                                  col_types = "text")

  station_name <- head_rows[[1, 1]] %>%
    sub("\\s+Virginia.*", "", x = .) %>%
    trimws()

  # Legacy schema (Wise): 4 header rows, single combined YYMMDDhhmm
  # datetime column, no separate date/time columns
  is_legacy <- grepl("YYMMDDhhmm", head_rows[[4, 1]] %||% "", fixed = TRUE)

  if (is_legacy) {
    df <- readxl::read_excel(path, skip = 4, col_names = FALSE,
                             col_types = "text")
    colnames(df) <- c("datetime_raw", "precip", "wind_speed", "wind_dir",
                      "temp_air_avg", "temp_fuel", "rh", "fuel_moisture",
                      "temp_air_max", "temp_air_min", "rh_max", "rh_min")
    df <- df %>%
      transmute(
        datetime      = as.POSIXct(datetime_raw, format = "%y%m%d%H%M", tz = "UTC"),
        precip, wind_speed, wind_dir, temp_air_avg, temp_fuel, rh,
        fuel_moisture, temp_air_max, temp_air_min, rh_max, rh_min
      ) %>%
      mutate(date = as.Date(datetime),
             time = format(datetime, "%H:%M")) %>%
      select(-datetime)

    return(
      df %>%
        mutate(
          station = as.character(station_name),
          across(-c(date, time, station), as.numeric),
          time = as.character(time),
          .before = 1
        )
    )
  }

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

readr::write_csv(all_raws_clean, "data/RAWS/RAWS_all_stations_092126.csv")
message("Done. Wrote ", nrow(all_raws_clean), " rows across ",
        n_distinct(all_raws_clean$station), " stations.")

