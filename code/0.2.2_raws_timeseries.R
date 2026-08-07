## ---- generating dataframe of all station/fire combos with all observations for week prior to fire

# load RAWS data

all.raws <- read.csv("data/RAWS/RAWS_all_stations.csv")
names(all.raws)

# load fire data

#unit.hist <- read.csv("data/rx_unit_history.csv", header=T, stringsAsFactors = F)
#unit.hist <- read.csv("data/rx_unit_hist_PAreconciled.csv", header=T, stringsAsFactors = F)
unit.hist <- read.csv("data/Rxb_unit_history_080726.csv", header=T, stringsAsFactors = F)

unit.hist <- unit.hist %>% 
  mutate(keep. = trimws(keep.),
         burn_entry = trimws(burn_entry)) %>% 
  # only retain rows the source data marks as usable
  filter(keep. == "yes",
         # multiday entries only have a real date for the first burn_entry;
         # later entries currently carry a placeholder date_fixed
         burn_entry == "1") %>% 
  mutate(date.clean = as.Date(date_fixed,
                              format = "%m/%d/%Y"),
         year = year(date.clean),
         month = month(date.clean),
         day = day(date.clean),
         RAWS_station = toupper(RAWS_station),
         Fire_name = Fire_name_fixed) %>% 
  filter(!is.na(date.clean),
         year>=RAWS1styr)

rawsMeta <- read.csv("data/RAWS/fw13/RAWSfw13list.csv") %>% 
  mutate(Name = sub("\\s+$", "", Name))

unit.hist <- unit.hist %>% 
  left_join(rawsMeta %>% 
              rename(RAWS_lat = LatDegrees,
                     RAWS_lon = LonDegrees,
                     RAWS_elev = Elevation),
            by = c("RAWS_station" = "Name"))


# station-level daily timeseries (full history, needed for drought metrics below)
raws.daily <- all.raws %>% 
  mutate(station = toupper(station),
         station = ifelse(station=="UPPER TRACT  WEST",
                          "UPPER TRACT",
                          station),
         date = as.Date(date)) %>% 
  group_by(station) %>% 
  mutate(hr.precip = precip - lag(precip),
         hr.precip = ifelse(hr.precip<0,0,hr.precip)) %>% 
  filter(!is.na(date)) %>% 
  # calculate days since last precipitation
  group_by(station, date) %>%
  mutate(daily_precip = sum(hr.precip, na.rm=T)) %>%
  group_by(station) %>%
  arrange(station, date, time) %>%
  mutate(last_precip_date = if_else(daily_precip > 0, date, NA_Date_)) %>%
  fill(last_precip_date, .direction = "down") %>%
  mutate(days_since_precip = as.integer(date - last_precip_date)) %>%
  # days since 0.1" precip
  mutate(last_precip_date_01 = if_else(daily_precip > 0.1, date, NA_Date_)) %>%
  fill(last_precip_date_01, .direction = "down") %>%
  mutate(days_since_precip_01 = as.integer(date - last_precip_date_01)) %>%
  
  # days since 0.25" precip
  mutate(last_precip_date_025 = if_else(daily_precip > 0.25, date, NA_Date_)) %>%
  fill(last_precip_date_025, .direction = "down") %>%
  mutate(days_since_precip_025 = as.integer(date - last_precip_date_025)) %>%
  
  ungroup()

# long-term mean annual precip per station, computed once from the full
# station history (independent of any fire's burn window)
station_map <- raws.daily %>% 
  distinct(station, date, daily_precip) %>% 
  mutate(yr = year(date)) %>% 
  group_by(station, yr) %>% 
  summarise(annual.precip = sum(daily_precip, na.rm=T), .groups = "drop_last") %>% 
  summarise(ts_MAP = mean(annual.precip, na.rm = T), .groups = "drop")

# join unit data, then restrict to each fire's burn window BEFORE summarizing
# so daily stats (and long-term MAP join below) aren't computed/warned over
# irrelevant dates outside the window of interest
raws.burn.ts <- raws.daily %>% 
  full_join(unit.hist %>% 
              mutate(Fire_key = paste0(Fire_name,"_",burn_entry)) %>% 
              select(Fire_name,
                     Fire_key, 
                     burn.date = date.clean,
                     burn_entry,
                     RAWS_station),
            relationship = "many-to-many",
            by = c("station" = "RAWS_station")) %>%
  filter(!is.na(Fire_name),
         date > burn.date - 7,
         date <= burn.date) %>% 
  # summarizing to daily values
  group_by(Fire_key, date) %>% 
  summarise(station = first(station),
            daily_mean_temp = mean(temp_air_avg, na.rm=T),
            daily_max_temp = max(temp_air_avg, na.rm=T),
            daily_mean_temp_aft = mean(temp_air_avg[time>="12:00" & time <= "18:00"], na.rm=T),
            daily_meanRH = mean(rh, na.rm=T),
            daily_minRH = min(rh, na.rm=T),
            daily_meanRH_aft = mean(rh[time>="12:00" & time <= "18:00"], na.rm=T),
            daily_minRH_aft = min(rh[time>="12:00" & time <= "18:00"], na.rm=T),
            daily_meanRH_lateaft = mean(rh[time>="15:00" & time <= "18:00"], na.rm=T),
            daily_minRH_lateaft = min(rh[time>="15:00" & time <= "18:00"], na.rm=T),
            daily_precip = first(daily_precip),
            days_since_precip = first(days_since_precip),
            days_since_precip_01 = first(days_since_precip_01),
            days_since_precip_025 = first(days_since_precip_025),
            burn.date = first(burn.date),
            burn_entry = first(burn_entry),
            fire_name = first(Fire_name),
            .groups = "drop") %>% 
  # attach pre-computed long-term MAP for each station
  left_join(station_map, by = "station")
  # kbdi = calc_kbdi_vector(temp_vec = daily_max_temp,
  #                         precip_vec = daily_precip,
  #                         map_vec = ts_MAP)

# join with new daily fire weather vars

fire.vars <- read.csv("data/RAWS/fire_danger_vars.csv", header = T, stringsAsFactors = F) %>% 
  filter(!is.na(KBDI)) %>% 
  mutate(date = as.Date(ObservationTime))

raws.burn.ts <- raws.burn.ts %>% 
  left_join(., fire.vars,
            by = c("station" = "StationName",
                   "date"))
# write and load

write.csv(raws.burn.ts, "data/raws_burn_ts_080726.csv", row.names = F)

raws.burn.ts <- read.csv("data/raws_burn_ts_080726.csv", header = T, stringsAsFactors = F)


# calculate trends

# classify a variable's within-window trajectory over the burn window
trend_class <- function(x, date) {
  # cor() errors (rather than warns) when fewer than 2 complete pairs remain,
  # e.g. fire.vars only covers 2005+, so pre-2005 fires have all-NA fire-weather vars
  if (sum(complete.cases(x, date)) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(NA_character_)
  }
  r <- suppressWarnings(cor(as.integer(date), x, use = "complete.obs"))
  case_when(
    is.na(r) ~ NA_character_,
    r >  0.5 ~ "increasing",
    r < -0.5 ~ "decreasing",
    TRUE ~ "stable"
  )
}

# fire-weather vars (from fire.vars) to compute trends for
fwvars_for_trend <- c("X1HrFM", "X10HrFM", "X100HrFM", "X1000HrFM",
                      "KBDI", "WoodyFM", "HerbFM", "ERC")

raws.burn.trends <- raws.burn.ts %>% 
  mutate(date=as.Date(date)) %>% 
  # calculate trends
  group_by(Fire_key) %>%
  summarise(temp_trend = trend_class(daily_mean_temp, date),
    rh_trend = trend_class(daily_meanRH, date),
    across(all_of(fwvars_for_trend),
           ~trend_class(.x, date),
           .names = "{sub('^X', '', .col)}_trend"),
    burn.date = first(burn.date),
    fire_name = first(fire_name),
    burn_entry = first(burn_entry))


# final fire weather dat

final.weather <- raws.burn.ts %>% 
  filter(date == burn.date) %>% 
  left_join(raws.burn.trends)

write.csv(final.weather,
          "data/fireweather_summary_Rx1b.csv",
          row.names = F)
