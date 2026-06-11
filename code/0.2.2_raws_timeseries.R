## ---- generating dataframe of all station/fire combos with all observations for week prior to fire

# load RAWS data

all.raws <- read.csv("data/RAWS/RAWS_all_stations.csv")
names(all.raws)

# load fire data

unit.hist <- read.csv("data/rx_unit_history.csv", header=T, stringsAsFactors = F)

unit.hist <- unit.hist %>% 
  mutate(date.clean = as.Date(date,
                              format = "%m/%d/%Y"),
         year = year(date.clean),
         month = month(date.clean),
         day = day(date.clean),
         RAWS_station = toupper(RAWS_station)) %>% 
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


# clean and calculate timeseries
raws.burn.ts <- all.raws %>% 
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
  
  ungroup() %>%
  # join unit data
  full_join(unit.hist %>% 
              mutate(Fire_key = paste0(Fire_name,"_",burn_entry)) %>% 
              select(Fire_name,
                     Fire_key, 
                     burn.date = date.clean,
                     burn_entry,
                     RAWS_station),
            relationship = "many-to-many",
            by = c("station" = "RAWS_station")) %>%
  filter(!is.na(Fire_name)) %>% 
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
            fire_name = first(Fire_name)) %>% 
  # calculate long-term MAP
  group_by(Fire_key, year(date)) %>% 
  mutate(annual.precip = sum(daily_precip, na.rm=T)) %>% 
  group_by(Fire_key) %>% 
  mutate(ts_MAP = mean(annual.precip,na.rm=T),
         kbdi = calc_kbdi_vector(temp_vec = daily_max_temp,
                                 precip_vec = daily_precip,
                                 map_vec = ts_MAP)) %>% 
  
  #filter out to relevant dates
  filter(date > burn.date-7,
         date <= burn.date)

raws.burn.ts <- read.csv("data/raws_burn_ts.csv", header = T, stringsAsFactors = F)

raws.burn.trends <- raws.burn.ts %>% 
  mutate(date=as.Date(date)) %>% 
  # calculate trends
  group_by(Fire_key) %>%
  summarise(temp_trend = case_when(
    cor(as.integer(date), daily_mean_temp, use = "complete.obs") >  0.5 ~ "increasing",
    cor(as.integer(date), daily_mean_temp, use = "complete.obs") < -0.5 ~ "decreasing",
    TRUE ~ "stable"),
    rh_trend = case_when(
      cor(as.integer(date), daily_meanRH, use = "complete.obs") >  0.5 ~ "increasing",
      cor(as.integer(date), daily_meanRH, use = "complete.obs") < -0.5 ~ "decreasing",
      TRUE ~ "stable"),
    kbdi_trend = case_when(
      cor(as.integer(date), kbdi, use = "complete.obs") >  0.5 ~ "increasing",
      cor(as.integer(date), kbdi, use = "complete.obs") < -0.5 ~ "decreasing",
      TRUE ~ "stable"),
    burn.date = first(burn.date),
    fire_name = first(fire_name),
    burn_entry = first(burn_entry))
