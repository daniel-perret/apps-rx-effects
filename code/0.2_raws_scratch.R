#### This script begins assembling RAWS data from stations proximate to prescribed burns on GWJ NF
#### The goal is to grab RAWS data from the relevant stations, and calculate temp and RH trends in the week prior to the burn event


## ---- read in Jean's burn data

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

## ---- exploring RAWS data 

rawsMeta <- read.csv("data/RAWS/RAWSfw13list.csv") %>% 
  mutate(Name = sub("\\s+$", "", Name))

unit.hist <- unit.hist %>% 
  left_join(rawsMeta %>% 
              rename(RAWS_lat = LatDegrees,
                     RAWS_lon = LonDegrees,
                     RAWS_elev = Elevation),
            by = c("RAWS_station" = "Name"))

rawslist <- list.files("data/RAWS",pattern = ".xlsx", full.names = T)

all.raws <- map_df(rawslist, function(x){readxl::read_excel(x)})


## ---- generating dataframe of all station/fire combos with all observations for week prior to fire
clean.raws <- all.raws %>% 
  mutate(date.clean = as.Date(observationDate,
                              format = "%Y%m%d")) %>% 
  select(nwsID, date.clean, observationTime, 
         airTemp = dryBulbTemp,
         RH = atmosMoisture) %>%
  full_join(unit.hist %>% 
              mutate(Fire_key = paste0(Fire_name,"_",burn_entry)) %>% 
              select(Fire_name,
                     Fire_key, 
                     burn.date = date.clean,
                     StationID) %>% 
              mutate(StationID = as.character(StationID)),
            relationship = "many-to-many",
            by = c("nwsID" = "StationID")) %>% filter(Fire_name == "Catback Mountain") %>% view()
  filter(date.clean > burn.date-7,
         date.clean <= burn.date) %>%
  # summarizing to a couple daily values to look at
  group_by(Fire_key, date.clean) %>% 
  mutate(meanTemp = mean(airTemp, na.rm=T),
         meanRH = mean(RH, na.rm=T),
         minRH = min(RH, na.rm=T)) %>%
  filter(observationTime == 1400) %>% 
  select(observationTime,
         Temp_1400 = airTemp,
         RH_1400 = RH,
         meanTemp,
         meanRH, 
         minRH,
         nwsID,
         burn.date) %>% 
  ungroup()
