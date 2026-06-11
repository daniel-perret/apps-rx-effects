#### Doing some initial climate timeseries and Rx outcome analyses 
#### Goal is to get an exploratory sense of the data and understand relationships between climate and Rx outcomes
#### 
#### 


# load Rx outcome polygons

rx.shp <- sf::read_sf("data/from_jean/GWJEFF post-burn forest structure/1burn_cnpy.shp")
names(rx.shp)
rx.shp$burn_name %>% unique()

raws.burn.trends$Fire_key %>% unique()

# summarize and calculate proportion in each severity class

rx.summary <- rx.shp %>% 
  sf::st_drop_geometry() %>% 
  filter(Severity_1 %in% 0:2) %>% 
  group_by(burn_name) %>% 
  mutate(tot_ac = sum(gisAcres)) %>% 
  group_by(burn_name, Severity_1) %>% 
  summarise(sev.prop = sum(gisAcres)/first(tot_ac))

# link up with RAWS climate trends and data

rx.cnpy.raws <- raws.burn.trends %>% 
  filter(burn_entry == 1) %>% 
  left_join(rx.summary %>% 
              pivot_wider(names_from = Severity_1, 
                          values_from = sev.prop, 
                          names_prefix = "sev_"),
            by = c("fire_name" = "burn_name")) %>% 
  pivot_longer(cols = starts_with("sev_"),
               names_to = "severity_class",
               values_to = "proportion") %>% 
  mutate(severity_class = case_when(
    severity_class == "sev_0" ~ "Closed",
    severity_class == "sev_1" ~ "Open",
    severity_class == "sev_2" ~ "Early Seral",
    TRUE ~ NA_character_)) %>% 
  left_join(raws.burn.ts %>% 
              filter(date == burn.date))

  
# box plots of climate trends by severity class

rx.cnpy.raws %>% 
  filter(!severity_class == "Closed") %>%
  ggplot(.,
         aes(x = temp_trend,
             y = proportion, 
             fill = severity_class)) +
  geom_boxplot()+
  labs(x = "Temperature trend",
       y = "Burn unit proportion") +
  scale_fill_manual(name = "Severity class",
                    values = c("Closed" = "darkgreen",
                               "Open" = "goldenrod",
                               "Early Seral" = "orangered"))

rx.cnpy.raws %>% 
  ggplot(.,
         aes(x = factor(severity_class,
                        levels = c("Closed",
                                   "Open",
                                   "Early Seral")),
             y = proportion, 
             fill = factor(temp_trend,
                           levels = c("decreasing",
                                      "stable",
                                      "increasing"))
             ))+
  geom_boxplot()+
  labs(x = "Canopy outcome",
       y = "Burn unit proportion",
       title = "Week-prior temperature trend") +
  scale_fill_manual(name = "Temperature trend",
                    values = c("decreasing" = "dodgerblue3",
                               "stable" = "goldenrod",
                               "increasing" = "orangered"))
rx.cnpy.raws %>% 
  ggplot(.,
         aes(x = factor(severity_class,
                        levels = c("Closed",
                                   "Open",
                                   "Early Seral")),
             y = proportion, 
             fill = factor(rh_trend,
                           levels = c("decreasing",
                                      "stable",
                                      "increasing"))
             ))+
  geom_boxplot()+
  labs(x = "Canopy outcome",
       y = "Burn unit proportion",
       title = "Week-prior RH trend") +
  scale_fill_manual(name = "RH trend",
                    values = c("decreasing" = "dodgerblue3",
                               "stable" = "goldenrod",
                               "increasing" = "orangered"))
rx.cnpy.raws %>% 
  ggplot(.,
         aes(x = factor(severity_class,
                        levels = c("Closed",
                                   "Open",
                                   "Early Seral")),
             y = proportion, 
             fill = factor(kbdi_trend,
                           levels = c("decreasing",
                                      "stable",
                                      "increasing"))
             ))+
  geom_boxplot()+
  labs(x = "Canopy outcome",
       y = "Burn unit proportion",
       title = "Week-prior KBDI trend") +
  scale_fill_manual(name = "KBDI trend",
                    values = c("decreasing" = "dodgerblue3",
                               "stable" = "goldenrod",
                               "increasing" = "orangered"))


# scatter plots of day of burn weather X severity class ----

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_minRH_aft,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Minimum afternoon RH on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_max_temp,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Max temp on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = days_since_precip,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Days since any precip",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = days_since_precip_01,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Days since precip > 0.01 in",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = days_since_precip_025,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Days since precip > 0.25 in",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_mean_temp,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Mean temp on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_meanRH,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Mean RH on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_meanRH_aft,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Mean afternoon RH on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_minRH,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Minimum RH on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_minRH_lateaft,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Minimum late afternoon RH on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = daily_precip,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "Precipitation on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")

rx.cnpy.raws %>% 
  ggplot(., 
         aes(x = kbdi,
             y = proportion,
             color = factor(severity_class,
                            levels = c("Closed",
                                       "Open",
                                       "Early Seral")))) +
  geom_point(pch = 19, size = 3, alpha = 0.6) +
  labs(x = "KBDI on burn day",
       y = "Burn unit proportion",
       color = "Canopy outcome") +
  scale_color_manual(values = c("Closed" = "forestgreen",
                                "Open" = "goldenrod3",
                                "Early Seral" = "firebrick3")) +
  facet_wrap(facets = ~severity_class,
             nrow = 2) +
  theme(legend.position = "none")
