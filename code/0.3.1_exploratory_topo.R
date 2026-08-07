#### Some exploratory analyses re: topographic variables


# load Rx outcome polygons

rx.shp <- sf::read_sf("data/from_jean/GWJEFF post-burn forest structure/1burn_cnpy.shp")

# load topographic/terrain variables

topo_stack <- rast("data/topography/topo_stack.tif")

# get a balanced spatial sample

set.seed(7984)

# Number of random points to generate within each severity class
n_per_class <- 100

# Generate random points within polygons for each severity class
balanced_sample <- rx.shp %>%
  nest(data = -Severity_1) %>%
  mutate(
    points = map(data, ~st_sample(.x, size = n_per_class, type = "random") %>%
                   st_as_sf())
  ) %>%
  select(-data) %>%
  unnest(points)

# join topo variables

balanced_sample <- balanced_sample %>%
  st_as_sf() %>% 
  st_join(rx.shp %>% 
            select(burn_name, Burn_Yea_1), 
          join = st_within, left = TRUE) %>% 
  bind_cols(
    terra::extract(topo_stack, vect(.), ID=F)
  )

# quick plot

balanced_sample %>% 
  sf::st_drop_geometry() %>% 
  filter(Severity_1 %in% 0:2) %>% 
  mutate(Severity = case_when(Severity_1 == 0 ~ "Closed",
                              Severity_1 == 1 ~ "Open",
                              Severity_1 == 2 ~ "Early seral")) %>% 
  ggplot(.,
         aes(x = chili)) +
  geom_density(aes(fill = (Severity)),
               alpha = 0.3) +
  scale_fill_manual(name = "Severity",
                    values = c("Closed" = "forestgreen",
                               "Open" = "goldenrod3",
                               "Early seral" = "firebrick2")) +
  labs(x = "CHILI")

balanced_sample %>% 
  sf::st_drop_geometry() %>% 
  filter(Severity_1 %in% 0:2) %>% 
  mutate(Severity = case_when(Severity_1 == 0 ~ "Closed",
                              Severity_1 == 1 ~ "Open",
                              Severity_1 == 2 ~ "Early seral")) %>% 
  ggplot(.,
         aes(x = tpi)) +
  geom_density(aes(fill = (Severity)),
               alpha = 0.3) +
  scale_fill_manual(name = "Severity",
                    values = c("Closed" = "forestgreen",
                               "Open" = "goldenrod3",
                               "Early seral" = "firebrick2")) +
  labs(x = "TPI")

balanced_sample %>% 
  sf::st_drop_geometry() %>% 
  filter(Severity_1 %in% 0:2) %>% 
  mutate(Severity = case_when(Severity_1 == 0 ~ "Closed",
                              Severity_1 == 1 ~ "Open",
                              Severity_1 == 2 ~ "Early seral")) %>% 
  ggplot(.,
         aes(x = slp.dg)) +
  geom_density(aes(fill = (Severity)),
               alpha = 0.3) +
  scale_fill_manual(name = "Severity",
                    values = c("Closed" = "forestgreen",
                               "Open" = "goldenrod3",
                               "Early seral" = "firebrick2")) + 
  labs(x = "Slope angle") 
