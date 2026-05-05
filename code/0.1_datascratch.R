rx1 <- sf::read_sf("data/from_jean/GWJEFF_Rxb_history/GWJEFF_Rxb_history.shp") %>% 
  sf::st_make_valid()

rx1 %>% 
  sf::st_drop_geometry() %>% 
  group_by(Severity_1) %>% 
  summarise(area_total = sum(gisAcres)) %>% 
  ungroup() %>% 
  mutate(prop = area_total/sum(area_total)) %>% view()

rx1 %>% 
  sf::st_drop_geometry() %>% 
  filter(Burn_Yea_1<8888) %>% 
  count(Burn_Yea_1) %>%
  ggplot(.,
         aes(x = Burn_Yea_1,
             y = n)) +
  geom_col()


rx1 %>% 
  sf::st_drop_geometry() %>% 
  filter(Burn_Yea_1<9999,
         Severity_1 %in% 0:2) %>%
  ggplot(.,
         aes(x = Burn_Yea_1,
             y = gisAcres))+
  geom_col(aes(fill=factor(Severity_1)))


ecosections <- sf::read_sf("../../SHARED_DATA/base_spatialdata/cleland_usfs_ecoregions/S_USA.EcomapSections.shp") %>% 
  sf::st_transform(st_crs(rx1)) %>% 
  sf::st_make_valid()

aoi <- sf::st_filter(ecosections, rx1, .predicate = st_intersects) %>% 
  sf::st_union()

sf::write_sf(aoi, "data/aoi.shp")

aoi <- sf::read_sf("data/aoi.shp")
plot(aoi,fill="red")


bps <- terra::rast("../../SHARED_DATA/")

tm <- terra::rast("../../SHARED_DATA/TREEMAP/RDS-2025-0031/Data/TreeMap2020_CONUS.tif")

aoi2 <- sf::st_transform(aoi, crs(tm))

tm <- tm %>% 
  terra::crop(aoi2,mask=T)

tm.rat <- terra::cats(tm) %>% 
  as.data.frame()

terra::writeRaster(tm, "data/TreeMap2020_aoi.tif")

tm.can <- tm %>% 
  terra::classify(., 
                  rcl = tm.rat %>% 
                    select(Value, CANOPYPCT) %>% 
                    as.matrix())

activeCat(tm) <- 1
#try this it might be faster
tm.can <- tm %>% 
  terra::subst(., from = tm.rat$Value, to = tm.rat$CANOPYPCT)






