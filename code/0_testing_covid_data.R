library(sf)
library(amt)
library(tidyverse)

#move1h <- readRDS("./data/Tucker_7704108/Tucker_1h_Move_Spatial.rds") 
move10d <- readRDS("./data/Tucker_7704108/Tucker_10d_Spatial.rds")
#road <- readRDS("./data/Tucker_7704108/Tucker_Road_Spatial.rds")

move10d <- move10d %>% arrange(ID, TimestampUTC)

move10d.tr <- mt_as_move2(move10d, coords = c("Longitude", "Latitude"), 
                          time_column = "TimestampUTC", track_id_column = "ID") |> sf::st_set_crs(4326) %>% mt_track_lines()
ggplot() +   geom_sf(data = move10d.tr[1:5,], aes(color = `ID`)) + theme_minimal()

move10d.sf <- st_as_sf(move10d.tr)

st_write(move10d.sf, "./data/Tucker_7704108/Tucker_10d_Tracks.shp")
