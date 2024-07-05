## this file turns the movement data used in the COVID analysis into tracks 

library(sf)
#library(amt)
library(tidyverse)
library(move2)

#move1h <- readRDS("./data/Tucker_7704108/Tucker_1h_Move_Spatial.rds") 
move10d <- readRDS("./data/Tucker_7704108/Tucker_10d_Spatial.rds")
#road <- readRDS("./data/Tucker_7704108/Tucker_Road_Spatial.rds")


move10d <- move10d %>% arrange(ID, TimestampUTC)

# plot global 10d displacement
move10d.subset <- move10d %>% filter(Displacement_km <= quantile(Displacement_km, 0.95))
quantile(move10d$Displacement_km)
hist(move10d.subset$Displacement_km, breaks = 100)
abline(v = 0.756, col = "red")
abline(v = 1.756, col = "red")
abline(v = 5.152, col = "red")

hist(log(move10d$Displacement_km), breaks = 100)


# plot mule deer 
move10d.muledeer <- move10d %>% filter(Species == "Odocoileus hemionus", Displacement_km <= quantile(Displacement_km, 0.95)) 
quantile(move10d.muledeer$Displacement_km)
hist(move10d.muledeer$Displacement_km, breaks = 100)
abline(v = 0.590, col = "red")
abline(v = 1.160, col = "red")
abline(v = 2.799, col = "red")

move10d.mongolia <- move10d %>% filter(StudyName == "Mongolian Gazelle" | StudyName == "Mongolian Khulan", Displacement_km <= quantile(Displacement_km, 0.95)) 
hist(move10d.mongolia$Displacement_km)
quantile(move10d.mongolia$Displacement_km)


# make into tracks #
move10d.move2 <- mt_as_move2(move10d, coords = c("Longitude", "Latitude"), 
                          time_column = "TimestampUTC", track_id_column = "ID")
move10d.tr <- move10d.move2 |> sf::st_set_crs(4326) %>% mt_track_lines()
move10d.sf <- st_as_sf(move10d.tr)
st_write(move10d.sf, "./data/Tucker_7704108/Tucker_10d_Tracks.shp")
