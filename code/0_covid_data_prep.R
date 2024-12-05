## this file generate 1d data from the 1h data. 
### also add HMI info and incorperate the species information into 10d and 1d covid data for later analysis. 

library(sf)
library(tidyverse)
library(geosphere)

# move1h <- readRDS("./data/movement/Tucker_7704108/Tucker_1h_Move_Spatial.rds")  %>%
#   dplyr::select(-"StringencyIndex", - "Lockdown(Y/N)", -"Country", - "RBrain_Size", -"Activity")

move10d <- readRDS("./data/movement/Tucker_7704108/Tucker_10d_Spatial.rds") %>%
  dplyr::select(-"StringencyIndex", - "Lockdown(Y/N)", -"Country", - "RBrain_Size", -"Activity")

move = move10d
target_day_scale = 10 # <- or 1

# ---------------- include species data ------------------- #
# pantheria data base
pant <- read_csv("./data/PanTHERIA_1-0_WR05_Aug2008.csv")
pant <- pant[,1:5] %>% 
  rename (Order = MSW05_Order,
          Family = MSW05_Family,
          Genus = MSW05_Genus,
          Species = MSW05_Species,
          Binomial = MSW05_Binomial) 

# find out which species are not listed in panethria database 
move  %>% filter(!Species %in% pant$Binomial) %>% dplyr::select(Species) %>% distinct() # elk and desert bighorn sheep is not included 
# add elk, cape bushbuck 
pant = rbind(pant,
             data.frame(Order = c("Artiodactyla", "Artiodactyla"), 
                        Family = c("Cervidae", "Bovidae"), 
                        Genus = c("Cervus", "Tragelaphus"), 
                        Species = c("canadensis", "sylvaticus"), 
                        Binomial = c("Cervus canadensis", "Tragelaphus sylvaticus")))

# treat Ovis canadensis nelsoni the same as Ovis canadensis
move <- move %>% mutate(Species = case_when(Species == "Ovis canadensis nelsoni" ~ "Ovis canadensis",
                                                            Species == "Ovis canadensis canadensis" ~ "Ovis canadensis",
                                                            Species == "Ovis canadensis californiana" ~ "Ovis canadensis",
                                                            Species == "Cervus_canadensis" ~ "Cervus canadensis",
                                                          .default = Species))

move <- move %>% rename (Binomial = Species ) %>% left_join(pant, by = "Binomial")
# move  %>% filter(!Binomial %in% pant$Binomial) %>% select(Binomial) %>% distinct() # elk and desert bighorn sheep is not included 
length(unique(move $Binomial)) # 21 species left for 10d, 15 for 1h

# ####################################
# # ----- turn 1h data to 1d----------
# ####################################
# move1d <- move %>%
#   dplyr::select(-Displacement_km) %>%
#   mutate(date = as.Date(TimestampUTC)) %>%
#   group_by(ID, date) %>%
#   summarize(
#     n = n(), # should be 24 to have full dataset 
#     # Get the first and last point of the day within the valid period
#     first_point = cbind(first(Longitude), first(Latitude)),
#     last_point = cbind(last(Longitude), last(Latitude)),
#     first_timestamp = first(TimestampUTC),  # Capture the first timestamp for reference
#     last_timestamp = last(TimestampUTC),    # Capture the last timestamp for reference
#     
#     # # Check if the first and last points match the intended time period boundaries
#     # is_first_point_start = hour(first_timestamp) %in% c(0,1),
#     # is_last_point_end = hour(last_timestamp) %in% c(22,23),
#     
#     # Calculate 1-day displacement (in meters)
#     one_day_displacement_m = distHaversine(first_point, last_point),
#     
#     # Convert to kilometers
#     one_day_displacement_km = one_day_displacement_m / 1000
#     
#   ) %>% # this calculates at least time stamp t, the movement distance prior to getting to this time stamp. i.e. the first row has NA displacement.
#   filter(n >= 0.75 * 24) %>%
#   dplyr::select(ID, first_timestamp, one_day_displacement_km) %>%
#   rename(TimestampUTC = first_timestamp,
#          Displacement_km = one_day_displacement_km) %>%
#   left_join(move %>% dplyr::select(-Displacement_km), 
#             by = c("ID", "TimestampUTC"))  # 493348 -> 20440
# 
# move <- move1d
####################################
# ----- get HMI data ----------
####################################
library(terra)
HMI <- rast("./data/human_modification_index/gHM.tif")
HMI <- project(HMI, "EPSG: 4326")

move_reshape <- move %>%
  group_by(ID) %>%
  arrange(TimestampUTC) %>%
  mutate(time_diff = as.numeric(difftime(lead(TimestampUTC), TimestampUTC, units = "days"))) %>%
   filter(time_diff < target_day_scale*1.25 & time_diff > target_day_scale*0.75) %>% # 14,753  --> 11,523 for 10d, 20,440 --> 18,071 for 1 d
  # get start and end point for averaging HMI values
  mutate(
    ID,
    start_long = Longitude,
    start_lat = Latitude,
    end_long = lead(Longitude),
    end_lat = lead(Latitude),
    start_time = TimestampUTC,
    end_time = lead(TimestampUTC),
    .keep = "none"
  ) %>%
  ungroup() %>%
  drop_na(end_long)  # 11523 -> 9999 for 10 day; 18,071 -> 17,881 for 1 day

# get midpoint for percolation distance extraction
move_reshape <- move_reshape %>% mutate(
    first_point = matrix(c(start_long, start_lat), nrow = nrow(move_reshape)),
    last_point = matrix(c(end_long, end_lat), nrow = nrow(move_reshape)),
    mid_point = midPoint(first_point, last_point),  # get middle point
    midLongitude = mid_point[,1], midLatitude = mid_point[,2]
  ) 
  
  
move.start.sf <-  move_reshape %>%
  st_as_sf(., coords = c("start_long", "start_lat"), crs = st_crs(4326)) 
move.end.sf <-  move_reshape %>%
  st_as_sf(., coords = c("end_long", "end_lat"), crs = st_crs(4326)) 

move.start.sf$HMI.start <- terra::extract(HMI, move.start.sf)[,-1]
move.end.sf$HMI.end <- terra::extract(HMI, move.end.sf)[,-1]

HMI.df <- tibble(
  ID = move_reshape$ID,
  TimestampUTC = move_reshape$start_time,
  Longitude = move_reshape$start_long,
  Latitude = move_reshape$start_lat,
  HMI1 = move.start.sf$HMI.start, 
  HMI2 = move.end.sf$HMI.end,
  midLongitude = move_reshape$midLongitude, 
  midLatitude = move_reshape$midLatitude) %>%
  mutate(HMI = (HMI1+HMI2)/2) %>%
  dplyr::select(-HMI1, -HMI2)

move <- move %>% left_join(HMI.df) # hmm lost almost 5000 row of data for this. no need to lose them if marlee can share raw data
#write_rds(move, "./data/movement/move1d_covid.rds")
# write_rds(move, "./data/movement/move10d_covid.rds")
