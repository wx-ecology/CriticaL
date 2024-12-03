# this script synthesize all movement data. 

# first, download move bank data, organized them into 10 day displacement, and combine them with Marlee's data
# second, add other data found from the web (Morato et al 2018 Jaguar)
# third, data collected from the barrier group emailing 

library(move2)
library(tidyverse)
library(geosphere)

setwd("/Users/wenjing/Senckenberg Dropbox/Wenjing Xu/CriticaL/")
movebank_store_credentials("wenjingxu", "Move_1877!0")

################################
##  download movebank data -----
################################

## 1. get all open study info ## 
# open_studies <- read_csv("./data/movement/open_movebank/filtered_study_id.csv") #study ID of open mammal data from movebank
# # check out the studies, pick ones to include 
# library(taxize)
# scientific_names <- open_studies$Binomial
# # Get id from databases
# uids <- get_uid(scientific_names)
# # keep only uids which you have in the database
# uids <- as.uid(uids)
# # get common names 
# Common_names <- sci2comm(uids, db = 'ncbi')
# 
# open_studies$Common_names <- unlist(lapply(Common_names, function(x) {
#   if (length(x) == 0) x = NA else x
# }))
##  write_csv(open_studies, "./data/movement/open_movebank/filtered_study_id.csv")
#
#open_studies <- read_csv("./data/movement/open_movebank/filtered_study_id.csv") %>% filter(select == 'y') #study ID of open mammal data from movebank
#
## 2. download and compile ## 
# each study might ned to mannually click through to accept license
#
#df <- tibble()
# for (id in unique(open_studies$study_id)) {
#   mv.i <- movebank_download_study(id, sensor_type_id = "gps", attributes = NULL, 
#                                   remove_movebank_outliers = T,'license-md5'='2f3d932ad96f9b2f49096f7ecf5e2559')
#   mv.i <- mt_as_event_attribute(mv.i, c(deployment_id, individual_id)) %>%
#     dplyr::mutate(coords_x=sf::st_coordinates(.)[,1], coords_y=sf::st_coordinates(.)[,2], crs = sf::st_crs(mv.i)$input) %>%
#     sf::st_drop_geometry(.)
#   df.i <- tibble(study_id = id, 
#                  study_name = movebank_download_study_info(study_id = id)$name,
#                  mv.i, 
#                  contact = movebank_download_study_info(study_id = id)$contact_person_name) %>%
#     left_join(tibble(deployment_id = mt_track_data(mv.i)$deployment_id,
#                      species =  mt_track_data(mv.i)$taxon_canonical_name), by = "deployment_id")
#   df <- rbind(df,df.i)
# }
# 
## 3. address special cases ## 
# # study 7023252 has 1 sec data. cut that down first to save space 
# # some individuals are tracked for over 1 month, some for less than a day. 
# study_7023252_summary <- df %>% filter(study_id == 7023252)  %>%
#   mutate(timestamp = ymd_hms(timestamp)) %>%
#   group_by(individual_id) %>%
#   summarise(duration = max(timestamp) - min(timestamp)) %>% 
#   filter(duration >= 60*60*24*10) # only keep individuals with at least 10 day data
# study_7023252  <- df %>% filter(study_id == 7023252, individual_id %in% study_7023252_summary$individual_id) %>%
#   mutate(month = month(timestamp), day = day(timestamp), hour = hour(timestamp), minute = minute(timestamp)) %>%
#   filter(minute == 0) %>%
#   group_by(individual_id, month, day, hour) %>%
#   mutate(coords_x = mean(coords_x), coords_y = mean(coords_y), timestamp = first(timestamp)) %>%
#   distinct() # turn the study into hourly coordinates
# df <- rbind((df %>% filter(study_id != 7023252)), 
#             (study_7023252 %>% ungroup() %>% dplyr::select(-month, -day, -hour, -minute)))
# 
# # check whether location data are in lat lon. If not, delete study. 
# 
# unique(df$species[which(!df$species %in% open_studies$Binomial)]) # it is all birds
# df <- df %>% filter(species %in% open_studies$Binomial)

# write_csv(df, "/Users/wenjing/Senckenberg Dropbox/Wenjing Xu/CriticaL/data/movement/open_movebank/midproduct_raw_movebank_data.csv")

###########################################################
#### data checks and cleaning/filtering  ------------------
###########################################################

df <- read_csv("./data/movement/open_movebank/midproduct_raw_movebank_data.csv") # n = 3902306
# df.raw = df

# target_time_scale_days = 1 # set up the target time scale
target_time_scale_days = 10 # set up the target time scale

# --- filter NA locations -- #
df <- df %>% filter(!is.na(coords_x), !is.na(coords_y)) # n = 3826440

# --- filter based on data duration -- #
df.summary <- df %>%
  group_by(study_id, individual_id, deployment_id) %>%
  summarise( temp_interval = median(diff(timestamp)),
             duration = max(timestamp) - min(timestamp)) %>%
  mutate( temp_interval_hours = as.numeric(temp_interval)/3600,
          duration_days = as.numeric(duration)/(3600*24))
df.keep <- df.summary %>% filter(duration_days >= target_time_scale_days)  # at least tracked for 10 days

# reduce dataset to individuals monitored for at least xx days at once
df.filtered <- df %>%
  filter(study_id %in% df.keep$study_id,
         individual_id %in% df.keep$individual_id,
         deployment_id %in% df.keep$deployment_id)

# delete study not in lat lon 
df.filtered <- df.filtered %>% filter(study_id != 4732684199 ) # coordinates of study id not in lat/lon; n = 3792215 points left

# generate remaining data info df, create unique ID for each individual (4 digit ID + species)
df.filtered.info <- df.filtered %>% ungroup() %>%
  dplyr::select(study_id, individual_id, species, study_name, contact) %>% distinct() %>% 
  mutate(
    individual_id = as.character(individual_id),
    ID = paste0(str_sub(individual_id, -4,-1), "_", species)
) 
nrow(df.filtered.info) #n_ind = 473
any(duplicated(df.filtered.info$ID))

###############################################################################
####  filter data into 10/1 day intervals and calculate disp -----------------   
###############################################################################
TenDayDisp_df <- df.filtered %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(study_id, individual_id, deployment_id) %>%
  arrange(timestamp) %>%
  mutate(day_group = (as.numeric(date) - min(as.numeric(date))) %/% target_time_scale_days) %>%
  group_by(study_id, individual_id, deployment_id, day_group) %>%
  summarize(
    # Get the first and last point of each 10-day period within the valid period
    first_point = cbind(first(coords_x), first(coords_y)),
    last_point = cbind(last(coords_x), last(coords_y)),
    first_timestamp = first(timestamp),  # Capture the first timestamp for reference
    last_timestamp = last(timestamp),    # Capture the last timestamp for reference
    
    # Check if the first point and the last point are 10 days apart
    timediffs = difftime(last_timestamp, first_timestamp, units = "hours"))  # n = 21218

# delete data if not enough recorded locations in each time segment
TenDayDisp_df  <- TenDayDisp_df %>% 
  filter(as.numeric(timediffs) > 24*0.9*target_time_scale_days) %>% # give it a 10% tolerance # n = 18081
  mutate(ten_day_displacement_km = distHaversine(first_point, last_point) / 1000) # calculate displacement in km 

# calculate middle point for percolation distance extraction 
TenDayDisp_df <- TenDayDisp_df %>% ungroup() %>%
  mutate(mid_point = midPoint(first_point, last_point),  # get middle point
         midLongitude = mid_point[,1], midLatitude = mid_point[,2]) %>% # only  # get middle position for extracting percolation distance 
  rename(Displacement_km = ten_day_displacement_km) %>%
  dplyr:: select(-deployment_id, -day_group, -timediffs)

move10d_movebank <- TenDayDisp_df %>% mutate(individual_id = as.character(individual_id)) %>%
  left_join(df.filtered.info, by = c('study_id', 'individual_id')) %>%
  dplyr::select(-individual_id) # n = 18,081

write_rds(move10d_movebank, "./data/movement/open_movebank/midproduct_move10d_movebank.rds")

###############################################################################
####  Environmental covariates extraction  ------------------------------------   
###############################################################################
#### next need to extract NDVI to each start and end point and take an average for each displacement segment
library(rgee)
library(sf)
library(geojsonio)
ee_Initialize()
ee_Authenticate() 

move10d_movebank <- read_rds("./data/movement/open_movebank/midproduct_move10d_movebank.rds") %>% mutate(unique_serial = row_number(.))

env.sf <- move10d_movebank %>% 
  dplyr::select(unique_serial, first_point, last_point, first_timestamp, last_timestamp) %>%
  pivot_longer(cols = c(first_point, last_point), names_to = "point", values_to = "coordinates" ) %>%
  mutate(Longitude = coordinates[,1], Latitude = coordinates[,2],
         TimestampUTC = ifelse (point == "first_point", first_timestamp, last_timestamp) ) %>%
  dplyr::select(-coordinates, -first_timestamp, -last_timestamp) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = st_crs(4326))

# Load MODIS NDVI collection
modis_ndvi <- ee$ImageCollection("MODIS/006/MOD13Q1")$select("NDVI")

# extract image dates from the image collection to later combine with points
modis_dates = ee_get_date_ic(modis_ndvi)$time_start

# match modis date with point date
env.sf <- env.sf %>% 
  mutate(modis_date = as.character(modis_dates[findInterval(TimestampUTC, modis_dates)]))

# Loop through each group and extract NDVI
ndvi_list <- tibble()
for (i in as.character(unique(env.sf $modis_date))[1:length(unique(env.sf $modis_date))] ) {
  
  points_group = env.sf %>% filter(modis_date == i) 
  
  filtered_modis <- modis_ndvi$
    filterDate(as.character(i))$ 
    median()  # reduce from IC to image 
  
  ndvi_values <- ee_extract(
    x = filtered_modis,
    y = points_group["unique_serial"],
    scale = 250,  # MODIS resolution
    fun = ee$Reducer$mean()    # Aggregation function (not critical for points)
  )
  
  ndvi_list <- rbind(ndvi_list, ndvi_values)
  
} # run time around 7 min
move10d_movebank <- move10d_movebank %>% 
  left_join(ndvi_list %>% group_by(unique_serial) %>% 
              summarise (NDVI = mean(NDVI)), # all environmental covariates for each segment are the mean of the start and the end point values
            by = "unique_serial")


# extract HFI and HMI ------------
library(terra)

HFI <- rast("./data/Human_footprint_index/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif")
HFI <- project(HFI, "EPSG: 4326")
HMI <- rast("./data/human_modification_index/gHM.tif")
HMI <- project(HMI, "EPSG: 4326")

env.sf$HFI <- terra::extract(HFI, env.sf)[,-1]
env.sf$HMI <- terra::extract(HMI, env.sf)[,-1]

# merge NDVI, HFI, and HMI back with movebank data
move10d_movebank <- move10d_movebank %>% 
  left_join(env.sf %>% st_drop_geometry(.) %>% 
              group_by(unique_serial) %>% 
              summarise (HFI = mean(HFI),
                      HMI = mean(HMI)), # all environmental covariates for each segment are the mean of the start and the end point values
            by = "unique_serial")

# write_rds(move10d_movebank, 
#           "./data/movement/open_movebank/midproduct_movebank_10d_disp_w_envinfo.rds")

# extract species, Diet, BodyMass_kg info -----------------   
diet <- read_csv("./data/EltonTraits/MamFuncDat.csv") %>% filter(Scientific %in% unique(move10d_movebank.sf$species))
sum(unique(move10d_movebank.sf$species) %in% unique(diet$Scientific))
unique(move10d_movebank.sf$species) [which(!unique(move10d_movebank.sf$species) %in% unique(diet$Scientific))]
# "Sapajus macrocephalus" is a subspecies of ""Sapajus apella". but both not in the database.

diet <- diet %>% mutate(Diet = 
                          case_when(
                            `Diet-Inv` + `Diet-Vend` + `Diet-Vect` + `Diet-Vfish` + `Diet-Vunk` + `Diet-Scav` == 0  ~ "Herbivore",
                            `Diet-Fruit` + `Diet-Nect` + `Diet-Seed` + `Diet-PlantO` == 0  ~ "Carnivore",
                            .default = "Omnivore"
                          ),
                        BodyMass_kg = `BodyMass-Value`/1000) %>%
  dplyr::select(Scientific, Diet, BodyMass_kg) %>%
  rename(species = Scientific)

move10d_movebank.sf <- move10d_movebank.sf %>% left_join(diet, by = "species")

# extract taxonomy info from PANTHERIA -----------------   
pant <- read_csv("./data/PanTHERIA_1-0_WR05_Aug2008.csv")
pant <- pant[,1:5] %>% 
  rename (Order = MSW05_Order,
          Family = MSW05_Family,
          Genus = MSW05_Genus,
          Species = MSW05_Species,
          Binomial = MSW05_Binomial) 

# find out which species are not listed in panethria database 
move10d_movebank.sf  %>% filter(!species %in% pant$Binomial) %>% dplyr::select(species) %>% distinct() 
#Sapajus macrocephalus is also not included in the pantheria databse. 

# add Sapajus apella (of which Sapajus macrocephalus is a subspecies of)
pant = rbind(pant,
             data.frame(Order = c("Primates"), 
                        Family = c("Cebidae"), 
                        Genus = c("Sapajus"), 
                        Species = c("apella macrocephalus"), 
                        Binomial = c("Sapajus apella macrocephalus")))

move10d_movebank.sf  <- move10d_movebank.sf  %>% rename (Binomial = species ) %>% left_join(pant, by = "Binomial")
# length(unique(move10d_movebank.sf $Binomial)) # 22 species

move10d_movebank <- move10d_movebank.sf %>% 
  st_drop_geometry() %>%                    
  cbind(st_coordinates(move10d_movebank.sf)) %>%
  rename(Longitude = X, Latitude = Y) %>% dplyr::select(-modis_date, -unique_serial)

write_rds(move10d_movebank, "./data/movement/open_movebank/move10d_movebank.rds")
