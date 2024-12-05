# this file organize other data from direct data sharing or downloading from open deposit other than movebank


library(tidyverse)
library(geosphere)

df <- read_csv("./data/movement/other_open_data/Morato_2018_Jaguar/jaguar_movement_data.csv")

df <- df %>% rename(StudyName = study.name, 
              ID = 'individual.local.identifier (ID)',
              TimestampUTC = timestamp,
              Longitude = location.long,
              Latitude = location.lat,
              Binomial = individual.taxon.canonical.name) %>%
  mutate(TimestampUTC = mdy_hm(TimestampUTC),
         ID = paste0(str_sub(StudyName, 1,2), ID, "_", Binomial),
         ContactPerson = "Ronaldo G Morato") %>%
  filter(!is.na(Longitude), !is.na(Latitude))  # 107 individuals jaguar

source = "MoratoJaguar"

###############################################################
# organize data into 1-day and 10-day displacement -----------
###############################################################

target_time_scale_days = 1 # set up the target time scale
# target_time_scale_days = 10 # set up the target time scale

# --- filter based on data duration -- #
# at least tracked for the targeted time scale (e.g. 10 days) #

df.summary <- df %>%
  group_by(StudyName, ID) %>%
  summarise( temp_interval = median(diff(TimestampUTC)),
             duration = max(TimestampUTC) - min(TimestampUTC)) %>%
  mutate( temp_interval_hours = as.numeric(temp_interval, "hours"),
          duration_days = as.numeric(duration, "days"))
df.keep <- df.summary %>% filter(duration_days >= target_time_scale_days) 

# reduce dataset to individuals monitored for at least xx days at once
df <- df %>%
  filter(StudyName %in% df.keep$StudyName,
         ID %in% df.keep$ID)
rm(df.summary, df.keep)

# --- filter based on data location -- #
bad_id = unique(df[which(abs(df$Longitude) > 180 ),]$StudyName)
df <- df %>% filter(!StudyName %in% bad_id ) # coordinates of study id not in lat/lon; n = 3792215 points left (10d) n = 3825479 (1d)

# generate remaining data info df, create unique ID for each individual (4 digit ID + species)
df.info <- df %>% ungroup() %>%
  dplyr::select(StudyName, ID, Binomial, StudyName, ContactPerson) %>% distinct() 
nrow(df.info) 
any(duplicated(df.info$ID))


###############################################################################
####  filter data into 10/1 day intervals and calculate disp -----------------   
###############################################################################
xDayDisp_df <- df %>%
  mutate(date = as.Date(TimestampUTC)) %>%
  group_by(StudyName, ID) %>%
  arrange(TimestampUTC) %>%
  mutate(day_group = (as.numeric(date) - min(as.numeric(date))) %/% target_time_scale_days) %>%
  group_by(StudyName, ID, day_group) %>%
  summarize(
    # Get the first and last point of each 10-day period within the valid period
    first_point = cbind(first(Longitude), first(Latitude)),
    last_point = cbind(last(Longitude), last(Latitude)),
    first_timestamp = first(TimestampUTC),  # Capture the first TimestampUTC for reference
    last_timestamp = last(TimestampUTC),    # Capture the last TimestampUTC for reference
    
    # Check if the first point and the last point are 10 days apart
    timediffs = difftime(last_timestamp, first_timestamp, units = "hours"))  # n = 21,218 (10d) n = 195,859) 

# delete data if not enough recorded locations in each time segment
xDayDisp_df  <- xDayDisp_df %>% 
  filter(as.numeric(timediffs) > 24*0.9*target_time_scale_days) %>% # give it a 10% tolerance
  mutate(x_day_displacement_km = distHaversine(first_point, last_point) / 1000) # calculate displacement in km 

# calculate middle point for percolation distance extraction 
xDayDisp_df <- xDayDisp_df %>% ungroup() %>%
  mutate(mid_point = midPoint(first_point, last_point),  # get middle point
         midLongitude = mid_point[,1], midLatitude = mid_point[,2]) %>% # only  # get middle position for extracting percolation distance 
  rename(Displacement_km = x_day_displacement_km) %>%
  dplyr:: select(-day_group, -timediffs)

move1d <- xDayDisp_df %>%
  left_join(df.info, by = c('StudyName', 'ID')) 
rm(xDayDisp_df)

# write_rds(move1d, paste0("./data/movement/other_open_data/move", target_time_scale_days, "d_oodata_",source,".rds")) #otheropendata

###############################################################################
####  Environmental covariates extraction  ------------------------------------   
###############################################################################
#### next need to extract NDVI to each start and end point and take an average for each displacement segment
library(rgee)
library(sf)
library(geojsonio)
ee_Initialize()
ee_Authenticate() 

oodata <- read_rds( paste0("./data/movement/other_open_data/move", target_time_scale_days, "d_oodata_",source,".rds")) %>%
  mutate(unique_serial = row_number(.))

env.sf <- oodata %>% 
  dplyr::select(unique_serial, first_point, last_point, first_timestamp, last_timestamp) %>%
  pivot_longer(cols = c(first_point, last_point), names_to = "point", values_to = "coordinates" ) %>%
  mutate(Longitude = coordinates[,1], Latitude = coordinates[,2],
         TimestampUTC = ifelse (point == "first_point", first_timestamp, last_timestamp) ) %>%
  dplyr::select(-coordinates, -first_timestamp, -last_timestamp) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = st_crs(4326))

# Load MODIS NDVI collection
modis_ndvi <- ee$ImageCollection("MODIS/006/MOD13Q1")$select("NDVI")

# extract image dates from the image collection to later combine with points
modis_dates = ee_get_date_ic(modis_ndvi)$time_start # n = 529

# match modis date with point date
env.sf <- env.sf %>% 
  mutate(modis_date = as.character(modis_dates[findInterval(TimestampUTC, modis_dates)]))

# Loop through each group and extract NDVI
ndvi_list <- tibble()
for (i in as.character(unique(env.sf $modis_date))[1:length(unique(env.sf $modis_date))] ) {
  
  points_group = env.sf %>% filter(modis_date == i) 
  
  # above limit when there are more than 5000 points to be extracted at once.
  n = ceiling(nrow(points_group)/4999)
  if (nrow(points_group) >= 5000) {
    points_group$Group <- rep(1:n, each = 4999, length.out = nrow(points_group))
  } else ( points_group$Group <- 1)
  
  for (ii in 1:n) {
    
    points_group_ii <- points_group %>% filter(Group == ii)
    
    filtered_modis <- modis_ndvi$
      filterDate(as.character(i))$ 
      median()  # reduce from IC to image 
    
    ndvi_values <- ee_extract(
      x = filtered_modis,
      y = points_group_ii["unique_serial"],
      scale = 250,  # MODIS resolution
      fun = ee$Reducer$mean()    # Aggregation function (not critical for points)
    )
    
    ndvi_list <- rbind(ndvi_list, ndvi_values)
  }
  
} # 

# merge NDVI with other displacement info 
oodata <- oodata %>% 
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

# merge NDVI, HFI, and HMI back with displacement df
oodata <- oodata %>% 
  left_join(env.sf %>% st_drop_geometry(.) %>% 
              group_by(unique_serial) %>% 
              summarise (HFI = mean(HFI),
                         HMI = mean(HMI)), # all environmental covariates for each segment are the mean of the start and the end point values
            by = "unique_serial")

write_rds(oodata %>% distinct(),
          paste0("./data/movement/other_open_data/midproduct_",source,"_",
                 target_time_scale_days, "d_disp_w_envinfo.rds"))

###############################################################################
##########  Species covariates extraction  ------------------------------------   
###############################################################################
 # oodata <- read_rds( paste0("./data/movement/other_open_data/midproduct_",source,"_",
 #                            target_time_scale_days, "d_disp_w_envinfo.rds"))

# extract species, Diet, BodyMass_kg info -----------------   
diet <- read_csv("./data/EltonTraits/MamFuncDat.csv") %>% filter(Scientific %in% unique(oodata$Binomial))
sum(unique(oodata$Binomial) %in% unique(diet$Scientific))
unique(oodata$Binomial) [which(!unique(oodata$Binomial) %in% unique(diet$Scientific))]

diet <- diet %>% mutate(Diet = 
                          case_when(
                            `Diet-Inv` + `Diet-Vend` + `Diet-Vect` + `Diet-Vfish` + `Diet-Vunk` + `Diet-Scav` == 0  ~ "Herbivore",
                            `Diet-Fruit` + `Diet-Nect` + `Diet-Seed` + `Diet-PlantO` == 0  ~ "Carnivore",
                            .default = "Omnivore"
                          ),
                        BodyMass_kg = `BodyMass-Value`/1000) %>%
  dplyr::select(Scientific, Diet, BodyMass_kg) %>%
  rename(Binomial = Scientific)

oodata <- oodata %>% left_join(diet, by = "Binomial")

# extract taxonomy info from PANTHERIA -----------------   
pant <- read_csv("./data/PanTHERIA_1-0_WR05_Aug2008.csv")
pant <- pant[,1:5] %>% 
  rename (Order = MSW05_Order,
          Family = MSW05_Family,
          Genus = MSW05_Genus,
          Species = MSW05_Species,
          Binomial = MSW05_Binomial) 

# find out which species are not listed in panethria database 
oodata  %>% filter(!Binomial %in% pant$Binomial) %>% dplyr::select(Binomial) %>% distinct() 

oodata <- oodata %>% left_join(pant, by = "Binomial")

oodata <- oodata %>%                    
  mutate(Longitude = last_point[,1], Latitude = last_point[,2],
         TimestampUTC = last_timestamp) %>%
  dplyr::select(-unique_serial, - first_point, -mid_point, -first_timestamp, -last_timestamp, -last_point)

#write_rds(oodata, "./data/movement/move10d_oodata.rds")
write_rds(oodata, paste0("./data/movement/move",target_time_scale_days,"d_oodata_", source,".rds"))
