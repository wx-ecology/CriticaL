# this script download prep all raw movement data into 10d or 1d dataframe with displacement calculated
# type 1: movebank data


library(move2)
library(tidyverse)
library(geosphere)

setwd("/Users/wenjing/Senckenberg Dropbox/Wenjing Xu/CriticaL/")

target_time_scale_days = 10 # set up the target time scale
# target_time_scale_days = 10 # set up the target time scale

################################
##  download movebank data -----
################################
# movebank_store_credentials("wenjingxu", "Move_1877!0")
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
#                  StudyName = movebank_download_study_info(study_id = id)$name,
#                  mv.i, 
#                  ContactPerson = movebank_download_study_info(study_id = id)$contact_person_name) %>%
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

# --- filter NA locations -- #
df <- df %>% filter(!is.na(coords_x), !is.na(coords_y)) # n = 3826440

# --- filter based on data duration -- #
# at least tracked for the targeted time scale (e.g. 10 days) #

df.summary <- df %>%
  group_by(study_id, individual_id, deployment_id) %>%
  summarise( temp_interval = median(diff(timestamp)),
             duration = max(timestamp) - min(timestamp)) %>%
  mutate( temp_interval_hours = as.numeric(temp_interval, "hours"),
          duration_days = as.numeric(duration, "days"))
df.keep <- df.summary %>% filter(duration_days >= target_time_scale_days) 

# reduce dataset to individuals monitored for at least xx days at once
df <- df %>%
  filter(study_id %in% df.keep$study_id,
         individual_id %in% df.keep$individual_id,
         deployment_id %in% df.keep$deployment_id)
rm(df.summary, df.keep)

# --- filter based on data duration -- #
bad_id = unique(df[which(abs(df$coords_x) > 180 ),]$study_id)
df <- df %>% filter(!study_id %in% bad_id ) # coordinates of study id not in lat/lon; n = 3792215 points left (10d) n = 3825479 (1d)

# generate remaining data info df, create unique ID for each individual (4 digit ID + species)
df.info <- df %>% ungroup() %>%
  dplyr::select(study_id, individual_id, species, StudyName, ContactPerson) %>% distinct() %>% 
  mutate(
    individual_id = as.character(individual_id),
    ID = paste0(str_sub(individual_id, -4,-1), "_", species)
  ) 
nrow(df.info) #n_ind = 473 (10d) n_ind = 536 (1d)
any(duplicated(df.info$ID))

###############################################################################
####  filter data into 10/1 day intervals and calculate disp -----------------   
###############################################################################
xDayDisp_df <- df %>%
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
    timediffs = difftime(last_timestamp, first_timestamp, units = "hours"))  # n = 21,218 (10d) n = 195,859) 

# delete data if not enough recorded locations in each time segment
xDayDisp_df  <- xDayDisp_df %>% 
  filter(as.numeric(timediffs) > 24*0.9*target_time_scale_days) %>% # give it a 10% tolerance # n = 18,081, n = 81,965
  mutate(Displacement_km = distHaversine(first_point, last_point) / 1000) # calculate displacement in km 

# calculate middle point for percolation distance extraction 
xDayDisp_df <- xDayDisp_df %>% ungroup() %>%
  mutate(mid_point = midPoint(first_point, last_point)) %>% # get middle point  
  dplyr:: select(-deployment_id, -day_group, -timediffs)

move_movebank <- xDayDisp_df %>% mutate(individual_id = as.character(individual_id)) %>%
  left_join(df.info, by = c('study_id', 'individual_id')) %>%
  dplyr::select(-individual_id, -midLongitude, -midLatitude) %>% # n = 18,081 - 10d; n = 81,965 - 1d
  rename(movebank_study_id = study_id, Binomial = species) %>%
  relocate(ID, first_point, mid_point, last_point, first_timestamp, last_timestamp, 
           Displacement_km, Binomial, movebank_study_id, study_name, contact)

rm(xDayDisp_df)

#write_rds(move_movebank, paste0("./data/movement/midproduct/midproduct_move", target_time_scale_days , "d_movebank.rds"))

