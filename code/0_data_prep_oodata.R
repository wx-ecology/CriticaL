# this script download prep all raw movement data into 10d or 1d dataframe with displacement calculated
# type 3: other open movement data

library(tidyverse)
library(geosphere)

df <- read_csv("./data/movement/other_open_data/Morato_2018_Jaguar/jaguar_movement_data.csv")

df <- df %>% rename(study_name = study.name, 
              ID = 'individual.local.identifier (ID)',
              TimestampUTC = timestamp,
              Longitude = location.long,
              Latitude = location.lat,
              Binomial = individual.taxon.canonical.name) %>%
  mutate(TimestampUTC = mdy_hm(TimestampUTC),
         ID = paste0(str_sub(study_name, 1,2), ID, "_", Binomial),
         contact = "Ronaldo G Morato") %>%
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
  group_by(study_name, ID) %>%
  summarise( temp_interval = median(diff(TimestampUTC)),
             duration = max(TimestampUTC) - min(TimestampUTC)) %>%
  mutate( temp_interval_hours = as.numeric(temp_interval, "hours"),
          duration_days = as.numeric(duration, "days"))
df.keep <- df.summary %>% filter(duration_days >= target_time_scale_days) 

# reduce dataset to individuals monitored for at least xx days at once
df <- df %>%
  filter(study_name %in% df.keep$study_name,
         ID %in% df.keep$ID)
rm(df.summary, df.keep)

# --- filter based on data location -- #
bad_id = unique(df[which(abs(df$Longitude) > 180 ),]$study_name)
df <- df %>% filter(!study_name %in% bad_id ) # coordinates of study id not in lat/lon; n = 3792215 points left (10d) n = 3825479 (1d)

# generate remaining data info df, create unique ID for each individual (4 digit ID + species)
df.info <- df %>% ungroup() %>%
  dplyr::select(study_name, ID, Binomial, contact) %>% distinct() 
nrow(df.info) 
any(duplicated(df.info$ID))


###############################################################################
####  filter data into 10/1 day intervals and calculate disp -----------------   
###############################################################################
xDayDisp_df <- df %>%
  mutate(date = as.Date(TimestampUTC)) %>%
  group_by(study_name, ID) %>%
  arrange(TimestampUTC) %>%
  mutate(day_group = (as.numeric(date) - min(as.numeric(date))) %/% target_time_scale_days) %>%
  group_by(study_name, ID, day_group) %>%
  summarize(
    # Get the first and last point of each 10-day period within the valid period
    first_point = cbind(first(Longitude), first(Latitude)),
    last_point = cbind(last(Longitude), last(Latitude)),
    first_timestamp = first(TimestampUTC),  # Capture the first TimestampUTC for reference
    last_timestamp = last(TimestampUTC),    # Capture the last TimestampUTC for reference
    
    # Check if the first point and the last point are 10 days apart
    timediffs = difftime(last_timestamp, first_timestamp, units = "hours"))  # n = 195,859) n = 134670

# delete data if not enough recorded locations in each time segment
xDayDisp_df  <- xDayDisp_df %>% 
  filter(as.numeric(timediffs) > 24*0.9*target_time_scale_days) %>% # give it a 10% tolerance
  mutate(mid_point = midPoint(first_point, last_point), # get middle point
         Displacement_km = distHaversine(first_point, last_point) / 1000) %>% # calculate displacement in km 
  dplyr:: select(-day_group, -timediffs) %>%
  ungroup() %>%
  left_join(df.info, by = c('study_name', 'ID')) %>%
  relocate(ID, first_point, mid_point, last_point, first_timestamp, last_timestamp, 
           Displacement_km, Binomial, study_name, contact)

write_rds(xDayDisp_df , paste0("./data/movement/midproduct/midproduct_move", target_time_scale_days, "d_oodata_",source,".rds")) #otheropendata
