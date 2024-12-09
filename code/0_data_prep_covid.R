# this script download prep all raw movement data into 10d or 1d dataframe with displacement calculated
# type 2: covid data

# because the open covid data only recorded starting location of the 10-day displacement, we reshaped the data to identify
# the starting and the ending locations and recalculated a roughly-10-day displacement, 
# with all other environmental covariates freshly extracted

library(sf)
library(tidyverse)
library(geosphere)

target_time_scale_days = 1

if (target_time_scale_days == 10) {
  move <- readRDS("./data/movement/Tucker_7704108/Tucker_10d_Spatial.rds") %>%
    dplyr::select(-"StringencyIndex", - "Lockdown(Y/N)", -"Country", - "RBrain_Size", -"Activity")
} else (
  move <- readRDS("./data/movement/Tucker_7704108/Tucker_1h_Move_Spatial.rds")  %>%
    dplyr::select(-"StringencyIndex", - "Lockdown(Y/N)", -"Country", - "RBrain_Size", -"Activity")
)


####################################
## reshape data frame ------ #######
# ----- identify start/end points and mid point, recalculate displacement ----------
####################################

if (target_time_scale_days == 1) {
  # ----- turn 1h data to 1d----------
  move_reshape <- move %>%
    dplyr::select(-Displacement_km) %>%
    mutate(date = as.Date(TimestampUTC)) %>%
    dplyr::arrange(ID, TimestampUTC) %>%
    group_by(ID, date) %>%
    summarize(
      
      # Get the first and last point of the day within the valid period
      first_point = cbind(first(Longitude), first(Latitude)),
      last_point = cbind(last(Longitude), last(Latitude)),
      first_timestamp = first(TimestampUTC),  # Capture the first timestamp for reference
      last_timestamp = last(TimestampUTC),    # Capture the last timestamp for reference
      
      # calculate time diff between the first and the last point of each day
      time_diff = as.numeric(difftime(last_timestamp, first_timestamp, units = "hours"))) %>% 
    
    # allow for 10% time diff
    filter(time_diff > target_time_scale_days*0.9*24) %>%
    
    # merge with other info 
    dplyr::select(-date, - time_diff) %>%
    left_join(move %>% dplyr::select(-Longitude, -Latitude) %>% rename(first_timestamp = TimestampUTC),
              by = c("ID", "first_timestamp")) %>%  # 493348 -> 20065
    
    ungroup () %>%
    mutate( mid_point = midPoint(first_point, last_point),
            Displacement_km = distHaversine(first_point, last_point) / 1000 ) %>%
    rename(study_name = StudyName, Binomial = Species, contact = ContactPerson) %>%
    dplyr::select (ID, first_point, mid_point, last_point,first_timestamp, last_timestamp, 
                   Displacement_km, Binomial,study_name, contact)
  
} else {
  
move_reshape <- move %>%
  arrange(ID, TimestampUTC) %>%
  group_by(ID) %>%
  mutate(time_diff = as.numeric(difftime(lead(TimestampUTC), TimestampUTC, units = "days")), # calculate consecutive points time difference
         # get start and end point 
         start_long = Longitude, 
         start_lat = Latitude,
         end_long = lead(Longitude),
         end_lat = lead(Latitude),
         first_timestamp = TimestampUTC,
         last_timestamp = lead(TimestampUTC)
         ) %>%
  filter(!is.na(time_diff), #end point of each individual has no next point
         time_diff < target_day_scale*1.1 & time_diff > target_day_scale*0.9) %>% # displacement for 10+- 10% days
  ungroup()  %>% # 11523 -> 12422 for 10 day; 18,071 -> 17,881 for 1 day
  
  # calculate displacement and get mid point.
  mutate(
  first_point = matrix(c(start_long, start_lat), nrow = nrow(move_reshape)),
  last_point = matrix(c(end_long, end_lat), nrow = nrow(move_reshape)),
  mid_point = midPoint(first_point, last_point),
  Displacement_km = distHaversine(first_point, last_point) / 1000 #recalculate displacement. checked already. almost identical as before.
)  %>%
  rename(study_name = StudyName, Binomial = Species, contact = ContactPerson) %>%
  dplyr::select (ID, first_point, mid_point, last_point,first_timestamp, last_timestamp, 
                 Displacement_km, Binomial,study_name, contact)
}

write_rds(move_reshape, paste0("./data/movement/midproduct/midproduct_move", target_time_scale_days , "d_covid.rds"))
