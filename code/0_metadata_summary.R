# 0_metadata_summary 

# this scripe summarize the data used for the critical analysis 


library(tidyverse)
target_time_scale_days = 1

move_covid <- read_rds( paste0("./data/movement/ready_data/move", target_time_scale_days,"d_covid.rds"))
move_movebank <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_movebank.rds"))
move_oomove <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_oodata_MoratoJaguar.rds"))

move <- rbind(move_covid %>%
                dplyr::select(study_name, ID, Binomial, TimestampUTC) %>% mutate(source = "covid"),
              move_movebank %>%
                dplyr::select(study_name, ID, Binomial, TimestampUTC) %>% mutate(source = "open_movebank"),
              move_oomove %>%
                dplyr::select(study_name, ID, Binomial, TimestampUTC) %>% mutate(source = "MoratoJaguar"))

meta <- move %>% group_by(study_name, Binomial) %>% summarise(
  n_ind = length(unique(ID)),
  study_start_date = date(min(TimestampUTC)),
  study_end_date = date(max(TimestampUTC)),
  study_span = difftime(study_end_date,study_start_date, "days")
)

meta2 <- move %>%
  group_by(study_name, Binomial, ID) %>%
  summarise(duration = max(TimestampUTC) - min(TimestampUTC)) %>%
  mutate( duration_days = as.numeric(duration, "days")) %>%
  group_by(study_name, Binomial) %>%
  summarise(mean_ind_duration = mean(duration_days))

meta <- meta %>% left_join(meta2)

write_csv(meta, paste0("./data/movement/ready_data/metadata_", target_time_scale_days,"d_non_globalbarrier.csv"))


# ------- 
