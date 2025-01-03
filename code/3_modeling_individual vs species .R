# this script runs models at the interspecies level, 
# but investigate whether it is the species level response to percolation distance 
# (species occurence effect, where certain species exhibit long-range movement simply do not occur in areas of low percolation distance)
# or individual behavior effect ( individual alter their movement to percolation distances)

library(tidyverse)
library(AICcmodavg)
source("./code/HaversineLMEfunctions.R")
library(nlme)
library(ggplot2)
#library(MuMIn)

###################################################
# ----- read data and prep data for modeling ------ 
###################################################
target_time_scale_days = 10
tartget_spatial_scale = "me" 

move_covid <- read_rds( paste0("./data/movement/ready_data/move", target_time_scale_days,"d_covid.rds"))
move_movebank <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_movebank.rds"))
move_oomove <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_oodata_MoratoJaguar.rds"))

# #library(GGally)
# ggpairs(move_covid, columns = c(2, 6:9,11:12, 14:15, 21))
# colSums(is.na(move_covid)) # some NAs in NDVI and HMI 

if (tartget_spatial_scale == "hi") {
  move <- rbind(move_covid %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi) ,
                move_movebank %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi) ,
                move_oomove %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi)) %>%
    rename(pd_adpt = pd_adpt_hi,
           bd_adpt = bd_adpt_hi,
           pd_scale = pd_scale_hi)
} else {
  move <- rbind(move_covid %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me),
                move_movebank %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me),
                move_oomove %>%
                  dplyr::select(ID, Order, Family, Genus, Species, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, Diet, NDVI, HFI, HMI, BodyMass_kg,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me)) %>%
    rename(pd_adpt = pd_adpt_me,
           bd_adpt = bd_adpt_me,
           pd_scale = pd_scale_me)
}

## prep other variables for modeling ## ------
# summary for each individual 
move.sum <- move %>%
  group_by(ID) %>%
  summarise(
    # displacement 
    Displacement_km = ifelse (tartget_spatial_scale == "hi",
                              quantile(Displacement_km, 0.95,  na.rm = T), 
                              median(Displacement_km, na.rm = T)),
    
    # percolation distances 
    pd_adpt_km = mean(pd_adpt[pd_adpt >=0], na.rm = T)/1000, # only summarize percolation distance when there is a value 
    
    #building densities 
    bd_adpt = mean(bd_adpt, na.rm = T),
    
    # other covairates  
    dist_2_build = mean(dist_2_build, na.rm = T), # distance to building 
    NDVI = mean(NDVI, na.rm=T), 
    HFI = mean(HFI, na.rm=T),
    HMI = mean(HMI, na.rm=T),
    BodyMass_kg = mean(BodyMass_kg, na.rm = T),
    Longitude = mean(Longitude, na.rm = T), 
    Latitude = mean(Latitude, na.rm = T)) %>% 
  
  # add spp info 
  left_join( distinct(
    move %>% dplyr::select(ID, Order, Family, Genus, Species, Binomial, Diet, pd_scale)
  ), by = "ID")  #n ind = 2075 for 10d // n ind = 650 for 1d

# remove ID with duplicated locations
if (any(duplicated(c(move.sum$Longitude, move.sum$Latitude)))) {
  move.sum <-move.sum[- which(duplicated(c(move.sum$Longitude, move.sum$Latitude))), ]}  

## remove species that have less than 3 individuals  
table(move.sum$Binomial)
move.sum <- move.sum %>% filter(Binomial %in% (move.sum %>% group_by(Binomial) %>% 
                                                 summarize (n = length(unique(ID))) %>% 
                                                 filter(n >= 3) %>% pull(Binomial))) # 39 -- > 26 spp for 10 d, 19 for 1d
length(unique(move.sum$Binomial))
nrow(move.sum) # 2044 or 625 for 1day

# ggpairs(move.sum, columns = c(2:8))
move.sum <- move.sum %>%
  mutate(log_BodyMass_kg = log(BodyMass_kg), 
         scale_NDVI = scale(NDVI),
         log_HFI = ifelse(HFI==0, log(HFI+0.001), log(HFI)),
         log_HMI = ifelse(HMI==0, log(HMI+0.001), log(HMI)),
         log_dist_2_build = log(dist_2_build),
         log_pd = log(pd_adpt_km),
         log_bd = ifelse(bd_adpt==0, log(bd_adpt+0.001), log(bd_adpt)),
         log_Displacement_km = log(Displacement_km)) %>%
  rename(lon = Longitude,
         lat = Latitude)

# get individual pd and species average pd 
pd.spp <- move.sum %>% group_by(Binomial) %>% summarise(log_pd_spp = mean(log_pd, na.rm = T))
move.sum <- move.sum %>% left_join(pd.spp) %>% mutate(log_pd_ind = log_pd - log_pd_spp)


colSums(is.na(move.sum))
sapply(move.sum, function(col) sum(is.nan(col)))
sapply(move.sum, function(col) sum(is.infinite(col)))

write_rds(move.sum, paste0("./data/movement/ready_data/move_",target_time_scale_days, "d_",
                           tartget_spatial_scale,
                           ".rds"))

###################################################
# -------------------- run models -----------------
###################################################
dat = move.sum %>% drop_na() 
dat <- dat[!apply(dat, 1, function(row) any(is.nan(row))), ] 
# n = 2026 for 10d // 597 for 1d

target_time_scale_days
tartget_spatial_scale 

results = tibble()
for (i in c("log_HFI", "log_HMI", "log_bd") ) {
  
  if (i == "log_bd") {
    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", i, " + log_dist_2_build + log_pd_spp + log_pd_ind"))
  } else {
    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", i, " + log_pd_spp + log_pd_ind"))
  }
  
  print (paste0("calculating for ", i, " SPP vs IND model..."))
  m = lme(formula,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
           random= ~1|Order/Family/Genus/Species,
           control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
           data = dat)
  saveRDS(m, paste0("./results/models/sppVSind_",target_time_scale_days, "d_",tartget_spatial_scale, "_", i, ".rds"))
  
  
  sum.i = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                 spatial_scale = paste0("disp_", tartget_spatial_scale),
                 amount_variable = i,
                 beta_pd_spp = (summary(m))$tTable["log_pd_spp", "Value"],
                 upper_pd_spp = intervals(m, which = "fixed")$fixed["log_pd_spp","upper"],
                 lower_pd_spp = intervals(m, which = "fixed")$fixed["log_pd_spp","lower"],
                 p_value_pd_spp = (summary(m))$tTable["log_pd_spp", "p-value"],
                 beta_pd_ind = (summary(m))$tTable["log_pd_ind", "Value"],
                 upper_pd_ind = intervals(m, which = "fixed")$fixed["log_pd_ind","upper"],
                 lower_pd_ind = intervals(m, which = "fixed")$fixed["log_pd_ind","lower"],
                 p_value_pd_ind = (summary(m))$tTable["log_pd_ind", "p-value"])
  
  results = rbind(results, sum.i)
} 


write_csv(results, paste0("./results/ModResults_SPPvsIND_move",target_time_scale_days, "d_",tartget_spatial_scale, ".csv"))

