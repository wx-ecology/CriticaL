# this script runs models at the interspecies level 

library(tidyverse)
library(AICcmodavg)
source("./code/HaversineLMEfunctions.R")
library(nlme)
library(ggplot2)
#library(MuMIn)

# function to identify outliers based on standardized residuals
get_clean_dat <- function (m, threshold) {
  dat$std_resid <- residuals(m, type = "normalized")
  dat_clean <- dat[abs(dat$std_resid) <= threshold, ]
}

###################################################
# ----- read data and prep data for modeling ------ 
###################################################
target_time_scale_days = 1
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
    pd_adpt = mean(pd_adpt[pd_adpt >=0], na.rm = T), # only summarize percolation distance when there is a value 
    
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

## remove species that have *****less than 3***** individuals  
table(move.sum$Binomial)
move.sum <- move.sum %>% filter(Binomial %in% (move.sum %>% group_by(Binomial) %>% 
                                                 summarize (n = length(unique(ID))) %>% 
                                                 filter(n >= 3) %>% pull(Binomial))) 
length(unique(move.sum$Binomial)) # 39 -- > 26 spp for 10 d, 19 for 1d
nrow(move.sum) # 2044 for 10 day or 625 for 1day

# ggpairs(move.sum, columns = c(2:8))
# transformation for modeling 
move.sum <- move.sum %>%
  mutate(log_BodyMass_kg = log(BodyMass_kg), 
         scale_NDVI = scale(NDVI),
         log_HFI = ifelse(HFI==0, log(HFI+0.001), log(HFI)),
         log_HMI = ifelse(HMI==0, log(HMI+0.001), log(HMI)),
         log_dist_2_build = log(dist_2_build),
         log_pd = log(pd_adpt),
         log_bd = ifelse(bd_adpt==0, log(bd_adpt+0.001), log(bd_adpt)),
         log_Displacement_km = log(Displacement_km)) %>%
  rename(lon = Longitude,
         lat = Latitude)

colSums(is.na(move.sum))
sapply(move.sum, function(col) sum(is.nan(col)))
sapply(move.sum, function(col) sum(is.infinite(col)))

###################################################
# -------------------- run models -----------------
###################################################
# remove na 
dat = move.sum %>% drop_na() 
dat <- dat[!apply(dat, 1, function(row) any(is.nan(row))), ] 
# n = 2026 for 10d // 597 for 1d

target_time_scale_days
tartget_spatial_scale 

results = tibble()
for (i in c("log_HFI", "log_HMI", "log_bd") ) {
  
  # null model for HFI and HMI 
  formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", i))
  
  print (paste0("calculating for ", i, " null model..."))
  m = lme(formula,
            correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
            data = dat)
  
  sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                  spatial_scale = paste0("disp_", tartget_spatial_scale),
                  amount_variable = i,
                  config_variable = NA,
                  AIC = AIC(m),
                  beta_pd = NA,
                  p_value_pd = NA)
  
  saveRDS(m, paste0("./results/models/allspp_",target_time_scale_days, "d_",tartget_spatial_scale, "_", i, "_null.rds"))
  results = rbind(results, sum)
  
  # plus pd 
  if (i != "log_bd") {

    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", 
                                 i, " + log_pd"))
    
    print (paste0("calculating for ", i, " pd model..."))
    m = lme(formula,
              correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
              random= ~1|Order/Family/Genus/Species,
              control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
              data = dat)
    # dat_clean <- get_clean_dat(m,3)
    # m = lme(formula,
    #         correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
    #         random= ~1|Order/Family/Genus/Species,
    #         control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
    #         data = dat_clean)
    
    sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                    spatial_scale = paste0("disp_", tartget_spatial_scale),
                    amount_variable = i,
                    config_variable = "log_pd",
                    AIC = AIC(m),
                    beta_pd = (summary(m))$tTable["log_pd", "Value"],
                    p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    
    saveRDS(m, paste0("./results/models/allspp_",target_time_scale_days, "d_",tartget_spatial_scale, "_",i, "_w_pd.rds"))
    results = rbind(results, sum)
    
  } else {
    
    print (paste0("calculating for ", i, " dist2build model..."))
    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", 
                                 i, "+ log_dist_2_build"))
    m <- lme(formula,
                   correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
                   random= ~1|Order/Family/Genus/Species,
                   control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
                   data = dat)
    
    sum.1 = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                   spatial_scale = paste0("disp_", tartget_spatial_scale),
                   amount_variable = i,
                   config_variable = "log_dist_2_build",
                   AIC = AIC(m),
                   beta_pd = NA,
                   p_value_pd = NA)
    saveRDS(m, paste0("./results/models/allspp_",target_time_scale_days, "d_",tartget_spatial_scale, "_", i, "_w_dist2build.rds"))
    
    print (paste0("calculating for ", i, " pd model..."))
    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", 
                                   i, "+ log_pd"))
    m <- lme(formula,
               correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
               random= ~1|Order/Family/Genus/Species,
               control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
               data = dat)
    
    sum.2 = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                   spatial_scale = paste0("disp_", tartget_spatial_scale),
                   amount_variable = i,
                   config_variable = "log_pd",
                   AIC = AIC(m),
                   beta_pd = (summary(m))$tTable["log_pd", "Value"],
                   p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    saveRDS(m, paste0("./results/models/allspp_",target_time_scale_days, "d_",tartget_spatial_scale, "_", i, "_w_pd.rds"))
    
    print (paste0("calculating for ", i, " dist2build + pd  model..."))
    formula <- as.formula(paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ", 
                                 i, "+ log_dist_2_build + log_pd"))
    m <- lme(formula,
                  correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
                  random= ~1|Order/Family/Genus/Species,
                  control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
                  data = dat)
    
    sum.3 = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                   spatial_scale = paste0("disp_", tartget_spatial_scale),
                   amount_variable = i,
                   config_variable = "log_dist_2_build + log_pd",
                   AIC = AIC(m),
                   beta_pd = (summary(m))$tTable["log_pd", "Value"],
                   p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    saveRDS(m, paste0("./results/models/allspp_",target_time_scale_days, "d_",tartget_spatial_scale, "_", i, "_w_dist2build_and_pd.rds"))
    
    results = rbind(results, sum.1, sum.2, sum.3)
  }
}

write_csv(results, paste0("./results/ModResults_move",target_time_scale_days, "d_",tartget_spatial_scale, "outlier_deleted.csv"))

