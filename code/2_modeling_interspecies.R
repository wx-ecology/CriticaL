# this script runs models at the interspecies level 
# update Feb 17 -> update the spatial autocorrelation structure
# and also extract individual marginal R2 

library(tidyverse)
library(AICcmodavg)
#source("./code/HaversineLMEfunctions.R")
library(nlme)
library(ggplot2)
library(glmm.hp)

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

## remove species that have *****less than 5***** individuals  
table(move.sum$Binomial)
move.sum <- move.sum %>% filter(Binomial %in% (move.sum %>% group_by(Binomial) %>% 
                                                 summarize (n = length(unique(ID))) %>% 
                                                 filter(n > 5) %>% pull(Binomial))) 
length(unique(move.sum$Binomial)) # 39 -- > 26 spp for 10 d, 19 for 1d
nrow(move.sum) # 2044 for 10 day or 625 for 1day

# ggpairs(move.sum, columns = c(2:8))
# transformation for modeling. scale and log so the estimates are at comparative scale.
move.sum <- move.sum %>%
  mutate(log_BodyMass_kg = log(BodyMass_kg), 
         scale_NDVI = as.vector(scale(NDVI)),
         scale_HFI = scale(HFI),
         scale_HMI = scale(HMI),
         log_dist_2_build = log(dist_2_build),
         log_pd = scale(log(pd_adpt_km)), # very similar scaled or not 
         slog_bd = scale(ifelse(bd_adpt==0, log(bd_adpt+0.001), log(bd_adpt))),
         log_Displacement_km = log(Displacement_km)) %>%
  rename(lon = Longitude,
         lat = Latitude)

# GGally::ggpairs(move.sum, columns = 19:26)

colSums(is.na(move.sum))
sapply(move.sum, function(col) sum(is.nan(col)))
sapply(move.sum, function(col) sum(is.infinite(col)))

###################################################
# -------------------- run models -----------------
###################################################
# remove na 
dat = move.sum %>% drop_na() 
dat <- dat[!apply(dat, 1, function(row) any(is.nan(row))), ] 
dat <- dat %>% filter(Binomial %in% (dat %>% group_by(Binomial) %>% 
                                      summarize (n = length(unique(ID))) %>% 
                                      filter(n > 5) %>% pull(Binomial))) 
# n = 2026 - 2018 for 10d // 597 for 1d

target_time_scale_days
tartget_spatial_scale 

models = tibble()
results = tibble()
mR2 = tibble()
for (i in c("scale_HFI", "scale_HMI", "slog_bd") ) {
  
  # ------------------------- HFI -----------------------------------------------------
  if(i == "scale_HFI") {
   
  # null model  
  print (paste0("calculating for ", i, " null model..."))
  m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + scale_HFI,
          correlation = corSpatial(form=~lon+lat, type = "exponential"),
          random= ~1|Order/Family/Genus/Species,
          control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000), method = "ML",
          data = dat)

  sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                  spatial_scale = paste0("disp_", tartget_spatial_scale),
                  amount_variable = i,
                  config_variable = NA,
                  AIC = AIC(m),
                  beta_pd = NA,
                upper_pd = NA,
                lower_pd = NA,
                  p_value_pd = NA)
  
  models <- rbind(models, tibble_row(comp_var = i, conf_var = NA, model = m, mR2 = NA))
  results = rbind(results, sum)
  
  # model plus pd 
  print (paste0("calculating ", i, " pd model..."))
  
  m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + scale_HFI + log_pd,
          correlation = corSpatial(form=~lon+lat, type = "exponential"),
          random= ~1|Order/Family/Genus/Species,
          control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, tolerance = 1e-4, niterEM = 50),
          method = "ML",
          data = dat)
  
  # trying opt = "nlminb" - worked 
  
  print (paste0("calculating mR2 for the ", i, " pd model..."))
  mR2.i <- glmm.hp(m)
  models <- rbind(models, tibble_row(comp_var = i, conf_var = "log_pd", model = m, mR2 = mR2.i ))
  
  sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                spatial_scale = paste0("disp_", tartget_spatial_scale),
                amount_variable = i,
                config_variable = "log_pd",
                AIC = AIC(m),
                beta_pd = (summary(m))$tTable["log_pd", "Value"],
                upper_pd = intervals(m, which = "fixed")$fixed["log_pd","upper"],
                lower_pd = intervals(m, which = "fixed")$fixed["log_pd","lower"],
                p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
  
  # summarizing results 
  mR2.i <- mR2.i$hierarchical.partitioning 
  mR2 <- rbind(mR2,
               tibble (comp_var = i, conf_var = "log_pd", fixed_effect = rownames(mR2.i), 
                       Unique = mR2.i[, "Unique"], Individual = mR2.i[, "Individual"], 
                       I.perc = mR2.i[, "I.perc(%)"] ))
  results = rbind(results, sum) 
   }
  
  # ------------------------- HMI -----------------------------------------------------
  if(i == "scale_HMI") {
    
    # null model  
    print (paste0("calculating ", i, " null model..."))
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + scale_HMI,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                  spatial_scale = paste0("disp_", tartget_spatial_scale),
                  amount_variable = i,
                  config_variable = NA,
                  AIC = AIC(m),
                  beta_pd = NA,
                  upper_pd = NA,
                  lower_pd = NA,
                  p_value_pd = NA)
    
    models <- rbind(models, tibble_row(comp_var = i, conf_var = NA, model = m, mR2 = NA))
    results = rbind(results, sum)
    
    # model plus pd 
    print (paste0("calculating ", i, " pd model..."))
    
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + scale_HMI + log_pd,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    print (paste0("calculating mR2 for the ", i, " pd model..."))
    mR2.i <- glmm.hp(m)
    models <- rbind(models, tibble_row(comp_var = i, conf_var = "log_pd", model = m, mR2 = mR2.i ))
    
    sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                  spatial_scale = paste0("disp_", tartget_spatial_scale),
                  amount_variable = i,
                  config_variable = "log_pd",
                  AIC = AIC(m),
                  beta_pd = (summary(m))$tTable["log_pd", "Value"],
                  upper_pd = intervals(m, which = "fixed")$fixed["log_pd","upper"],
                  lower_pd = intervals(m, which = "fixed")$fixed["log_pd","lower"],
                  p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    
    # summarizing results 
    mR2.i <- mR2.i$hierarchical.partitioning 
    mR2 <- rbind(mR2,
                 tibble (comp_var = i, conf_var = "log_pd", fixed_effect = rownames(mR2.i), 
                         Unique = mR2.i[, "Unique"], Individual = mR2.i[, "Individual"], 
                         I.perc = mR2.i[, "I.perc(%)"] ))
    results = rbind(results, sum) 
  }

  # ------------------------- building density -----------------------------------------------------
  if (i == "slog_bd") {
    
    # null model  
    print (paste0("calculating ", i, " null model..."))
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + slog_bd,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    sum =  tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                  spatial_scale = paste0("disp_", tartget_spatial_scale),
                  amount_variable = i,
                  config_variable = NA,
                  AIC = AIC(m),
                  beta_pd = NA,
                  upper_pd = NA,
                  lower_pd = NA,
                  p_value_pd = NA)
    
    models <- rbind(models, tibble_row(comp_var = i, conf_var = NA, model = m, mR2 = NA))
    results = rbind(results, sum)
    
    # model plus distance 2 building
    print (paste0("calculating ", i, " dist2build model..."))
    
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + slog_bd + log_dist_2_build,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    print (paste0("calculating mR2 for the ", i, " dist2build model..."))
    mR2.i <- glmm.hp(m)
    
    models <- rbind(models, tibble_row(comp_var = i, conf_var = "log_dist_2_build", model = m, mR2 = mR2.i ))
    
    sum = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                 spatial_scale = paste0("disp_", tartget_spatial_scale),
                 amount_variable = i,
                 config_variable = "log_dist_2_build",
                 AIC = AIC(m),
                 beta_pd = NA,
                 upper_pd = NA,
                 lower_pd = NA,
                 p_value_pd = NA)
    
    mR2.i <- mR2.i$hierarchical.partitioning 
    mR2 <- rbind(mR2,
                 tibble (comp_var = i, conf_var = "log_dist_2_build", fixed_effect = rownames(mR2.i), 
                         Unique = mR2.i[, "Unique"], Individual = mR2.i[, "Individual"], 
                         I.perc = mR2.i[, "I.perc(%)"] ))
    
    results = rbind(results, sum) 
    
    # model plus log_pd
    print (paste0("calculating ", i, " log_pd model..."))
    
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + slog_bd + log_pd,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    print (paste0("calculating mR2 for the ", i, " log_pd model..."))
    mR2.i <- glmm.hp(m)
    
    models <- rbind(models, tibble_row(comp_var = i, conf_var = "log_pd", model = m, mR2 = mR2.i ))
    
    sum = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                 spatial_scale = paste0("disp_", tartget_spatial_scale),
                 amount_variable = i,
                 config_variable = "log_pd",
                 AIC = AIC(m),
                 beta_pd = (summary(m))$tTable["log_pd", "Value"],
                 upper_pd = intervals(m, which = "fixed")$fixed["log_pd","upper"],
                 lower_pd = intervals(m, which = "fixed")$fixed["log_pd","lower"],
                 p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    
    mR2.i <- mR2.i$hierarchical.partitioning 
    mR2 <- rbind(mR2,
                 tibble (comp_var = i, conf_var = "log_pd", fixed_effect = rownames(mR2.i), 
                         Unique = mR2.i[, "Unique"], Individual = mR2.i[, "Individual"], 
                         I.perc = mR2.i[, "I.perc(%)"] ))
    results = rbind(results, sum) 
    
    # model dist2build + log_pd
    print (paste0("calculating ", i, " dist2build + pd model..."))
    
    m = lme(log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + slog_bd + log_dist_2_build + log_pd,
            correlation = corSpatial(form=~lon+lat, type = "exponential"),
            random= ~1|Order/Family/Genus/Species,
            control = list(opt = "nlminb", msMaxIter = 1000, msMaxEval = 1000, 
                           tolerance = 1e-4, niterEM = 50), 
            method = "ML",
            data = dat)
    
    print (paste0("calculating mR2 for the ", i, " dist2build + pd model..."))
    mR2.i <- glmm.hp(m)
    
    models <- rbind(models, tibble_row(comp_var = i, conf_var = "log_dist_2_build + log_pd", model = m, mR2 = mR2.i ))
    
    sum = tibble(temporal_scale = paste0(target_time_scale_days, "d"),
                 spatial_scale = paste0("disp_", tartget_spatial_scale),
                 amount_variable = i,
                 config_variable = "log_dist_2_build + log_pd",
                 AIC = AIC(m),
                 beta_pd = (summary(m))$tTable["log_pd", "Value"],
                 upper_pd = intervals(m, which = "fixed")$fixed["log_pd","upper"],
                 lower_pd = intervals(m, which = "fixed")$fixed["log_pd","lower"],
                 p_value_pd = (summary(m))$tTable["log_pd", "p-value"])
    
    mR2.i <- mR2.i$hierarchical.partitioning 
    mR2 <- rbind(mR2,
                 tibble (comp_var = i, conf_var = "log_dist_2_build + log_pd", fixed_effect = rownames(mR2.i), 
                         Unique = mR2.i[, "Unique"], Individual = mR2.i[, "Individual"], 
                         I.perc = mR2.i[, "I.perc(%)"] ))
    results = rbind(results, sum) 
  }
}

saveRDS(models, paste0("./results/models/zscored/Mods_allspp_",target_time_scale_days, "d_",tartget_spatial_scale,".rds"))
write_csv(results, paste0("./results/ModResults_allSpp/scale_zscored/ModResults_move",target_time_scale_days, "d_",tartget_spatial_scale, ".csv"))
write_csv(mR2, paste0("./results/ModResults_allSpp/scale_zscored/mR2Results_move",target_time_scale_days, "d_",tartget_spatial_scale, ".csv"))


# # when mR2 has nan issue
# mR2 <- rbind(mR2,
#              tibble (comp_var = i, conf_var = "XXXX", fixed_effect = names(fixed.effects(m))[2:length(names(fixed.effects(m)))], 
#                      Unique = NA, Individual = NA, 
#                      I.perc = NA ))
