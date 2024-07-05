### this script annotate each 10 day and 1 h steps with percolation distance 
### (the average value between each sequential pair of locations according to tucker 2018)

library(terra)
library(sf)
library(tidyverse)
library(amt)

#################################################################
###  data 
#################################################################
# read movement data
# move1h <- readRDS("./data/Tucker_7704108/Tucker_1h_Move_Spatial.rds") 
move10d <- readRDS("./data/Tucker_7704108/Tucker_10d_Spatial.rds")

# pantheria data base
pant <- read_csv("./data/PanTHERIA_1-0_WR05_Aug2008.csv")
pant <- pant[,1:5] %>% 
  rename (Order = MSW05_Order,
          Family = MSW05_Family,
          Genus = MSW05_Genus,
          Species = MSW05_Species,
          Binomial = MSW05_Binomial)

# read perc dist layers
perc_dist_1km <- rast("./data/perc_dists_may2024/radius_1km/d1km_r1km_results_perc_dist_fin.tif")
perc_dist_5km <- rast("./data/perc_dists_may2024/radius_5km/d1km_r5km_results_perc_dist_fin.tif")
perc_dist_10km <- rast("./data/perc_dists_may2024/radius_10km/d1km_r10km_results_perc_dist_fin.tif")

# read building amount layers
bd_area_1km <- rast("./data/perc_dists_may2024/radius_1km/d1km_r1km_results_global_overlap.tif") 
bd_area_5km <- rast("./data/perc_dists_may2024/radius_5km/d1km_r5km_results_global_overlap.tif") 
bd_area_10km <- rast("./data/perc_dists_may2024/radius_10km/d1km_r10km_results_global_overlap.tif") 

#################################################################
############ calculate repsonse variables #############
#################################################################
# calculate 0.5 and 0.95 quantile displacement distance 
move10d.track <- make_track(move10d, .x = Longitude, .y = Latitude, .t = TimestampUTC, crs = 4326, all_cols = T )
move10d.steps.all <- move10d.track %>% 
  transform_coords(., 3395) %>% # to world mercator
  nest(data = -"ID") %>%
  # calculate step length and step time diff 
  mutate(steps = map(data, function(x) {steps(x)})) %>%
  # unnest back to the big df 
  select(ID, steps) %>%
  unnest(cols = steps) 

# filter out the steps that are not 10 day 
move10d.steps <- move10d.steps.all %>%
  # filter out the steps that are not 10 day 
  mutate(dt_ = as.numeric(dt_, units = "days")) %>%
  filter(dt_ < 11 & dt_ > 9)

# summarize 50% and 95% percentile displacements 
move10d.summary <- move10d.steps %>% 
  group_by(ID) %>%
  summarise (dist_50 = median(sl_),
          dist_95 = quantile(sl_, 0.95))

#################################################################
############ calculate percolation distance variables #############
#################################################################
# construct a df with all movement locations that are actually used 
movement.df <- unique(
  rbind (
  tibble(ID = move10d.steps$ID, x = move10d.steps$x1_, y = move10d.steps$y1_, t = move10d.steps$t1_), 
  tibble(ID = move10d.steps$ID, x = move10d.steps$x2_, y = move10d.steps$y2_, t = move10d.steps$t2_)
  )
  ) # less points than the original data frame

# turn data frame to spatial 
movement.sf <- st_as_sf(movement.df, coords = c("x", "y"), crs = 3395) %>%
  st_transform(crs = 4326)

###  extract value to points
movement.sf$perc_dist_10km <- terra::extract(perc_dist_10km, movement.sf)[,-1]
movement.sf$perc_dist_5km <- terra::extract(perc_dist_5km, movement.sf)[,-1]
movement.sf$perc_dist_1km <- terra::extract(perc_dist_1km, movement.sf)[,-1]

movement.sf$bd_area_10km <- terra::extract(bd_area_10km, movement.sf)[,-1]
movement.sf$bd_area_5km <- terra::extract(bd_area_5km, movement.sf)[,-1]
movement.sf$bd_area_1km <- terra::extract(bd_area_1km, movement.sf)[,-1]

# if building area is 0, means no building. Turn NA to 0
movement.sf <- movement.sf %>% mutate(bd_area_10km = ifelse(is.na(bd_area_10km), 0, bd_area_10km),
                                      bd_area_5km = ifelse(is.na(bd_area_5km), 0, bd_area_5km),
                                      bd_area_1km = ifelse(is.na(bd_area_1km), 0, bd_area_1km))

# for the test run we filter out the points that with negative perc distance 
# -100 -> If even the minimum distance of 100m leads to the largest cluster
# -200 -> completely uninhabited terrestrial investigation areas
# -300 -> water areas

movement.pd10 <- movement.sf %>% filter(perc_dist_10km >= 0) %>% # 14449 -> 12218
  group_by(ID) %>%
  summarise(perc_dist_10km = mean(perc_dist_10km),
            bd_area_10km = mean(bd_area_10km)) #total ind 1511, after filtering perc dist 1370 remains
movement.pd5 <- movement.sf %>% filter(perc_dist_5km >= 0) %>% # 14449 -> 10461
  group_by(ID) %>%
  summarise(perc_dist_5km = mean(perc_dist_5km),
            bd_area_5km = mean(bd_area_5km)) #total ind 1511, after filtering perc dist 1299 remains
movement.pd1 <- movement.sf %>% filter(perc_dist_1km >= 0) %>% # 14449 -> 2485
  group_by(ID) %>%
  summarise(perc_dist_1km = mean(perc_dist_1km),
            bd_area_1km = mean(bd_area_1km)) #total ind 1511, after filtering perc dist 639 remains

#################################################################
############ calculate other predictor variables #############
#################################################################
# average for each individual 
move10d.var <- move10d %>% group_by(ID) %>%
  summarise(NDVI = mean(NDVI, na.rm=T), 
         HFI = mean(HFI, na.rm=T),
         BodyMass_kg = mean(BodyMass_kg, na.rm = T),
         Longitude = mean(Longitude, na.rm = T), 
         Latitude = mean(Latitude, na.rm = T)) %>%
  rename(lon = Longitude,
         lat = Latitude)

move10d.var <- distinct(move10d %>% dplyr::select(ID, Species, Diet, RBrain_Size, Activity)) %>%
  left_join(move10d.var, by = "ID")

### add species info 
# first manually change some of the names to remove the subspecies
move10d.var <- move10d.var %>% 
  mutate ( Species = if_else(Species %in% c("Ovis canadensis nelsoni", "Ovis canadensis canadensis", "Ovis canadensis californiana"), 
                             "Ovis canadensis", Species)) %>% 
  rename(Binomial = Species) %>%
  left_join(pant, by = "Binomial" ) 

# add info that pentheria did not cover 
a1 <- move10d.var %>%
  filter(Binomial %in% c("Cervus canadensis", "Cervus_canadensis")) %>%
  mutate (Order = "Artiodactyla", Family = "Cervidae", Genus = "Cervus", Species = "canadensis") 

a2 <- move10d.var %>%
  filter(Binomial == "Tragelaphus sylvaticus") %>%
  mutate (Order = "Artiodactyla", Family = "Bovidae", Genus = "Tragelaphus", Species = "sylvaticus") 
                  
move10d.var <- move10d.var %>% 
  filter(!Binomial %in% c("Cervus canadensis", "Cervus_canadensis", "Tragelaphus sylvaticus")) %>%
  rows_append(a1) %>%
  rows_append(a2)
  
### combine all data 
# we prep three dataframe for modeling 
pd10.df <- movement.pd10  %>%
  left_join(move10d.summary, by = "ID") %>%
  left_join(move10d.var, by = "ID") %>%
  st_drop_geometry()
  
pd5.df <- movement.pd5 %>% 
  left_join(move10d.summary, by = "ID") %>%
  left_join(move10d.var, by = "ID")%>%
  st_drop_geometry()

pd1.df <- movement.pd1 %>% 
  left_join(move10d.summary, by = "ID") %>%
  left_join(move10d.var, by = "ID")%>%
  st_drop_geometry()

#################################################################
###  model fitting ##############################################
#################################################################
source("./code/HaversineLMEfunctions.R")
library(nlme)

# function to identify outliers based on standardized residuals
get_clean_dat <- function(m, threshold) {
  dat$std_resid <- residuals(m, type = "normalized")
  dat_clean <- dat[abs(dat$std_resid) <= threshold, ]
}

######################### r = 10km, 95% displacement #############################

dat = pd10.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_10km = log(perc_dist_10km),
                      bd_area_10km = log(bd_area_10km), 
                      dist_95 = log(dist_95))  

# without perc dist or bd area
m10 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
             correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
             random= ~1|Order/Family/Genus/Species,
             control =list(sing.tol=1e-20),
             data = dat)  

dat_clean <- get_clean_dat(m10,3)
m1 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control =list(sing.tol=1e-20, msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)  

# m20 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km,
#                     correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
#                     random= ~1|Order/Family/Genus/Species,
#                     control =list(sing.tol=1e-20),
#                     data = dat)
# 
# dat_clean <- get_clean_dat(m20,3)
# m2 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km,
#           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
#           random= ~1|Order/Family/Genus/Species,
#           data = dat_clean)

m30 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
    correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
    random= ~1|Order/Family/Genus/Species,
    data = dat)

dat_clean <- get_clean_dat(m30,3)
m3 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control = list(msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)

######################### r = 5km, 95% displacement #############################
dat = pd5.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_5km = log(perc_dist_5km),
                      bd_area_5km = log(bd_area_5km),
                      dist_95 = log(dist_95))  

# without perc dist or bd area
m40 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat)  

dat_clean <- get_clean_dat(m40,3)
m4 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control = list(msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)  

m50 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat)

dat_clean <- get_clean_dat(m50,3)
m5 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

m60 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km + bd_area_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat)

dat_clean <- get_clean_dat(m60,3)
m6 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km + bd_area_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

######################### r = 1km, 95% displacement #############################
dat = pd1.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_10km = log(perc_dist_1km),
                      bd_area_10km = log(bd_area_1km),
                      dist_95 = log(dist_95))  

# without perc dist or bd area
m70 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control =list(sing.tol=1e-20),
          data = dat)  

dat_clean <- get_clean_dat(m70,3)
m7 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control =list(sing.tol=1e-20),
          data = dat_clean)  

m80 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat)

dat_clean <- get_clean_dat(m80,3)
m8 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

m90 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km + bd_area_1km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat)

dat_clean <- get_clean_dat(m90,3)
m9 <- lme(dist_95 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

######################### r = 10km, 50% displacement #############################

dat = pd10.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_10km = log(perc_dist_10km),
                      bd_area_10km = log(bd_area_10km), 
                      dist_50 = log(dist_50))  

# without perc dist or bd area
m100 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           control =list(sing.tol=1e-20, msMaxIter = 1000, msMaxEval = 1000),
           data = dat)  
dat_clean <- get_clean_dat(m100,3)
m10 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control =list(sing.tol=1e-20, msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)  

m110 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km,
                    correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
                    random= ~1|Order/Family/Genus/Species,
                    control =list(sing.tol=1e-20),
                    data = dat)
dat_clean <- get_clean_dat(m110,3)
m11 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"),
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

m120 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)
dat_clean <- get_clean_dat(m120,3)
m12 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control = list(msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)

######################### r = 5km, 50% displacement #############################
dat = pd5.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_5km = log(perc_dist_5km),
                      bd_area_5km = log(bd_area_5km),
                      dist_50 = log(dist_50))  

# without perc dist or bd area
m130 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)  
dat_clean <- get_clean_dat(m130,3)
m13 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control = list(msMaxIter = 1000, msMaxEval = 1000),
          data = dat_clean)  

m140 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)
dat_clean <- get_clean_dat(m140,3)
m14 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

m150 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km + bd_area_5km,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)
dat_clean <- get_clean_dat(m150,3)
m15 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_5km + bd_area_5km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

######################### r = 1km, 50% displacement #############################
dat = pd1.df # start with one perc distance scale 

any(duplicated(c(dat$lon, dat$lat)))
dat <- dat[-which(duplicated(c(dat$lon, dat$lat))),] # remove ID with duplicated locations 

# log transform body mass and scale NDVI
dat <- dat %>% mutate(BodyMass_kg = log(BodyMass_kg), 
                      NDVI = scale(NDVI),
                      perc_dist_10km = log(perc_dist_1km),
                      bd_area_10km = log(bd_area_1km),
                      dist_50 = log(dist_50))  

# without perc dist or bd area
m160 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           control =list(sing.tol=1e-20),
           data = dat)  
dat_clean <- get_clean_dat(m160,3)
m16 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          control =list(sing.tol=1e-20),
          data = dat_clean)  

m170 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)
dat_clean <- get_clean_dat(m170,3)
m17 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)

m180 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_1km + bd_area_1km,
           correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
           random= ~1|Order/Family/Genus/Species,
           data = dat)
dat_clean <- get_clean_dat(m180,3)
m18 <- lme(dist_50 ~ BodyMass_kg + NDVI + HFI + Diet + perc_dist_10km + bd_area_10km,
          correlation = corHaversine(form=~lon+lat, mimic = "corGaus"), 
          random= ~1|Order/Family/Genus/Species,
          data = dat_clean)
#################################################################
###  model comparison and summary ###############################
#################################################################
Cand.models.r10km <- list ("original" = m1, "pd + bd" = m3)
aictab(Cand.models.r10km)

Cand.models.r5km <- list ("original" =m4, "pd" = m5, "pd + bd" = m6)
aictab(Cand.models.r5km)

Cand.models.r1km <- list ("original" =m7, "pd" = m8, "pd + bd" = m9)
aictab(Cand.models.r1km)
