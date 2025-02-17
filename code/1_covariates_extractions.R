# this script extrat environmental covariates form movement data. 
# previous movement data has columns: ID, first_point, mid_point, last_point,first_timestamp, last_timestamp, 
# Displacement_km, Binomial,study_name, contact

library(tidyverse)
library(sf)
library(geojsonio)
library(rgee)
library(terra)

#ee_Initialize()
# ee_Authenticate() 
#ee_check()

#setwd("/Users/wenjing/Senckenberg Dropbox/Wenjing Xu/CriticaL/")

target_time_scale_days = 10
target_set = "oodata" # "covid", "movebank", "oodata" 

move = read_rds(paste0(
  "./data/movement/midproduct/",
  list.files("./data/movement/midproduct", 
             pattern = paste0("move", target_time_scale_days,"d.+",target_set)))) %>% 
  mutate(unique_serial = row_number(.))

# source = "MoratoJaguar"
###############################################################################
####  Environmental covariates extraction  ------------------------------------   
###############################################################################


#### step 1: extract NDVI from GEE ----------------------------

# to each start and end point and take an average for each displacement segment
env.sf <- move  %>%
  dplyr::select(unique_serial, first_point, last_point, first_timestamp, last_timestamp) %>%
  pivot_longer(cols = c(first_point, last_point), names_to = "point", values_to = "coordinates" ) %>%
  mutate(Longitude = coordinates[,1], Latitude = coordinates[,2],
         TimestampUTC = ifelse (point == "first_point", first_timestamp, last_timestamp) ) %>%
  dplyr::select(-coordinates, -first_timestamp, -last_timestamp) %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = st_crs(4326))

# # Load MODIS NDVI collection
# modis_ndvi <- ee$ImageCollection("MODIS/006/MOD13Q1")$select("NDVI")
#
# # extract image dates from the image collection to later combine with points
# modis_dates = ee_get_date_ic(modis_ndvi)$time_start # n = 529
# 
# # match modis date with point date
# env.sf <- env.sf %>%
#   mutate(modis_date = as.character(modis_dates[findInterval(TimestampUTC, modis_dates)]))
# 
# # Loop through each group and extract NDVI
# ndvi_list <- tibble()
# milestones <- seq(10, 90, by = 10)
# last_reported_milestone <- 0
# 
# for (i in as.character(unique(env.sf $modis_date))[1:length(unique(env.sf $modis_date))] ) {
# 
#   points_group = env.sf %>% filter(modis_date == i)
# 
#   n = ceiling(nrow(points_group)/4999)
# 
#   # above limit when there are more than 5000 points to be extracted at once.
#   if (n > 1) {
#     points_group$Group <- rep(1:n, each = 4999, length.out = nrow(points_group))
#   } else {points_group$Group <- 1}
# 
#   for (ii in 1:n) {
# 
#     points_group_ii <- points_group %>% filter(Group == ii)
# 
#     filtered_modis <- modis_ndvi$
#     filterDate(as.character(i))$
#     median()  # reduce from IC to image
# 
#   ndvi_values <- ee_extract(
#     x = filtered_modis,
#     y = points_group_ii["unique_serial"],
#     scale = 250,  # MODIS resolution
#     fun = ee$Reducer$mean()    # Aggregation function (not critical for points)
#   )
# 
#   ndvi_list <- rbind(ndvi_list, ndvi_values)
# 
#   # progress tracking
#   percentage_completed  = which(as.character(unique(env.sf$modis_date)) == i)/length(unique(env.sf $modis_date))*100
# 
#     # Check if the percentage has reached the next milestone
#     for (milestone in milestones) {
#       if (percentage_completed >= milestone && milestone > last_reported_milestone) {
#         message(sprintf("%.0f%% is completed", milestone))
#         last_reported_milestone <- milestone  # Update the last reported milestone
#       }
#     }
#   }}
# 
# # merge NDVI with other displacement info
# move <- move %>%
#   left_join(ndvi_list %>% group_by(unique_serial) %>%
#               summarise (NDVI = mean(NDVI)), # all environmental covariates for each segment are the mean of the start and the end point values
#             by = "unique_serial")
# 
# ## step 2: extract HFI and HMI and distance to building --------------------------------------
# HFI <- rast("./data/Human_footprint_index/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint-4326.tif")
# HMI <- rast("./data/human_modification_index/gHM-4326.tif")
# build_dist <- rast("./data/Percolation_distances/geodetic_distance/window_distance_d1km.tif")
# 
# env.sf$HFI <- terra::extract(HFI, env.sf)[,-1]
# env.sf$HMI <- terra::extract(HMI, env.sf)[,-1]
# env.sf$dist_2_build <- terra::extract(build_dist, env.sf)[,-1]
# 
# # merge HFI, and HMI back with displacement df
# move <- move %>%
#   left_join(env.sf %>% st_drop_geometry(.) %>%
#               group_by(unique_serial) %>%
#               summarise (HFI = mean(HFI),
#                          HMI = mean(HMI),
#                          dist_2_build = mean(dist_2_build)), # all environmental covariates for each segment are the mean of the start and the end point values
#             by = "unique_serial")

####### step 3: extract percolation distance, building density ---------------------
# because both derived from a moving window approach (search radii), using mid point of the displacement to extract 

# read perc dist layers
perc_dist_1km <- rast("./data/Percolation_distances/perc_dists_202410/radius_1km/d1km_r1km_results_global_fin.tif")
perc_dist_2p5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_2p5km/d1km_r2p5km_results_global_fin.tif")
perc_dist_5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_5km/d1km_r5km_results_global_fin.tif")
perc_dist_7p5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_7p5km/d1km_r7p5km_results_global_fin.tif")
perc_dist_10km <- rast("./data/Percolation_distances/perc_dists_202410/radius_10km/d1km_r10km_results_global_fin.tif")
perc_dist_15km <- rast("./data/Percolation_distances/perc_dists_202410/radius_15km/d1km_r15km_results_global_fin.tif")

# read building density layers
build_dens_1km <- rast("./data/Percolation_distances/perc_dists_202410/radius_1km/d1km_r1km_results_global_overlap.tif")
build_dens_2p5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_2p5km/d1km_r2p5km_results_global_overlap.tif")
build_dens_5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_5km/d1km_r5km_results_global_overlap.tif")
build_dens_7p5km <- rast("./data/Percolation_distances/perc_dists_202410/radius_7p5km/d1km_r7p5km_results_global_overlap.tif")
build_dens_10km <- rast("./data/Percolation_distances/perc_dists_202410/radius_10km/d1km_r10km_results_global_overlap.tif")
build_dens_15km <- rast("./data/Percolation_distances/perc_dists_202410/radius_15km/d1km_r15km_results_global_overlap.tif")

###  extract adaptive pd values ##########################################

# for each species, summarize median and 95% displacement
mobility_summary <- 
  move %>%
  group_by(Binomial) %>% summarise(
    displacement_me = median (Displacement_km),
    displacement_hi = quantile(Displacement_km, 0.95)
  ) 

# turn mid point data frame to spatial 
move.mid.sf <- move %>% mutate(mid_lon = move$mid_point[,1], 
                               mid_lat = move$mid_point[,2]) %>%
  st_as_sf(., coords = c("mid_lon", "mid_lat"), crs = 4326)

# adaptive pd based on displacement scale 
move.pd = data.frame()
for (i in unique(move.mid.sf$Binomial)) {
  move.mid.sf.spp.i <- move.mid.sf %>% filter(Binomial == i)
  disp_range_hi <- (mobility_summary %>% filter(Binomial == i))$displacement_hi
  disp_range_me <- (mobility_summary %>% filter(Binomial == i))$displacement_me
  
  # first extract pd that matches long movements 
  if (disp_range_hi < 2 ) {
    move.mid.sf.spp.i$pd_scale_hi = 1
    move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_1km, move.mid.sf.spp.i)[,-1]
    move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_1km, move.mid.sf.spp.i)[,-1]
  } else {
    if (disp_range_hi < 5 ) {
      move.mid.sf.spp.i$pd_scale_hi = 2.5
      move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_2p5km, move.mid.sf.spp.i)[,-1]
      move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_2p5km, move.mid.sf.spp.i)[,-1]
    } else {
      if (disp_range_hi < 10 ) {
        move.mid.sf.spp.i$pd_scale_hi = 5
        move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_5km, move.mid.sf.spp.i)[,-1]
        move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_5km, move.mid.sf.spp.i)[,-1]
      } else {
        if (disp_range_hi < 15 ) {
          move.mid.sf.spp.i$pd_scale_hi = 7.5
          move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_7p5km, move.mid.sf.spp.i)[,-1]
          move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_7p5km, move.mid.sf.spp.i)[,-1]
        } else {
          if (disp_range_hi < 20) {
            move.mid.sf.spp.i$pd_scale_hi = 10
            move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_10km, move.mid.sf.spp.i)[,-1]
            move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_10km, move.mid.sf.spp.i)[,-1]
          } else {
            move.mid.sf.spp.i$pd_scale_hi = 15 
            move.mid.sf.spp.i$pd_adpt_hi <- terra::extract(perc_dist_15km, move.mid.sf.spp.i)[,-1]
            move.mid.sf.spp.i$bd_adpt_hi <- terra::extract(build_dens_15km, move.mid.sf.spp.i)[,-1]
          }
        }
      }
    }
  }
  
  # then extract pd that matches median movements 
  if (disp_range_me < 2 ) {
    move.mid.sf.spp.i$pd_scale_me = 1
    move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_1km, move.mid.sf.spp.i)[,-1]
    move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_1km, move.mid.sf.spp.i)[,-1]
  } else {
    if (disp_range_me < 5 ) {
      move.mid.sf.spp.i$pd_scale_me = 2.5
      move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_2p5km, move.mid.sf.spp.i)[,-1]
      move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_2p5km, move.mid.sf.spp.i)[,-1]
    } else {
      if (disp_range_me < 10 ) {
        move.mid.sf.spp.i$pd_scale_me = 5
        move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_5km, move.mid.sf.spp.i)[,-1]
        move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_5km, move.mid.sf.spp.i)[,-1]
      } else {
        if (disp_range_me < 15 ) {
          move.mid.sf.spp.i$pd_scale_me = 7.5
          move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_7p5km, move.mid.sf.spp.i)[,-1]
          move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_7p5km, move.mid.sf.spp.i)[,-1]
        } else {
          if (disp_range_me < 20) {
            move.mid.sf.spp.i$pd_scale_me = 10
            move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_10km, move.mid.sf.spp.i)[,-1]
            move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_10km, move.mid.sf.spp.i)[,-1]
          } else {
            move.mid.sf.spp.i$pd_scale_me = 15 
            move.mid.sf.spp.i$pd_adpt_me <- terra::extract(perc_dist_15km, move.mid.sf.spp.i)[,-1]
            move.mid.sf.spp.i$bd_adpt_me <- terra::extract(build_dens_15km, move.mid.sf.spp.i)[,-1]
          }
        }
      }
    }
  }
  
  move.pd <- rbind (move.pd , move.mid.sf.spp.i)
}

## ---- fill NAs -------------------------
## fill NAs using the diameter values 
move.mid.sf <-  move.pd  %>% mutate(
  pd_adpt_hi = ifelse(pd_adpt_hi == -200, pd_scale_hi*1000*2,pd_adpt_hi),
  pd_adpt_me = ifelse(pd_adpt_me == -200, pd_scale_me*1000*2,pd_adpt_me)
)

## turn NA building density to 0 ----- ##
move.mid.sf <-  move.mid.sf  %>% mutate(
  bd_adpt_hi = ifelse(is.na(bd_adpt_hi), 0, bd_adpt_hi),
  bd_adpt_me = ifelse(is.na(bd_adpt_me), 0, bd_adpt_me)
)

move <- move.mid.sf %>% st_drop_geometry()

rm(move.pd, move.mid.sf)


# ###### step 4: extract taxonomy info from PANTHERIA -----------------   
# pant <- read_csv("./data/PanTHERIA_1-0_WR05_Aug2008.csv")
# pant <- pant[,1:5] %>% 
#   rename (Order = MSW05_Order,
#           Family = MSW05_Family,
#           Genus = MSW05_Genus,
#           Species = MSW05_Species,
#           Binomial = MSW05_Binomial) 
# 
# # find out which species are not listed in panethria database 
# move  %>% filter(!Binomial %in% pant$Binomial) %>% dplyr::select(Binomial) %>% distinct() 
# #Sapajus macrocephalus is also not included in the pantheria databse. 
# 
# # # for movebank, add Sapajus apella (of which Sapajus macrocephalus is a subspecies of)
# if  (target_set == "movebank") {
#   pant = rbind(pant,
#                data.frame(Order = c("Primates"),
#                           Family = c("Cebidae"),
#                           Genus = c("Sapajus"),
#                           Species = c("apella macrocephalus"),
#                           Binomial = c("Sapajus macrocephalus")))
# } 
# 
# # # for covid dataset - add elk, cape bushbuck --------
# if (target_set == "covid") {
#   pant = rbind(pant,
#                data.frame(Order = c("Artiodactyla", "Artiodactyla"),
#                           Family = c("Cervidae", "Bovidae"),
#                           Genus = c("Cervus", "Tragelaphus"),
#                           Species = c("canadensis", "sylvaticus"),
#                           Binomial = c("Cervus canadensis", "Tragelaphus sylvaticus")))
#   
#   # treat Ovis canadensis nelsoni the same as Ovis canadensis
#   move <- move %>% mutate(Binomial = case_when(Binomial == "Ovis canadensis nelsoni" ~ "Ovis canadensis",
#                                                Binomial == "Ovis canadensis canadensis" ~ "Ovis canadensis",
#                                                Binomial == "Ovis canadensis californiana" ~ "Ovis canadensis",
#                                                Binomial == "Cervus_canadensis" ~ "Cervus canadensis",
#                                               .default = Binomial))
# }
# 
# move  <- move %>% left_join(pant, by = "Binomial")
# 
# # check again
# move  %>% filter(!Binomial %in% pant$Binomial) %>% dplyr::select(Binomial) %>% distinct() 
# 
# 
# ###### step 5: extract species, Diet, BodyMass_kg info -----------------   
# diet <- read_csv("./data/EltonTraits/MamFuncDat.csv") %>% filter(Scientific %in% unique(move$Binomial))
# # check which species not covered in the database
# length(unique(move$Binomial)) == sum(unique(move$Binomial) %in% unique(diet$Scientific))
# # unique(move$Binomial) [which(!unique(move$Binomial) %in% unique(diet$Scientific))]
# 
# diet <- diet %>% mutate(Diet = 
#                           case_when(
#                             `Diet-Inv` + `Diet-Vend` + `Diet-Vect` + `Diet-Vfish` + `Diet-Vunk` + `Diet-Scav` == 0  ~ "Herbivore",
#                             `Diet-Fruit` + `Diet-Nect` + `Diet-Seed` + `Diet-PlantO` == 0  ~ "Carnivore",
#                             .default = "Omnivore"
#                           ),
#                         BodyMass_kg = `BodyMass-Value`/1000) %>%
#   dplyr::select(Scientific, Diet, BodyMass_kg) %>%
#   rename(Binomial = Scientific)
# 
# move <- move %>% left_join(diet, by = "Binomial")
# 
# # which one has no diet info 
# unique(move %>% filter(is.na(Diet)) %>% pull(Binomial) )
# # mannually adding 
# if (target_set == "covid") {
#   move <- move %>% mutate(Diet = case_when(Binomial %in% c("Cervus canadensis", "Tragelaphus sylvaticus") ~ "Herbivore",
#                                                .default = Diet),
#                           BodyMass_kg = case_when(Binomial == "Cervus canadensis" ~ 300.51, # number from Tucker covid dataset
#                                                   Binomial == "Tragelaphus sylvaticus" ~ 37.50, # number from Tucker covid dataset 
#                                                   .default = BodyMass_kg))
# }
# 
# # check again
# unique(move %>% filter(is.na(Diet)) %>% pull(Binomial) )

###### step 6: final data organization  -----------------   
# move <- move %>%                    
#   mutate(Longitude = first_point[,1], Latitude = first_point[,2], # follow previous practice, only keep first location
#          TimestampUTC = first_timestamp) %>%
#   dplyr::select(-unique_serial, - first_point, -mid_point, -first_timestamp, -last_timestamp, -last_point)

#write_rds(movebank, "./data/movement/move10d_movebank.rds")
if (target_set == "oodata") {
  write_rds(move, paste0("./data/movement/midproduct/move",target_time_scale_days,"d_", target_set,"_",source,"firstmidlastpoints.rds"))
} else {
  write_rds(move, paste0("./data/movement/midproduct/move",target_time_scale_days,"d_", target_set,"firstmidlastpoints.rds"))
}

