# this script runs models at the interspecies level 

library(tidyverse)
library(mapview)
library(ggplot2)
library(sf)
library(geosphere)
library(rgee)
library(terra)

ee_Initialize()
ee_check()

# Function to calculate UTM zone EPSG code
get_utm_crs <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1  # Calculate UTM zone
  if (lat >= 0) {
    epsg <- 32600 + zone  # Northern Hemisphere
  } else {
    epsg <- 32700 + zone  # Southern Hemisphere
  }
  return(epsg)
}


###################################################
# ----- read data and prep data for viz ------ 
###################################################
target_time_scale_days = 10
tartget_spatial_scale = "hi" 

move_covid <- read_rds( paste0("./data/movement/ready_data/move", target_time_scale_days,"d_covid.rds")) %>% mutate(source = "covid")
move_movebank <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_movebank.rds")) %>% mutate(source = "movebank")
move_oomove <- read_rds(paste0("./data/movement/ready_data/move", target_time_scale_days,"d_oodata_MoratoJaguar.rds")) %>% mutate(source = "oodata")

if (tartget_spatial_scale == "hi") {
  move <- rbind(move_covid %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi, source) ,
                move_movebank %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi, source) ,
                move_oomove %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_hi, bd_adpt_hi, dist_2_build, pd_scale_hi, source)) %>%
    rename(pd_adpt = pd_adpt_hi,
           bd_adpt = bd_adpt_hi,
           pd_scale = pd_scale_hi)
} else {
  move <- rbind(move_covid %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me, source),
                move_movebank %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me, source),
                move_oomove %>%
                  dplyr::select(ID, Binomial, Longitude, Latitude, 
                                Displacement_km, TimestampUTC, HFI, HMI,
                                pd_adpt_me, bd_adpt_me, dist_2_build, pd_scale_me, source)) %>%
    rename(pd_adpt = pd_adpt_me,
           bd_adpt = bd_adpt_me,
           pd_scale = pd_scale_me)
}

  
# raw mule deer movement 
# raw_move <- readRDS("./data/movement/Tucker_7704108/Tucker_1h_Move_Spatial.rds")  
# raw_move <- read_csv("./data/movement/other_open_data/Morato_2018_Jaguar/jaguar_movement_data.csv") %>%
#   rename(ID = `individual.local.identifier (ID)`, Longitude = location.long, Latitude = location.lat, 
#          Species = individual.taxon.canonical.name, TimestampUTC = timestamp) %>%
#   mutate(TimestampUTC = mdy_hm(TimestampUTC), ID = paste0(str_sub(study.name, 1,2), ID, "_", Species))
raw_move <- read_csv("./data/movement/open_movebank/midproduct_raw_movebank_data.csv") %>%
  rename( TimestampUTC = timestamp, Longitude = coords_x, Latitude = coords_y) %>%
  mutate( individual_id = as.character(individual_id),
          ID = paste0(str_sub(individual_id, -4,-1), "_", species) )

# we want to visualize animals experience different percolation distances under the same building density 
move.summary <- move %>%
  mutate(bd_adpt = round(bd_adpt,0)) %>%
  filter(ID %in% raw_move$ID) %>%
  group_by(ID, bd_adpt) %>%
  summarise(min_pd = min(pd_adpt),
            max_pd = max(pd_adpt),
            sd_pd = sd(pd_adpt),
            n = n()) %>%
  filter(sd_pd != 0, 
         bd_adpt !=0, 
         !is.na(sd_pd))

ID.list = move.summary$ID

raw <- raw_move %>% filter(ID %in% ID.list) 

# Load the World Settlement Footprint dataset
wsf <- ee$Image('DLR/WSF/WSF2015/v1')

for (id in ID.list[18:length(ID.list)]) {
  
  id = as.character(id)
  
  id.minpd <- move.summary %>% filter(ID == id) %>% pull (min_pd)
  id.maxpd <- move.summary %>% filter(ID == id) %>% pull (max_pd)

  id.disp <- move %>% filter(ID == id) %>% filter(pd_adpt == id.minpd | pd_adpt == id.maxpd)
  
  for (i in (1:nrow(id.disp))) {
    
    disp.id.i <- id.disp[i,]
    
    # get raw movement for the corresponding displacement 
    traj.id.i <- raw %>% filter( ID == id,
                                 TimestampUTC >= disp.id.i$TimestampUTC,
                                 TimestampUTC <= (disp.id.i$TimestampUTC + days(10))) # since the location/time represent the beginning location of each displacement
    
    # convert trajectory to spatial obj
    id.i.pts.sp <- st_as_sf(traj.id.i, coords = c("Longitude", "Latitude"), crs = 4326)
    traj.id.i.sp <- id.i.pts.sp %>%
      arrange(TimestampUTC) %>%
      summarise(do_union = FALSE) %>%
      st_cast("LINESTRING")
    
    disp.id.i.pts <- rbind(
      id.i.pts.sp %>% arrange(TimestampUTC) %>% slice_head(n = 1),
      id.i.pts.sp %>% arrange(TimestampUTC) %>% slice_tail(n = 1)
      )
    disp.id.i.sp <- disp.id.i.pts %>% 
      summarise(do_union = FALSE) %>%
      st_cast("LINESTRING")
    
    
    # get mid point and make a circle that matches the pd scale 
    midpoint <- disp.id.i.pts %>%
      summarise(
        geometry = st_sfc(
          st_point(c(
            mean(st_coordinates(geometry)[, 1]),  # Average longitude
            mean(st_coordinates(geometry)[, 2])   # Average latitude
          ))
        )
      )
    st_crs(midpoint) = st_crs(id.i.pts.sp)
    
    # get the best projected EPSG that fits so that the mapped circle is a circle
    EPSG <- get_utm_crs(st_coordinates(midpoint)[,"X"], st_coordinates(midpoint)[,"Y"])
    
    midpoint_projected <- st_transform(midpoint, crs = EPSG)
    buffer_pd_scale <- st_buffer(midpoint_projected, dist = disp.id.i$pd_scale*1000)
   
    # Clip the dataset to the polygon
    polygon_ee <- sf_as_ee(buffer_pd_scale)
    clipped_wsf <- wsf$clip(polygon_ee)
    # Download the clipped image as a GeoTIFF
    clipped_wsf_tif <- ee_as_rast(
      image = clipped_wsf,
      region = polygon_ee$geometry(),
      scale = 250,  # Set resolution (in meters)
      via = "drive"  # Export via Google Drive
    ) 
    clipped_wsf_tif <- ifel(clipped_wsf_tif == 255, 1, 0)
    writeRaster(clipped_wsf_tif, filename = paste0("./visualization/", disp.id.i$Binomial, "/", id, 
                                              "_bd", round(disp.id.i$bd_adpt,0), 
                                              "_pd",disp.id.i$pd_scale,"km_", disp.id.i$pd_adpt, "_",i, ".tif"),
                overwrite=TRUE)
    
    # visualization 
    p <- ggplot() +
      geom_spatraster(data = clipped_wsf_tif) +
      scale_fill_gradient(low = NA, high = "red") +
      geom_sf(data = buffer_pd_scale, color = "grey", fill = NA, linewidth = 1.2) +
      geom_sf(data = traj.id.i.sp, color = "#404040") +
      geom_sf(data = disp.id.i.sp, color = "orange", linewidth = 1.5) + 
      theme_void() +
      ggtitle(paste0(disp.id.i$ID, " ", target_time_scale_days, "-day", 
                     "\n building density = ", round(disp.id.i$bd_adpt,1), 
                     "\n percolation distance = ", disp.id.i$pd_adpt/1000, 
                     " km (radius = ", disp.id.i$pd_scale, " km)")) +
      theme(legend.position = "none",
            plot.title = element_text(size=10)) 
    
  ggsave(paste0("./visualization/", disp.id.i$Binomial, "/", id,
                "_bd", round(disp.id.i$bd_adpt,0), 
                "_pd",disp.id.i$pd_scale,"km_", disp.id.i$pd_adpt, "_",i, ".png"))
    
  st_write(buffer_pd_scale, paste0("./visualization/", disp.id.i$Binomial, "/", id,"_buffer_", i, ".shp"))
  st_write(traj.id.i.sp, paste0("./visualization/", disp.id.i$Binomial, "/", id,"_traj_", i, ".shp"))
  st_write(disp.id.i.sp, paste0("./visualization/", disp.id.i$Binomial, "/", id,"_disp_", i, ".shp"))
  
    # mapviewOptions(basemaps = c("OpenStreetMap"))
    # mapview(clipped_wsf_tif, na.color = NA, col.regions = c(NA, "red"))+
    #   mapview(disp.id.i.sp, color = "orange", legend = FALSE) +  #displacement
    #   mapview(traj.id.i.sp, col.regions = "grey", legend = FALSE)  + #trajectory
    #   mapview(disp.id.i.pts, col.regions = "red", cex = 4, legend = FALSE) +  #displacement pts
    #   mapview(buffer_pd_scale, legend = FALSE)  #buffer circle
  }
  
}







