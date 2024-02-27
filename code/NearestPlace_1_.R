source("./code/gdal_polygonizeR.R")
library(terra)
library(tidyverse)

# read geotiff 
WSF2015 <- rast("./data/WSF2019/WSF2015_v1_EPSG4326_PercentSettlementArea_100m.tif"). ## <<-- hmm maybe change to 2015 and just use GEE 
names(WSF2015) <- "WSF2015" 

# # read land shapefile 
# Land <- vect("./data/land-polygons-complete-4326/land_polygons.shp") #open street map land data

WSF2015 <- subst(WSF2015, from = 0, to = NA) #make all 0 to NA
hist(WSF2015)

# polygonize rasters 

# generate buffers and use the buffer as mask 

# extracting masked value 