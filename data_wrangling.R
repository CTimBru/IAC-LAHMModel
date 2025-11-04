rm(list=ls())
# Plan: To load & process the data for our select study areas: https://egis-lacounty.hub.arcgis.com/maps/ad51845ea5fb4eb483bc2a7c38b2370c/about
# Additional Plan: Prepare to extend to the whole LA County.

# Libraries
require(data.table) #CSV & data processing: https://www.rdocumentation.org/packages/data.table/versions/1.8.4/topics/data.table
require(raster) #tif & spatial processing: https://rspatial.org/raster/pkg/index.html
require(sf) #For geospatial encoding: https://r-spatial.github.io/sf/reference/index.html
require(terra) # For spatial analysis: https://rspatial.org/pkg/index.html
require(fs) # File tools: https://www.rdocumentation.org/packages/fs/versions/1.6.6
#require(ff) # File tools: https://www.rdocumentation.org/packages/ff/versions/1.0-1
require(geodata) # Grab bioclim data for USA only: https://cran.r-project.org/web/packages/geodata/refman/geodata.html
require(dplyr) # Easier shapefile manipulation: https://dplyr.tidyverse.org/reference/index.html


# Load users working directory
source("rvar/var.R")

# Set random number string
set.seed(1)

# Set working directory
# RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

####
# 1) Initial unzipping & sorting
####

# I) Study area(s): https://egis-lacounty.hub.arcgis.com/maps/ad51845ea5fb4eb483bc2a7c38b2370c/about
# Define study area boundary
# NOTE: The fulldatasets need to be sourced and placed in the .fullDataset folder
Eaton_dir <- ".fullDatasets/Eaton_Perimeter"
Palisades_dir <- ".fullDatasets/Palisades_Perimeter"
# Unzip the shapefiles
unzip(".fullDatasets/Eaton_Perimeter_20250121.zip", overwrite=TRUE, exdir=Eaton_dir)
unzip(".fullDatasets/Palisades_Perimeter_20250121.zip", overwrite=TRUE, exdir=Palisades_dir)

# II) California area: https://data.ca.gov/dataset/ca-geographic-boundaries
# Additionally load a California shapefile for some preprocessing on US-sized datasets
California_county_dir <- ".fullDatasets/CA_Counties"
# Unzip the shapefiles
unzip(".fullDatasets/ca_counties.zip", overwrite=TRUE, exdir=California_county_dir)

# III) Wildland-Urban Interface (WUI): https://silvis.forest.wisc.edu/globalwui/
WUI_dir <- ".fullDatasets/WUI_LA"
# Unzip the 
unzip(".fullDatasets/NA.zip", overwrite=TRUE, exdir=WUI_dir)
# zip has a folder NA inside
WUI_NA_dir <- ".fullDatasets/WUI_LA/NA"

# WUI is in Equi7 coordinates, convert: https://epsg.io/transform#s_srs=4326&t_srs=27705&x=NaN&y=NaN
# LA County bounding Box -> Equi7 + floor() ceiling(): https://anthonylouisdagostino.com/bounding-boxes-for-all-us-counties/
# -118.951721	32.75004	-117.646374	34.823302
# -119, 32, -117, 35
# 6238766.765700199  3004898.593003723, 6412629.719152724, 3198038.1350843543
# -> X0062,Y0030, X0065,Y0032
x_coord <- 62:65
y_coord <- 66:69

#We need to find all of the tifs relevant to these
k <- 1
wui_folder_list <- array()
for(i in x_coord){
  x_component <- paste0("X00",i)
  for(j in y_coord){
    wui_folder <- paste0(x_component,"_Y00",j)
    wui_folder_list[k] <- wui_folder
    k <- k+1
  }
}
print(wui_folder_list)

# Move & rename the tifs within the bounding box
for(folder in wui_folder_list){
  file_move(paste0(WUI_NA_dir,"/",folder,"/WUI.tif"),paste0(WUI_dir,"/WUI-",folder,".tif"))
}

# Delete the unnecessary WUI files
unlink(WUI_NA_dir,recursive = TRUE)

# var cleanup
rm(WUI_NA_dir,x_coord,y_coord,i,j,k,folder,wui_folder_list)

# IV) Wildfire Hazard Potential (WHP): https://usfs-public.box.com/shared/static/et4mghz8sq0kxsk2ag6fpey5uec4f8fn.zip
WHP_dir <- ".fullDatasets/WHP_CA"
unzip(".fullDatasets/RDS-2020-0016-2_California.zip", overwrite=TRUE, exdir=WHP_dir)
# ovir files are image pyramids. Efficient for switching between scales.
WHP_files <- list.files(path=paste0(WHP_dir,"/CA"),pattern="WHP_*")
for(WHP_file in WHP_files){
  file_move(paste0(WHP_dir,"/CA/",WHP_file),paste0(WHP_dir,"/",WHP_file))
}

#Delete the unnecessary RDS files
unlink(paste0(WHP_dir,"/CA"),recursive = TRUE)
unlink(paste0(WHP_dir,"/WHP_CA.tif.ovr"))

# var cleanup
rm(WHP_file,WHP_files)

# V) Bioclim BClim: https://www.worldclim.org/data/worldclim21.html
# USA only bioclim: https://geodata.ucdavis.edu/climate/worldclim/2_1/tiles/iso/USA_wc2.1_30s_bio.tif
# worldclim_country(var="bio",country="USA",res=0.5,path=".fullDatasets/") path gathered from geodata library
BioClim_dir <- ".fullDatasets/Bioclim_USA"

#File was sourced manually and is in uncompressed format, no further actions

# VI) LA County assessor data
# Data: https://data.lacounty.gov/datasets/lacounty::assessor-parcel-data-rolls-2021-present/about
assessor_dir <- ".fullDatasets/Assessor"

#File was sourced manually and is in uncompressed format, no further actions

####
# 2) Combined/cut down data to LA County bounding box
####

# Set output dir
data_dir <- "LA_Data/"

# I) Study area(s)
# Get the shapefile names
Eaton_shape <- list.files(path=Eaton_dir,pattern="*\\.shp$")
Eaton_shape <- paste0(Eaton_dir,"/",Eaton_shape)
Palisades_shape <- list.files(path=Palisades_dir,pattern="*\\.shp$")
Palisades_shape <- paste0(Palisades_dir,"/",Palisades_shape)

# Load the shapefiles
Eaton_Perimeter <- read_sf(Eaton_shape)
Palisades_Perimeter<- read_sf(Palisades_shape)
# Reproject to WSG84.
# 4326 is WGS84: https://spatialreference.org/ref/epsg/4326/
Eaton_Perimeter <- st_transform(Eaton_Perimeter,crs=st_crs(4326))
Palisades_Perimeter <- st_transform(Palisades_Perimeter,crs=st_crs(4326))

# Output Eaton & Palisades shapefile to data dir
st_write(Eaton_Perimeter,paste0(data_dir,"Eaton_WF.shp"))
st_write(Palisades_Perimeter,paste0(data_dir,"Palisades_WF.shp"))

# var cleanup
rm(Eaton_shape,Palisades_shape)

# II) County Shapefile
# Get the shapefile names
CA_county_shape <- list.files(path=California_county_dir,pattern="*\\.shp$")
CA_county_shape <- paste0(California_county_dir,"/",CA_county_shape)

# Load the shapefile
LA_County_Perimeter <- read_sf(CA_county_shape) %>% dplyr::filter(NAME == "Los Angeles")
LA_County_Perimeter <- st_transform(LA_County_Perimeter,crs=st_crs(4326))

# Output LA County shapefile to data dir
st_write(LA_County_Perimeter,paste0(data_dir,"LA_County.shp"))
# save bounding box
LA_County_bbox <- extent(LA_County_Perimeter)
LA_County_bbox <- st_as_sfc(st_bbox(LA_County_bbox))
st_crs(LA_County_bbox) = 4326

# var cleanup
rm(CA_county_shape)

# III) Wildlands-Urban-Interface
# Get WUI files as a list including dir
WUI_files <- paste0(WUI_dir,"/",list.files(path=WUI_dir,pattern="*\\.tif$"))

# Make of raster set from all of the separate tifs
combined_WUI <- vrt(WUI_files)

# Project to WSG84: ESPG4326
## Warning projecting takes ~20mins. Multithread target
combined_WUI <- project(combined_WUI,"epsg:4326")

# Crop to LA County
LA_County_WUI <- crop(combined_WUI,LA_County_Perimeter)
LA_County_WUI <- mask(combined_WUI,LA_County_Perimeter)

#Crop to less nuts bounding box
LA_County_WUI <- crop(LA_County_WUI,LA_County_bbox)

# Write the raster
raster::writeRaster(LA_County_WUI,paste0(data_dir,"WUI.tif"),overwrite=TRUE)

# var cleanup
rm(WUI_files,combined_WUI)

# IV) WHP
WHP_CA <- raster(paste0(WHP_dir,"/WHP_CA.tif"))

# Large tif, lets cut this down
# Convoluted way to get the bbox from WSG84 -> AEA
AEA_bbox <- st_transform(LA_County_bbox,crs=crs(WHP_CA))
AEA_df <- cbind(st_sf(geometry=AEA_bbox),do.call(rbind,lapply(AEA_bbox, st_bbox)))

#Crop & Mask WHP_CA before reprojecting
LA_County_WHP <- crop(WHP_CA,st_bbox(c(xmin=AEA_df[[1]],xmax=AEA_df[[3]],ymin=AEA_df[[2]],ymax=AEA_df[[4]])))

# Reproject to WSG84: EPSG4326
LA_County_WHP <- projectRaster(LA_County_WHP,crs="epsg:4326")

# Crop to LA County
LA_County_WHP <- crop(LA_County_WHP,LA_County_Perimeter)
LA_County_WHP <- mask(LA_County_WHP,LA_County_Perimeter)

# Write the raster
raster::writeRaster(LA_County_WHP,paste0(data_dir,"WHP.tif"),overwrite=TRUE)


# var cleanup
rm(AEA_bbox,AEA_df,WHP_CA)

# V) Bioclim
BioClim_USA <- brick(paste0(BioClim_dir,"/USA_wc2.1_30s_bio.tif"))
BioClim_USA <- unstack(BioClim_USA)


#Cali Bounding Box + floor() ceiling()
#-124.41060660766607,32.5342307609976,-114.13445790587905,42.00965914828148: https://observablehq.com/@rdmurphy/u-s-state-bounding-boxes
#-125,32,-114,42