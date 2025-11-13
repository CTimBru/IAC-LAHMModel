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
require(geodata) # Grab elevation data: https://cran.r-project.org/web/packages/geodata/refman/geodata.html
require(future.apply) # Multithreading, not great on Windows 11: https://cran.r-project.org/web/packages/future.apply/refman/future.apply.html
require(tidyverse) # Dataframes: https://www.rdocumentation.org/packages/tidyverse/versions/2.0.0
require(corrplot) # Correlation plots: https://www.rdocumentation.org/packages/corrplot/versions/0.95
library(caret) # Find correlations about |p| > 0.7: https://www.rdocumentation.org/packages/caret/versions/7.0-1

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
BioClim_dir <- ".fullDatasets/Bioclim"

unzip(".fullDatasets/wc2.1_30s_bio.zip", overwrite=TRUE, exdir=BioClim_dir)

#Delete after Cropping

WorldClim_dir <- ".fullDatasets/Worldclim"

unzip(".fullDatasets/wc2.1_30s_prec.zip", overwrite=TRUE, exdir=WorldClim_dir)
unzip(".fullDatasets/wc2.1_30s_srad.zip", overwrite=TRUE, exdir=WorldClim_dir)
unzip(".fullDatasets/wc2.1_30s_wind.zip", overwrite=TRUE, exdir=WorldClim_dir)

#Delete after Cropping

#VII) Terrain data
elevation_dir <- ".fullDatasets/Elevation"
#Download Elevation data in approximate LA Area.
#lon=mean(st_bbox(LA_County_Perimeter)[c(1,3)]),lat=mean(st_bbox(LA_County_Perimeter)[c(2,4)])
elevation_3s(lon=-118.2991,lat=33.78668,path='.fullDatasets')
#https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_13_06.zip wget if the download fails, server close to the internet backbone works.
unzip(".fullDatasets/srtm_13_06.zip", overwrite=TRUE, exdir=elevation_dir)

# VIII) LA County assessor data
# Data: https://data.lacounty.gov/datasets/lacounty::assessor-parcel-data-rolls-2021-present/about
assessor_dir <- ".fullDatasets/Assessor"

#File was sourced manually and is in uncompressed format, no further actions

####
# 2) Combined/cut down data to LA County bounding box
####

# Set output dir
data_dir <- "LA_Data/"

# I) Study area(s)
if(file.exists(paste0(data_dir,"Eaton_WF.shp")) && file.exists(paste0(data_dir,"Palisades_WF.shp"))){
    Eaton_Perimeter <- read_sf(paste0(data_dir,"Eaton_WF.shp"))
    Palisades_Perimeter <- read_sf(paste0(data_dir,"Palisades_WF.shp"))
  } else  {
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
  }

if(file.exists(paste0(data_dir,"LA_County.shp"))){
  LA_County_Perimeter <- read_sf(paste0(data_dir,"LA_County.shp"))
} else  {
# II) County Shapefile
  # Get the shapefile names
  CA_county_shape <- list.files(path=California_county_dir,pattern="*\\.shp$")
  CA_county_shape <- paste0(California_county_dir,"/",CA_county_shape)
  
  # Load the shapefile
  LA_County_Perimeter <- read_sf(CA_county_shape) %>% dplyr::filter(NAME == "Los Angeles")
  LA_County_Perimeter <- st_transform(LA_County_Perimeter,crs=st_crs(4326))
  
  # Output LA County shapefile to data dir
  st_write(LA_County_Perimeter,paste0(data_dir,"LA_County.shp"))

  # var cleanup
  rm(CA_county_shape)
}
# save bounding box
LA_County_bbox <- extent(LA_County_Perimeter)
LA_County_bbox <- st_as_sfc(st_bbox(LA_County_bbox))
st_crs(LA_County_bbox) = 4326

# III) Wildlands-Urban-Interface
if(file.exists(paste0(data_dir,"WUI.tif"))){
  LA_County_WUI <- raster(paste0(data_dir,"WUI.tif"))
} else  {
  # Choose singlethread or multithread
  # ISSUE: Multithread causes tif artifacts (gaps in data, projection+vrt artifact?)
  WUI_multithread <- FALSE
  if (WUI_multithread == TRUE){
    
    WUI_files <- paste0(WUI_dir,"/",list.files(path=WUI_dir,pattern="*\\.tif$"))
    
    # Define function to run in parallel
    WGS84_raster_reprojection <- function(i){
      rasterFile <- WUI_files[i]
      print(paste("Beginning run:",rasterFile))
      print(Sys.time())
      rasterReproj <- projectRaster(raster(rasterFile),crs ="+init=epsg:4326")
      rasterRename <- gsub("WUI_LA/","WUI_LA/WGS84/",rasterFile)
      raster::writeRaster(rasterReproj,rasterRename,overwrite=TRUE)
      print(paste("Ending run:",rasterFile))
      rm(rasterRename,rasterReproj,rasterFile)
      print(Sys.time())
    }
    
    #Set plan for cores (use half available). >2 led to memory issues!
    plan(multisession, workers=2)
    
    future_lapply(1:length(WUI_files),FUN=WGS84_raster_reprojection,future.seed=1)
    
    #Back to singlecore#Back to singleFUN = core
    plan(sequential)
    
    WUI_files <- paste0(WUI_dir,"/WGS84/",list.files(path=paste0(WUI_dir,"/WGS84"),pattern="*\\.tif$"))
    
    combined_WUI <- vrt(WUI_files)
    LA_County_WUI <- crop(combined_WUI,LA_County_Perimeter)
    LA_County_WUI <- mask(LA_County_WUI,LA_County_Perimeter)
    
    # Write the raster
    raster::writeRaster(LA_County_WUI,paste0(data_dir,"WUI.tif"),overwrite=TRUE)
    
  } else {
    # Get WUI files as a list including dir
    WUI_files <- paste0(WUI_dir,"/",list.files(path=WUI_dir,pattern="*\\.tif$"))
    
    # Make of raster set from all of the separate tifs
    combined_WUI <- vrt(WUI_files)
    
    # Project to WSG84: ESPG4326
    ## Warning projecting takes ~20mins. Multithread target
    combined_WUI <- project(combined_WUI,"epsg:4326")
    
    # Crop to LA County
    LA_County_WUI <- crop(combined_WUI,LA_County_Perimeter)
    LA_County_WUI <- mask(LA_County_WUI,LA_County_Perimeter)
    
    # Write the raster
    raster::writeRaster(LA_County_WUI,paste0(data_dir,"WUI.tif"),overwrite=TRUE)
    
    # var cleanup
    rm(WUI_files,combined_WUI)
  }
}

# IV) WHP
if(file.exists(paste0(data_dir,"WHP.tif"))){
  LA_County_WHP <- raster(paste0(data_dir,"WHP.tif"))
} else  {
  # Get WHP file
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
}
# V) Bioclim
# Check the data dir if bioclims are already processed
if(file.exists(paste0(data_dir,"wc2.1_30s_bio_9.tif"))){
  bioclim_rasters <- paste0(data_dir,list.files(path=data_dir,pattern="wc2.1_30s_bio_.*\\.tif$")) %>%
    stack()
} else  {
  bioclim_list <- list.files(path=data_dir,pattern="\\bwc2.1_30s_bio.*")
  if(length(bioclim_list)<19){
    # Set list to source dir
    bioclim_list <- list.files(path=BioClim_dir,pattern="\\bwc2.1_30s_bio.*")
    for(bioclim in bioclim_list){
      #Read in bioclim layers
      bioclim_layer <- raster(paste0(BioClim_dir,"/",bioclim))
      #printing datum, if projection is not EPSG4326 | WSG84 then recode for reprojection!
      print(crs(bioclim_layer, asText = TRUE))
      #n.b. suggested to build a bounding box in that co-ords & then crop before reprojection, due to massive files
      if(grepl("WGS84", crs(bioclim_layer, asText = TRUE))){
        # Crop & Mask to LA
        bioclim_layer <- crop(bioclim_layer,LA_County_Perimeter)
        bioclim_layer <- mask(bioclim_layer,LA_County_Perimeter)
        # Crop bioclim to bbox around LA County
        # Output to LA dir
        raster::writeRaster(bioclim_layer,paste0(data_dir,bioclim),overwrite=TRUE)
      } else {print('Coordinate system requires reprojection to continue!')}
    }
  }
  bioclim_rasters <- paste0(data_dir,list.files(path=data_dir,pattern="wc2.1_30s_bio_.*\\.tif$")) %>%
    stack()
}

#VI) Worldclim
# Check the data dir if worlclims are already processed
if(file.exists(paste0(data_dir,"wc2.1_30s_wind_12.tif"))){
  worldclim_rasters <- paste0(data_dir,c(list.files(path=data_dir,pattern="\\bwc2.1_30s_prec.*"),
                                         list.files(path=data_dir,pattern="\\bwc2.1_30s_srad.*"),
                                         list.files(path=data_dir,pattern="\\bwc2.1_30s_wind.*"))) %>%
    stack()
} else  {
  worldclim_list <- c(list.files(path=data_dir,pattern="\\bwc2.1_30s_prec.*"),
                      list.files(path=data_dir,pattern="\\bwc2.1_30s_srad.*"),
                      list.files(path=data_dir,pattern="\\bwc2.1_30s_wind.*"))
  if(length(worldclim_list)<37){
    # Set list to source dir
    worldclim_list <- c(list.files(path=WorldClim_dir,pattern="\\bwc2.1_30s_prec.*"),
                        list.files(path=WorldClim_dir,pattern="\\bwc2.1_30s_srad.*"),
                        list.files(path=WorldClim_dir,pattern="\\bwc2.1_30s_wind.*"))
    for(worldclim in worldclim_list){
      #Read in bioclim layers
      worldclim_layer <- raster(paste0(WorldClim_dir,"/",worldclim))
      #printing datum, if projection is not EPSG4326 | WSG84 then recode for reprojection!
      print(crs(worldclim_layer, asText = TRUE))
      #n.b. suggested to build a bounding box in that co-ords & then crop before reprojection, due to massive files
      if(grepl("WGS84", crs(worldclim_layer, asText = TRUE))){
        # Crop & Mask to LA
        worldclim_layer <- crop(worldclim_layer,LA_County_Perimeter)
        worldclim_layer <- mask(worldclim_layer,LA_County_Perimeter)
        # Crop bioclim to bbox around LA County
        # Output to LA dir
        raster::writeRaster(worldclim_layer,paste0(data_dir,worldclim),overwrite=TRUE)
      } else {print('Coordinate system requires reprojection to continue!')}
    }
  }
  worldclim_rasters <- paste0(data_dir,c(list.files(path=data_dir,pattern="\\bwc2.1_30s_prec.*"),
                                         list.files(path=data_dir,pattern="\\bwc2.1_30s_srad.*"),
                                         list.files(path=data_dir,pattern="\\bwc2.1_30s_wind.*"))) %>%
    stack()
}

#VII) Terrain data
if(file.exists(paste0(data_dir,"Elevation.tif"))){
  LA_County_Elevation <- raster(paste0(data_dir,"Elevation.tif"))
} else  {
  # Get Elevation file
  Elevation_raster <- raster(paste0(elevation_dir,"/srtm_13_06.tif"))
  # Crop & Mask to LA
  Elevation_raster <- crop(Elevation_raster,LA_County_Perimeter)
  Elevation_raster <- mask(Elevation_raster,LA_County_Perimeter)
  # Output to LA dir
  raster::writeRaster(Elevation_raster,paste0(data_dir,'Elevation.tif'),overwrite=TRUE)
  LA_County_Elevation <- raster(paste0(data_dir,"Elevation.tif"))
}

#Calculate slope & aspect
LA_County_Slope <- terra::terrain(LA_County_Elevation,'slope',units='radians')
LA_County_Aspect <- terra::terrain(LA_County_Elevation,'aspect',units='radians')

#VIII) Assessor
# Check the data dir if Parcel Data is processed
if(file.exists(paste0(data_dir,"ParcelData_2021_12.csv"))){
  LA_County_Parcels <- paste0(data_dir,list.files(path=data_dir,pattern="ParcelData_2021_.*\\.csv$")) %>%
    lapply(read_csv) %>% bind_rows()
  LA_County_Parcels <- LA_County_Parcels[complete.cases(LA_County_Parcels),]
} else  {
  LA_County_Parcels <- read_csv(paste0(assessor_dir,"/Parcel_Data_2021_Table_8468414499611436475.csv"))
  #Property Use Type: SFR- Single Family Residence, NA, VAC- Vacant, C/I- Commercial/Industrial, R-I- ?Multi-family?, CND- Condominium
  #Year Built: Year 50%+ of the original structure was complete
  #Effective Year: Last year of any 'substantial' new construction or major rehabilitation.
  #Total Value: Tax assessed value, this doesn't represent market value. Assessed at sale & capped at 2% per year increase
  #Number of Buildings:
  #Number of Bathrooms:
  #Number of Units:
  #Square Footage:
  LA_County_Parcels <- LA_County_Parcels[c('Property Use Type','Year Built','Effective Year','Total Value','Number of Buildings','Number of Bathrooms','Square Footage','Number of Units','Location Latitude','Location Longitude')] %>% as_tibble(.name_repair = 'universal')
  
  LA_County_Parcels$Property.Use.Type <- as.factor(LA_County_Parcels$Property.Use.Type)
  
  #Split into 5 csv
  rowCount <- nrow(LA_County_Parcels)
  splits <- 12
  rowGroups <- cut(
    seq_len(rowCount),
    breaks = splits,
    labels = FALSE
  )
  split_list <- split(LA_County_Parcels,rowGroups)
  
  i <- 1
  while(i < splits+1){
    write_csv(split_list[[i]],paste0("LA_Data/ParcelData_2021_",i,".csv"))
    i <- i + 1
  }
  LA_County_Parcels <- LA_County_Parcels[complete.cases(LA_County_Parcels),]
  rm(split_list)
}

#IX) Extract data from the datasets based on the assessor points

#convert assessor points into spatial features
assessor_points <- st_as_sf(LA_County_Parcels, coords=c("Location.Longitude","Location.Latitude"), crs = 4326)
rm(LA_County_Parcels)

#Set Palisades & Eaton Fire Perimeter to have a value of 1
Eaton_Perimeter$Eaton <- 1
Palisades_Perimeter$Palisades <- 1

#Extract WUI
extracted_WUI <- raster::extract(LA_County_WUI,assessor_points)
rm(LA_County_WUI)
#dataframe
extracted_WUI <- as.data.frame(extracted_WUI)

#Extract WHP
extracted_WHP <- raster::extract(LA_County_WHP,assessor_points)
rm(LA_County_WHP)
#dataframe
extracted_WHP <- as.data.frame(extracted_WHP)

#Extract bioclim
extracted_bioclim <- raster::extract(bioclim_rasters,assessor_points)
rm(bioclim_rasters)
#dataframe
extracted_bioclim <- as.data.frame(extracted_bioclim)

#Extract worldclim
extracted_worldclim <- raster::extract(worldclim_rasters,assessor_points)
rm(worldclim_rasters)
#dataframe
extracted_worldclim <- as.data.frame(extracted_worldclim)

#Extract Elevation, Slope & Aspect
extracted_elevation <- raster::extract(LA_County_Elevation,assessor_points)
rm(LA_County_Elevation)
extracted_elevation <- as.data.frame(extracted_elevation)

extracted_slope <- raster::extract(LA_County_Slope,assessor_points)
rm(LA_County_Slope)
extracted_slope <- as.data.frame(extracted_slope)

extracted_aspect <- raster::extract(LA_County_Aspect,assessor_points)
rm(LA_County_Aspect)
extracted_aspect <- as.data.frame(extracted_aspect)

#Load in data
#Change N/As in Eaton & Palisades to 0s
assessor_points <- st_join(assessor_points,Eaton_Perimeter['Eaton'])
assessor_points$Eaton[is.na(assessor_points$Eaton)] <- 0
assessor_points <- st_join(assessor_points,Palisades_Perimeter['Palisades'])
assessor_points$Palisades[is.na(assessor_points$Palisades)] <- 0

#Load in WHP & WUI
assessor_points$WUI <- extracted_WUI$extracted_WUI
rm(extracted_WUI)
assessor_points$WHP <- extracted_WHP$extracted_WHP
rm(extracted_WHP)


#load in the bioclim 1 by 1, joins of this size are not advisable.
for(col in names(extracted_bioclim)){
  assessor_points[[col]] <- extracted_bioclim[[col]]
}
rm(col,extracted_bioclim)


#load in the worlclim 1 by 1, joins of this size are not advisable.
for(col in names(extracted_worldclim)){
  assessor_points[[col]] <- extracted_worldclim[[col]]
}
rm(col,extracted_worldclim)

#Factor Eaton & Palisades
assessor_points$Eaton <- as.factor(assessor_points$Eaton)
assessor_points$Palisades <- as.factor(assessor_points$Palisades)



#Keep complete cases & Remove cases where year.built is zero
assessor_points <- na.omit(assessor_points)
assessor_points <- assessor_points[assessor_points$Year.Built != 0,]

#If effective year is zero, set equal to year built
assessor_points <- assessor_points %>%
  mutate(
    Effective.Year = if_else(Effective.Year == 0, Year.Built, Effective.Year)
  )

#VII) Lead Data
#nb Lead data is not a part of this git
sensitive_dir <- RVar_sensitive
sensitive_filename <- RVar_sensitive_filename
read_csv(paste0(RVar_sensitive,RVar_sensitive_filename))


# Correlations check:
temp <- assessor_points[sample(nrow(assessor_points),10000), ] %>% as.data.frame
is_numeric_col <- sapply(temp,is.numeric)
temp <- temp[,is_numeric_col] %>% as.data.frame
temp2 <- cor(temp,method='spearman')
collinear_data <- findCorrelation(abs(temp2),cutoff=0.7,verbose=TRUE)
collinear_data[[length(collinear_data)]] <- 2
corrplot(temp2)


