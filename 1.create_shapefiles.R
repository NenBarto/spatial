#australia-wide basins


library(ggplot2)
library(tidyverse)
library(plyr)
library(gtools)
library(galah)
library(remotes)
library(rgeos)
library(geosphere)
library(geojsonsf)
library(sf)
library(rmapshaper)
#rm(list=ls())

#doing all of australia
#start with geofabric
#subsect it by subdivisions
#however watch out for those over the borders

#borders
states<-st_read("../shapefiles/state_borders/STE_2021_AUST_GDA94.shp",type=1)
nc_geom<-st_geometry(states)

#geofabric
load("../shapefiles/geofabric/_Level2_basin_WKT.RData")
basins<-basin_wkt
sf_use_s2(FALSE)
geofabric <- bind_rows(lapply(basins, function(x) x %>% st_as_sf(crs = "wgs84"))) %>%
  mutate(basin = names(basin_wkt)) %>%
  st_crop(xmin = 110, xmax = 155, ymin = -45, ymax = - 5) %>%
  st_transform(crs = 3112)

  
####### 1. subsect geofabric with Vic basins
sf <- st_read("../shapefiles/vic/BASIN100.shp") %>% st_transform(crs = 3112)
map <- sf::st_as_sf(sf)
map<-map[,"geometry"]
mapSimp<-ms_simplify(map)
geofabricSimp<-ms_simplify(geofabric)

#find intersect
basin_intersect<-st_intersection(geofabric,map)

#join with geofabric
geofabric_vic<-bind_rows(geofabric,basin_intersect)
#geofabric_vicSimp<-ms_simplify(geofabric_vic,weighting = 2,snap_interval = 2000)

####### 2. subsect geofabric with MW basins

sf <- st_read("../shapefiles/mw/HWS2018_Subcatchment_Boundaries.shp") %>% st_transform(crs = 3112)
map <- sf::st_as_sf(sf)
map<-map[,"geometry"]
mapSimp<-ms_simplify(map)
basin_intersect<-st_intersection(geofabric_vic,map)

#join with geofabric
geofabric_vic_MW<-bind_rows(geofabric_vic,basin_intersect)
geofabric_vic_MW_simp<-ms_simplify(geofabric_vic_MW,snap_interval = 500,keep=0.01)
geofabric_vic_MW_simp %>% ggplot() + geom_sf() + coord_sf(crs = "wgs84", xlim = c(143, 147), ylim = c(-37, -39))

#clean from small polygons
geofabric_vic_MW_clean <-geofabric_vic_MW_simp %>% dplyr::mutate(area=st_area(.)) %>% dplyr::filter(area>=set_units(1e+06,m^2))
basins<-geofabric_vic_MW_clean$basin
for(basin in unique(basins)){cat(".");basins[grep(basin,basins,fixed=T)]<-paste0(basins[grep(basin,basins,fixed=T)],".",1:length(grep(basin,basins,fixed=T)))}
geofabric_vic_MW_clean$basin<-basins
st_write(geofabric_vic_MW_clean,"../shapefiles/subcatchments_v1.shp",delete_layer=TRUE)

MWfile<-paste0(sourceDir,"MW_surveys_bySubcatchment.RData")
MWcorrectedfile<-gsub(".RData",".corrected.RData",MWfile)

if(!file.exists(MWcorrectedfile)){
  catchments<-load(paste0(sourceDir,"MW_surveys_bySubcatchment.RData"))
  fishc$Longitude<-gsub(".*\\(","",fishc$geometry)
  fishc$Longitude<-as.numeric(gsub(", .*","",fishc$Longitude))
  fishc$Latitude<-gsub(".* ","",fishc$geometry)
  fishc$Latitude<-as.numeric(gsub("\\)","",fishc$Latitude))
  
  catchments<-data.frame(species=c(birdc$species,fishc$species,frogc$species),
                         catch=c(birdc$catch,fishc$catch,frogc$catch),
                         longitude=c(birdc$Longitude,fishc$Longitude,frogc$Longitude),
                         latitude= c(birdc$Latitude,fishc$Latitude,frogc$Latitude))
  dsf <- sf::st_as_sf(catchments, coords=c("longitude","latitude"), crs=4326)
  map <- sf::st_as_sf(sf)
  #sf_use_s2(FALSE)
  int <- sf::st_intersects(dsf, map)
  
  lengthInt<-sapply(int,length)
  isEmpty<-unlist(lengthInt)
  
  #clean up the missing data
  #make the column to combine data
  catchments$longlat<-paste0(catchments$longitude,catchments$latitude)
  
  #subselect and make unique
  catchmentsShort<-unique(catchments[isEmpty==0,c("longitude","latitude","longlat")])
  
  closest <- list()
  for(i in seq_len(nrow(catchmentsShort))){
    cat(".")
    catchSF<-sf::st_as_sf(catchmentsShort, coords=c("longitude","latitude"), crs=4326)
    closest[[i]] <- map[which.min(
      st_distance(map, catchSF[i,])),]
  }
  catchmentsShort$subcatch_id<-sapply(closest,function(x){x$BASIN_NO})
  catchmentsShort$subcatch<-sapply(closest,function(x){x$BASIN_NAME})
  
  #now fix the non-overlapping points in two goes - first those that overlap
  catchments$subcatch_id[isEmpty==1] <- as.character(map$BASIN_NO[unlist(int)])
  catchments$subcatch[isEmpty==1] <- as.character(map$BASIN_NAME[unlist(int)])
  
  #then those that do not
  catchments$subcatch_id[isEmpty==0] <- plyr::mapvalues(catchments$longlat[isEmpty==0],catchmentsShort$longlat,catchmentsShort$subcatch_id)
  catchments$subcatch[isEmpty==0] <- plyr::mapvalues(catchments$longlat[isEmpty==0],catchmentsShort$longlat,catchmentsShort$subcatch)

  write.csv(catchments,MWcorrectedfile,row.names=F)
} else {
  catchments<-read.csv(MWcorrectedfile)
}

#continue with the regular workflow
catchments<-unique(catchments)
corrected<-read.csv(paste0(sourceDir,"fishfrogbird.corrections.csv"),header=F)
colnames(corrected)<-c("MW","correctedSpecies")
catchments<-merge(catchments,corrected,by.x="species",by.y="MW")



load(paste0("../sources/VBA_by_basin.RData"))
vba$correctedSpecies<-gsub(" \\(.*","",vba$SCIENTIFIC_NME)
vba$correctedSpecies<-gsub("(.*) (.*) (.*)","\\1 \\2",vba$correctedSpecies)
vba$correctedSpecies<-gsub(" \\(.*","",vba$correctedSpecies)

inVBAnotInDB<-unique(vba$correctedSpecies[!vba$correctedSpecies %in% c(currentDB$scientific_name,multiDF$species)])
countNames<-sapply(inVBAnotInDB,function(x){grep(" ",x)})==1
doubleNameMissingVBA<-inVBAnotInDB[countNames]
doubleNameMissingVBA<-doubleNameMissingVBA[!is.na(doubleNameMissingVBA)]
doubleNameMissingVBA<-doubleNameMissingVBA[!grepl("sp.",doubleNameMissingVBA)]
write.csv(doubleNameMissingVBA,"../tables/VicBasinVBAspeciesMissingFromVertDB.csv",row.names=F)

vba<-as.data.frame(vba)
vba$longitude<-vba$LONG_DD94
vba$latitude<-vba$LAT_DD94

VBAfile<-paste0(sourceDir,"VBA_surveys_bySubcatchment.RData")
VBAcorrectedfile<-gsub(".RData",".corrected.RData",VBAfile)


if(!file.exists(VBAcorrectedfile)){
  
  dsf <- sf::st_as_sf(vba, coords=c("longitude","latitude"), crs=4326)
  map <- sf::st_as_sf(sf)
  #sf_use_s2(FALSE)
  int <- sf::st_intersects(dsf, map)
  
  lengthInt<-sapply(int,length)
  isEmpty<-unlist(lengthInt)
  
  #clean up the missing data
  #make the column to combine data
  vba$longlat<-paste0(vba$longitude,vba$latitude)
  
  #subselect and make unique
  catchmentsShort<-unique(vba[isEmpty==0,c("longitude","latitude","longlat")])
  
  closest <- list()
  for(i in seq_len(nrow(catchmentsShort))){
    cat(".")
    catchSF<-sf::st_as_sf(catchmentsShort, coords=c("longitude","latitude"), crs=4326)
    #if the distance is less than a KM from the border of victoria
    if(as.numeric(min(st_distance(map, catchSF[i,])))<2000){
      closest[[i]] <- map[which.min(
        st_distance(map, catchSF[i,])),]
    } else {
      cat(paste0("\ntoo far: ",catchmentsShort$latitude[i],"\t",catchmentsShort$longitude[i],"\n"))
      closest[[i]] <- data.frame(NULL)
    }
    
  }
  isNull<-sapply(closest,function(x){sum(nrow(x))==0})
  #eliminate sites outside Vic
  closest<-closest[!isNull]
  catchmentsShort<-catchmentsShort[!isNull,]
  catchmentsShort$subcatch_id<-sapply(closest,function(x){x$BASIN_NO})
  catchmentsShort$subcatch<-sapply(closest,function(x){x$BASIN_NAME})
  
  #now fix the non-overlapping points in two goes - first those that overlap
  vba$subcatch_id[isEmpty==1] <- as.character(map$BASIN_NO[unlist(int)])
  vba$subcatch[isEmpty==1] <- as.character(map$BASIN_NAME[unlist(int)])
  
  #then those that do not
  vba$subcatch_id[isEmpty==0] <- plyr::mapvalues(vba$longlat[isEmpty==0],catchmentsShort$longlat,catchmentsShort$subcatch_id)
  vba$subcatch[isEmpty==0] <- plyr::mapvalues(vba$longlat[isEmpty==0],catchmentsShort$longlat,catchmentsShort$subcatch)
  
  write.csv(vba,VBAcorrectedfile,row.names=F)
} else {
  vba<-read.csv(MWcorrectedfile)
}










