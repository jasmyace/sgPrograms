

# deal with the raster
require(raster)
require(rgdal)
require(rgeos)      # dissolve in rasterToPolygons

dataDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data"
rastDir <- paste0(dataDir,"/Spatial/Density Rasters")
polyDir <- paste0(dataDir,"/Spatial/Lek Density Polygons")
zoneDir <- paste0(dataDir,"/Spatial/Management Zones/SG_MgmtZones_ver2_20061018")

scale <- 10

# i need spatial extent of leks.

# i need spatial extent of management zones.

# https://www.sciencebase.gov/catalog/item/4fc68959e4b0f02c1d6a8151

zoneShp <- readOGR(zoneDir,"SG_MgmtZones_ver2_20061018")

# i need each of the rasters. 

mZonesNum <- c(1:7)
mZonesRom <- c("I","II","III","IV","V","VI","VII")

zonesDat <- data.frame(mZonesNum=mZonesNum,mZonesRom=mZonesRom)

rasters <- list("vector",7)
polys <- list("vector",7)
for(i in 1:7){
  thisZone <- as.character(droplevels(zonesDat[zonesDat$mZonesNum == i,]$mZonesRom))
  
  setwd("//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Density Rasters")
  ugh <- paste0("60to85 by 5% Relative Percent Population_MZ ",thisZone,".img")
  rasters[[i]] <- raster(ugh)
  
  yExt <- rasters[[i]]@ncols
  xExt <- rasters[[i]]@nrows
  xScale <- c(rep( floor(xExt / scale), scale - 1), xExt - sum(rep( floor(xExt / scale), scale - 1)))    # make size of x blocks
  yScale <- c(rep( floor(yExt / scale), scale - 1), yExt - sum(rep( floor(yExt / scale), scale - 1)))    # make size of y blocks
  
  tmpRast <- list("vector",scale^2)
  tmpPoly <- list("vector",scale^2)
  for(j in 1:scale){
    
    for(k in 1:scale){
    
      if(j == 1){
        xMin <- 1
        xMax <- xScale[j]    
      } else {
        xMin <- sum(xScale[1:(j - 1)]) + 1
        xMax <- sum(xScale[1:j])     
      }
      
      if(k == 1){
        yMin <- 1
        yMax <- yScale[k]      
      } else {
        yMin <- sum(yScale[1:(k - 1)]) + 1
        yMax <- sum(yScale[1:k])     
      }

      putHere <- (j - 1)*scale + k                                                                         # j is x-dir, k is y-dir
      tmpRast[[putHere]] <- crop(rasters[[i]],extent(rasters[[i]],xMin,xMax,yMin,yMax))
      if(!(tmpRast[[putHere]]@data@min %in% c(Inf,-Inf) | tmpRast[[putHere]]@data@min %in% c(Inf,-Inf))){  # check if rast has other than NA
        tmpPoly[[putHere]] <- rasterToPolygons(tmpRast[[putHere]],dissolve=TRUE)
      }
      
      if(j == 1 & k == 1){
        thePoly <- tmpPoly[[putHere]]
      } else {
        thePoly <- gUnion(thePoly,tmpPoly[[putHere]])
      }
    }
    
    # here!
    
  }
  
      plot(tmpPoly[[putHere]])
      
  
  
  writeOGR(polys[[i]],polyDir,paste0("Density Polygons - Zone - ",thisZone),driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  #   writeOGR(zoneShp,polyDir,paste0("Density Polygons - Zone ",thisZone),driver="ESRI Shapefile",overwrite_layer=TRUE)
}


