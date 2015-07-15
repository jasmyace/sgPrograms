makeDensityPolyBits <- function(dat,scale)
  
# dataDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data"
# rastDir <- paste0(dataDir,"/Spatial/Density Rasters")
# polyDir <- paste0(dataDir,"/Spatial/Density Polygons")
# tempDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Density Polygons/tmp"
# zoneDir <- paste0(dataDir,"/Spatial/Management Zones/SG_MgmtZones_ver2_20061018")

  dat <- datList[[1]]
# scale <- 10

  rasters <- list("vector",7)
  polys <- list("vector",7)
  the0 <- list("vector",7)
  the1 <- list("vector",7)
  the2 <- list("vector",7)
  the3 <- list("vector",7)
  the4 <- list("vector",7)
  the5 <- list("vector",7)
  the6 <- list("vector",7)
  
  mZones <- unique(dat$mZoneNum)
  mZones <- mZones[order(mZones)]

  for(i in 1:7){
    if(i <= 7){
      thisZone <- mZones[i]
      
      thisZoneRoman <- as.character(droplevels(dat[dat$mZoneNum == thisZone,][1,]$Mgmt_zone))
      
      setwd("//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Density Rasters")
      ugh <- paste0("60to85 by 5% Relative Percent Population_",thisZoneRoman,".img")
      rasters[[i]] <- raster(ugh)
      
      xExt <- rasters[[i]]@ncols
      yExt <- rasters[[i]]@nrows
      xScale <- c(rep( floor(xExt / scale), scale - 1), xExt - sum(rep( floor(xExt / scale), scale - 1)))    # make size of x blocks
      yScale <- c(rep( floor(yExt / scale), scale - 1), yExt - sum(rep( floor(yExt / scale), scale - 1)))    # make size of y blocks
      
      tmpRast <- list("vector",scale^2)
      tmpPoly <- list("vector",scale^2)
      thePoly <- NULL
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
          tmpRast[[putHere]] <- crop(rasters[[i]],extent(rasters[[i]],yMin,yMax,xMin,xMax))
          if(!(tmpRast[[putHere]]@data@min %in% c(Inf,-Inf) | tmpRast[[putHere]]@data@min %in% c(Inf,-Inf))){  # check if rast has other than NA
            tmpPoly[[putHere]] <- rasterToPolygons(tmpRast[[putHere]],dissolve=TRUE)
            names(tmpPoly[[putHere]]@data)[1] <- "Class"
    
            #writeOGR(tmpPoly[[putHere]],tempDir,paste0("shp",i,j,k),driver="ESRI Shapefile",overwrite_layer=TRUE)#gUnion(thePoly,tmpPoly[[putHere]],byid=TRUE)
       
          } else {
            tmpPoly[[putHere]] <- -999
          }
          cat(paste0("Completed loop for file ",i," - row ",j," - column ",k," on a ",scale,"x",scale," scale.\n"))
        }    
      }
      
      # put all the shapefile bits into one file.  
      j <- 1            
      theFirst <- 0
      for(j in 1:(scale*scale)){
        
        if(!is.numeric(tmpPoly[[j]])){
        
          theFirst <- ifelse(theFirst == 0,1,2)                  # j == 1 may not be the first truly shapefile.  need to find it
          
          nR             <- length(slot(tmpPoly[[j]],"polygons"))   
          if(theFirst == 1){                                                                                                                                                                               # for 1st shp file, do this
            uidR         <- 1;                                                                                                                                                               # make a unique id
            polys[[i]]   <- spChFIDs(tmpPoly[[j]], as.character(uidR:(uidR + nR - 1)))                                                                                              # make feature id of polygons unique
            uidR         <- uidR + nR                                                                                                                                                       # make unique id for all polys   
            theFirst     <- 2                                         # turn theFirst away from its magic value
          } else {                                                                                                                                                                                  # for other than 1st shp file, do this
            tmpPoly[[j]] <- spChFIDs(tmpPoly[[j]], as.character(uidR:(uidR + nR - 1)))                                                                                               # make feature id of polygons unique
            uidR         <- uidR + nR                                                                                                                                                       # make unique id for all polys
            polys[[i]]   <- spRbind(polys[[i]],tmpPoly[[j]])                                                                                                                               # union ith shp file with all previous ones
          }                                                                                                                                                                                         # close out if
        }
      }
    } 
  }

  
  save(polys,file=paste0(polyDir,"/polys.RData"))
  
  load(paste0(polyDir,"/polys.RData"))

  # make the true summed-up files i really need.
  for(h in 1:7){          # management zone
    for(i in 1:7){        # controller
      for(j in 1:i){      # 5% cutoff
        if(j == 1){
          assign(paste0("MZone_",h,"_the",5*(j - 1) + 60),gUnaryUnion(polys[[h]][polys[[h]]@data$Class == (7 - j),]))
        } else if (j > 1 & j < 7){
          assign(paste0("MZone_",h,"_the",5*(j - 1) + 60),gUnion(get(paste0("MZone_",h,"_the",5*(j - 2) + 60)),gUnaryUnion(polys[[h]][polys[[h]]@data$Class == (7 - j),])))
        } else {
          assign(paste0("MZone_",h,"_the100"),gUnion(get(paste0("MZone_",h,"_the",5*(j - 2) + 60)),gUnaryUnion(polys[[h]][polys[[h]]@data$Class == (7 - j),])))
        }
      }
      if(i <= 6){
        writeOGR(as(get(paste0("MZone_",h,"_the",5*(j - 1) + 60)),"SpatialPolygonsDataFrame"),polyDir,paste0("MZone_",h,"_the",5*(j - 1) + 60),driver="ESRI Shapefile",overwrite_layer=TRUE)          
      } else {
        writeOGR(as(get(paste0("MZone_",h,"_the100")),"SpatialPolygonsDataFrame"),polyDir,paste0("MZone_",h,"_the100"),driver="ESRI Shapefile",overwrite_layer=TRUE)              
      }
    }
  }

  # make shapefiles for mzone 8 = mzone 2 + mzone 7
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the60") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the60") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the60" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the65") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the65") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the65" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the70") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the70") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the70" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the75") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the75") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the75" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the80") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the80") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the80" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the85") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the85") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the85" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the100"),width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_7_the100"),width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame"),polyDir,"MZone_8_the100",driver="ESRI Shapefile",overwrite_layer=TRUE)


  z60 <- z65 <- z70 <- z75 <- z80 <- z85 <- z100 <- list("vector",6)
  # make shapefiles for mzone 9 = mzone 1 + mzone 2 + mzone 3 + mzone 4 + mzone 5 + mzone 6 + mzone 7
  for(i in 1:6){
    if(i == 1){
      z60[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the60") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the60") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z65[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the65") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the65") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z70[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the70") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the70") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z75[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the75") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the75") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z80[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the80") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the80") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z85[[i]]  <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the85") ,width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the85") ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z100[[i]] <- as(gUnion(as(gBuffer(readOGR(polyDir,"MZone_2_the100"),width=0.01),"SpatialPolygonsDataFrame"),as(gBuffer(readOGR(polyDir,"MZone_1_the100"),width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
    } else {
      z60[[i]]  <- as(gUnion(z60[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the60")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z65[[i]]  <- as(gUnion(z65[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the65")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z70[[i]]  <- as(gUnion(z70[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the70")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z75[[i]]  <- as(gUnion(z75[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the75")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z80[[i]]  <- as(gUnion(z80[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the80")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z85[[i]]  <- as(gUnion(z85[[i - 1]] ,as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the85")) ,width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
      z100[[i]] <- as(gUnion(z100[[i - 1]],as(gBuffer(readOGR(polyDir,paste0("MZone_",i + 1,"_the100")),width=0.01),"SpatialPolygonsDataFrame")),"SpatialPolygonsDataFrame")
    } 
  }

  # output the zone 9 = national shapefiles!
  writeOGR(z60[[6]] ,polyDir,"MZone_9_the60" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z65[[6]] ,polyDir,"MZone_9_the65" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z70[[6]] ,polyDir,"MZone_9_the70" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z75[[6]] ,polyDir,"MZone_9_the75" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z80[[6]] ,polyDir,"MZone_9_the80" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z85[[6]] ,polyDir,"MZone_9_the85" ,driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(z100[[6]],polyDir,"MZone_9_the100",driver="ESRI Shapefile",overwrite_layer=TRUE)









}    # close out the function