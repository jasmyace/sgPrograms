assignCore <- function(dat,shpfile){
  
  
  dat <- datList[[1]]#preCountDat
  
  
  coords = cbind(dat$Long,dat$Lat)
  sp = SpatialPoints(coords,proj4string=CRS(PROJlat))
  
  sp2 <- as(sp, "SpatialPointsDataFrame")
  sp2@data <- dat
  
  sp2 <- spTransform(sp2,CRS(PROJaea))
  
  leks <- list("vector",9)
  core <- list("vector",9)
  buff0 <- list("vector",9)
  diff0 <- rep(NA,9)
  perc0 <- rep(NA,9)
  buff1 <- list("vector",9)
  diff1 <- rep(NA,9)
  perc1 <- rep(NA,9)
  keep <- list("vector",9)
  
  allLek <- list("vector",9)
  corLek <- list("vector",9)
  inCore <- list("vector",9)
  
  Fleks1stZero <- list("vector",9)
  Fcore1stZero <- list("vector",9)
  Fnoco1stZero <- list("vector",9)
  
  FleksAllZero <- list("vector",9)
  FcoreAllZero <- list("vector",9)
  FnocoAllZero <- list("vector",9)
  
  SumFleks1stZero <- NULL
  SumFcore1stZero <- NULL
  SumFnoco1stZero <- NULL
  
  SumFleksAllZero <- NULL
  SumFcoreAllZero <- NULL
  SumFnocoAllZero <- NULL
  
  for(i in 1:7){
    
    leks[[i]] <- sp2[sp2@data$mZoneNum == i,]
    core[[i]] <- readOGR(polyDir,paste0('MZone_',i,'_the75'))    # fix self-intersections   
    
    if(gIsValid(core[[i]]) == TRUE){
      buff0[[i]] <- core[[i]]
      buff1[[i]] <- core[[i]]
    } else {
      buff0[[i]] <- gBuffer(core[[i]],width=0) 
      buff1[[i]] <- gBuffer(core[[i]],width=0.01) 
    }
    
    diff0[i] <- ( gArea(buff0[[i]]) - gArea(core[[i]]) ) / 120^2
    perc0[i] <- round( 100* (diff0[i] / ( gArea(core[[i]]) / 120^2 )),4)
    
    diff1[i] <- ( gArea(buff1[[i]]) - gArea(core[[i]]) ) / 120^2
    perc1[i] <- round( 100* (diff1[i] / ( gArea(core[[i]]) / 120^2 )),4)
    
    ugh <- gWithin(leks[[i]],buff0[[i]],byid=TRUE)
    ughNames <- colnames(ugh)[ugh]
    
    keep[[i]] <- leks[[i]][rownames(leks[[i]]@data) %in% ughNames,]
    
    allLek[[i]] <- leks[[i]]@data
    allLek[[i]]$Rsort <- c(1:nrow(allLek[[i]]))
    corLek[[i]] <- keep[[i]]@data[,c('Mgmt_zone','Lek_ID','Year')]
    corLek[[i]]$inCore <- 1
    inCore[[i]] <- merge(allLek[[i]],corLek[[i]],by=c('Mgmt_zone','Lek_ID','Year'),all.x=TRUE)    # mzone not needed but eh
    inCore[[i]]$inCore[is.na(inCore[[i]]$inCore)] <- 0
    inCore[[i]] <- inCore[[i]][order(inCore[[i]]$Rsort),]
    
    FleksAllZero[[i]] <- leks[[i]]
    
    FleksAllZero[[i]]@data <- inCore[[i]]                                              # reduce to core and non-core, all zeros
    FcoreAllZero[[i]] <- FleksAllZero[[i]][FleksAllZero[[i]]@data$inCore == 1,]        # reduce to core             , all zeros
    FnocoAllZero[[i]] <- FleksAllZero[[i]][FleksAllZero[[i]]@data$inCore == 0,]        # reduce to          non-core, all zeros
    
    Fleks1stZero[[i]] <- FleksAllZero[[i]][FleksAllZero[[i]]@data$DupZero == 0,]       # reduce to core and non-core, 1st zeros
    Fcore1stZero[[i]] <- Fleks1stZero[[i]][Fleks1stZero[[i]]@data$inCore == 1,]        # reduce to core             , all zeros
    Fnoco1stZero[[i]] <- Fleks1stZero[[i]][Fleks1stZero[[i]]@data$inCore == 0,]        # reduce to          non-core, all zeros    
    
    writeOGR(FleksAllZero[[i]],analDir,paste0('Zone ',i,' Both-75 - All Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)     # core & non-core, all zero 
    writeOGR(FcoreAllZero[[i]],analDir,paste0('Zone ',i,' Core-75 - All Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)                # core           , all zero 
    writeOGR(FnocoAllZero[[i]],analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)            #        non-core, all zero
    
    writeOGR(Fleks1stZero[[i]],analDir,paste0('Zone ',i,' Both-75 - 1st Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)     # core & non-core, 1st zero 
    writeOGR(Fcore1stZero[[i]],analDir,paste0('Zone ',i,' Core-75 - 1st Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)                # core           , 1st zero 
    writeOGR(Fnoco1stZero[[i]],analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)            #        non-core, 1st zero
    
  }

  # deal with zones 8 and 9
  for(i in 8:9){
    if(i == 8){       # combine zones 2 and 7
      for(j in 1:3){
        if(j == 1){
          doit <- 'Both'
        } else if(j == 2){
          doit <- 'Core-75'
        } else {
          doit <- 'Non-Core-75'
        }
        pts8.A1 <- readOGR(analDir,paste0('Zone 2 ',doit,' - All Zero'))
        pts8.A2 <- readOGR(analDir,paste0('Zone 7 ',doit,' - All Zero'))
        pts8.A <- spRbind(pts8.A1,pts8.A2)
        pts8.A@data$Mgmt_zone <- "MZ VIII"
        pts8.A@data$mZoneNum <- 8
        writeOGR(pts8.A,analDir,paste0('Zone ',i,' ',doit,' - All Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)    
        
        pts8.B1 <- readOGR(analDir,paste0('Zone 2 ',doit,' - 1st Zero'))
        pts8.B2 <- readOGR(analDir,paste0('Zone 7 ',doit,' - 1st Zero'))
        pts8.B <- spRbind(pts8.B1,pts8.B2)
        pts8.B@data$Mgmt_zone <- "MZ VIII"
        pts8.B@data$mZoneNum <- 8
        writeOGR(pts8.B,analDir,paste0('Zone ',i,' ',doit,' - 1st Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE) 
        
        rm(pts8.A1,pts8.A2,pts8.A,pts8.B1,pts8.B2,pts8.B)
      }
    } else {          # make big file
      for(j in 1:3){
        if(j == 1){
          doit <- 'Both'
        } else if(j == 2){
          doit <- 'Core-75'
        } else {
          doit <- 'Non-Core-75'
        }
        for(k in 1:7){
          if(k == 1){            
            pts9.Aa <- readOGR(analDir,paste0('Zone ',k,' ',doit,' - All Zero'))
            pts9.Bb <- readOGR(analDir,paste0('Zone ',k,' ',doit,' - 1st Zero'))
          } else {
            pts9.Aa <- spRbind(pts9.Aa,readOGR(analDir,paste0('Zone ',k,' ',doit,' - All Zero')))
            pts9.Bb <- spRbind(pts9.Bb,readOGR(analDir,paste0('Zone ',k,' ',doit,' - 1st Zero')))
          }
        }
        #pts9.Aa@data$Mgmt_zone <- "MZ IX"
        pts9.Aa@data$mZoneNum <- 9
        writeOGR(pts9.Aa,analDir,paste0('Zone 9 ',doit,' - All Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)    
        
        #pts9.Bb@data$Mgmt_zone <- "MZ IX"
        pts9.Bb@data$mZoneNum <- 9
        writeOGR(pts9.Bb,analDir,paste0('Zone 9 ',doit,' - 1st Zero'),driver="ESRI Shapefile",overwrite_layer=TRUE)             
        
        rm(pts9.Aa,pts9.Bb)
      }  
    }  
  }
  
}
