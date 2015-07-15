

require(R2WinBUGS)      # run WinBUGS
require(lme4)           # glmm for poisson
require(sp)             # read in spatial points
require(rgdal)          # read in and write shpfiles
require(rgeos)          # geometry processing
require(raster)         # read in rasters
require(maptools)       #spRbind

dataDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data"
origDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Original"
progDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Programs"
# wBugDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Programs"
outpDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Output"
polyDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Density Polygons"
analDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Analysis Sets"

PROJaea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
PROJlat <- "+init=epsg:4326"

source(paste0(progDir,"/helperFunctions.R"))                # helpful functions available for use, when needed
source(paste0(progDir,"/prepareForAnalysis.R"))             # prep the data; deal with zeros
#source(paste0(progDir,"/makeDensityPolyBits.R"))            # convert rasters to shapefiles
#source(paste0(progDir,"/assignCore.R"))                     # assign core
#source(paste0(progDir,"/pointsSummary.R"))                  # summarize results of assigning by core

# source(paste0(progDir,"/runMA.R"))
# source(paste0(progDir,"/runMB.R"))
source(paste0(progDir,"/runMC.R"))
source(paste0(progDir,"/runMD.R"))
# source(paste0(progDir,"/runME.R"))
source(paste0(progDir,"/runMF.R"))
# source(paste0(progDir,"/runMG.R"))
# source(paste0(progDir,"/runMH.R"))
# source(paste0(progDir,"/runMI.R"))
source(paste0(progDir,"/runMJ.R"))
# source(paste0(progDir,"/runMK.R"))
source(paste0(progDir,"/runML.R"))
# source(paste0(progDir,"/runMM.R"))
# source(paste0(progDir,"/runMN.R"))
# source(paste0(progDir,"/runMO.R"))
# source(paste0(progDir,"/runMP.R"))
# source(paste0(progDir,"/runMQ.R"))


sg <- read.csv(paste0(origDir,"/allStatesFinal.csv"))

datList <- prepareForAnalysis(sg)

#makeDensityPolyBits(datList[[1]])                           # make shapefiles from rasters



# datList returns the sage grouse data in four different ways
#      data cut       nrows (all zeros)     nrows (special zeros)
# 1. all the data               92,031                    75,726
# 2.  1965 - 1989               xx,xxx                     3,180
# 3.  1990 - 2014               xx,xxx                     6,349
# 4. ??

#assignCore(datList[[1]])   # data gets output as shapefiles
#pointsSummary()             # summarize the points shapefiles




# do the analysis

mZones <- unique(datList[[1]]$mZoneNum)
mZones <- mZones[order(mZones)]

datAllZerosCore <- vector("list",9)
datAllZerosNoco <- vector("list",9) 
datAllZerosLeks <- vector("list",9) 
dat1stZerosCore <- vector("list",9) 
dat1stZerosNoco <- vector("list",9) 
dat1stZerosLeks <- vector("list",9)
#resultsA <- vector("list",8)
#resultsB <- vector("list",8)
resultsC <- vector("list",8)
resultsD <- vector("list",8)
#resultsE <- vector("list",8)
resultsF <- vector("list",8)

#resultsG <- vector("list",1)
#resultsH <- vector("list",1)
#resultsI <- vector("list",1)
resultsJ <- vector("list",1)
#resultsK <- vector("list",1)
resultsL <- vector("list",1)

for(i in 5:7){

#   datAllZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
#   datAllZerosNoco[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
#   datAllZerosLeks[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data        # read in all zeros, all data, ith mzone
 
  dat1stZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
#   dat1stZerosNoco[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
#   dat1stZerosLeks[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
  
  if(i <= 8){   
    
    BUGSDir <- 'C:/Users/jmitchell/WinBUGS14'
    
   #resultsA[[i]] <- runMA(datAllZerosCore[[i]],progDir,BUGSDir,i)    # all zeros, core,     no mzone effect
   #resultsB[[i]] <- runMB(datAllZerosNoco[[i]],progDir,BUGSDir,i)    # all zeros, non-core, no mzone effect
   #resultsC[[i]] <- runMC(datAllZerosLeks[[i]],progDir,BUGSDir,i)    # all zeros, all leks, no mzone effect
  
    resultsD[[i]] <- runMD(dat1stZerosCore[[i]],progDir,BUGSDir,i)    # 1st zeros, core,     no mzone effect
   #resultsE[[i]] <- runME(dat1stZerosNoco[[i]],progDir,BUGSDir,i)    # 1st zeros, non-core, no mzone effect
   #resultsF[[i]] <- runMF(dat1stZerosLeks[[i]],progDir,BUGSDir,i)    # 1st zeros, all leks, no mzone effect     <--- already done ---
    
  } else if(i == 9){
    
   # mZoneNum == 9, but Mgmt_zone == different things.
    datAllZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosCore[[9]]$Mgmt_zone))))
    datAllZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosNoco[[9]]$Mgmt_zone))))
    datAllZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosLeks[[9]]$Mgmt_zone))))
    dat1stZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosCore[[9]]$Mgmt_zone))))
    dat1stZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosNoco[[9]]$Mgmt_zone))))
    dat1stZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosLeks[[9]]$Mgmt_zone))))

   #resultsG[[1]] <- runMG(datAllZerosCore[[9]],progDir)    # all zeros, core,     fixed mzone int + slope
   #resultsH[[1]] <- runMH(datAllZerosNoco[[9]],progDir)    # all zeros, non-core, fixed mzone int + slope
   #resultsI[[1]] <- runMI(datAllZerosLeks[[9]],progDir)    # all zeros, all leks, fixed mzone int + slope
  
    resultsJ[[1]] <- runMJ(dat1stZerosCore[[9]],progDir)    # 1st zeros, core,     fixed mzone int + slope
   #resultsK[[1]] <- runMK(dat1stZerosNoco[[9]],progDir)    # 1st zeros, non-core, fixed mzone int + slope
    resultsL[[1]] <- runML(dat1stZerosLeks[[9]],progDir)    # 1st zeros, all leks, fixed mzone int + slope 
    
  }
}

