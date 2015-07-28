

require(R2WinBUGS)      # run WinBUGS
require(lme4)           # glmm for poisson
require(sp)             # read in spatial points
require(rgdal)          # read in and write shpfiles
require(rgeos)          # geometry processing
require(raster)         # read in rasters
require(maptools)       #spRbind
require(Cairo)          # nice output
require(RColorBrewer)   # nice colors

dataDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data"
origDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Original"
progDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/sgPrograms"
# wBugDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Programs"
outpDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Output"
rsltDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Results"
polyDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Density Polygons"
analDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Analysis Sets"

PROJaea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
PROJlat <- "+init=epsg:4326"

source(paste0(progDir,"/helperFunctions.R"))                # helpful functions available for use, when needed
source(paste0(progDir,"/prepareForAnalysis.R"))             # prep the data; deal with zeros
#source(paste0(progDir,"/makeDensityPolyBits.R"))            # convert rasters to shapefiles
#source(paste0(progDir,"/assignCore.R"))                     # assign core
#source(paste0(progDir,"/pointsSummary.R"))                  # summarize results of assigning by core



#source(paste0(progDir,"/runMA.R"))
#source(paste0(progDir,"/runMB.R"))
#source(paste0(progDir,"/runMC.R"))
#source(paste0(progDir,"/runMD.R"))
#source(paste0(progDir,"/runME.R"))
#source(paste0(progDir,"/runMF.R"))
#source(paste0(progDir,"/runMG.R"))
#source(paste0(progDir,"/runMH.R"))
#source(paste0(progDir,"/runMI.R"))
#source(paste0(progDir,"/runMJ.R"))
#source(paste0(progDir,"/runMK.R"))
#source(paste0(progDir,"/runML.R"))
#source(paste0(progDir,"/runMM.R"))
#source(paste0(progDir,"/runMN.R"))
#source(paste0(progDir,"/runMO.R"))
source(paste0(progDir,"/runMP.R"))
#source(paste0(progDir,"/runMQ.R"))
source(paste0(progDir,"/runMR.R"))

sg <- read.csv(paste0(origDir,"/allStatesFinal2015.csv"))

datList <- prepareForAnalysis(sg)

rm(sg)
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

# datAllZerosCore <- vector("list",9)
# datAllZerosNoco <- vector("list",9) 
# datAllZerosLeks <- vector("list",9) 
dat1stZerosCore <- vector("list",9) 
# dat1stZerosNoco <- vector("list",9) 
# dat1stZerosLeks <- vector("list",9)

#resultsA <- vector("list",8)
#resultsB <- vector("list",8)
#resultsC <- vector("list",8)
#resultsD <- vector("list",8)
#resultsE <- vector("list",8)
#resultsF <- vector("list",8)

#resultsG <- vector("list",1)
#resultsH <- vector("list",1)
#resultsI <- vector("list",1)
#resultsJ <- vector("list",1)
#resultsK <- vector("list",1)
#resultsL <- vector("list",1)
#resultsM <- vector("list",1)
#resultsN <- vector("list",1)
#resultsO <- vector("list",1)
resultsP <- vector("list",1)
resultsQ <- vector("list",1)
resultsR <- vector("list",1)

for(i in 9:9){

#   datAllZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
#   datAllZerosNoco[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
#   datAllZerosLeks[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data        # read in all zeros, all data, ith mzone
 
  dat1stZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
#   dat1stZerosNoco[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
#   dat1stZerosLeks[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
  
  BUGSDir <- 'C:/Program Files/WinBUGS14'                         # my machine  
  BUGSDir <- 'C:/Users/jmitchell/WinBUGS_14/winbugs14/WinBUGS14'  # sonic 1
  BUGSDir <- 'C:/Program Files/WinBUGS14'  ?                      # sonic 2 & others

  if(i <= 8){   
    
   #resultsA[[i]] <- runMA(datAllZerosCore[[i]],progDir,BUGSDir,i)    # all zeros, core,     no mzone effect
   #resultsB[[i]] <- runMB(datAllZerosNoco[[i]],progDir,BUGSDir,i)    # all zeros, non-core, no mzone effect
   #resultsC[[i]] <- runMC(datAllZerosLeks[[i]],progDir,BUGSDir,i)    # all zeros, all leks, no mzone effect
  
   #resultsD[[i]] <- runMD(dat1stZerosCore[[i]],progDir,BUGSDir,i)    # 1st zeros, core,     no mzone effect
   #resultsE[[i]] <- runME(dat1stZerosNoco[[i]],progDir,BUGSDir,i)    # 1st zeros, non-core, no mzone effect
   #resultsF[[i]] <- runMF(dat1stZerosLeks[[i]],progDir,BUGSDir,i)    # 1st zeros, all leks, no mzone effect     <--- already done ---
    
  } else if(i == 9){
    
   # mZoneNum == 9, but Mgmt_zone == different things.
    
   #datAllZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosCore[[9]]$Mgmt_zone))))
   #datAllZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosNoco[[9]]$Mgmt_zone))))
   #datAllZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosLeks[[9]]$Mgmt_zone))))
   #dat1stZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosCore[[9]]$Mgmt_zone))))
   #dat1stZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosNoco[[9]]$Mgmt_zone))))
   #dat1stZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosLeks[[9]]$Mgmt_zone))))
    
   #resultsG[[1]] <- runMG(datAllZerosCore[[9]],progDir,BUGSDir,9)    # all zeros, core,     fixed mzone int + slope
   #resultsH[[1]] <- runMH(datAllZerosNoco[[9]],progDir,BUGSDir,9)    # all zeros, non-core, fixed mzone int + slope
   #resultsI[[1]] <- runMI(datAllZerosLeks[[9]],progDir,BUGSDir,9)    # all zeros, all leks, fixed mzone int + slope
  
   #resultsJ[[1]] <- runMJ(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,     fixed mzone int + slope
   #resultsK[[1]] <- runMK(dat1stZerosNoco[[9]],progDir,BUGSDir,9)    # 1st zeros, non-core, fixed mzone int + slope
   #resultsL[[1]] <- runML(dat1stZerosLeks[[9]],progDir,BUGSDir,9)    # 1st zeros, all leks, fixed mzone int + slope 

    Core1stZeros15p <- readOGR(analDir,'Core1stZeros15p')@data
    source(paste0(progDir,"/makeSmallCoreSamp.R"))
    smallCoreSamp <- makeSmallCoreSamp()
    
   #resultsM[[1]] <- runMM(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,                  triple-level random mzone -- k lek vars
   #resultsN[[1]] <- runMN(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,                  triple-level fixed mzone -- 1 lek var
   #resultsN[[1]] <- runMN(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level fixed mzone -- 1 lek var
   #resultsO[[1]] <- runMO(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone -- 1 lek var
    resultsP[[1]] <- runMP(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level fixed mzone, collapsed model-- 1 lek var
    resultsQ[[1]] <- runMQ(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone, collapsed model-- 1 lek var
    resultsR[[1]] <- runMR(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone, collapsed model, exp val monitoring
   
  }
}
#write.csv(resultsN[[1]][[4]]$summary,'//LAR-FILE-SRV/Data/Jason/sage grouse/checkit.csv')
# resultsN[[1]][[4]]$summary


source(paste0(progDir,"/summarizeModels.R"))
source(paste0(progDir,"/compileResults.R"))



