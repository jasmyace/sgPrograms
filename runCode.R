
# packDir <- "C:/Program Files/R/R-3.2.1/library"
# install.packages(c("R2WinBUGS","lme4","sp","rgdal","rgeos","raster","maptools","Cairo","RColorBrewer"),packDir,repos="http://cran.us.r-project.org")

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
#source(paste0(progDir,"/makeDensityPolyBits.R"))           # convert rasters to shapefiles
#source(paste0(progDir,"/assignCore.R"))                    # assign core
#source(paste0(progDir,"/pointsSummary.R"))                 # summarize results of assigning by core
source(paste0(progDir,"/readInAnalyticFiles.R"))            # read in analytic shapefiles
source(paste0(progDir,"/makeTracePlots.R"))                 # make trace plots
source(paste0(progDir,"/makePosteriorPlots.R"))             # make posterior plots
source(paste0(progDir,"/makeTrendPlots.R"))                 # make trend plots
source(paste0(progDir,"/makeHistogramPlots.R"))             # make histogram plots

#source(paste0(progDir,"/runMA.R"))
#source(paste0(progDir,"/runMB.R"))
#source(paste0(progDir,"/runMC.R"))
source(paste0(progDir,"/runMD.R"))
source(paste0(progDir,"/runME.R"))
source(paste0(progDir,"/runMF.R"))
#source(paste0(progDir,"/runMG.R"))
#source(paste0(progDir,"/runMH.R"))
#source(paste0(progDir,"/runMI.R"))
#source(paste0(progDir,"/runMJ.R"))
#source(paste0(progDir,"/runMK.R"))
#source(paste0(progDir,"/runML.R"))
#source(paste0(progDir,"/runMM.R"))
#source(paste0(progDir,"/runMN.R"))
#source(paste0(progDir,"/runMO.R"))
#source(paste0(progDir,"/runMP.R"))
#source(paste0(progDir,"/runMQ.R"))
#source(paste0(progDir,"/runMR.R"))
#source(paste0(progDir,"/runMS.R"))
#source(paste0(progDir,"/runMT.R"))
#source(paste0(progDir,"/runMU.R"))
source(paste0(progDir,"/runMV.R"))
source(paste0(progDir,"/runMZ.R"))

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

datAllZerosCore <- vector("list",9)
#datAllZerosNoco <- vector("list",9) 
#datAllZerosLeks <- vector("list",9) 
dat1stZerosCore <- vector("list",9) 
dat1stZerosNoco <- vector("list",9) 
dat1stZerosLeks <- vector("list",9)

#resultsA <- vector("list",8)
#resultsB <- vector("list",8)
#resultsC <- vector("list",8)
resultsD <- vector("list",8)
resultsE <- vector("list",8)
resultsF <- vector("list",8)
#resultsG <- vector("list",1)
#resultsH <- vector("list",1)
#resultsI <- vector("list",1)
#resultsJ <- vector("list",1)
#resultsK <- vector("list",1)
#resultsL <- vector("list",1)
#resultsM <- vector("list",1)
#resultsN <- vector("list",1)
#resultsO <- vector("list",1)
#resultsP <- vector("list",1)
#resultsQ <- vector("list",1)
#resultsR <- vector("list",1)
#resultsS <- vector("list",1)
#resultsT <- vector("list",1)
#resultsU <- vector("list",1)
resultsV <- vector("list",8)
resultsZ <- vector("list",1)

for(i in 1:6){

  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  if(i == 7){mZone <- 9}
  
  datAllZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
 #datAllZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
 #datAllZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - All Zero'))@data        # read in all zeros, all data, ith mzone
 
  dat1stZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
  dat1stZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
  dat1stZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
  
 #BUGSDir <- 'C:/Program Files/WinBUGS14'                                                              # my machine  
 #BUGSDir <- 'C:/Users/jmitchell/WinBUGS_14/winbugs14/WinBUGS14'                                       # sonic 1
  BUGSDir <- 'C:/Users/jmitchell/WinBUGS14'                                                            # sonic 2 & others

  if(mZone <= 8){   
    
   #resultsA[[mZone]] <- runMA(datAllZerosCore[[mZone]],progDir,BUGSDir,mZone)                         # all zeros, core,     no mzone effect
   #resultsB[[mZone]] <- runMB(datAllZerosNoco[[mZone]],progDir,BUGSDir,mZone)                         # all zeros, non-core, no mzone effect
   #resultsC[[mZone]] <- runMC(datAllZerosLeks[[mZone]],progDir,BUGSDir,mZone)                         # all zeros, all leks, no mzone effect
  
    resultsD[[mZone]] <- runMD(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Try 1',mZone)                 # 1st zeros, core,     ind mzone, B mat
    resultsE[[mZone]] <- runME(dat1stZerosNoco[[mZone]],progDir,BUGSDir,'Try 1',mZone)                 # 1st zeros, non-core, ind mzone, B mat
    resultsF[[mZone]] <- runMF(dat1stZerosLeks[[mZone]],progDir,BUGSDir,'Try 1',mZone)                 # 1st zeros, all leks, no mzone effect     
    
    resultsV[[mZone]] <- runMV(datAllZerosCore[[mZone]],progDir,BUGSDir,'ZInf Try 2',mZone)            # ALL zeros, core,     ind mzone, B mat, zero inf
    
   } #else if(mZone == 9){
#     
#     # run models based on cutting up mZone 9, or all the data, into constituent states
#     nStates <- length(unique(dat1stZerosLeks[[mZone]]$State))
#     states <- as.character(droplevels(unique(dat1stZerosLeks[[mZone]]$State)))
#     
#    #dat1stZerosCoreStates <- vector("list",nStates)
#    #dat1stZerosNocoStates <- vector("list",nStates)
#     dat1stZerosLeksStates <- vector("list",nStates)
#     
#    #resultsD.st <- vector("list",nStates)
#    #resultsE.st <- vector("list",nStates)
#     resultsF.st <- vector("list",nStates)
#     
#     for(j in 1:nStates){
#       thisState <- states[j]
#      #dat1stZerosCoreStates[[j]] <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == thisState,]
#      #dat1stZerosNocoStates[[j]] <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == thisState,]
#       dat1stZerosLeksStates[[j]] <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == thisState,]
#      #resultsD.st[[j]] <- runMD(dat1stZerosCoreStates[[j]],progDir,BUGSDir,'Try 1',thisState)       # 1 state, 1st zeros, core    
#      #resultsE.st[[j]] <- runME(dat1stZerosNocoStates[[j]],progDir,BUGSDir,'Try 1',thisState)       # 1 state, 1st zeros, non-core 
#       resultsF.st[[j]] <- runMF(dat1stZerosLeksStates[[j]],progDir,BUGSDir,'Try 1',thisState)       # 1 state, 1st zeros, all leks
#     }
#      
#    #datAllZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosCore[[9]]$Mgmt_zone))))
#    #datAllZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosNoco[[9]]$Mgmt_zone))))
#    #datAllZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(datAllZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(datAllZerosLeks[[9]]$Mgmt_zone))))
#    #dat1stZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosCore[[9]]$Mgmt_zone))))
#    #dat1stZerosNoco[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosNoco[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosNoco[[9]]$Mgmt_zone))))
#    #dat1stZerosLeks[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosLeks[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosLeks[[9]]$Mgmt_zone))))
#     
#    #resultsG[[1]] <- runMG(datAllZerosCore[[9]],progDir,BUGSDir,9)    # all zeros, core,     fixed mzone int + slope
#    #resultsH[[1]] <- runMH(datAllZerosNoco[[9]],progDir,BUGSDir,9)    # all zeros, non-core, fixed mzone int + slope
#    #resultsI[[1]] <- runMI(datAllZerosLeks[[9]],progDir,BUGSDir,9)    # all zeros, all leks, fixed mzone int + slope
#   
#    #resultsJ[[1]] <- runMJ(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,     fixed mzone int + slope
#    #resultsK[[1]] <- runMK(dat1stZerosNoco[[9]],progDir,BUGSDir,9)    # 1st zeros, non-core, fixed mzone int + slope
#    #resultsL[[1]] <- runML(dat1stZerosLeks[[9]],progDir,BUGSDir,9)    # 1st zeros, all leks, fixed mzone int + slope 
# 
#    #Core1stZeros15p <- readOGR(analDir,'Core1stZeros15p')@data
#    #source(paste0(progDir,"/makeSmallCoreSamp.R"))
#    #smallCoreSamp <- makeSmallCoreSamp()
#     
#    #smallCoreSamp2 <- smallCoreSamp[smallCoreSamp$Mgmt_zone %in% c('MZ I','MZ III','MZ V'),]
#    #save(smallCoreSamp2,file=paste0(outpDir,"/smallCoreSamp2.Sonic1.RData"))
#     
#    #resultsM[[1]] <- runMM(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,                  triple-level random mzone -- k lek vars
#    #resultsN[[1]] <- runMN(dat1stZerosCore[[9]],progDir,BUGSDir,9)    # 1st zeros, core,                  triple-level fixed mzone -- 1 lek var
#    #resultsN[[1]] <- runMN(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level fixed mzone -- 1 lek var
#    #resultsO[[1]] <- runMO(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone -- 1 lek var
#    #resultsP[[1]] <- runMP(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level fixed mzone, collapsed model-- 1 lek var
#    #resultsQ[[1]] <- runMQ(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone, collapsed model-- 1 lek var
#    #resultsR[[1]] <- runMR(Core1stZeros15p     ,progDir,BUGSDir,9)    # 1st zeros, core, 15% BAS sample,  triple-level random mzone, collapsed model, exp val monitoring
#    #resultsS[[1]] <- runMS(dat1stZerosCore[[9]],progDir,BUGSDir,'All MZones Core, Try 5',9)          # 1st zeros,     core, all 2015 data, Sauer-Link model with B-matrices
#    #resultsS[[1]] <- runMS(dat1stZerosNoco[[9]],progDir,BUGSDir,'All MZones Non-Core, Try 5',9)      # 1st zeros, non-core, all 2015 data, Sauer-Link model with B-matrices
# 
#    #resultsS[[1]] <- runMS(smallCoreSamp2      ,progDir,BUGSDir,'MZones 1,3,5 Core, Try 1',9)        # 1st zeros,     core, all 2015 data, Sauer-Link model with B-matrices, exp adj M   
#    #resultsT[[1]] <- runMT(smallCoreSamp2      ,progDir,BUGSDir,'MZones 1,3,5 Core, Adj M,Try 1',9)  # 1st zeros,     core, all 2015 data, Sauer-Link model with B-matrices, model adj M
# 
#    #resultsU[[1]] <- runMU(smallCoreSamp2      ,progDir,BUGSDir,'MZones 1,3,5 Core, Adj M,Try 1',9)  # 1st zeros,     core, all 2015 data, Sauer-Link model with B-matrices, model adj M
#    
#    #resultsU[[1]] <- runMU(dat1stZerosCore[[9]],progDir,BUGSDir,'All MZones Core, Adj M,Try 1',9)    # 1st zeros,     core, all 2015 data, Sauer-Link model with B-matrices, model adj M
#    #resultsU[[1]] <- runMU(dat1stZerosNoco[[9]],progDir,BUGSDir,'All MZones Non-Core, Adj M,Try 1',9)# 1st zeros, non-core, all 2015 data, Sauer-Link model with B-matrices, model adj M 
#   }
}

resultsZ[[1]] <- runMZ(fakeData,progDir,BUGSDir,'Fake Data Test',99)            # fake data test




source(paste0(progDir,"/summarizeModels.R"))
source(paste0(progDir,"/compileResults.R"))



