
# packDir <- "C:/Program Files/R/R-3.2.1/library"
# install.packages(c("R2WinBUGS","lme4","sp","rgdal","rgeos","raster","maptools","Cairo","RColorBrewer","spsurvey","RGtk2","plyr","reshape2","AICcmodavg","car","R2wd","jpeg","matrixcalc","countreg"),packDir,repos="http://cran.us.r-project.org")

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
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.09.02'

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
source(paste0(progDir,"/make90pCredInt.R"))                 # calculate 5% 95% percentiles
source(paste0(progDir,"/makeLekTrendPlots.R"))              # make lek plots
source(paste0(progDir,"/makeRandomLekPlots.R"))             # make random lek plots
source(paste0(progDir,"/makeRangeWide.R"))                  # make rangewide


datAllZerosCore <- vector("list",9)
datAllZerosNoco <- vector("list",9) 
datAllZerosLeks <- vector("list",9) 
dat1stZerosCore <- vector("list",9) 
dat1stZerosNoco <- vector("list",9) 
dat1stZerosLeks <- vector("list",9)

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

# testing-for-discontinuity programs

source(paste0(progDir,"/runTestAn.R"))
source(paste0(progDir,"/runTestBn.R"))
source(paste0(progDir,"/runTestCn.R"))
source(paste0(progDir,"/runTestDn.R"))
source(paste0(progDir,"/runTestEn.R"))
source(paste0(progDir,"/runTestFn.R"))
source(paste0(progDir,"/runTestGn.R"))

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

resultsAllAn <- vector("list",8)
resultsAllBn <- vector("list",8)
resultsAllCn <- vector("list",8)
resultsAllDn <- vector("list",8)
resultsAllEn <- vector("list",8)
resultsAllFn <- vector("list",8)

results1stAn <- vector("list",8)
results1stBn <- vector("list",8)
results1stCn <- vector("list",8)
results1stDn <- vector("list",8)
results1stEn <- vector("list",8)
results1stFn <- vector("list",8)

for(i in 1:6){

  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  if(i == 7){mZone <- 9}
  
  datAllZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
  datAllZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
  datAllZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - All Zero'))@data        # read in all zeros, all data, ith mzone
 
  #dat1stZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
  #dat1stZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
  #dat1stZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
  
 #BUGSDir <- 'C:/Program Files/WinBUGS14'                                                              # my machine  
 #BUGSDir <- 'C:/Users/jmitchell/WinBUGS_14/winbugs14/WinBUGS14'                                       # sonic 1
  BUGSDir <- 'C:/Users/jmitchell/WinBUGS14'                                                            # sonic 2 & others

  if(mZone <= 8){   
    
   #resultsD[[mZone]] <- runMD(datAllZerosCore[[mZone]],progDir,BUGSDir,'All Zeros 1',mZone)                         # all zeros, core,     no mzone effect
    resultsE[[mZone]] <- runME(datAllZerosNoco[[mZone]],progDir,BUGSDir,'All Zeros 1',mZone)                         # all zeros, non-core, no mzone effect
    resultsF[[mZone]] <- runMF(datAllZerosLeks[[mZone]],progDir,BUGSDir,'All Zeros 1',mZone)                         # all zeros, all leks, no mzone effect
  
   #resultsD[[mZone]] <- runMD(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Try constant mu.a+mu.b',mZone)                 # 1st zeros, core,     ind mzone, B mat
   #resultsE[[mZone]] <- runME(dat1stZerosNoco[[mZone]],progDir,BUGSDir,'Try 1',mZone)                 # 1st zeros, non-core, ind mzone, B mat
   #resultsF[[mZone]] <- runMF(dat1stZerosLeks[[mZone]],progDir,BUGSDir,'Try 1',mZone)                 # 1st zeros, all leks, no mzone effect     
    
   #resultsV[[mZone]] <- runMV(datAllZerosCore[[mZone]],progDir,BUGSDir,'ZInf Try 2',mZone)            # ALL zeros, core,     ind mzone, B mat, zero inf
   
   # TESTING TESTING TESTING
    
    
    
#    runTestAn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-An-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Bn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_a
#    runTestCn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Cn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_b
#    runTestDn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Dn-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop mu_a, mu_b correlation
#    runTestEn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-En-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Fn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-An-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Bn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_a
#    runTestCn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Cn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_b
#    runTestDn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Dn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-En-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Fn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    
#    runTestAn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-An-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Bn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_a
#    runTestCn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Cn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_b
#    runTestDn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Dn-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop mu_a, mu_b correlation
#    runTestEn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-En-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Fn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-An-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Bn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_a
#    runTestCn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Cn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_b
#    runTestDn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Dn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-En-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Fn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#     
#    runTestAn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-An-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Bn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_a
#    runTestCn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Cn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_b
#    runTestDn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Dn-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop mu_a, mu_b correlation
#    runTestEn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-En-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(datAllZerosCore[[mZone]],progDir,BUGSDir,'Test All-Fn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-An-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Bn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_a
#    runTestCn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Cn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_b
#    runTestDn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Dn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-En-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(dat1stZerosCore[[mZone]],progDir,BUGSDir,'Test 1st-Fn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    
#    
#    
#   
#    
# 
#    
#    runTestAn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-An-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Bn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_a
#    runTestCn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Cn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_b
#    runTestEn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-En-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Fn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    runTestGn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Gn-4k-76k-1',mZone,76000,80000,1)   # all zeros - original
#    runTestGn(simMZ1PoisRandCoefCorrData877Leks,progDir,BUGSDir,'Test simCorr877Lek-Gn-4k-76k-1',mZone,76000,80000,1)   # all zeros - original   
#    
# #    runTestAn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-An-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
# #    runTestBn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Bn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_a
# #    runTestCn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Cn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_b
# #    runTestEn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-En-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
# #    runTestFn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Fn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
# #    runTestGn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Gn-4k-76k-1',mZone,76000,80000,1)   # all zeros - original
#    
# #    runTestAn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-An-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
# #    runTestBn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Bn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_a
# #    runTestCn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Cn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix mu_b
# #    runTestDn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Dn-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop mu_a, mu_b correlation
# #    runTestEn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-En-4k-76k-1',mZone,76000,80000,1)   # all zeros - drop mu_a altogether (regression through the origin)
# #    runTestFn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Fn-4k-76k-1',mZone,76000,80000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
# #    runTestGn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Gn-4k-76k-1',mZone,76000,80000,1)   # all zeros - original
# 
#    runTestAn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-An-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Bn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_a
#    runTestCn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Cn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix mu_b
#    runTestDn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Dn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-En-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Fn-4k-76k-1',mZone,76000,80000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    runTestGn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Gn-4k-76k-1',mZone,76000,80000,1)   # all zeros - original   
#    
#    
#    
# 
#    
#    
#    
#    
#    runTestAn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-An-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Bn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_a
#    runTestCn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Cn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_b
#    runTestEn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-En-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Fn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-An-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Bn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_a
#    runTestCn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Cn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_b
#    runTestEn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-En-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Fn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    
#    runTestAn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-An-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Bn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_a
#    runTestCn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Cn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix mu_b
#    runTestDn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Dn-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop mu_a, mu_b correlation
#    runTestEn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-En-4k-304k-4',mZone,304000,320000,4)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Fn-4k-304k-4',mZone,304000,320000,4)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-An-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Bn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_a
#    runTestCn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Cn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix mu_b
#    runTestDn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Dn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-En-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Fn-4k-304k-4',mZone,304000,320000,4)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    runTestAn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-An-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Bn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_a
#    runTestCn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Cn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_b
#    runTestEn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-En-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefCorrData31Leks,progDir,BUGSDir,'Test simCorr31Lek-Fn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-An-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Bn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_a
#    runTestCn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Cn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_b
#    runTestEn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-En-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefCorrData5Leks,progDir,BUGSDir,'Test simCorr5Lek-Fn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
#    
#    runTestAn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-An-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Bn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_a
#    runTestCn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Cn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix mu_b
#    runTestDn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Dn-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop mu_a, mu_b correlation
#    runTestEn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-En-4k-316k-1',mZone,316000,320000,1)   # all zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefUncorrData31Leks,progDir,BUGSDir,'Test simUncorr31Lek-Fn-4k-316k-1',mZone,316000,320000,1)   # all zeros - fix 0.5*sigma^2_{\epsilon}
#    
#    runTestAn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-An-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop 0.5*sigma^2_{\epsilon}
#    runTestBn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Bn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_a
#    runTestCn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Cn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix mu_b
#    runTestDn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Dn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop mu_a, mu_b correlation
#    runTestEn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-En-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - drop mu_a altogether (regression through the origin)
#    runTestFn(simPoisRandCoefUncorrData5Leks,progDir,BUGSDir,'Test simUncorr5Lek-Fn-4k-316k-1',mZone,316000,320000,1)   # 1st zeros - fix 0.5*sigma^2_{\epsilon} 
   
   
   
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




mz1datAllZerosCore10 <- vector("list",41)         # read in all zeros, core data, ith mzone
mz1datAllZerosNoco10 <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz1datAllZerosLeks10 <- vector("list",41)         # read in all zeros, all data, ith mzone
mz1dat1stZerosCore10 <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz1dat1stZerosNoco10 <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz1dat1stZerosLeks10 <- vector("list",41)         # read in 1st zeros, all data, ith mzone

mz3datAllZerosCore10[[mZone]] <- vector("list",41)         # read in all zeros, core data, ith mzone
mz3datAllZerosNoco10[[mZone]] <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz3atAllZerosLeks10[[mZone]] <- vector("list",41)         # read in all zeros, all data, ith mzone
mz31dat1stZerosCore10[[mZone]] <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz3dat1stZerosNoco10[[mZone]] <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz3dat1stZerosLeks10[[mZone]] <- vector("list",41)         # read in 1st zeros, all data, ith mzone

mz4datAllZerosCore10[[mZone]] <- vector("list",41)         # read in all zeros, core data, ith mzone
mz4datAllZerosNoco10[[mZone]] <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz4datAllZerosLeks10[[mZone]] <- vector("list",41)         # read in all zeros, all data, ith mzone
mz4dat1stZerosCore10[[mZone]] <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz4dat1stZerosNoco10[[mZone]] <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz4dat1stZerosLeks10[[mZone]] <- vector("list",41)         # read in 1st zeros, all data, ith mzone

mz5datAllZerosCore10[[mZone]] <- vector("list",41)         # read in all zeros, core data, ith mzone
mz5datAllZerosNoco10[[mZone]] <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz5datAllZerosLeks10[[mZone]] <- vector("list",41)         # read in all zeros, all data, ith mzone
mz5dat1stZerosCore10[[mZone]] <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz5dat1stZerosNoco10[[mZone]] <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz5dat1stZerosLeks10[[mZone]] <- vector("list",41)         # read in 1st zeros, all data, ith mzone

mz6datAllZerosCore10[[mZone]] <- vector("list",41)         # read in all zeros, core data, ith mzone
mz6datAllZerosNoco10[[mZone]] <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz6datAllZerosLeks10[[mZone]] <- vector("list",41)         # read in all zeros, all data, ith mzone
mz6dat1stZerosCore10[[mZone]] <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz6dat1stZerosNoco10[[mZone]] <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz6dat1stZerosLeks10[[mZone]] <- vector("list",41)         # read in 1st zeros, all data, ith mzone

mz8datAllZerosCore10[[mZone]] <- vector("list",41)         # read in all zeros, core data, ith mzone
mz8datAllZerosNoco10[[mZone]] <- vector("list",41)         # read in all zeros, non-core data, ith mzone
mz8datAllZerosLeks10[[mZone]] <- vector("list",41)         # read in all zeros, all data, ith mzone
mz8dat1stZerosCore10[[mZone]] <- vector("list",41)         # read in 1st zeros, core data, ith mzone
mz8dat1stZerosNoco10[[mZone]] <- vector("list",41)         # read in 1st zeros, non-core data, ith mzone
mz8dat1stZerosLeks10[[mZone]] <- vector("list",41)         # read in 1st zeros, all data, ith mzone

# MAKE 10 YEAR LOOK BACKS
for(i in 1:7){
  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  if(i == 7){mZone <- 9}
    
  # build correct data sets.
  datAllZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
  datAllZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
  datAllZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - All Zero'))@data        # read in all zeros, all data, ith mzone
    
  dat1stZerosCore[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
  dat1stZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
  dat1stZerosLeks[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
    
 #BUGSDir <- 'C:/Program Files/WinBUGS14'                                                              # my machine  
 #BUGSDir <- 'C:/Users/jmitchell/WinBUGS_14/winbugs14/WinBUGS14'                                       # sonic 1
  BUGSDir <- 'C:/Users/jmitchell/WinBUGS14'                                                            # sonic 2 & others
    
  for(j in 1:41){
    
    tenAllZerosCore <- datAllZerosCore[[mZone]][datAllZerosCore[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]
    tenAllZerosNoco <- datAllZerosNoco[[mZone]][datAllZerosNoco[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]    
    tenAllZerosLeks <- datAllZerosLeks[[mZone]][datAllZerosLeks[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]    
      
    ten1stZerosCore <- dat1stZerosCore[[mZone]][dat1stZerosCore[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]      
    ten1stZerosNoco <- dat1stZerosNoco[[mZone]][dat1stZerosNoco[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]      
    ten1stZerosLeks <- dat1stZerosLeks[[mZone]][dat1stZerosLeks[[mZone]]$Year %in% seq(1964 + j,1964 + j + 10),]      
      
    if(mZone <= 8){   
        
      runMD10(tenAllZerosCore,progDir,BUGSDir,paste0('All Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1)            # all zeros, core,     
      runME10(tenAllZerosNoco,progDir,BUGSDir,paste0('All Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1)            # all zeros, non-core
      runMF10(tenAllZerosLeks,progDir,BUGSDir,paste0('All Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1)            # all zeros, all leks
        
      runMD10(ten1stZerosCore,progDir,BUGSDir,paste0('1st Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1,)           # 1st zeros, core
      runME10(ten1stZerosNoco,progDir,BUGSDir,paste0('1st Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1,)           # 1st zeros, non-core
      runMF10(ten1stZerosLeks,progDir,BUGSDir,paste0('1st Zeros ',1964 + j,'-',1964 + j + 10),mZone,j,76000,80000,1,)           # 1st zeros, all leks    
    } 
    rm(tenAllZerosCore,tenAllZerosNoco,tenAllZerosLeks,ten1stZerosCore,ten1stZerosNoco,ten1stZerosLeks)
  }
}