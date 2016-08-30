
require(R2WinBUGS)      # run WinBUGS
require(lme4)           # glmm for poisson
require(sp)             # read in spatial points
require(rgdal)          # read in and write shpfiles
require(rgeos)          # geometry processing
require(raster)         # read in rasters
require(maptools)       # spRbind
require(Cairo)          # nice output
require(RColorBrewer)   # nice colors

# you'll need to change these to whatever they need to be.  
# just dump the attachments in my email, in the folder structure as zipped,
# wherever you want them to go.  then just change the prefix of these file stems.  
# you can keep them all the same to keep things easy.  
progDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/sgPrograms/sentToPaul--2016.05.05"
outpDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/sgPrograms/sentToPaul--2016.05.05" 
analDir <- "//LAR-FILE-SRV/Data/Jason/sage grouse/sgPrograms/sentToPaul--2016.05.05/Spatial Analysis Sets"
#BUGSDir <- 'C:/WHEREVER THE WinBUGS14.exe FILE RESIDES.'
BUGSDir <- 'C:/WinBUGS14'   # jason's machine  

# read in the program that runs the first-zeros, periphery (non-core) model.  
source(paste0(progDir,"/runME.R"))

# constants for this small application.
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZone <- 6

# list of data and/or results.
dat1stZerosNoco <- vector("list",9)       # holds the shapefile
resultsE <- vector("list",8)              # holds the results

dat1stZerosNoco[[mZone]] <- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - 1st Zero'))@data     # read in 1st zeros, non-core data, 6th mzone (WA)
resultsE[[mZone]] <- runME(dat1stZerosNoco[[mZone]],progDir,BUGSDir,'Try 1',mZone)                    # 1st zeros, non-core, ind mzone, B mat



# get 90% credible intervals.  
# this is extracted from program make90pCredInt.R

bayes <- resultsE[[mZone]][[4]]
blist <- bayes$sims.list
bsums <- bayes$summary

# helper function
get5.95 <- function(metric){
  if(is.matrix(blist[[metric]])){         # matrix
    the.5.95 <- data.frame(X5.=apply(blist[[metric]],2,function(x) quantile(x,0.05)),X95.=apply(blist[[metric]],2,function(x) quantile(x,0.95)))
  } else if(is.array(blist[[metric]])){   # array
    the.5.95 <- data.frame(X5.=apply(blist[[metric]][,1,],2,function(x) quantile(x,0.05)),X95.=apply(blist[[metric]][,1,],2,function(x) quantile(x,0.95)))
  } else {                                # vector
    the.5.95 <- data.frame(X5.=quantile(blist[[metric]],0.05),X95.=quantile(blist[[metric]],0.95))      
  }
  the.5.95
}

ans5.95 <- NULL
for(x in 1:length(bayes$sims.list)){
  metric <- attributes(bayes$sims.list[x])$names[1]
  ans <- get5.95(metric)
  ans5.95 <- rbind(ans5.95,ans)
}

newAns <- cbind(bsums,ans5.95)
newAns <- newAns[,c('mean','sd','2.5%','X5.','25%','50%','75%','X95.','97.5%')]
names(newAns)[names(newAns) == 'X5.'] <- '5%'
names(newAns)[names(newAns) == 'X95.'] <- '95%'

newAns