


# loadThis <- 'Model S MZone 9 All MZones Core, Try 2'
# loadThis <- 'Model S MZone 9 All MZones Non-Core, Try 2'
# 
# loadThis <- 'Model S MZone 9 MZones 1,3,5 Core, Try 1'
# loadThis <- 'Model T MZone 9 MZones 1,3,5 Core, Adj M,Try 1'
# 
# loadThis <- 'Model U MZone 9 All MZones Non-Core, Adj M,Try 1'
# loadThis <- 'Model U MZone 9 All MZones Core, Adj M,Try 1'

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("RData", theFiles) == 1]
colVec <- brewer.pal(8,"Set1")


states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")

nStates <- length(states)
nMZones <- length(mZones)


readInAnalyticFiles(analDir)





for(h in 1:2){
  if(h == 1){
    units <- mZones
    nUnits <- nMZones
    start <- 1
    string <- 'Management Zone'
  } else {
    units <- states
    nUnits <- nStates
    start <- 3
    string <- 'State'
  }

  # loop over mZone files in theFiles
  for(i in 1:nUnits){
  
    theUnit <- units[i]
    numUnit <- as.numeric(substr(theUnit,7,7))                                               # applicable only to mZone?
    unitFile <- theFiles[grepl(units[i],theFiles)]
  
    for(j in start:3){
      
      # get data for analysis
      if(h == 1){
        load(paste0(outpDir,'/',unitFile[j]))                                                # get bayesian stuff -- mzone
        #load(paste0(outpDir,"/Model V MZone 6 ZInf All Leks Try 2",".RData"))
      } else {
        load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
          dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type      
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
          dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type      
      }
      
      dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
      
      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else {
        file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      }

      ifelse(!dir.exists(file.path(outpDir,file)), dir.create(file.path(outpDir,file)), FALSE)      # make new folder
      file.copy(paste0(outpDir,"/",file,".RData"),paste0(outpDir,"/",file))                         # copy bayes output to new folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), FALSE)        # make new trace plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), FALSE)# make new posterior plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), FALSE)        # make new mzone trend plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), FALSE)        # make new mzone trend plots folder
      
      # make trace, posterior, histogram plots
      nParms <- dim(bayes$sims.array)[3]
      parmList <- dimnames(bayes$sims.array)[[3]]
#       makeTracePlots(nParms,parmList,paste0(outpDir,'/',file,'/Trace Plots'),file,bayes)    
#       makePosteriorPlots(nParms,parmList,paste0(outpDir,'/',file,'/Posterior Plots'),file,bayes)
#       makeHistogramPlots(dat,paste0(" - ",string," ",theUnit),paste0(outpDir,'/',file,'/Zeros Plots'))
      
      
      bsums90 <- make90pCredInt(bayes)
      
      
      # make bayes summary file of estimates
      write.csv(bsums90,paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
      
      # make trend plots
#       makeTrendPlots(dat,runType,paste0(" - ",string," ",theUnit),1,paste0(outpDir,'/',file,'/Trend Plots'),file,bayes)
  
      rm(nParms,parmList)
      
      
#       if(j == 3){
#         ifelse(!dir.exists(file.path(outpDir,'Rangewide')), dir.create(file.path(outpDir,'Rangewide')), FALSE)      # make new folder
#         ifelse(!dir.exists(file.path(paste0(outpDir,'/','Rangewide','/Trend Plots'))), dir.create(file.path(paste0(outpDir,'/','Rangewide','/Trend Plots'))), FALSE)        # make new mzone trend plots folder
#         
#         # make rangewide estimates
#         makeRangeWide("Model D",units,dat1stZerosCore,"Core",paste0(outpDir,'/Rangewide/Trend Plots'))
#         makeRangeWide("Model E",units,dat1stZerosNoco,"Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
#         makeRangeWide("Model F",units,dat1stZerosLeks,"All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
#         # make core + non-core + all lek plot
#       } 
    }
  }
}




# make a big table of results

theSumms <- NULL
# make table
for(h in 1:1){
  if(h == 1){
    units <- mZones
    nUnits <- nMZones
    start <- 1
    string <- 'Management Zone'
  } else {
    units <- states
    nUnits <- nStates
    start <- 3
    string <- 'State'
  }
  
  # loop over mZone files in theFiles
  for(i in 1:nUnits){
    
    theUnit <- units[i]
    numUnit <- as.numeric(substr(theUnit,7,7))                                               # applicable only to mZone?
    unitFile <- theFiles[grepl(units[i],theFiles)]
    
    for(j in start:3){
      
      # get data for analysis
      if(h == 1){
        load(paste0(outpDir,'/',unitFile[j]))                                                # get bayesian stuff -- mzone
        #load(paste0(outpDir,"/Model V MZone 6 ZInf All Leks Try 2",".RData"))
      } else {
        load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
          dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type      
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
          dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type      
      }
      
      dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
      
      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else {
        file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      }
      
      thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
      mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.a'   ,mean=thisOne[thisOne$X == 'mu.a'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a'   ,]$X95.)
      mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.b'   ,mean=thisOne[thisOne$X == 'mu.b'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b'   ,]$X95.)
      beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='beta'   ,mean=thisOne[thisOne$X == 'beta[1]',]$mean,X5.=thisOne[thisOne$X == 'beta[1]',]$X5.,X95.=thisOne[thisOne$X == 'beta[1]',]$X95.)
      sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sdnoise',mean=thisOne[thisOne$X == 'sdnoise',]$mean,X5.=thisOne[thisOne$X == 'sdnoise',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise',]$X95.)
      rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='rho'    ,mean=thisOne[thisOne$X == 'rho'    ,]$mean,X5.=thisOne[thisOne$X == 'rho'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho'    ,]$X95.)
      sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.a',mean=thisOne[thisOne$X == 'sigma.a',]$mean,X5.=thisOne[thisOne$X == 'sigma.a',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a',]$X95.)
      sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.b',mean=thisOne[thisOne$X == 'sigma.b',]$mean,X5.=thisOne[thisOne$X == 'sigma.b',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b',]$X95.)
      
      summs <- rbind(mu.a,mu.b,beta,sdnoise,rho,sigma.a,sigma.b)
      
      theSumms <- rbind(theSumms,summs)
       
    }
  }
}

theSumms <- theSumms[order(theSumms$focus,theSumms$stat,theSumms$group),]
theSumms$mean <- ifelse(theSumms$stat %in% c('mu.a','mu.b'),print(formatC(signif(exp(theSumms$mean),digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(theSumms$mean,digits=3), digits=3,format="fg", flag="#")))
theSumms$X5.  <- ifelse(theSumms$stat %in% c('mu.a','mu.b'),print(formatC(signif(exp(theSumms$X5.) ,digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(theSumms$X5. ,digits=3), digits=3,format="fg", flag="#")))
theSumms$X95. <- ifelse(theSumms$stat %in% c('mu.a','mu.b'),print(formatC(signif(exp(theSumms$X95.),digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(theSumms$X95.,digits=3), digits=3,format="fg", flag="#")))


coreSumms <- theSumms[theSumms$cut == 'Core',]
names(coreSumms)[names(coreSumms) == 'mean'] <- 'Coremean'
names(coreSumms)[names(coreSumms) == 'X5.'] <- 'CoreX5.'
names(coreSumms)[names(coreSumms) == 'X95.'] <- 'CoreX95.'
names(coreSumms)[names(coreSumms) == 'Nleks'] <- 'CoreNleks'
coreSumms$cut <- NULL

nocoSumms <- theSumms[theSumms$cut == 'Non-Core',]
names(nocoSumms)[names(nocoSumms) == 'mean'] <- 'Nocomean'
names(nocoSumms)[names(nocoSumms) == 'X5.'] <- 'NocoX5.'
names(nocoSumms)[names(nocoSumms) == 'X95.'] <- 'NocoX95.'
names(nocoSumms)[names(nocoSumms) == 'Nleks'] <- 'NocoNleks'
nocoSumms$stat <- NULL
nocoSumms$cut <- NULL
nocoSumms$focus <- NULL
nocoSumms$group <- NULL

leksSumms <- theSumms[theSumms$cut == 'All Leks',]
names(leksSumms)[names(leksSumms) == 'mean'] <- 'Leksmean'
names(leksSumms)[names(leksSumms) == 'X5.'] <- 'LeksX5.'
names(leksSumms)[names(leksSumms) == 'X95.'] <- 'LeksX95.'
names(leksSumms)[names(leksSumms) == 'Nleks'] <- 'LeksNleks'
leksSumms$stat <- NULL
leksSumms$cut <- NULL
leksSumms$focus <- NULL
leksSumms$group <- NULL

table1 <- cbind(coreSumms,nocoSumms,leksSumms)
table1 <- table1[,c('focus','group','stat','CoreNleks','Coremean','CoreX5.','CoreX95.','NocoNleks','Nocomean','NocoX5.','NocoX95.','LeksNleks','Leksmean','LeksX5.','LeksX95.')]


# make plot of B51 trend with 90% cred int (core, non-core, all leks) for each mzone & state
# make plot of all states together
# make plot of all mzones together (core, non-core, all leks)
# make B10 plots

# compare of B51s to betas
# time elapsed analysis
















# load in data for sample means 
dat <- dat1stZerosCore[[9]]    # loop over cores?
#dat <- dat1stZerosNoco[[9]]  

# load(paste0(outpDir,'/smallCoreSamp2.Sonic2.RData'))
# dat <- smallCoreSamp2



# nZones <- length(unique(dat$mZone_num))


# load(paste0(outpDir,'/',loadThis,'.RData')) 









# try2 <- bayes$summary
# try3 <- bayes$summary


# load model 
# these <- grep(paste0("Model ",modelLetter), theFiles , ignore.case=FALSE, fixed=TRUE)
# zoneString <- unlist(lapply(strsplit(theFiles,".",fixed=TRUE),function(x) strsplit(x,".",fixed=TRUE)[[1]][1]))
# zones <- as.numeric(substr(zoneString,nchar(zoneString),nchar(zoneString)))
# loadThese <- paste0(outpDir,"/",theFiles[these])
# assign(paste0("nModel",modelLetter),length(loadThese))
# 




# these <- 6
# load(paste0(outpDir,"/",theFiles[these]))  



Year <- data.frame(Year=seq(1965,2015,1))

indx <- data.frame(bayes$summary)
indx$vars <- rownames(indx)
indx$x <- ifelse(substr(indx$vars,1,1) == 'N',substr(indx$vars,3,3),-99)
indx$y <- ifelse(substr(indx$vars,1,1) == 'N',ifelse(nchar(indx$vars) == 6,substr(indx$vars,5,5),substr(indx$vars,5,6)),-99)
indx$x0 <- ifelse(substr(indx$vars,1,2) == 'N0',ifelse(nchar(indx$vars) == 5,substr(indx$vars,4,4),substr(indx$vars,4,5)),-99)
  
# build the true B-matrix
beta.mzone <- matrix(NA,nrow=nZones,ncol=51)
list.beta.mzone <- indx[substr(indx$vars,1,2) == 'N[',]
for(i in 1:nZones){
  for(j in 1:51){
    beta.mzone[i,j] <- list.beta.mzone[list.beta.mzone$x == i & list.beta.mzone$y == j,]$mean
  }
}





if(grepl("Non-Core", loadThis) == 1){
  coreText <- "Non-Core"
} else {
  coreText <- "Core"
}

dev.off()
par(mfrow=c(3,2))

for(i in 1:nZones){
  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  
  theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])
  
  thisOne <- dat[dat$mZone_num == mZone,]
  obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
  obsMeans$Year <- rownames(obsMeans)
  
  plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
  
  x  <- plotYears$Year
  y1 <- plotYears$N
  y2 <- plotYears$MeanPMales
  
  yMax <- 100
  
  plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",coreText," - Management Zone ",mZone))
  par(new=TRUE)
  plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  axis(side=1,labels=TRUE,seq(1965,2015,5))
  axis(side=2,labels=TRUE,seq(0,yMax,10))
  
  rm(plotYears,x,y1,y2)
}








  
  
dev.off()
# investigate regionwide trends
# build the true B-matrix
beta.range <- rep(NA,51)
list.beta.range <- indx[substr(indx$vars,1,3) == 'N0[',]
for(k in 1:51){
  beta.range[k] <- list.beta.range[list.beta.range$x0 == k,]$mean
}

theN <- data.frame(Year=seq(1965,2015,1),N=beta.range)

thisOne <- dat
obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
obsMeans$Year <- rownames(obsMeans)

plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)

x  <- plotYears$Year
y1 <- plotYears$N
y2 <- plotYears$MeanPMales
i <- 1

yMax <- 100

plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\n75% ",coreText," Core - Rangewide"))
par(new=TRUE)
plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
axis(side=1,labels=TRUE,seq(1965,2015,5))
axis(side=2,labels=TRUE,seq(0,yMax,10))

rm(plotYears,x,y1,y2)  











B10 <- data.frame(bayes$summary[substr(rownames(bayes$summary),1,3) == 'B10',])

B10$mZone <- as.numeric(substr(rownames(B10),11,11))

for(i in 1:nZones){
  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  
  theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])
  
  thisOne <- dat[dat$mZone_num == mZone,]
  obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
  obsMeans$Year <- rownames(obsMeans)
  
  plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
  
  x  <- plotYears$Year
  y1 <- plotYears$N
  y2 <- plotYears$MeanPMales
  
  yMax <- 100
  
  plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% Core - Management Zone ",mZone))
  par(new=TRUE)
  plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  axis(side=1,labels=TRUE,seq(1965,2015,5))
  axis(side=2,labels=TRUE,seq(0,yMax,10))
  
  rm(plotYears,x,y1,y2)  
}


g <- smallCoreSamp2$Peak_Males[smallCoreSamp2$Peak_Males <= 100]
h <- hist(g, breaks=100, density=10, col="lightgray", xlab="Accuracy", main="Overall") 
xfit <- seq(min(g),max(g),length=40) 
yfit <- dpois(xfit,lambda=mean(g)) 
yfit <- yfit*diff(h$mids[1:2])*length(g) 
lines(xfit, yfit, col="black", lwd=2)

























allResults <- rbind(estsC,estsD,estsF)
rm(estsC,estsD,estsF)
names(allResults)[names(allResults) == 'Zone'] <- 'mZone'
allResults <- data.frame(allResults,ParmType=unlist(lapply(strsplit(as.character(droplevels(allResults$Parameter)),'[',fixed=TRUE), function(x) x[[1]][1])))

allResults <- allResults[order(allResults$mZone,allResults$modCode),]

mu.aAll <- allResults[allResults$Parameter == 'mu.a',]
mu.bAll <- allResults[allResults$Parameter == 'mu.b',]
rhosAll <- allResults[allResults$Parameter == 'rho',]
taonoiseAll <- allResults[allResults$Parameter == 'taonoise',]

SumData <- read.csv(paste0(analDir,'/SumData.csv'))

yMax <- max(datList[[1]]$Peak_Males)

helper <- allResults[,c('Model','Zeros','Cut','mZone','modCode')]
helper <- unique(helper)
rownames(helper) <- NULL





doThese <- c('D')
nDoThese <- length(doThese)
table1 <- NULL
colVec <- brewer.pal(8,"Set1")

for(i in 1:8){
  
  #A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
  #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
  #C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
  
  D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
  D <- smallCoreSamp
  #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
  #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
  
  for(j in 1:nDoThese){
    
    thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]
    
    if(doThese[j] %in% c('D','E','F')){
      thisOne <- get(doThese[j])
      #thisOne <- thisOne[thisOne$DupZero == 0,]  # dont think this is necessary when doing core
    } else {
      thisOne <- get(doThese[j])   
    }
    leks <- unique(thisOne$Lek_IDnum)
    thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
    thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID))                     
    nLeks <- length(leks)
  
    Model <- as.character(droplevels(thisHelper$Model))
    Zeros <- as.character(droplevels(thisHelper$Zeros))
    Cut   <- as.character(droplevels(thisHelper$Cut))
    mZone <- i
    
    # make individual lek trends
#     for(k in 1:nLeks){    
#       # do this for winbugs numbering (so on factor)
#       xlek <- thisOne[thisOne$lekCls == k,]$Year - 1990
#       ylek <- thisOne[thisOne$lekCls == k,]$Peak_Males
#       the.lek <- thisOne[thisOne$lekCls == k,]$Lek_ID[1]
#       lek.mu.a <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('a[',k,']'),]$mean
#       lek.mu.b <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('b[',k,']'),]$mean
#       ypred <- exp(lek.mu.a + lek.mu.b*xlek)
#       
#       CairoPNG(filename=paste0(rsltDir,'/Lek Trends/MZone ',i,'/',the.lek,' - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',mZone,'.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)
#       
#       plot(xlek + 1990,ylek,col='gray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5,pch=19)
#       par(new=TRUE)
#       plot(xlek + 1990,ypred,col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Lek ",the.lek," in Management Zone ",i,"\nYears 1965-2014"))
#       axis(side=1,labels=TRUE,seq(1965,2015,5))
#       axis(side=2,labels=TRUE,seq(0,yMax,10))
#       
#       dev.off()
#     }
    
   CairoPNG(filename=paste0(rsltDir,'/Trends + Yes Tao - No Leks - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',mZone,'.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)

    # make gray individual lek trends
    for(k in 1:nLeks){
      
      # do this on k numbering
      x <- thisOne[thisOne$Lek_IDnum == leks[k],]$Year
      y <- thisOne[thisOne$Lek_IDnum == leks[k],]$Peak_Males
            
      yMax <- 100#max(thisOne$Peak_Males)
      
#       if(k == 1){
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
#       } else if(k > 1 & k < nLeks){
#         par(new=TRUE)
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
#       } else if(k == nLeks){
#         par(new=TRUE)
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='Year',ylab='Peak Males',lwd=0.5)        
#       }
    }
     
#     get mean trend info for this mzone and model
    mu.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$mean
    lo.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X5.
    hi.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X95.
    mu.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$mean
    lo.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X5.
    hi.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X95. 
    tao     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$mean
    lo.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X5.
    hi.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X95.
    rho     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$mean
    lo.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X5.
    hi.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X95.
    sigma.a <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.a',]$mean
    sigma.b <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.b',]$mean

    thisRow <- data.frame(thisHelper,mu.a=mu.a,lo.a,hi.a,mu.b,lo.b,hi.b,tao,lo.tao,hi.tao,rho,lo.rho,hi.rho,sigma.a,sigma.b)
    table1 <- rbind(table1,thisRow)

    obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
    obsMeans$Year <- rownames(obsMeans)

    year <- seq(-25,25,1)
    est <- exp(mu.a + mu.b*year + 0.5*tao)
    se <- sqrt(sigma.a^2 + sigma.b^2 + 2*rho*sigma.a*sigma.b)
    lo <- exp(mu.a + mu.b*year + 0.5*tao) - 1.645*se 
    hi <- exp(mu.a + mu.b*year + 0.5*tao) + 1.645*se 
    x <- year + 1990
    
    #par(new=TRUE)
    plot(obsMeans$Year,obsMeans$MeanPMales,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\n75% Core - Management Zone ",i))
    par(new=TRUE)
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
    par(new=TRUE)
    plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
    par(new=TRUE)
    plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)    
    axis(side=1,labels=TRUE,seq(1965,2015,5))
    axis(side=2,labels=TRUE,seq(0,yMax,10))
    legend("topleft",bty = "n",legend=c("Annual Mean","Temporal Trend & 90% Credible Intervals","Individual Lek Counts"),pch=c(19,NA,NA),lwd=c(NA,2,1),col=c('darkgray',colVec[i],"gray95"))
    dev.off()
  }
}

write.csv(table1,paste0(rsltDir,'/table1 + Yes Tao.csv'))

table1 <- read.csv(paste0(rsltDir,'/table1 + No Tao.csv'))


CairoPNG(filename=paste0(rsltDir,'/Trends + No Tao - All MZones - 1st Zeros - Core Leks.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)

for(i in 1:8){
  thisRow <- table1[i,]
  year <- seq(-25,25,1)
  est <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao)
  se <- sqrt(thisRow$sigma.a^2 + thisRow$sigma.b^2 + 2*thisRow$rho*thisRow$sigma.a*thisRow$sigma.b)
  lo <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) - 1.645*se
  hi <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) + 1.645*se
  x <- year + 1990
  if(i <= 7){
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     par(new=TRUE)
#     plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     par(new=TRUE)
#     plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)   
    par(new=TRUE)
  } else if(i == 8){
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\nAll Management Zones"))
    axis(side=1,labels=TRUE,seq(1965,2015,5))
    axis(side=2,labels=TRUE,seq(0,yMax,10))
    legend("topright",bty="n",legend=c('I','II','III','IV','V','VI','VII','VIII'),col=colVec,lwd=rep(2,8))
  }

}
dev.off()





doThese <- c('A','C')
nDoThese <- length(doThese)
the.zeros <- NULL


# read in all the data
for(i in 9:9){
  
   A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
  #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
   C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
  
  #D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
  #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
  #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
}
                                                                                                                                                                             # close out if
states <- as.character(droplevels(unique(A$State)))
nStates <- length(states)

for(i in 1:nStates){

  for(j in 1:nDoThese){
    
#   thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]

    thisOne <- get(doThese[j])
    thisOne <- thisOne[thisOne$State == states[i],]
    
    thisOne$zeroPresent <- ifelse(thisOne$Peak_Males == 0,1,0)

    leks <- length(unique(thisOne$Lek_IDnum))
    leks.zero <- length(unique(thisOne[thisOne$zeroPresent == 1,]$Lek_IDnum))
    leks.zero.dupped <- length(unique(thisOne[thisOne$zeroPresent == 1 & thisOne$DupZero == 1,]$Lek_IDnum))
    
#     thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
#     thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID))                     
#     nLeks <- length(leks)
#     
#     Model <- as.character(droplevels(thisHelper$Model))
#     Zeros <- as.character(droplevels(thisHelper$Zeros))
#     Cut   <- as.character(droplevels(thisHelper$Cut))
#     mZone <- i
    
    thisOne$Lek_ID <- as.character(droplevels(thisOne$Lek_ID))

    zeroInv <- data.frame(NDuppedZeros=tapply(thisOne$DupZero, thisOne$Lek_ID, FUN=sum))
    zeroInv$Lek <- rownames(zeroInv) 
    the.zero.mean <- data.frame(Model=doThese[j],NLeks=leks,NLeksZero=leks.zero,NLeksZeroDupped=leks.zero.dupped,AvgDuppedZeroMeanLek=mean(zeroInv[zeroInv$NDuppedZeros > 0,]$NDuppedZeros))
    the.zero.mean$State <- states[i]
    the.zero.mean$PercentLeksWithDup <- round(100*the.zero.mean$NLeksZeroDupped / the.zero.mean$NLeksZero,2)
    the.zeros <- rbind(the.zeros,the.zero.mean)
    write.csv(the.zeros[the.zeros$Model == 'C',],paste0(rsltDir,'/the.zeros.modelC.csv'))
  }
}




