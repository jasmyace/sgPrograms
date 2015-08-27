makeTrendPlots <- function(dat,runType,theUnit,nZones,tracDir,files,bayes){

#   runType <- 'Core'
#   theUnit <- 'MZone 1'
#   nZones <- 1
#   tracDir <- paste0(outpDir,'/',file,'/Trend Plots')
#   
#   file                                             ,dat                            ,string           ,theUnit  ,runType
#   "Model D MZone 6 Test simUncorr31Lek-Dn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core"
#   
#   
#   
#   file <- "Model D MZone 6 Test simUncorr31Lek-Dn-4k-76k-1"
#   dat <- simPoisRandCoefUncorrData31Leks
#   string <- 'Management Zone'
#   theUnit <- 'mZone 6'
#   runType <- "Core"
#   
#   
#   
#   dat <- dat
#   runType <- runType
#   theUnit <- paste0(" - ",string," ",theUnit)
#   nZones <- 1
#   tracDir <- paste0(outpDir,'/',file,'/Trend Plots')
#   files <- file
#   bayes <- bayes
#   
  
  
  
  
  if(min(dat$Year) == 1){
    dat$Year <- dat$Year + 1964
  }
  
  
  bsums <- bayes$summary
  Year <- data.frame(Year=seq(1965,2015,1))
  mu.a <- bayes$summary[rownames(data.frame(bayes$summary)) == 'mu.a',][1]
  mu.b <- bayes$summary[rownames(data.frame(bayes$summary)) == 'mu.b',][1]
  sdnoise <- bayes$summary[rownames(data.frame(bayes$summary)) == 'sdnoise',][1]
  
  indx <- data.frame(bsums)
  indx$vars <- rownames(indx)
  indx$x <- ifelse(substr(indx$vars,1,1) == 'N',substr(indx$vars,3,3),-99)
  indx$y <- ifelse(substr(indx$vars,1,1) == 'N',ifelse(nchar(indx$vars) == 6,substr(indx$vars,5,5),substr(indx$vars,5,6)),-99)
  indx$x0 <- ifelse(substr(indx$vars,1,2) == 'N0',ifelse(nchar(indx$vars) == 5,substr(indx$vars,4,4),substr(indx$vars,4,5)),-99)
  
  # build the true B-matrix
  nTimes <- length(unique(dat$Year))
  beta.mzone <- matrix(NA,nrow=nZones,ncol=nTimes)
  list.beta.mzone <- indx[substr(indx$vars,1,2) == 'N[',]
  for(i in 1:nZones){
    for(j in 1:nTimes){
      beta.mzone[i,j] <- list.beta.mzone[list.beta.mzone$x == i & list.beta.mzone$y == j,]$mean
    }
  }
  
  coreText <- 'ugh'

  if(length(unique(dat$Year)) == 11){
    minYear <- min(dat$Year)
    medYear <- quantile(dat$Year,0.5) - minYear
    maxYear <- max(dat$Year)
  } else {
    minYear <- 1965
    medYear <- 26
    maxYear <- 2015
  }
  
  theN <- data.frame(Year=seq(minYear,maxYear,1),N=beta.mzone[i,])  
  obsMeans <- data.frame(MeanPMales=tapply(dat$Peak_Males, list(dat$Year), mean))
  obsMeans$Year <- rownames(obsMeans)
  
  plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
  plotYears$runType <- ''
  plotYears$numMZone <- theUnit 
  plotYears$dataCut <- runType
  plotYears$Nhat <- exp(mu.a + mu.b*(as.numeric(plotYears$Year) - (minYear - 1) - medYear) + 0.5*sdnoise*sdnoise)
    
  x  <- plotYears$Year
  y1 <- plotYears$N
  y2 <- plotYears$MeanPMales
  y3 <- plotYears$Nhat
    
  yMax <- 100
    
  png(filename=paste0(tracDir,'/Trend Plot - ',files,' - ',runType,'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years ",minYear,"-",maxYear,"\n75% ",runType,theUnit))
  par(new=TRUE)
  plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  par(new=TRUE)
  plot(x,y3,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  axis(side=1,labels=TRUE,seq(1965,2015,5))
  axis(side=2,labels=TRUE,seq(0,yMax,10))
    
  dev.off()
  
  write.csv(plotYears,paste0(tracDir,'/plotYears - ',files,' - ',runType,'.csv'))
  rm(plotYears,x,y1,y2)
}