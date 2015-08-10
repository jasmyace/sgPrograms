makeRangeWide <- function(runType,units,datList,runType2,tracDir){
  
#   runType <- "Model D"
#   datList <- dat1stZerosCore
#   runType2 <- "Core"
#   tracDir <- paste0(outpDir,'/','Rangewide','/Trend Plots')
  
  
  
  
  # build lek weight matrix for each mzone
  dat <- datList[[9]]
  dat$Mgmt_zone <- as.factor(ifelse(dat$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat$Mgmt_zone))))

  dat$zoneCls = as.numeric(as.factor(droplevels(dat$Mgmt_zone)))                   # management zone
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  zoneCls <- sort(rep(unique(dat$zoneCls),length(sort(unique(dat$yearCls)))))
  yearCls <- rep(sort(unique(dat$yearCls)),length(unique(dat$zoneCls)))
  BackBone <- as.data.frame(cbind(zoneCls,yearCls))
  tempAgg <- aggregate(data=dat, lekCls ~ zoneCls + yearCls, function(x) length(unique(x)))
  preLekMat <- merge(BackBone,tempAgg,by=c('zoneCls','yearCls'),all.x=TRUE)
  preLekMat <- preLekMat[order(preLekMat$zoneCls,preLekMat$yearCls),]
  preLekMat[is.na(preLekMat)] <- 0
  
  lekMat <- matrix(NA,nrow=length(unique(dat$yearCls)),ncol=length(unique(dat$zoneCls)))
  mZoneClasses <- sort(unique(dat$zoneCls))
  for(i in 1:length(unique(zoneCls))){
    if(i == 1){
      getMZone <-mZoneClasses[1] 
      lekMat[,1] <- preLekMat[preLekMat$zoneCls == getMZone,]$lekCls
    } else {
      getMZone <-mZoneClasses[i] 
      lekMat[,i] <- preLekMat[preLekMat$zoneCls == getMZone,]$lekCls
    }
  }
  nLeksMat <- t(lekMat)
  nLeksVec <- colSums(nLeksMat)
  
  
  
  
  
  
  
  
  # get bsims for each
  bsimsList <-  vector("list",9)
  allTimes <- NULL
  for(x in 1:6){
  
    theUnit <- units[x]
    numUnit <- as.numeric(substr(theUnit,7,7))                                               # applicable only to mZone?
    
    load(paste0(outpDir,'/',runType,' MZone ',numUnit,' Try ',1,".RData")) # generalize the '1' in 'Try 1' later
    
    rTime <- data.frame(time=time,group=theUnit)
    allTimes <- rbind(allTimes,rTime)
    bsimsList[[numUnit]]  <- bayes$sims.list
    
    rm(data,inits,parameters,bayes,CYear,freqs,time)
    
  }
  
  
  
  
  
  
  getSampleN <- function(metric,bsimsList){
    NmatList <- vector("list",51)
    wgtMeanNum <- vector("list",51)
    wgtMeanDen <- vector("list",51)
    wgtMean <- vector("list",51)
    allMeans <- NULL
    for(x in 1:51){  
      NmatList[[x]] <- cbind(bsimsList[[1]][metric][[1]][,1,x],bsimsList[[3]][metric][[1]][,1,x],bsimsList[[4]][metric][[1]][,1,x],bsimsList[[5]][metric][[1]][,1,x],bsimsList[[6]][metric][[1]][,1,x],bsimsList[[8]][metric][[1]][,1,x]) 
      wgtMeanNum[[x]] <- NmatList[[x]] %*% nLeksMat[,x]
      wgtMeanDen[[x]] <- sum(nLeksMat[,x])
      wgtMean[[x]] <- wgtMeanNum[[x]] / wgtMeanDen[[x]]
      thisMean <- data.frame(mean=mean(wgtMean[[x]]),sd=sd(wgtMean[[x]]),x2.5=quantile(wgtMean[[x]],0.025),x5=quantile(wgtMean[[x]],0.05),x25=quantile(wgtMean[[x]],0.25),x50=quantile(wgtMean[[x]],0.5),x75=quantile(wgtMean[[x]],0.75),x95=quantile(wgtMean[[x]],0.95),x97.5=quantile(wgtMean[[x]],0.975))
      allMeans <- rbind(allMeans,thisMean)
    }
    colnames(allMeans) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')
    rownames(allMeans) <- paste0(metric,'0[1,',seq(1,51,1),']') 
    allMeans <- list(allMeans,wgtMean)
  }
    
  getSampleNotN <- function(metric,bsimsList){
    NmatList <- vector("list",1)
    wgtMeanNum <- vector("list",1)
    wgtMeanDen <- vector("list",1)
    wgtMean <- vector("list",1)
    allMeans <- NULL

    NmatList[[1]] <- cbind(bsimsList[[1]][metric][[1]],bsimsList[[3]][metric][[1]],bsimsList[[4]][metric][[1]],bsimsList[[5]][metric][[1]],bsimsList[[6]][metric][[1]],bsimsList[[8]][metric][[1]])         
    wgtMeanNum[[1]] <- NmatList[[1]] %*% nLeksMat[,51]
    wgtMeanDen[[1]] <- sum(nLeksMat[,51])
    wgtMean[[1]] <- wgtMeanNum[[1]] / wgtMeanDen[[1]]
    allMeans <- data.frame(mean=mean(wgtMean[[1]]),sd=sd(wgtMean[[1]]),x2.5=quantile(wgtMean[[1]],0.025),x5=quantile(wgtMean[[1]],0.05),x25=quantile(wgtMean[[1]],0.25),x50=quantile(wgtMean[[1]],0.5),x75=quantile(wgtMean[[1]],0.75),x95=quantile(wgtMean[[1]],0.95),x97.5=quantile(wgtMean[[1]],0.975))

    names(allMeans) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')
    rownames(allMeans) <- paste0(metric,'0') 
    allMeans
  }
  
  N0 <- getSampleN("N",bsimsList)[[1]]
  mu.a0 <- getSampleNotN("mu.a",bsimsList)
  mu.b0 <- getSampleNotN("mu.b",bsimsList)
  sdnoise0 <- getSampleNotN("sdnoise",bsimsList)
  rho0 <- getSampleNotN("rho",bsimsList)
  sigma.a0 <- getSampleNotN("sigma.a",bsimsList)
  sigma.b0 <- getSampleNotN("sigma.b",bsimsList)  
  theNforB <- getSampleN("N",bsimsList)[[2]]
  
  # function for making posterior of Bs
  makeB <- function(theNforB,tb,ta,label){
    preBeta0 <- (theNforB[[tb]] / theNforB[[ta]])^(1/(tb-ta))    # note that it's 50 because it's 51 - 1.
    beta0 <- data.frame(mean=mean(preBeta0),sd=sd(preBeta0),x2.5=quantile(preBeta0,0.025),x5=quantile(preBeta0,0.05),x25=quantile(preBeta0,0.25),x50=quantile(preBeta0,0.5),x75=quantile(preBeta0,0.75),x95=quantile(preBeta0,0.95),x97.5=quantile(preBeta0,0.975))
    names(beta0) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')  
    rownames(beta0) <- label   
    beta0
  }
  
  beta0 <- makeB(theNforB,51,1,'beta0')
  B0.10.15.05 <- makeB(theNforB,51,41,'B0.10.15.05')
  B0.10.14.04 <- makeB(theNforB,50,40,'B0.10.14.04')
  B0.10.13.03 <- makeB(theNforB,49,39,'B0.10.13.03')
  B0.10.12.02 <- makeB(theNforB,48,38,'B0.10.12.02')
  B0.10.11.01 <- makeB(theNforB,47,37,'B0.10.11.01')
  B0.10.10.00 <- makeB(theNforB,46,36,'B0.10.10.00')
  B0.10.09.99 <- makeB(theNforB,45,35,'B0.10.09.99')
  B0.10.08.98 <- makeB(theNforB,44,34,'B0.10.08.98')
  B0.10.07.97 <- makeB(theNforB,43,33,'B0.10.07.97')
  B0.10.06.96 <- makeB(theNforB,42,32,'B0.10.06.96')
  B0.10.05.95 <- makeB(theNforB,41,31,'B0.10.05.95')
  B0.10.04.94 <- makeB(theNforB,40,30,'B0.10.04.94')
  B0.10.03.93 <- makeB(theNforB,39,29,'B0.10.03.93')
  B0.10.02.92 <- makeB(theNforB,38,28,'B0.10.02.92')
  B0.10.01.91 <- makeB(theNforB,37,27,'B0.10.01.91')
  B0.10.00.90 <- makeB(theNforB,36,26,'B0.10.00.90')
  B0.10.99.89 <- makeB(theNforB,35,25,'B0.10.99.89')
  B0.10.98.88 <- makeB(theNforB,34,24,'B0.10.98.88')
  B0.10.97.87 <- makeB(theNforB,33,23,'B0.10.97.87')
  B0.10.96.86 <- makeB(theNforB,32,22,'B0.10.96.86')
  B0.10.95.85 <- makeB(theNforB,31,21,'B0.10.95.85')
  B0.10.94.84 <- makeB(theNforB,30,20,'B0.10.94.84')
  B0.10.93.83 <- makeB(theNforB,29,19,'B0.10.93.83')
  B0.10.92.82 <- makeB(theNforB,28,18,'B0.10.92.82')
  B0.10.91.81 <- makeB(theNforB,27,17,'B0.10.91.81')
  B0.10.90.80 <- makeB(theNforB,26,16,'B0.10.90.80')
  B0.10.89.79 <- makeB(theNforB,25,15,'B0.10.89.79')
  B0.10.88.78 <- makeB(theNforB,24,14,'B0.10.88.78')
  B0.10.87.77 <- makeB(theNforB,23,13,'B0.10.87.77')
  B0.10.86.76 <- makeB(theNforB,22,12,'B0.10.86.76')
  B0.10.85.75 <- makeB(theNforB,21,11,'B0.10.85.75')
  B0.10.84.74 <- makeB(theNforB,20,10,'B0.10.84.74')
  B0.10.83.73 <- makeB(theNforB,19, 9,'B0.10.83.73')
  B0.10.82.72 <- makeB(theNforB,18, 8,'B0.10.82.72')
  B0.10.81.71 <- makeB(theNforB,17, 7,'B0.10.81.71')
  B0.10.80.70 <- makeB(theNforB,16, 6,'B0.10.80.70')
  B0.10.79.69 <- makeB(theNforB,15, 5,'B0.10.79.69')
  B0.10.78.68 <- makeB(theNforB,14, 4,'B0.10.78.68')
  B0.10.77.67 <- makeB(theNforB,13, 3,'B0.10.77.67')
  B0.10.76.66 <- makeB(theNforB,12, 2,'B0.10.76.66')
  B0.10.75.65 <- makeB(theNforB,11, 1,'B0.10.75.65')

  # put all the estimates together
  rangewide <- rbind(N0,mu.a0,mu.b0,beta0,sdnoise0,rho0,sigma.a0,sigma.b0,
                     B0.10.15.05,B0.10.14.04,B0.10.13.03,B0.10.12.02,B0.10.11.01,B0.10.10.00,B0.10.09.99,B0.10.08.98,B0.10.07.97,B0.10.06.96,
                     B0.10.05.95,B0.10.04.94,B0.10.03.93,B0.10.02.92,B0.10.01.91,B0.10.00.90,B0.10.99.89,B0.10.98.88,B0.10.97.87,B0.10.96.86,
                     B0.10.95.85,B0.10.94.84,B0.10.93.83,B0.10.92.82,B0.10.91.81,B0.10.90.80,B0.10.89.79,B0.10.88.78,B0.10.87.77,B0.10.86.76,
                     B0.10.85.75,B0.10.84.74,B0.10.83.73,B0.10.82.72,B0.10.81.71,B0.10.80.70,B0.10.79.69,B0.10.78.68,B0.10.77.67,B0.10.76.66,
                     B0.10.75.65)
  
 
  
  # make graph
  Year <- data.frame(Year=seq(1965,2015,1))
  mu.a <- mu.a0$mean
  mu.b <- mu.b0$mean
  sdnoise <- sdnoise0$mean
  
  indx <- data.frame(rangewide)
  indx$vars <- rownames(indx)
  indx$x <- ifelse(substr(indx$vars,1,2) == 'N0',substr(indx$vars,4,4),-99)
  indx$y <- ifelse(substr(indx$vars,1,2) == 'N0',ifelse(nchar(indx$vars) == 7,substr(indx$vars,6,6),substr(indx$vars,6,7)),-99)
  
  # build the true B-matrix
  beta.mzone <- matrix(NA,nrow=1,ncol=51)
  list.beta.mzone <- indx[substr(indx$vars,1,3) == 'N0[',]
  for(i in 1:1){
    for(j in 1:51){
      beta.mzone[i,j] <- list.beta.mzone[list.beta.mzone$x == i & list.beta.mzone$y == j,]$mean
    }
  }
  
  coreText <- 'ugh'
  
  theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])  
  obsMeans <- data.frame(MeanPMales=tapply(dat$Peak_Males, list(dat$Year), mean))
  obsMeans$Year <- rownames(obsMeans)
  
  plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
  plotYears$runType <- runType
  plotYears$numMZone <- 'Rangewide'
  plotYears$dataCut <- runType2
  plotYears$Nhat <- exp(mu.a + mu.b*(as.numeric(plotYears$Year) - 1964 - 26) + 0.5*sdnoise*sdnoise)
  
  x  <- plotYears$Year
  y1 <- plotYears$N
  y2 <- plotYears$MeanPMales
  y3 <- plotYears$Nhat
  
  yMax <- 100
  
  png(filename=paste0(tracDir,'/Trend Plot - Rangewide - ',runType2,'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Rangewide Temporal Trend of Observed Peak Males\nYears 1965-2015 - 75% ",runType2))
  par(new=TRUE)
  plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  par(new=TRUE)
  plot(x,y3,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  axis(side=1,labels=TRUE,seq(1965,2015,5))
  axis(side=2,labels=TRUE,seq(0,yMax,10))
  
  dev.off()
  
  write.csv(plotYears,paste0(tracDir,'/plotYears - Rangewide - ',runType2,'.csv'))
  write.csv(rangewide,paste0(tracDir,'/bayesSummary - ',runType,'.csv'))
  rm(plotYears,x,y1,y2)
  
}