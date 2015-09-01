makeRangeWide <- function(runType,zeros,units,y10,datList,runType2,tracDir){
  
#   runType <- "Model D"
#   zeros <- "1st"
#   y10 <- 'N'
#   datList <- dat1stZerosCore
#   runType2 <- "1st Zeros - Core"
#   tracDir <- paste0(outpDir,'/','Rangewide','/Trend Plots')

  

  
  # build lek weight matrix for each mzone
  if(y10 == 'N'){
    dat <- datList[[9]]
  } else {
    dat <- datList[[9]][datList[[9]]$Year %in% seq(1964 + 41,1964 + 41 + 10),]
  }
  dat$Mgmt_zone <- as.factor(ifelse(dat$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat$Mgmt_zone))))

  dat$zoneCls = as.numeric(as.factor(droplevels(dat$Mgmt_zone)))                   # management zone
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  nYears <- max(dat$yearCls)
  
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
    
    if(zeros == 'all' & y10 == 'N'){
      load(paste0(outpDir,'/',runType,' MZone ',numUnit,' All Zeros ',1,".RData")) # swap out " Try " with " All Zeros "
    } else if(zeros == '1st' & y10 == 'N'){
      load(paste0(outpDir,'/',runType,' MZone ',numUnit,' Try ',1,".RData")) # swap out " Try " with " All Zeros "
    } else if(zeros == 'all' & y10 == 'Y'){
      load(paste0(outpDir,'/',runType,' MZone ',numUnit,' All Zeros 2005-2015.RData')) # swap out " Try " with " All Zeros "
    } else if(zeros == '1st' & y10 == 'Y'){
      load(paste0(outpDir,'/',runType,' MZone ',numUnit,' 1st Zeros 2005-2015.RData')) # swap out " Try " with " All Zeros "
    }
    
    rTime <- data.frame(time=time,group=theUnit)
    allTimes <- rbind(allTimes,rTime)
    bsimsList[[numUnit]]  <- bayes$sims.list
    
    rm(data,inits,parameters,bayes,CYear,freqs,time)
    
  }
  
  
  
  
  
  
  getSampleN <- function(metric,bsimsList,nYears){
    NmatList <- vector("list",nYears)
    wgtMeanNum <- vector("list",nYears)
    wgtMeanDen <- vector("list",nYears)
    wgtMean <- vector("list",nYears)
    allMeans <- NULL
    for(x in 1:nYears){  
      NmatList[[x]] <- cbind(bsimsList[[1]][metric][[1]][,1,x],bsimsList[[3]][metric][[1]][,1,x],bsimsList[[4]][metric][[1]][,1,x],bsimsList[[5]][metric][[1]][,1,x],bsimsList[[6]][metric][[1]][,1,x],bsimsList[[8]][metric][[1]][,1,x]) 
      wgtMeanNum[[x]] <- NmatList[[x]] %*% nLeksMat[,x]
      wgtMeanDen[[x]] <- sum(nLeksMat[,x])
      wgtMean[[x]] <- wgtMeanNum[[x]] / wgtMeanDen[[x]]
      thisMean <- data.frame(mean=mean(wgtMean[[x]]),sd=sd(wgtMean[[x]]),x2.5=quantile(wgtMean[[x]],0.025),x5=quantile(wgtMean[[x]],0.05),x25=quantile(wgtMean[[x]],0.25),x50=quantile(wgtMean[[x]],0.5),x75=quantile(wgtMean[[x]],0.75),x95=quantile(wgtMean[[x]],0.95),x97.5=quantile(wgtMean[[x]],0.975))
      allMeans <- rbind(allMeans,thisMean)
    }
    colnames(allMeans) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')
    rownames(allMeans) <- paste0(metric,'0[1,',seq(1,nYears,1),']') 
    allMeans <- list(allMeans,wgtMean)
  }
    
  getSampleNotN <- function(metric,bsimsList,nYears){
    
#     metric <- "mu.b"
#     bsimsList <- bsimsList
#     nYears <- 51
    
    
    NmatList <- vector("list",1)
    wgtMeanNum <- vector("list",1)
    wgtMeanDen <- vector("list",1)
    wgtMean <- vector("list",1)
    allMeans <- NULL

    NmatList[[1]] <- cbind(bsimsList[[1]][metric][[1]],bsimsList[[3]][metric][[1]],bsimsList[[4]][metric][[1]],bsimsList[[5]][metric][[1]],bsimsList[[6]][metric][[1]],bsimsList[[8]][metric][[1]])         
    wgtMeanNum[[1]] <- NmatList[[1]] %*% nLeksMat[,nYears]
    wgtMeanDen[[1]] <- sum(nLeksMat[,nYears])
    wgtMean[[1]] <- wgtMeanNum[[1]] / wgtMeanDen[[1]]
    allMeans <- data.frame(mean=mean(wgtMean[[1]]),sd=sd(wgtMean[[1]]),x2.5=quantile(wgtMean[[1]],0.025),x5=quantile(wgtMean[[1]],0.05),x25=quantile(wgtMean[[1]],0.25),x50=quantile(wgtMean[[1]],0.5),x75=quantile(wgtMean[[1]],0.75),x95=quantile(wgtMean[[1]],0.95),x97.5=quantile(wgtMean[[1]],0.975))

    names(allMeans) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')
    rownames(allMeans) <- paste0(metric,'0') 
    allMeans
  }
  
  N0 <- getSampleN("N",bsimsList,nYears)[[1]]
  mu.a0 <- getSampleNotN("mu.a",bsimsList,nYears)
  mu.b0 <- getSampleNotN("mu.b",bsimsList,nYears)
  sdnoise0 <- getSampleNotN("sdnoise",bsimsList,nYears)
  rho0 <- getSampleNotN("rho",bsimsList,nYears)
  sigma.a0 <- getSampleNotN("sigma.a",bsimsList,nYears)
  sigma.b0 <- getSampleNotN("sigma.b",bsimsList,nYears)  
  theNforB <- getSampleN("N",bsimsList,nYears)[[2]]
  
  # function for making posterior of Bs
  makeB <- function(theNforB,tb,ta,label){
    preBeta0 <- (theNforB[[tb]] / theNforB[[ta]])^(1/(tb-ta))    # note that it's 50 because it's 51 - 1.
    beta0 <- data.frame(mean=mean(preBeta0),sd=sd(preBeta0),x2.5=quantile(preBeta0,0.025),x5=quantile(preBeta0,0.05),x25=quantile(preBeta0,0.25),x50=quantile(preBeta0,0.5),x75=quantile(preBeta0,0.75),x95=quantile(preBeta0,0.95),x97.5=quantile(preBeta0,0.975))
    names(beta0) <- c('mean','sd','2.50%','5%','25%','50%','75%','95%','97.50%')  
    rownames(beta0) <- label   
    beta0
  }
  
  if(nYears == 51){
    beta0 <- makeB(theNforB,51,1,'beta0')
    B0.10.05.15 <- makeB(theNforB,51,41,'B0.10.05.15')
    B0.10.04.14 <- makeB(theNforB,50,40,'B0.10.04.14')
    B0.10.03.13 <- makeB(theNforB,49,39,'B0.10.03.13')
    B0.10.02.12 <- makeB(theNforB,48,38,'B0.10.02.12')
    B0.10.01.11 <- makeB(theNforB,47,37,'B0.10.01.11')
    B0.10.00.10 <- makeB(theNforB,46,36,'B0.10.00.10')
    B0.10.99.09 <- makeB(theNforB,45,35,'B0.10.99.09')
    B0.10.98.08 <- makeB(theNforB,44,34,'B0.10.98.08')
    B0.10.97.07 <- makeB(theNforB,43,33,'B0.10.97.07')
    B0.10.96.06 <- makeB(theNforB,42,32,'B0.10.96.06')
    B0.10.95.05 <- makeB(theNforB,41,31,'B0.10.95.05')
    B0.10.94.04 <- makeB(theNforB,40,30,'B0.10.94.04')
    B0.10.93.03 <- makeB(theNforB,39,29,'B0.10.93.03')
    B0.10.92.02 <- makeB(theNforB,38,28,'B0.10.92.02')
    B0.10.91.01 <- makeB(theNforB,37,27,'B0.10.91.01')
    B0.10.90.00 <- makeB(theNforB,36,26,'B0.10.90.00')
    B0.10.89.99 <- makeB(theNforB,35,25,'B0.10.89.99')
    B0.10.88.98 <- makeB(theNforB,34,24,'B0.10.88.98')
    B0.10.87.97 <- makeB(theNforB,33,23,'B0.10.87.97')
    B0.10.86.96 <- makeB(theNforB,32,22,'B0.10.86.96')
    B0.10.85.95 <- makeB(theNforB,31,21,'B0.10.85.95')
    B0.10.84.94 <- makeB(theNforB,30,20,'B0.10.84.94')
    B0.10.83.93 <- makeB(theNforB,29,19,'B0.10.83.93')
    B0.10.82.92 <- makeB(theNforB,28,18,'B0.10.82.92')
    B0.10.81.91 <- makeB(theNforB,27,17,'B0.10.81.91')
    B0.10.80.90 <- makeB(theNforB,26,16,'B0.10.80.90')
    B0.10.79.89 <- makeB(theNforB,25,15,'B0.10.79.89')
    B0.10.78.88 <- makeB(theNforB,24,14,'B0.10.78.88')
    B0.10.77.87 <- makeB(theNforB,23,13,'B0.10.77.87')
    B0.10.76.86 <- makeB(theNforB,22,12,'B0.10.76.86')
    B0.10.75.85 <- makeB(theNforB,21,11,'B0.10.75.85')
    B0.10.74.84 <- makeB(theNforB,20,10,'B0.10.74.84')
    B0.10.73.83 <- makeB(theNforB,19, 9,'B0.10.73.83')
    B0.10.72.82 <- makeB(theNforB,18, 8,'B0.10.72.82')
    B0.10.71.81 <- makeB(theNforB,17, 7,'B0.10.71.81')
    B0.10.70.80 <- makeB(theNforB,16, 6,'B0.10.70.80')
    B0.10.69.79 <- makeB(theNforB,15, 5,'B0.10.69.79')
    B0.10.68.78 <- makeB(theNforB,14, 4,'B0.10.68.78')
    B0.10.67.77 <- makeB(theNforB,13, 3,'B0.10.67.77')
    B0.10.66.76 <- makeB(theNforB,12, 2,'B0.10.66.76')
    B0.10.65.75 <- makeB(theNforB,11, 1,'B0.10.65.75')
    
    # put all the estimates together
    rangewide <- rbind(N0,mu.a0,mu.b0,beta0,sdnoise0,rho0,sigma.a0,sigma.b0,
                       B0.10.05.15,B0.10.04.14,B0.10.03.13,B0.10.02.12,B0.10.01.11,B0.10.00.10,B0.10.99.09,B0.10.98.08,B0.10.97.07,B0.10.96.06,
                       B0.10.95.05,B0.10.94.04,B0.10.93.03,B0.10.92.02,B0.10.91.01,B0.10.90.00,B0.10.89.99,B0.10.88.98,B0.10.87.97,B0.10.86.96,
                       B0.10.85.95,B0.10.84.94,B0.10.83.93,B0.10.82.92,B0.10.81.91,B0.10.80.90,B0.10.79.89,B0.10.78.88,B0.10.77.87,B0.10.76.86,
                       B0.10.75.85,B0.10.74.84,B0.10.73.83,B0.10.72.82,B0.10.71.81,B0.10.70.80,B0.10.69.79,B0.10.68.78,B0.10.67.77,B0.10.66.76,
                       B0.10.65.75)
  } else {
    # put all the estimates together
    rangewide <- rbind(N0,mu.a0,mu.b0,beta0,sdnoise0,rho0,sigma.a0,sigma.b0)    
  }
  

    
  # make graph
  if(nYears == 51){
    Year <- data.frame(Year=seq(1965,2015,1))
  } else {
    Year <- data.frame(Year=seq(2005,2015,1))
  }
  mu.a <- mu.a0$mean
  mu.b <- mu.b0$mean
  sdnoise <- sdnoise0$mean
  
  indx <- data.frame(rangewide)
  indx$vars <- rownames(indx)
  indx$x <- ifelse(substr(indx$vars,1,2) == 'N0',substr(indx$vars,4,4),-99)
  indx$y <- ifelse(substr(indx$vars,1,2) == 'N0',ifelse(nchar(indx$vars) == 7,substr(indx$vars,6,6),substr(indx$vars,6,7)),-99)
  
  # build the true B-matrix
  beta.mzone <- matrix(NA,nrow=1,ncol=nYears)
  list.beta.mzone <- indx[substr(indx$vars,1,3) == 'N0[',]
  for(i in 1:1){
    for(j in 1:nYears){
      beta.mzone[i,j] <- list.beta.mzone[list.beta.mzone$x == i & list.beta.mzone$y == j,]$mean
    }
  }
  
  coreText <- 'ugh'
  
  if(nYears == 51){
    theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])  
  } else {
    theN <- data.frame(Year=seq(2005,2015,1),N=beta.mzone[i,])      
  }
  obsMeans <- data.frame(MeanPMales=tapply(dat$Peak_Males, list(dat$Year), mean))
  obsMeans$Year <- rownames(obsMeans)
  
  plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
  plotYears$runType <- runType
  plotYears$numMZone <- 'Rangewide'
  plotYears$dataCut <- runType2
  if(nYears == 51){
    plotYears$Nhat <- exp(mu.a + mu.b*(as.numeric(plotYears$Year) - 1964 - 26) + 0.5*sdnoise*sdnoise)
    xLim <- c(1965,2015)
  } else {
    plotYears$Nhat <- exp(mu.a + mu.b*(as.numeric(plotYears$Year) - 2004 - 6) + 0.5*sdnoise*sdnoise)
    xLim <- c(2005,2015)
  }
  
  x  <- plotYears$Year
  y1 <- plotYears$N
  y2 <- plotYears$MeanPMales
  y3 <- plotYears$Nhat
  
  yMax <- 100
  
  png(filename=paste0(tracDir,'/Trend Plot - Rangewide - ',runType2,'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=xLim,ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Rangewide Temporal Trend of Observed Peak Males\nYears 1965-2015 - 75% ",runType2))
  par(new=TRUE)
  plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=xLim,ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  par(new=TRUE)
  plot(x,y3,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=xLim,ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  if(nYears == 51){
    axis(side=1,labels=TRUE,seq(1965,2015,5))
  } else {
    axis(side=1,labels=TRUE,seq(2005,2015,1))    
  }
  axis(side=2,labels=TRUE,seq(0,yMax,10))
  
  dev.off()
  
  write.csv(plotYears,paste0(tracDir,'/plotYears - Rangewide - ',runType2,'.csv'))
  write.csv(rangewide,paste0(tracDir,'/bayesSummary - ',runType2,'.csv'))
  rm(plotYears,x,y1,y2)
  
}