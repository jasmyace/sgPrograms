makeRandomLekPlots <- function(dat,tracDir,file,bayes,bsums90){
  
     
  nLeksI <- length(rownames(bsums90)[substr(rownames(bsums90),1,2) == 'a['])
  nLeksS <- length(rownames(bsums90)[substr(rownames(bsums90),1,2) == 'b['])
  if(nLeksI != nLeksS){stop('[a does not equal [b.  Check it out.')}
    
  
  # make a few different graphs, so we have options.
  for(j in 1:10){
  
    # make a random sample
    bsums90$lekCls <- c(1:nLeksI,1:nLeksS,rep(NA,nrow(bsums90) - nLeksI - nLeksS))
    lekSamp <- data.frame(lekCls=sample(unique(bsums90$lekCls),16,replace=FALSE))   
    bsums90Samp <- bsums90[bsums90$lekCls %in% lekSamp$lekCls,]   # the model results we need
      
    dat$lekCls <- as.numeric(as.factor(droplevels(dat$Lek_ID)))
    df <- dat[,c('lekCls','Lek_ID','Year','Cluster','Long','Lat','Peak_Males')]
  
    # the merge puts the lekCls in order.  use this to make a key.
    theSamp <- merge(lekSamp,dat,by=c('lekCls'),all.x=TRUE)   # the raw data we need
    printTheSamp <- unique(theSamp[,c('lekCls','Mgmt_zone','Lek_ID','Cluster','Population','mZoneNum')])
    printTheSamp$key <- seq(1,16)
    printTheSamp$help <- 'Position #1 is the top left; position #16 is the bottom right;  sequencing first follows row 1, row 2, etc.'
    
    write.csv(printTheSamp,paste0(tracDir,'/Random Lek Plot Key ',j,' - ',file,'.csv'))
    png(filename=paste0(tracDir,'/Random Lek Plot ',j,' - ',file,'.png'),width=8,height=6,units="in",res=300,pointsize=12)
    
    lay <- layout(matrix(seq(1,25),5,5,byrow=TRUE),c(0.08,0.23,0.23,0.23,0.23,0.08),c(0.23,0.23,0.23,0.23,0.08))
    layout.show(lay)
    
    leks <- unique(bsums90Samp$lekCls)
    for(i in 1:20){
      
      # nothingness in left space
      if(i %% 5 == 1){
        par(mar = c(0,0,0,0)); plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE); u <- par("usr")
      } else {
       
        theLek <- leks[1]
        thisSamp <- theSamp[theSamp$lekCls == theLek,]
      
        backbone <- data.frame(Year=seq(1965,2015))
        forPlot <- merge(backbone,thisSamp,by=c('Year'),all.x=TRUE)
      
        mu.a <- bsums90Samp[bsums90Samp$lekCls == theLek,][1,1]
        mu.b <- bsums90Samp[bsums90Samp$lekCls == theLek,][2,1]
        yMin <- 0   #min(forPlot$Peak_Males)
        yMax <- 100 #max(forPlot$Peak_Males)
        xMin <- min(forPlot[!is.na(forPlot$Peak_Males),]$Year)
        xMax <- max(forPlot[!is.na(forPlot$Peak_Males),]$Year)
      
        forPlot$trend <- exp(mu.a + mu.b*(forPlot$Year - 1964 - 26))
        forPlot$trend <- ifelse(forPlot$Year < xMin | forPlot$Year > xMax,NA,forPlot$trend)
  
        X <- forPlot$Year
        yObs <- forPlot$Peak_Males
        yEst <- forPlot$trend
      
        par(mar = c(0.5,0.5,0.5,0.5))
        if(i == 17){
          plot(X,yObs,type='p',pch=19,col='gray50',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='')
          par(new=TRUE)
          plot(X,yEst,type='l',col='red',axes=FALSE,xpd=NA,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='Year',ylab='Peak Males / Lek',lwd=2)
          axis(side=2,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),las=2,cex.axis=0.8)
          axis(side=1,at=c(1965,1975,1985,1995,2005,2015),labels=c(1965,1975,1985,1995,2005,2015),cex.axis=0.8)
        } else if (i %% 5 == 2 & i != 17){
          plot(X,yObs,type='p',pch=19,col='gray50',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='')
          par(new=TRUE)
          plot(X,yEst,type='l',col='red',axes=FALSE,xpd=NA,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='Peak Males / Lek',lwd=2)
          axis(side=2,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),las=2,cex.axis=0.8)
        } else if(i > 17){
          plot(X,yObs,type='p',pch=19,col='gray50',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='')
          par(new=TRUE)
          plot(X,yEst,type='l',col='red',axes=FALSE,xpd=NA,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='Year',ylab='',lwd=2)      
          axis(side=1,at=c(1965,1975,1985,1995,2005,2015),labels=c(1965,1975,1985,1995,2005,2015),cex.axis=0.8)
        } else {
          plot(X,yObs,type='p',pch=19,cex=0.8,col='gray50',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='')
          par(new=TRUE)
          plot(X,yEst,type='l',col='red',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)              
        }
        leks <- leks[-1]
      } 
    }
    dev.off()
  }
}