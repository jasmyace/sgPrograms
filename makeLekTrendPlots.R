makeLekTrendPlots <- function(dat,tracDir,bayes){
  
#   dat <- dat
#   tracDir <- paste0(outpDir,'/',file,'/Lek Plots')
#   bayes <- bayes

  
  bsums <- bayes$summary  
  bsums.a <- bsums[substr(rownames(bsums),1,2) == 'a[',]
  bsums.b <- bsums[substr(rownames(bsums),1,2) == 'b[',]
  
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year    
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  nLeks <- max(dat$lekCls)
  
  if(min(dat$Year) == 1){
    dat$Year <- dat$Year + 1964
  }
  
  for(x in 1:nLeks){
    
    lek_ID <- gsub("/"," - ",dat[dat$lekCls == x,]$Lek_ID[1])
    cluster <- dat[dat$lekCls == x,]$Cluster[1]
    int.lek <- bsums.a[x,1]
    slp.lek <- bsums.b[x,1]
    
    if(length(unique(dat$Year)) == 11){
      minYear <- min(dat$Year)
      medYear <- quantile(dat$Year,0.5) - minYear
    } else {
      minYear <- 1965
      medYear <- 26
    }
    
    Nest <- data.frame(Nest=exp(int.lek + slp.lek*(seq(1965,2015,1) - (minYear - 1) - medYear)))
    Nest$Year <- as.numeric(rownames(Nest)) + 1964
    Nobs <- dat[dat$lekCls == x,c('Year','yearCls','Peak_Males')]
    
    trendLek <- merge(Nobs,Nest,by=c('Year'),all.x=TRUE)
  
    X <- trendLek$Year
    y1 <- trendLek$Peak_Males
    y2 <- trendLek$Nest
    
    yMax <- 100
    
    png(filename=paste0(tracDir,'/Lek Plot - ',lek_ID,'.png'),width=8,height=6,units="in",res=300,pointsize=12)
    
    plot(X,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years ",minYear,"-",maxYear,"\nLek ID: ",lek_ID,' - Cluster: ',cluster))
    par(new=TRUE)
    plot(X,y2,type='l',col=colVec[1],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
    axis(side=1,labels=TRUE,seq(1965,2015,5))
    axis(side=2,labels=TRUE,seq(0,yMax,10))
    
    dev.off()
  }
}