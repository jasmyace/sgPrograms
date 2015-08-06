makeHistogramPlots <- function(dat,theUnit,tracDir){
  
#   dat <- dat
#   theUnit <- 'WY'
#   tracDir <- paste0(outpDir,'/',file,'/Zeros Plots')
  
  
  
  dat <- dat[order(dat$Year),]
  
  nYears <- length(unique(dat$Year))
  years <- unique(dat$Year)
  
  png(filename=paste0(tracDir,'/Zeros Plots - 1st - ',file,'.png'),width=5,height=80,units="in",res=300,pointsize=12)
  
  lay <- layout(matrix(seq(1,51,1),51,1,byrow=TRUE),1,rep(1/51,51))
  layout.show(lay)
  
  for(z in 1:nYears){
    if(z <= nYears){
      thisYear <- years[z]
      yearDat <- dat[dat$Year == thisYear & dat$Peak_Males <= 100,]
      hist(yearDat$Peak_Males,breaks=100,main=paste0('Histogram',theUnit,' - First Zeros - Year ',thisYear))
    } else {
      plot.new()
    }
  }
  dev.off() 
}