makeHistogramPlots <- function(dat,theUnit,tracDir,file){
  
#   dat <- datAllZerosCore[[6]] simPoisRandCoefCorrData31Leks dat
#   theUnit <- paste0(" - ",string," ",theUnit)
#   tracDir <- paste0(outpDir,'/',file,'/Zeros Plots')
#   file <- file
  
  
  dat <- dat[order(dat$Year),]
  
  if(min(dat$Year) == 1){
    dat$Year <- dat$Year + 1964
  }
  
  nYears <- length(unique(dat$Year))
  years <- unique(dat$Year)
  
  png(filename=paste0(tracDir,'/Zeros Plots - All Zeros ',file,'.png'),width=5,height=80,units="in",res=300,pointsize=12)
  
  lay <- layout(matrix(seq(1,nYears + 1,1),nYears + 1,1,byrow=TRUE),1,c(rep(1/(nYears + 4),nYears),4/(nYears + 4)))
  layout.show(lay)
  
  for(z in 1:nYears){
    if(z <= nYears){
      thisYear <- years[z]
      yearDat <- dat[dat$Year == thisYear & dat$Peak_Males <= 100,]
      lambda <- mean(yearDat$Peak_Males)
      xpois <- trunc(rpois(1000,lambda))
      hist(yearDat$Peak_Males,breaks=100,main=paste0('Histogram',theUnit,' - ',file,'\nYear ',thisYear))
      lines(seq(0,100), length(yearDat$Peak_Males)*dpois(seq(0,100),lambda),col='red')
    } else {
      plot.new()
    }
  }
  allDat <- dat[dat$Peak_Males <= 100,]
  lambda <- mean(allDat$Peak_Males)
  xpois <- trunc(rpois(1000,lambda))
  hist(allDat$Peak_Males,breaks=100,main=paste0('Histogram',theUnit,' - ',file,'\n All Years'))
  lines(seq(0,100), length(allDat$Peak_Males)*dpois(seq(0,100),lambda),col='red')
  
  # zero-truncated poisson
#   n <- length(allDat$Peak_Males) # desired size of sample
#   T <- lambda                    # pre-truncation mean of Poisson
#   U <- runif(n)                  # the uniform sample
#   t <- -log(1 - U*(1 - exp(-T))) # the "first" event-times
#   T1 <- (T - t)                  # the set of (T-t)
#   X <- rpois(n,T1)+1
#   lines(seq(0,100), length(allDat$Peak_Males)*dpois(seq(0,100),mean(X)),col='blue')
  
  dev.off() 
}


# foo <- rpois(100, lambda=1)
# hist(foo, prob=TRUE)
# plot(table(foo)/length(foo)); points(0:8, dpois(0:8,lambda=1))
# curve(dpois(x, lambda=mean(foo)), add=TRUE)
# 
# set.seed(123)
# xpois <- trunc(rpois(100, 4))
# hist(xpois)
# lines(seq(0,10), 100*dpois(seq(0,10), 4))





