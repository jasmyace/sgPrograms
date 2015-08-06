makePosteriorPlots <- function(nParms,parmList,tracDir,file,bayes){
  
  #   nParms <- nParms
  #   parmList <- parmList
  #   tracDir <- paste0(outpDir,'/',file,'/Trace Plots')
  
  bsims <- bayes$sims.array
  bsums <- bayes$summary
  
  nPlots <- floor(nParms / 196) + 1
  nSamps <- nrow(bsims)
  
  for(k in 1:nPlots){
    
    png(filename=paste0(tracDir,'/Posterior Plots - ',file,' - ',196*(k - 1) + 1,'-',196*k,'.png'),width=40,height=40,units="in",res=300,pointsize=12)
    
    lay <- layout(matrix(seq(1,196*k,1),14,14,byrow=TRUE),rep(1/14,14),rep(1/14,14))
    layout.show(lay)
    for(l in 1:196){
      if(196*(k - 1) + l <= nParms){
        hist(bsims[1:nSamps,1,196*(k - 1) + l],col='blue',border='blue',main=dimnames(bsims)[[3]][196*(k - 1) + l],breaks=100)
        abline(v = mean(bsims[1:4000,1,196*(k - 1) + l]), lwd=2, col='red')
        abline(v = median(bsims[1:4000,1,196*(k - 1) + l]), lwd=2, col='green')
      } else {
        plot.new()
      }
    }
    dev.off()
  }
}