makeTracePlots <- function(nParms,parmList,tracDir,bayes$sims.array,bayes$summary){

  nParms <- nParms
  parmList <- parmList
  tracDir <- tracDir
  bsims <- bayes$sims.array
  bsums <- bayes$summary

  nPlots <- floor(nParms / 196) + 1
  nSamps <- nrow(bsims)
  
  for(k in 1:nPlots){
  
    png(filename=paste0(outpDir,tracDir,'/Trace Plots - ',196*(k - 1) + 1,'-',196*k,'.png'),width=40,height=40,units="in",res=300,pointsize=12)
  
    lay <- layout(matrix(seq(1,196*k,1),14,14,byrow=TRUE),rep(1/14,14),rep(1/14,14))
    layout.show(lay)
    for(i in 1:196*k){
      plot(c(1:nSamps),bsims[1:nSamps,1,i],type='l',lwd=0.5,col='red',main=dimnames(bsims)[[3]][i])
    }
    
    dev.off()

  }
}