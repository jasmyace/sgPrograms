
      
    readInAnalyticFiles(analDir)
    
    colVec <- brewer.pal(12,"Paired")
    colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]
    
    
    
    
  
    
    compileResults("Model D MZone 6 Test All-An-4k-76k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")
      
    
    
    
    
    
    
    
    
    
    
    compileResults <- function(file,dat,string,theUnit,runType){
      
#       file <- "Model D MZone 6 Test All-An-4k-76k-1"
#       dat <- datAllZerosCore[[6]]
#       string <- "Management Zone"
#       theUnit <- "mZone 6"
#       runType <- "Core"
      
      load(paste0(outpDir,'/',file,".RData"))  
      
      ifelse(!dir.exists(file.path(outpDir,file)), dir.create(file.path(outpDir,file)), FALSE)      # make new folder
      file.copy(paste0(outpDir,"/",file,".RData"),paste0(outpDir,"/",file))                         # copy bayes output to new folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), FALSE)        # make new trace plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), FALSE)# make new posterior plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), FALSE)        # make new mzone trend plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), FALSE)        # make new mzone trend plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Lek Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Lek Plots'))), FALSE)            # make new mzone lek plots folder
      
      # make trace, posterior, histogram plots
      nParms <- dim(bayes$sims.array)[3]
      parmList <- dimnames(bayes$sims.array)[[3]]
      makeTracePlots(nParms,parmList,paste0(outpDir,'/',file,'/Trace Plots'),file,bayes)    
      makePosteriorPlots(nParms,parmList,paste0(outpDir,'/',file,'/Posterior Plots'),file,bayes)
      makeHistogramPlots(dat,paste0(" - ",string," ",theUnit),paste0(outpDir,'/',file,'/Zeros Plots'),file)
      makeLekTrendPlots(dat,paste0(outpDir,'/',file,'/Lek Plots'),bayes)
      
      bsums90 <- make90pCredInt(bayes)
      
      
      # make bayes summary file of estimates
      write.csv(bsums90,paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
      
      # make trend plots
      makeTrendPlots(dat,runType,paste0(" - ",string," ",theUnit),1,paste0(outpDir,'/',file,'/Trend Plots'),file,bayes)
      
      rm(nParms,parmList)
    }  

