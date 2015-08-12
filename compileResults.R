


# loadThis <- 'Model S MZone 9 All MZones Core, Try 2'
# loadThis <- 'Model S MZone 9 All MZones Non-Core, Try 2'
# 
# loadThis <- 'Model S MZone 9 MZones 1,3,5 Core, Try 1'
# loadThis <- 'Model T MZone 9 MZones 1,3,5 Core, Adj M,Try 1'
# 
# loadThis <- 'Model U MZone 9 All MZones Non-Core, Adj M,Try 1'
# loadThis <- 'Model U MZone 9 All MZones Core, Adj M,Try 1'

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("RData", theFiles) == 1]
colVec <- brewer.pal(12,"Paired")

colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]

states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZones2 <- c("Zone 1","Zone 3","Zone 4","Zone 5","Zone 6","Zone 2 & 7")
textyZones <- c("Great Plains","Southern Great Basin","Snake River Plain","Northern Great Basin","Columbian Basin","Wyoming Basin & Colorado Plateau")

nStates <- length(states)
nMZones <- length(mZones)
textyStates <- c("California","Colorado","Idaho","Montana","North Dakota","Nevada","Oregon","South Dakota","Utah","Washington","Wyoming")

readInAnalyticFiles(analDir)





for(h in 1:2){
  if(h == 1){
    units <- mZones
    nUnits <- nMZones
    start <- 1
    string <- 'Management Zone'
  } else {
    units <- states
    nUnits <- nStates
    start <- 3
    string <- 'State'
  }

  # loop over mZone files in theFiles
  for(i in 1:nUnits){
  
    theUnit <- units[i]
    numUnit <- as.numeric(substr(theUnit,7,7))                                               # applicable only to mZone?
    unitFile <- theFiles[grepl(units[i],theFiles)]
  
    for(j in start:3){
      
      # get data for analysis
      if(h == 1){
        load(paste0(outpDir,'/',unitFile[j]))                                                # get bayesian stuff -- mzone
        #load(paste0(outpDir,"/Model V MZone 6 ZInf All Leks Try 2",".RData"))
      } else {
        load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
          dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type      
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
          dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type      
      }
      
      dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
      
      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else {
        file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      }

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
      makeHistogramPlots(dat,paste0(" - ",string," ",theUnit),paste0(outpDir,'/',file,'/Zeros Plots'))
      makeLekTrendPlots(dat,paste0(outpDir,'/',file,'/Lek Plots'),bayes)
      
      bsums90 <- make90pCredInt(bayes)
      
      
      # make bayes summary file of estimates
      write.csv(bsums90,paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
      
      # make trend plots
      makeTrendPlots(dat,runType,paste0(" - ",string," ",theUnit),1,paste0(outpDir,'/',file,'/Trend Plots'),file,bayes)
  
      rm(nParms,parmList)
      
      # make sure h == 1 and j == 3!  don't run over states -- won't work!
      if(h == 1 & j == 3){
        ifelse(!dir.exists(file.path(outpDir,'Rangewide')), dir.create(file.path(outpDir,'Rangewide')), FALSE)      # make new folder
        ifelse(!dir.exists(file.path(paste0(outpDir,'/','Rangewide','/Trend Plots'))), dir.create(file.path(paste0(outpDir,'/','Rangewide','/Trend Plots'))), FALSE)        # make new mzone trend plots folder
        
        # make rangewide estimates
        makeRangeWide("Model D",units,dat1stZerosCore,"Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model E",units,dat1stZerosNoco,"Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model F",units,dat1stZerosLeks,"All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
        # make core + non-core + all lek plot
      } 
    }
  }
}




# make a big table of results

theSumms <- NULL
# make table
for(h in 1:3){
  if(h == 1){
    units <- mZones
    nUnits <- nMZones
    start <- 1
    string <- 'Management Zone'
  } else if(h == 2){
    units <- states
    nUnits <- nStates
    start <- 3
    string <- 'State'
  } else {
    units <- 'Rangewide'
    nUnits <- 1
    start <- 1
    string <- units
  }
  
  # loop over mZone files in theFiles
  for(i in 1:nUnits){
    
    theUnit <- units[i]
    numUnit <- as.numeric(substr(theUnit,7,7))                                               # applicable only to mZone?
    unitFile <- theFiles[grepl(units[i],theFiles)]
    
    for(j in start:3){
      
      # get data for analysis
      if(h == 1){
        load(paste0(outpDir,'/',unitFile[j]))                                                # get bayesian stuff -- mzone
        #load(paste0(outpDir,"/Model V MZone 6 ZInf All Leks Try 2",".RData"))
      } else if(h == 2){
        load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
#       } else {
#         load(paste0(outpDir,'/',))
#       }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
          dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
        allNLeks <- length(unique(dat1stZerosCore[[9]]$Lek_ID))
        file2 <- 'D'
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type  
        allNLeks <- length(unique(dat1stZerosNoco[[9]]$Lek_ID))
        file2 <- 'E'
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
          dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type   
        allNLeks <- length(unique(dat1stZerosLeks[[9]]$Lek_ID))
        file2 <- 'F'
      }
      
      if(h == 1 | h == 2){
        dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
        
      }

      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else if(h == 2){
        file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      } else {
        file <- 'Rangewide/Trend Plots'                                                             # get name of new folder 
      }

      if(h == 1 | h == 2){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
        mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.a'   ,mean=thisOne[thisOne$X == 'mu.a'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a'   ,]$X95.)
        mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.b'   ,mean=thisOne[thisOne$X == 'mu.b'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b'   ,]$X95.)
        beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='beta'   ,mean=thisOne[thisOne$X == 'beta[1]',]$mean,X5.=thisOne[thisOne$X == 'beta[1]',]$X5.,X95.=thisOne[thisOne$X == 'beta[1]',]$X95.)
        sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sdnoise',mean=thisOne[thisOne$X == 'sdnoise',]$mean,X5.=thisOne[thisOne$X == 'sdnoise',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise',]$X95.)
        rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='rho'    ,mean=thisOne[thisOne$X == 'rho'    ,]$mean,X5.=thisOne[thisOne$X == 'rho'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho'    ,]$X95.)
        sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.a',mean=thisOne[thisOne$X == 'sigma.a',]$mean,X5.=thisOne[thisOne$X == 'sigma.a',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a',]$X95.)
        sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.b',mean=thisOne[thisOne$X == 'sigma.b',]$mean,X5.=thisOne[thisOne$X == 'sigma.b',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b',]$X95.)
      } else if (h == 3){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - Model ',file2,'.csv'))
        mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='mu.a0'   ,mean=thisOne[thisOne$X == 'mu.a0'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a0'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a0'   ,]$X95.)
        mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='mu.b0'   ,mean=thisOne[thisOne$X == 'mu.b0'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b0'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b0'   ,]$X95.)
        beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='beta0'   ,mean=thisOne[thisOne$X == 'beta0'   ,]$mean,X5.=thisOne[thisOne$X == 'beta0'   ,]$X5.,X95.=thisOne[thisOne$X == 'beta0'   ,]$X95.)
        sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sdnoise0',mean=thisOne[thisOne$X == 'sdnoise0',]$mean,X5.=thisOne[thisOne$X == 'sdnoise0',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise0',]$X95.)
        rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='rho0'    ,mean=thisOne[thisOne$X == 'rho0'    ,]$mean,X5.=thisOne[thisOne$X == 'rho0'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho0'    ,]$X95.)
        sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sigma.a0',mean=thisOne[thisOne$X == 'sigma.a0',]$mean,X5.=thisOne[thisOne$X == 'sigma.a0',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a0',]$X95.)
        sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sigma.b0',mean=thisOne[thisOne$X == 'sigma.b0',]$mean,X5.=thisOne[thisOne$X == 'sigma.b0',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b0',]$X95.)
        
      }

      summs <- rbind(mu.a,mu.b,beta,sdnoise,rho,sigma.a,sigma.b)
      
      theSumms <- rbind(theSumms,summs)
       
    }
  }
}

theSumms <- theSumms[order(theSumms$focus,theSumms$stat,theSumms$group),]
# theSumms$mean <- ifelse(theSumms$stat %in% c('mu.a','mu.b','mu.a0','mu.b0'),print(formatC(signif(exp(as.numeric(theSumms$mean)),digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(as.numeric(theSumms$mean),digits=3), digits=3,format="fg", flag="#")))
# theSumms$X5.  <- ifelse(theSumms$stat %in% c('mu.a','mu.b','mu.a0','mu.b0'),print(formatC(signif(exp(as.numeric(theSumms$X5.)) ,digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(as.numeric(theSumms$X5. ),digits=3), digits=3,format="fg", flag="#")))
# theSumms$X95. <- ifelse(theSumms$stat %in% c('mu.a','mu.b','mu.a0','mu.b0'),print(formatC(signif(exp(as.numeric(theSumms$X95.)),digits=3), digits=3,format="fg", flag="#")),print(formatC(signif(as.numeric(theSumms$X95.),digits=3), digits=3,format="fg", flag="#")))


coreSumms <- theSumms[theSumms$cut == 'Core',]
names(coreSumms)[names(coreSumms) == 'mean'] <- 'Coremean'
names(coreSumms)[names(coreSumms) == 'X5.'] <- 'CoreX5.'
names(coreSumms)[names(coreSumms) == 'X95.'] <- 'CoreX95.'
names(coreSumms)[names(coreSumms) == 'Nleks'] <- 'CoreNleks'
coreSumms$cut <- NULL

nocoSumms <- theSumms[theSumms$cut == 'Non-Core',]
names(nocoSumms)[names(nocoSumms) == 'mean'] <- 'Nocomean'
names(nocoSumms)[names(nocoSumms) == 'X5.'] <- 'NocoX5.'
names(nocoSumms)[names(nocoSumms) == 'X95.'] <- 'NocoX95.'
names(nocoSumms)[names(nocoSumms) == 'Nleks'] <- 'NocoNleks'
nocoSumms$cut <- NULL
nocoSumms$focus <- NULL

leksSumms <- theSumms[theSumms$cut == 'All Leks',]
names(leksSumms)[names(leksSumms) == 'mean'] <- 'Leksmean'
names(leksSumms)[names(leksSumms) == 'X5.'] <- 'LeksX5.'
names(leksSumms)[names(leksSumms) == 'X95.'] <- 'LeksX95.'
names(leksSumms)[names(leksSumms) == 'Nleks'] <- 'LeksNleks'

preTable1 <- merge(nocoSumms,coreSumms,by=c('group','stat'),all.x=TRUE)
preTable2 <- merge(preTable1,leksSumms,by=c('group','stat'),all.y=TRUE)
preTable2$focus.y <- NULL
preTable2$cut <- NULL
preTable2$focus.x <- ifelse(is.na(preTable2$focus.x),'State',preTable2$focus.x)   # should probably fix above.  eh.
preTable2$focus.x <- ifelse(preTable2$focus.x == 3,'Rangewide',preTable2$focus.x)
  
preTable3 <- preTable2[order(preTable2$focus.x,preTable2$stat,preTable2$group),]
table1 <- preTable3[,c('focus.x','group','stat','LeksNleks','Leksmean','LeksX5.','LeksX95.','CoreNleks','Coremean','CoreX5.','CoreX95.','NocoNleks','Nocomean','NocoX5.','NocoX95.')]
write.csv(table1,paste0(manuDir,'/table1.csv'))











# make table 2 of b10s.
Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))

names(Rcore.CSV) <- c('X','core.mean','core.sd','coreX2.5','coreX5','coreX25','coreX50','coreX75','coreX95','coreX97.5')
names(Rnoco.CSV) <- c('X','noco.mean','noco.sd','nocoX2.5','nocoX5','nocoX25','nocoX50','nocoX75','nocoX95','nocoX97.5')
names(Rleks.CSV) <- c('X','leks.mean','leks.sd','leksX2.5','leksX5','leksX25','leksX50','leksX75','leksX95','leksX97.5')

Rcore.CSV2 <- Rcore.CSV[,c('X','core.mean','coreX5','coreX95')]
Rnoco.CSV2 <- Rnoco.CSV[,c('X','noco.mean','nocoX5','nocoX95')]
Rleks.CSV2 <- Rleks.CSV[,c('X','leks.mean','leksX5','leksX95')]

preTable2   <- merge(Rleks.CSV2,Rcore.CSV2,by=c('X'))
preTable2.B <- merge(preTable2,Rnoco.CSV2,by=c('X'))
preTable2.C <- preTable2.B[1:41,]
table2 <- rbind(preTable2.C[7:41,],preTable2.C[1:6,] )
rownames(table2) <- NULL
write.csv(table2,paste0(manuDir,'/table2.csv'))

























# make plot, for each mzone, of core, non-core, all leks together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  

#   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
#   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
#   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
#   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
#   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
#   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
#   
#   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
#   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
#   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
#   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
#   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
#   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
  coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
  RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 

  coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
  RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
  coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
  RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.


  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
  Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                                                 Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                                                 Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                                                Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
                                                Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                                                Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC','Y.core','Y.core.X5','Y.core.X95',
                                       'Y.noco','Y.noco.X5','Y.noco.X95',
                                       'Y.leks','Y.leks.X5','Y.leks.X95',
                                       'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
                                       'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
                                       'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
  
  
    x <- plotYears$Year
  yA1 <- plotYears$Y.core
  yA2 <- plotYears$Y.core.X5
  yA3 <- plotYears$Y.core.X95
  yB1 <- plotYears$Y.noco
  yB2 <- plotYears$Y.noco.X5
  yB3 <- plotYears$Y.noco.X95
  yC1 <- plotYears$Y.leks
  yC2 <- plotYears$Y.leks.X5
  yC3 <- plotYears$Y.leks.X95

  yRA1 <- plotYears$Y.Rcore
  yRA2 <- plotYears$Y.Rcore.X5
  yRA3 <- plotYears$Y.Rcore.X95
  yRB1 <- plotYears$Y.Rnoco
  yRB2 <- plotYears$Y.Rnoco.X5
  yRB3 <- plotYears$Y.Rnoco.X95
  yRC1 <- plotYears$Y.Rleks
  yRC2 <- plotYears$Y.Rleks.X5
  yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  # management zone specific
  if(i == 6){
    
    plot(x,yA1,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=2,lwd=2,main=paste0("Management Zone 2 & 7: ",textyZones[i]))    
  } else {
    plot(x,yA1,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=2,lwd=2,main=paste0("Management Zone ",substr(units[i],7,7),": ",textyZones[i]))
  }
    #   par(new=TRUE)
#   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE)
#   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
   
  par(new=TRUE)
  plot(x,yB1,type='l',col='black'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=3,lwd=2)
#   par(new=TRUE)
#   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  par(new=TRUE)
  plot(x,yC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)
#   par(new=TRUE)
#   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)

  # rangewide specific
#   par(new=TRUE)
#   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE)
#   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   
#   par(new=TRUE)
#   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   
#   par(new=TRUE)
#   plot(x,yRC1,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  legend("topright",c("Core","Periphery","Combined"),lty=c(2,3,1),lwd=c(2,2,2),col=c('black','black','black'),bty="n")
  axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
  axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
  mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
  mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
  dev.off()
  
}




# make plot, for each state, of all leks together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:11){
  units <- paste0("MZone ",states)
  nUnits <- nStates
#   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
#   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
#   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
#   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
#   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
#   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
#   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
#   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
#   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
#   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
#   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
#   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
#   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
#   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
#   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
#   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
#   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
#   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
#   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
#   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
#   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
#   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
#   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
#   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
#   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
#   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
#   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
#   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
#   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
#   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
#   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
#   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
#                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
#                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
#                           Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
#                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                          Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
#                         'Y.core','Y.core.X5','Y.core.X95',
#                         'Y.noco','Y.noco.X5','Y.noco.X95',
                        'Y.leks','Y.leks.X5','Y.leks.X95',
#                         'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
#                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
                        'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
  
  
  x <- plotYears$Year
#   yA1 <- plotYears$Y.core
#   yA2 <- plotYears$Y.core.X5
#   yA3 <- plotYears$Y.core.X95
#   yB1 <- plotYears$Y.noco
#   yB2 <- plotYears$Y.noco.X5
#   yB3 <- plotYears$Y.noco.X95
  yC1 <- plotYears$Y.leks
  yC2 <- plotYears$Y.leks.X5
  yC3 <- plotYears$Y.leks.X95
  
#   yRA1 <- plotYears$Y.Rcore
#   yRA2 <- plotYears$Y.Rcore.X5
#   yRA3 <- plotYears$Y.Rcore.X95
#   yRB1 <- plotYears$Y.Rnoco
#   yRB2 <- plotYears$Y.Rnoco.X5
#   yRB3 <- plotYears$Y.Rnoco.X95
  yRC1 <- plotYears$Y.Rleks
  yRC2 <- plotYears$Y.Rleks.X5
  yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  CairoPNG(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  # state specific
#   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
#   par(new=TRUE)
#   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE)
#   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   
#   par(new=TRUE)
#   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  plot(x,yC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0(textyStates[i]))
#   par(new=TRUE)
#   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
#   par(new=TRUE)
#   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE)
#   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   
#   par(new=TRUE)
#   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
#   par(new=TRUE)
#   plot(x,yRC1,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   par(new=TRUE)
#   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  legend("topright",c("Combined"),lty=c(1),lwd=c(2),col=c('black'),bty="n")
  axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
  axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
  mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
  mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
  dev.off()
  
}












# make plot, for rangewide alone, that includes core, non-core, all leks
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:1){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
  #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  #leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
  RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
  #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  #leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
  RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
  #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  #leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
  RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
#   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
#   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
                          #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                          #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          #Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                                                     Y.Rcore=Y.Rcore,#Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
                                                     Y.Rnoco=Y.Rnoco,#Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                          Y.Rleks=Y.Rleks)#Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
                        #                         'Y.core','Y.core.X5','Y.core.X95',
                        #                         'Y.noco','Y.noco.X5','Y.noco.X95',
                        #'Y.leks','Y.leks.X5','Y.leks.X95',
                                                 'Y.Rcore',#'Y.Rcore.X5','Y.Rcore.X95',
                                                 'Y.Rnoco',#'Y.Rnoco.X5','Y.Rnoco.X95',
                        'Y.Rleks')#'Y.Rleks.X5','Y.Rleks.X95')
  
  
  x <- plotYears$Year
  #   yA1 <- plotYears$Y.core
  #   yA2 <- plotYears$Y.core.X5
  #   yA3 <- plotYears$Y.core.X95
  #   yB1 <- plotYears$Y.noco
  #   yB2 <- plotYears$Y.noco.X5
  #   yB3 <- plotYears$Y.noco.X95
#   yC1 <- plotYears$Y.leks
#   yC2 <- plotYears$Y.leks.X5
#   yC3 <- plotYears$Y.leks.X95
  
    yRA1 <- plotYears$Y.Rcore
  #   yRA2 <- plotYears$Y.Rcore.X5
  #   yRA3 <- plotYears$Y.Rcore.X95
    yRB1 <- plotYears$Y.Rnoco
  #   yRB2 <- plotYears$Y.Rnoco.X5
  #   yRB3 <- plotYears$Y.Rnoco.X95
  yRC1 <- plotYears$Y.Rleks
#   yRC2 <- plotYears$Y.Rleks.X5
#   yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - Rangewide.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
  # state specific
  #   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
  #   par(new=TRUE)
  #   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
#   plot(x,yC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0(textyStates[i]))
  #   par(new=TRUE)
  #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
  #   par(new=TRUE)
    plot(x,yRA1,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=2,lwd=2,main="Rangewide")
  #   par(new=TRUE)
  #   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
    par(new=TRUE)
    plot(x,yRB1,type='l',col='black'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=3,lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
    par(new=TRUE)
    plot(x,yRC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
legend("topright",c("Core","Periphery","Combined"),lty=c(2,3,1),lwd=c(2,2,2),col=c('black','black','black'),bty="n")
axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis

  dev.off()
  
}










# make the new B10 plots -- histograms
# Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
# Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
# Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
# 
# Rcore.CSV      <- Rcore.CSV[59:99,]
# Rnoco.CSV      <- Rnoco.CSV[59:99,]
# Rleks.CSV      <- Rleks.CSV[59:99,]

png(filename=paste0(manuDir,'/B10 Rangewide Histograms - C, NC, AL.png'),width=12.5,height=9,units="in",res=600,pointsize=12)
par(mfrow=c(2,2))
min <- min(table2$core.mean)
max <- max(table2$core.mean)
hist(table2$core.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="Core",xlab="Estimates of Annual % Change",col="gray",las=1)
abline(v=table2[41,]$core.mean,lwd=3)
axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
axis(2,at=c(0:10),labels=c(0:10),las=1)
text(table2[41,]$core.mean-0.00005,9.9,"2005-2015",cex=0.8,pos=2)
text(table2[41,]$core.mean-0.00005,9.5,"Annualized % Change",cex=0.8,pos=2)

min <- min(table2$noco.mean)
max <- max(table2$noco.mean)
hist(table2$noco.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="Periphery",xlab="Estimates of Annual % Change",col="gray",las=1)
abline(v=table2[41,]$noco.mean,lwd=3)
axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
axis(2,at=c(0:10),labels=c(0:10),las=1)
text(table2[41,]$noco.mean-0.0005,9.9,"2005-2015",cex=0.8,pos=2)
text(table2[41,]$noco.mean-0.0005,9.5,"Annualized % Change",cex=0.8,pos=2)

min <- min(table2$leks.mean)
max <- max(table2$leks.mean)
hist(table2$leks.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="All Leks Combined",xlab="Estimates of Annual % Change",col="gray",las=1)
abline(v=table2[41,]$leks.mean,lwd=3)
# ylab <- ifelse(seq(min,max,(max-min)/9) >= 1,seq(min,max,(max-min)/9)-1,seq(min,max,(max-min)/9)-1)
axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
axis(2,at=c(0:10),labels=c(0:10),las=1)
text(table2[41,]$leks.mean-0.0001,9.9,"2005-2015",cex=0.8,pos=2)
text(table2[41,]$leks.mean-0.0001,9.5,"Annualized % Change",cex=0.8,pos=2)

plot.new()
dev.off()

par(mfrow=c(1,1))



print(formatC(signif(round(seq(min,max,(max-min)/9),3))),digits=3)






# make the B10 plots -- old with all B10s on one graph
Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))

Year <- seq(1965,2015,1)
YearC <- Year - 1964 - 26

RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean

RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean

RYcore <- exp(RcoreMu.a)*RcoreBeta51^YearC
RYnoco <- exp(RnocoMu.a)*RnocoBeta51^YearC
RYleks <- exp(RleksMu.a)*RleksBeta51^YearC

yMin <- 0
yMax <- 20


# all leks
png(filename=paste0(manuDir,'/B10 Plot of All Leks - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
plot(Year,RYleks,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)

for(i in 1:41){
  x <- Year[i:(i + 10)]                      # restrict to years of interest
  B10 <- Rleks.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
  mid <- median(x)                           # calculate temporal median of this B10 stretch
  B10.mu.a <- RYleks[mid - 1964]             # this is the intercept to use for this 10-year stretch
  B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
   
  par(new=TRUE)
  plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
}

axis(side=1,labels=TRUE,seq(1965,2015,5))
axis(side=2,labels=TRUE,seq(yMin,yMax,1))
dev.off()


# core
png(filename=paste0(manuDir,'/B10 Plot of Core - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
plot(Year,RYcore,type='l',col='darkgreen' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)

for(i in 1:41){
  x <- Year[i:(i + 10)]                      # restrict to years of interest
  B10 <- Rcore.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
  mid <- median(x)                           # calculate temporal median of this B10 stretch
  B10.mu.a <- RYcore[mid - 1964]             # this is the intercept to use for this 10-year stretch
  B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
  
  par(new=TRUE)
  plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
}

axis(side=1,labels=TRUE,seq(1965,2015,5))
axis(side=2,labels=TRUE,seq(yMin,yMax,1))
dev.off()


# non-core
png(filename=paste0(manuDir,'/B10 Plot of Non-Core - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
plot(Year,RYnoco,type='l',col='darkred' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)

for(i in 1:41){
  x <- Year[i:(i + 10)]                      # restrict to years of interest
  B10 <- Rnoco.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
  mid <- median(x)                           # calculate temporal median of this B10 stretch
  B10.mu.a <- RYnoco[mid - 1964]             # this is the intercept to use for this 10-year stretch
  B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
  
  par(new=TRUE)
  plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
}

axis(side=1,labels=TRUE,seq(1965,2015,5))
axis(side=2,labels=TRUE,seq(yMin,yMax,1))
dev.off()













# make plot of all mzone-core, mzone-non-core, mzone-all leks, states-all leks















# make plot of all states - all leks together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - All States.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:11){
  if(i != 8){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
  #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
  #   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  #   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
  #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
  #   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  #   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
  #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
  #   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  #   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
                          #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                          #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                          #                           Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
                          #                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                          Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
                        #                         'Y.core','Y.core.X5','Y.core.X95',
                        #                         'Y.noco','Y.noco.X5','Y.noco.X95',
                        'Y.leks','Y.leks.X5','Y.leks.X95',
                        #                         'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
                        #                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
                        'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
  
  
  x <- plotYears$Year
  #   yA1 <- plotYears$Y.core
  #   yA2 <- plotYears$Y.core.X5
  #   yA3 <- plotYears$Y.core.X95
  #   yB1 <- plotYears$Y.noco
  #   yB2 <- plotYears$Y.noco.X5
  #   yB3 <- plotYears$Y.noco.X95
  yC1 <- plotYears$Y.leks
  yC2 <- plotYears$Y.leks.X5
  yC3 <- plotYears$Y.leks.X95
  
  #   yRA1 <- plotYears$Y.Rcore
  #   yRA2 <- plotYears$Y.Rcore.X5
  #   yRA3 <- plotYears$Y.Rcore.X95
  #   yRB1 <- plotYears$Y.Rnoco
  #   yRB2 <- plotYears$Y.Rnoco.X5
  #   yRB3 <- plotYears$Y.Rnoco.X95
  yRC1 <- plotYears$Y.Rleks
  yRC2 <- plotYears$Y.Rleks.X5
  yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
    
  # state specific
  #   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
  #   par(new=TRUE)
  #   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  if(i > 1){
    par(new=TRUE)
  }
  if(i == 1){
    plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0("All States - Combined Leks"))    
  } else {
    plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)     
  }
  if(i == 11){
    par(new=TRUE)
    plot(x,yRC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=5)
  }
#   par(new=TRUE)
#   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
  #   par(new=TRUE)
  #   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  

#   par(new=TRUE)
#   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   par(new=TRUE) 
#   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  }
}
legend("topright",c(states[-8],"Rangewide"),lty=c(rep(1,10),1),lwd=c(rep(2,10),5),col=c(colVec[-8],'black'),bty="n",ncol=2)
axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
dev.off()





















# make plot of all mzones - all leks together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
  #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
  #   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  #   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
  #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
  #   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  #   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
  #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
  #   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  #   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
                          #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                          #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                          #                           Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
                          #                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                          Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
                        #                         'Y.core','Y.core.X5','Y.core.X95',
                        #                         'Y.noco','Y.noco.X5','Y.noco.X95',
                        'Y.leks','Y.leks.X5','Y.leks.X95',
                        #                         'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
                        #                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
                        'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
  
  
  x <- plotYears$Year
  #   yA1 <- plotYears$Y.core
  #   yA2 <- plotYears$Y.core.X5
  #   yA3 <- plotYears$Y.core.X95
  #   yB1 <- plotYears$Y.noco
  #   yB2 <- plotYears$Y.noco.X5
  #   yB3 <- plotYears$Y.noco.X95
  yC1 <- plotYears$Y.leks
  yC2 <- plotYears$Y.leks.X5
  yC3 <- plotYears$Y.leks.X95
  
  #   yRA1 <- plotYears$Y.Rcore
  #   yRA2 <- plotYears$Y.Rcore.X5
  #   yRA3 <- plotYears$Y.Rcore.X95
  #   yRB1 <- plotYears$Y.Rnoco
  #   yRB2 <- plotYears$Y.Rnoco.X5
  #   yRB3 <- plotYears$Y.Rnoco.X95
  yRC1 <- plotYears$Y.Rleks
  yRC2 <- plotYears$Y.Rleks.X5
  yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  # state specific
  #   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
  #   par(new=TRUE)
  #   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  if(i > 1){
    par(new=TRUE)
  }
  if(i == 1){
    plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0("All Management Zones - Combined Leks"))    
  } else {
    plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)#,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% All States"))
  }
  if(i == 6){
    par(new=TRUE)
    plot(x,yRC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=5)
  }
  #   par(new=TRUE)
  #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
  #   par(new=TRUE)
  #   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE)
  #   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)

  #   par(new=TRUE)
  #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
}
legend("topright",c(c(mZones2[1],mZones2[6],mZones2[2:5]),"Rangewide"),lty=c(rep(1,6),1),lwd=c(rep(2,6),5),col=c(c(colVec[1],colVec[6],colVec[2:5]),'black'),bty="n",ncol=2)
axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
dev.off()
































# make plot of all mzones - core together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of C - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
    bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
    Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
    coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
    RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  #   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #     coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
    coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  #   leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
    RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  #   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  #   RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
    coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  #   leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
    RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  #   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  #   RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
    coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  #   leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
    RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  #   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  #   RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
    Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
    Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
    Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
    Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
    Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
    Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  #   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  #   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  #   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
                                                    Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                          #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          #                           Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                                                    Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95)
                          #                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
                          #                           Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
                                                'Y.core','Y.core.X5','Y.core.X95',
                        #                         'Y.noco','Y.noco.X5','Y.noco.X95',
                        #                       'Y.leks','Y.leks.X5','Y.leks.X95',
                                                'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95')
                        #                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
                        #                        'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95'
  
  
  x <- plotYears$Year
    yA1 <- plotYears$Y.core
    yA2 <- plotYears$Y.core.X5
    yA3 <- plotYears$Y.core.X95
  #   yB1 <- plotYears$Y.noco
  #   yB2 <- plotYears$Y.noco.X5
  #   yB3 <- plotYears$Y.noco.X95
  #   yC1 <- plotYears$Y.leks
  #   yC2 <- plotYears$Y.leks.X5
  #   yC3 <- plotYears$Y.leks.X95
  
    yRA1 <- plotYears$Y.Rcore
    yRA2 <- plotYears$Y.Rcore.X5
    yRA3 <- plotYears$Y.Rcore.X95
  #   yRB1 <- plotYears$Y.Rnoco
  #   yRB2 <- plotYears$Y.Rnoco.X5
  #   yRB3 <- plotYears$Y.Rnoco.X95
  #   yRC1 <- plotYears$Y.Rleks
  #   yRC2 <- plotYears$Y.Rleks.X5
  #   yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  if(i > 1){
    par(new=TRUE)
  }
  if(i == 1){
    plot(x,yA1,type='l',col=colVec[i],axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0("All Management Zones - Core Leks"))
  } else {
    plot(x,yA1,type='l',col=colVec[i],axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)    
  }
  if(i == 6){ 
    par(new=TRUE)
    plot(x,yRA1,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=5)
  }
  #     par(new=TRUE)
  #     plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #     par(new=TRUE)
  #     plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  

  #   plot(x,yC1,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)#,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% All States"))
  #   par(new=TRUE)
  #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
  #     par(new=TRUE)
  #     plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #     par(new=TRUE)
  #     plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
  #   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  #   par(new=TRUE)
  #   plot(x,yRC1,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
}
legend("topright",c(c(mZones2[1],mZones2[6],mZones2[2:5]),"Rangewide"),lty=c(rep(1,6),1),lwd=c(rep(2,6),5),col=c(c(colVec[1],colVec[6],colVec[2:5]),'black'),bty="n",ncol=2)
axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
dev.off()










# make plot of all mzones - non-core together
manuDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Results/Manuscript 2015.08.12'

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of NC - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  #   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
  #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
  #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
  #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
  #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
  #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
  #   
  #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
  #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
  #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
  #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
  #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
  #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
  
  #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
  nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
  #   leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
  #   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
  RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
  #   RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
  
  #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
  nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
  #   leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
  #   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
  RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
  #   RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
  
  #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
  nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
  #   leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
  #   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
  RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
  #   RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
  
  
  
  Year <- data.frame(Year=seq(1965,2015,1))
  YearC <- Year - 1964 - 26
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  #   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  #   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  #   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  plotYears <- data.frame(Year=Year,YearC=YearC, 
                          #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
                                                     Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
                          #                           Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
                          #   Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95)
                              Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95)
  #                           Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
  names(plotYears) <- c('Year','YearC',
                        #'Y.core','Y.core.X5','Y.core.X95',
                                                 'Y.noco','Y.noco.X5','Y.noco.X95',
                        #                       'Y.leks','Y.leks.X5','Y.leks.X95',
                        #'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95')
                           'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95')
  #                        'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95'
  
  
  x <- plotYears$Year
  #   yA1 <- plotYears$Y.core
  #   yA2 <- plotYears$Y.core.X5
  #   yA3 <- plotYears$Y.core.X95
  yB1 <- plotYears$Y.noco
  yB2 <- plotYears$Y.noco.X5
  yB3 <- plotYears$Y.noco.X95
  #   yC1 <- plotYears$Y.leks
  #   yC2 <- plotYears$Y.leks.X5
  #   yC3 <- plotYears$Y.leks.X95
  
  #   yRA1 <- plotYears$Y.Rcore
  #   yRA2 <- plotYears$Y.Rcore.X5
  #   yRA3 <- plotYears$Y.Rcore.X95
  yRB1 <- plotYears$Y.Rnoco
  yRB2 <- plotYears$Y.Rnoco.X5
  yRB3 <- plotYears$Y.Rnoco.X95
  #   yRC1 <- plotYears$Y.Rleks
  #   yRC2 <- plotYears$Y.Rleks.X5
  #   yRC3 <- plotYears$Y.Rleks.X95
  
  yMax <- 30
  
  if(i > 1){
    par(new=TRUE)
  }
  if(i == 1){
    plot(x,yB1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main="All Management Zones - Periphery Leks")
  } else {
    plot(x,yB1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)
  }
  if(i == 6){
    par(new=TRUE)
    plot(x,yRB1,type='l',col='black'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=5)
  }
  #     plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)#,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
  #     par(new=TRUE)
  #     plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #     par(new=TRUE)
  #     plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   
  #   par(new=TRUE)
 
#     plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     par(new=TRUE) 
#     plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  
  #   plot(x,yC1,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)#,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% All States"))
  #   par(new=TRUE)
  #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  # rangewide specific
  #   par(new=TRUE)
  #   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #     par(new=TRUE)
  #     plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #     par(new=TRUE)
  #     plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   

  #     par(new=TRUE)
#     plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     par(new=TRUE) 
#     plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
  #   par(new=TRUE)
  #   plot(x,yRC1,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
  #   par(new=TRUE)
  #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  #   par(new=TRUE) 
  #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
  
}
legend("topright",c(c(mZones2[1],mZones2[6],mZones2[2:5]),"Rangewide"),lty=c(rep(1,6),1),lwd=c(rep(2,6),5),col=c(c(colVec[1],colVec[6],colVec[2:5]),'black'),bty="n",ncol=2)
axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
dev.off()





































































# 
# 
# 
# 
# 
# 
# 
# 
# 
# # load in data for sample means 
# dat <- dat1stZerosCore[[9]]    # loop over cores?
# #dat <- dat1stZerosNoco[[9]]  
# 
# # load(paste0(outpDir,'/smallCoreSamp2.Sonic2.RData'))
# # dat <- smallCoreSamp2
# 
# 
# 
# # nZones <- length(unique(dat$mZone_num))
# 
# 
# # load(paste0(outpDir,'/',loadThis,'.RData')) 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # try2 <- bayes$summary
# # try3 <- bayes$summary
# 
# 
# # load model 
# # these <- grep(paste0("Model ",modelLetter), theFiles , ignore.case=FALSE, fixed=TRUE)
# # zoneString <- unlist(lapply(strsplit(theFiles,".",fixed=TRUE),function(x) strsplit(x,".",fixed=TRUE)[[1]][1]))
# # zones <- as.numeric(substr(zoneString,nchar(zoneString),nchar(zoneString)))
# # loadThese <- paste0(outpDir,"/",theFiles[these])
# # assign(paste0("nModel",modelLetter),length(loadThese))
# # 
# 
# 
# 
# 
# # these <- 6
# # load(paste0(outpDir,"/",theFiles[these]))  
# 
# 
# 
# Year <- data.frame(Year=seq(1965,2015,1))
# 
# indx <- data.frame(bayes$summary)
# indx$vars <- rownames(indx)
# indx$x <- ifelse(substr(indx$vars,1,1) == 'N',substr(indx$vars,3,3),-99)
# indx$y <- ifelse(substr(indx$vars,1,1) == 'N',ifelse(nchar(indx$vars) == 6,substr(indx$vars,5,5),substr(indx$vars,5,6)),-99)
# indx$x0 <- ifelse(substr(indx$vars,1,2) == 'N0',ifelse(nchar(indx$vars) == 5,substr(indx$vars,4,4),substr(indx$vars,4,5)),-99)
#   
# # build the true B-matrix
# beta.mzone <- matrix(NA,nrow=nZones,ncol=51)
# list.beta.mzone <- indx[substr(indx$vars,1,2) == 'N[',]
# for(i in 1:nZones){
#   for(j in 1:51){
#     beta.mzone[i,j] <- list.beta.mzone[list.beta.mzone$x == i & list.beta.mzone$y == j,]$mean
#   }
# }
# 
# 
# 
# 
# 
# if(grepl("Non-Core", loadThis) == 1){
#   coreText <- "Non-Core"
# } else {
#   coreText <- "Core"
# }
# 
# dev.off()
# par(mfrow=c(3,2))
# 
# for(i in 1:nZones){
#   if(i == 1){mZone <- 1}
#   if(i == 2){mZone <- 3}
#   if(i == 3){mZone <- 4}
#   if(i == 4){mZone <- 5}
#   if(i == 5){mZone <- 6}
#   if(i == 6){mZone <- 8}
#   
#   theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])
#   
#   thisOne <- dat[dat$mZone_num == mZone,]
#   obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
#   obsMeans$Year <- rownames(obsMeans)
#   
#   plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
#   
#   x  <- plotYears$Year
#   y1 <- plotYears$N
#   y2 <- plotYears$MeanPMales
#   
#   yMax <- 100
#   
#   plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",coreText," - Management Zone ",mZone))
#   par(new=TRUE)
#   plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   axis(side=1,labels=TRUE,seq(1965,2015,5))
#   axis(side=2,labels=TRUE,seq(0,yMax,10))
#   
#   rm(plotYears,x,y1,y2)
# }
# 
# 
# 
# 
# 
# 
# 
# 
#   
#   
# dev.off()
# # investigate regionwide trends
# # build the true B-matrix
# beta.range <- rep(NA,51)
# list.beta.range <- indx[substr(indx$vars,1,3) == 'N0[',]
# for(k in 1:51){
#   beta.range[k] <- list.beta.range[list.beta.range$x0 == k,]$mean
# }
# 
# theN <- data.frame(Year=seq(1965,2015,1),N=beta.range)
# 
# thisOne <- dat
# obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
# obsMeans$Year <- rownames(obsMeans)
# 
# plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
# 
# x  <- plotYears$Year
# y1 <- plotYears$N
# y2 <- plotYears$MeanPMales
# i <- 1
# 
# yMax <- 100
# 
# plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\n75% ",coreText," Core - Rangewide"))
# par(new=TRUE)
# plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
# axis(side=1,labels=TRUE,seq(1965,2015,5))
# axis(side=2,labels=TRUE,seq(0,yMax,10))
# 
# rm(plotYears,x,y1,y2)  
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# B10 <- data.frame(bayes$summary[substr(rownames(bayes$summary),1,3) == 'B10',])
# 
# B10$mZone <- as.numeric(substr(rownames(B10),11,11))
# 
# for(i in 1:nZones){
#   if(i == 1){mZone <- 1}
#   if(i == 2){mZone <- 3}
#   if(i == 3){mZone <- 4}
#   if(i == 4){mZone <- 5}
#   if(i == 5){mZone <- 6}
#   if(i == 6){mZone <- 8}
#   
#   theN <- data.frame(Year=seq(1965,2015,1),N=beta.mzone[i,])
#   
#   thisOne <- dat[dat$mZone_num == mZone,]
#   obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
#   obsMeans$Year <- rownames(obsMeans)
#   
#   plotYears <- merge(obsMeans,theN,by=c('Year'),all.x=TRUE)
#   
#   x  <- plotYears$Year
#   y1 <- plotYears$N
#   y2 <- plotYears$MeanPMales
#   
#   yMax <- 100
#   
#   plot(x,y1,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% Core - Management Zone ",mZone))
#   par(new=TRUE)
#   plot(x,y2,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#   axis(side=1,labels=TRUE,seq(1965,2015,5))
#   axis(side=2,labels=TRUE,seq(0,yMax,10))
#   
#   rm(plotYears,x,y1,y2)  
# }
# 
# 
# g <- smallCoreSamp2$Peak_Males[smallCoreSamp2$Peak_Males <= 100]
# h <- hist(g, breaks=100, density=10, col="lightgray", xlab="Accuracy", main="Overall") 
# xfit <- seq(min(g),max(g),length=40) 
# yfit <- dpois(xfit,lambda=mean(g)) 
# yfit <- yfit*diff(h$mids[1:2])*length(g) 
# lines(xfit, yfit, col="black", lwd=2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# allResults <- rbind(estsC,estsD,estsF)
# rm(estsC,estsD,estsF)
# names(allResults)[names(allResults) == 'Zone'] <- 'mZone'
# allResults <- data.frame(allResults,ParmType=unlist(lapply(strsplit(as.character(droplevels(allResults$Parameter)),'[',fixed=TRUE), function(x) x[[1]][1])))
# 
# allResults <- allResults[order(allResults$mZone,allResults$modCode),]
# 
# mu.aAll <- allResults[allResults$Parameter == 'mu.a',]
# mu.bAll <- allResults[allResults$Parameter == 'mu.b',]
# rhosAll <- allResults[allResults$Parameter == 'rho',]
# taonoiseAll <- allResults[allResults$Parameter == 'taonoise',]
# 
# SumData <- read.csv(paste0(analDir,'/SumData.csv'))
# 
# yMax <- max(datList[[1]]$Peak_Males)
# 
# helper <- allResults[,c('Model','Zeros','Cut','mZone','modCode')]
# helper <- unique(helper)
# rownames(helper) <- NULL
# 
# 
# 
# 
# 
# doThese <- c('D')
# nDoThese <- length(doThese)
# table1 <- NULL
# colVec <- brewer.pal(8,"Set1")
# 
# for(i in 1:8){
#   
#   #A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
#   #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
#   #C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
#   
#   D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
#   D <- smallCoreSamp
#   #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
#   #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
#   
#   for(j in 1:nDoThese){
#     
#     thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]
#     
#     if(doThese[j] %in% c('D','E','F')){
#       thisOne <- get(doThese[j])
#       #thisOne <- thisOne[thisOne$DupZero == 0,]  # dont think this is necessary when doing core
#     } else {
#       thisOne <- get(doThese[j])   
#     }
#     leks <- unique(thisOne$Lek_IDnum)
#     thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
#     thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID))                     
#     nLeks <- length(leks)
#   
#     Model <- as.character(droplevels(thisHelper$Model))
#     Zeros <- as.character(droplevels(thisHelper$Zeros))
#     Cut   <- as.character(droplevels(thisHelper$Cut))
#     mZone <- i
#     
#     # make individual lek trends
# #     for(k in 1:nLeks){    
# #       # do this for winbugs numbering (so on factor)
# #       xlek <- thisOne[thisOne$lekCls == k,]$Year - 1990
# #       ylek <- thisOne[thisOne$lekCls == k,]$Peak_Males
# #       the.lek <- thisOne[thisOne$lekCls == k,]$Lek_ID[1]
# #       lek.mu.a <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('a[',k,']'),]$mean
# #       lek.mu.b <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('b[',k,']'),]$mean
# #       ypred <- exp(lek.mu.a + lek.mu.b*xlek)
# #       
# #       CairoPNG(filename=paste0(rsltDir,'/Lek Trends/MZone ',i,'/',the.lek,' - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',mZone,'.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)
# #       
# #       plot(xlek + 1990,ylek,col='gray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5,pch=19)
# #       par(new=TRUE)
# #       plot(xlek + 1990,ypred,col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Lek ",the.lek," in Management Zone ",i,"\nYears 1965-2014"))
# #       axis(side=1,labels=TRUE,seq(1965,2015,5))
# #       axis(side=2,labels=TRUE,seq(0,yMax,10))
# #       
# #       dev.off()
# #     }
#     
#    CairoPNG(filename=paste0(rsltDir,'/Trends + Yes Tao - No Leks - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',mZone,'.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)
# 
#     # make gray individual lek trends
#     for(k in 1:nLeks){
#       
#       # do this on k numbering
#       x <- thisOne[thisOne$Lek_IDnum == leks[k],]$Year
#       y <- thisOne[thisOne$Lek_IDnum == leks[k],]$Peak_Males
#             
#       yMax <- 100#max(thisOne$Peak_Males)
#       
# #       if(k == 1){
# #         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
# #       } else if(k > 1 & k < nLeks){
# #         par(new=TRUE)
# #         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
# #       } else if(k == nLeks){
# #         par(new=TRUE)
# #         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='Year',ylab='Peak Males',lwd=0.5)        
# #       }
#     }
#      
# #     get mean trend info for this mzone and model
#     mu.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$mean
#     lo.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X5.
#     hi.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X95.
#     mu.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$mean
#     lo.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X5.
#     hi.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X95. 
#     tao     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$mean
#     lo.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X5.
#     hi.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X95.
#     rho     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$mean
#     lo.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X5.
#     hi.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X95.
#     sigma.a <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.a',]$mean
#     sigma.b <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.b',]$mean
# 
#     thisRow <- data.frame(thisHelper,mu.a=mu.a,lo.a,hi.a,mu.b,lo.b,hi.b,tao,lo.tao,hi.tao,rho,lo.rho,hi.rho,sigma.a,sigma.b)
#     table1 <- rbind(table1,thisRow)
# 
#     obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
#     obsMeans$Year <- rownames(obsMeans)
# 
#     year <- seq(-25,25,1)
#     est <- exp(mu.a + mu.b*year + 0.5*tao)
#     se <- sqrt(sigma.a^2 + sigma.b^2 + 2*rho*sigma.a*sigma.b)
#     lo <- exp(mu.a + mu.b*year + 0.5*tao) - 1.645*se 
#     hi <- exp(mu.a + mu.b*year + 0.5*tao) + 1.645*se 
#     x <- year + 1990
#     
#     #par(new=TRUE)
#     plot(obsMeans$Year,obsMeans$MeanPMales,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\n75% Core - Management Zone ",i))
#     par(new=TRUE)
#     plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     par(new=TRUE)
#     plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     par(new=TRUE)
#     plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)    
#     axis(side=1,labels=TRUE,seq(1965,2015,5))
#     axis(side=2,labels=TRUE,seq(0,yMax,10))
#     legend("topleft",bty = "n",legend=c("Annual Mean","Temporal Trend & 90% Credible Intervals","Individual Lek Counts"),pch=c(19,NA,NA),lwd=c(NA,2,1),col=c('darkgray',colVec[i],"gray95"))
#     dev.off()
#   }
# }
# 
# write.csv(table1,paste0(rsltDir,'/table1 + Yes Tao.csv'))
# 
# table1 <- read.csv(paste0(rsltDir,'/table1 + No Tao.csv'))
# 
# 
# CairoPNG(filename=paste0(rsltDir,'/Trends + No Tao - All MZones - 1st Zeros - Core Leks.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)
# 
# for(i in 1:8){
#   thisRow <- table1[i,]
#   year <- seq(-25,25,1)
#   est <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao)
#   se <- sqrt(thisRow$sigma.a^2 + thisRow$sigma.b^2 + 2*thisRow$rho*thisRow$sigma.a*thisRow$sigma.b)
#   lo <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) - 1.645*se
#   hi <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) + 1.645*se
#   x <- year + 1990
#   if(i <= 7){
#     plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
# #     par(new=TRUE)
# #     plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
# #     par(new=TRUE)
# #     plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)   
#     par(new=TRUE)
#   } else if(i == 8){
#     plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\nAll Management Zones"))
#     axis(side=1,labels=TRUE,seq(1965,2015,5))
#     axis(side=2,labels=TRUE,seq(0,yMax,10))
#     legend("topright",bty="n",legend=c('I','II','III','IV','V','VI','VII','VIII'),col=colVec,lwd=rep(2,8))
#   }
# 
# }
# dev.off()
# 
# 
# 
# 
# 
# doThese <- c('A','C')
# nDoThese <- length(doThese)
# the.zeros <- NULL
# 
# 
# # read in all the data
# for(i in 9:9){
#   
#    A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
#   #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
#    C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
#   
#   #D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
#   #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
#   #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
# }
#                                                                                                                                                                              # close out if
# states <- as.character(droplevels(unique(A$State)))
# nStates <- length(states)
# 
# for(i in 1:nStates){
# 
#   for(j in 1:nDoThese){
#     
# #   thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]
# 
#     thisOne <- get(doThese[j])
#     thisOne <- thisOne[thisOne$State == states[i],]
#     
#     thisOne$zeroPresent <- ifelse(thisOne$Peak_Males == 0,1,0)
# 
#     leks <- length(unique(thisOne$Lek_IDnum))
#     leks.zero <- length(unique(thisOne[thisOne$zeroPresent == 1,]$Lek_IDnum))
#     leks.zero.dupped <- length(unique(thisOne[thisOne$zeroPresent == 1 & thisOne$DupZero == 1,]$Lek_IDnum))
#     
# #     thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
# #     thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID))                     
# #     nLeks <- length(leks)
# #     
# #     Model <- as.character(droplevels(thisHelper$Model))
# #     Zeros <- as.character(droplevels(thisHelper$Zeros))
# #     Cut   <- as.character(droplevels(thisHelper$Cut))
# #     mZone <- i
#     
#     thisOne$Lek_ID <- as.character(droplevels(thisOne$Lek_ID))
# 
#     zeroInv <- data.frame(NDuppedZeros=tapply(thisOne$DupZero, thisOne$Lek_ID, FUN=sum))
#     zeroInv$Lek <- rownames(zeroInv) 
#     the.zero.mean <- data.frame(Model=doThese[j],NLeks=leks,NLeksZero=leks.zero,NLeksZeroDupped=leks.zero.dupped,AvgDuppedZeroMeanLek=mean(zeroInv[zeroInv$NDuppedZeros > 0,]$NDuppedZeros))
#     the.zero.mean$State <- states[i]
#     the.zero.mean$PercentLeksWithDup <- round(100*the.zero.mean$NLeksZeroDupped / the.zero.mean$NLeksZero,2)
#     the.zeros <- rbind(the.zeros,the.zero.mean)
#     write.csv(the.zeros[the.zeros$Model == 'C',],paste0(rsltDir,'/the.zeros.modelC.csv'))
#   }
# }
# 
# 
# 
# 
