


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
        makeRangeWide("Model D","1st",units,'N',dat1stZerosCore,"1st Zeros - Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model E","1st",units,'N',dat1stZerosNoco,"1st Zeros - Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model F","1st",units,'N',dat1stZerosLeks,"1st Zeros - All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
        
        makeRangeWide("Model D","all",units,'N',datAllZerosCore,"All Zeros - Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model E","all",units,'N',datAllZerosNoco,"All Zeros - Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model F","all",units,'N',datAllZerosLeks,"All Zeros - All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
        
        makeRangeWide("Model D","1st",units,'Y',dat1stZerosCore,"Y05-15 - 1st Zeros - Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model E","1st",units,'Y',dat1stZerosNoco,"Y05-15 - 1st Zeros - Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model F","1st",units,'Y',dat1stZerosLeks,"Y05-15 - 1st Zeros - All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
        
        makeRangeWide("Model D","all",units,'Y',datAllZerosCore,"Y05-15 - All Zeros - Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model E","all",units,'Y',datAllZerosNoco,"Y05-15 - All Zeros - Non-Core",paste0(outpDir,'/Rangewide/Trend Plots'))
        makeRangeWide("Model F","all",units,'Y',datAllZerosLeks,"Y05-15 - All Zeros - All Leks",paste0(outpDir,'/Rangewide/Trend Plots'))
        # make core + non-core + all lek plot
      } 
    }
  }
}












# build big file for first zeros, mzones and states, all years

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("RData", theFiles) == 1]
theFiles <- theFiles[grepl("Try 1", theFiles) == 1]
# theFiles <- theFiles[grepl(c("CA"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("CO"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("ID"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("MT"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("ND"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("NV"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("OR"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("SD"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("UT"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("WA"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("WY"), theFiles) == 0]
theFiles <- theFiles[grepl(c("ZInf"), theFiles) == 0]
colVec <- brewer.pal(12,"Paired")

colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]

states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZones2 <- c("Zone 1","Zone 3","Zone 4","Zone 5","Zone 6","Zone 2 & 7")
textyZones <- c("Great Plains","Southern Great Basin","Snake River Plain","Northern Great Basin","Columbian Basin","Wyoming Basin & Colorado Plateau")

nStates <- length(states)
nMZones <- length(mZones)
textyStates <- c("California","Colorado","Idaho","Montana","North Dakota","Nevada","Oregon","South Dakota","Utah","Washington","Wyoming")

# readInAnalyticFiles(analDir)


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
        file3 <- '1st Zeros - Core'
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type  
        allNLeks <- length(unique(dat1stZerosNoco[[9]]$Lek_ID))
        file2 <- 'E'
        file3 <- '1st Zeros - Non-Core'
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
        file3 <- '1st Zeros - All Leks'
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
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file3,'.csv'))
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
write.csv(table1,paste0(manuDir,'/table1 - 1st Zeros.csv'))









# build big file for all zeros, mzones and states, all years

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("RData", theFiles) == 1]
theFiles <- theFiles[grepl("All Zeros", theFiles) == 1]
theFiles <- theFiles[grepl("2005-2015", theFiles) == 0]
# theFiles <- theFiles[grepl("State", theFiles) == 0]
colVec <- brewer.pal(12,"Paired")

colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]

states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZones2 <- c("Zone 1","Zone 3","Zone 4","Zone 5","Zone 6","Zone 2 & 7")
textyZones <- c("Great Plains","Southern Great Basin","Snake River Plain","Northern Great Basin","Columbian Basin","Wyoming Basin & Colorado Plateau")

nStates <- length(states)
nMZones <- length(mZones)
textyStates <- c("California","Colorado","Idaho","Montana","North Dakota","Nevada","Oregon","South Dakota","Utah","Washington","Wyoming")

# readInAnalyticFiles(analDir)


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
        file3 <- 'All Zeros - Core'
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
          dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type  
        allNLeks <- length(unique(dat1stZerosNoco[[9]]$Lek_ID))
        file2 <- 'E'
        file3 <- 'All Zeros - Non-Core'
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
        file3 <- 'All Zeros - All Leks'
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
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file3,'.csv'))
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
write.csv(table1,paste0(manuDir,'/table1 - All Zeros.csv'))






























# build big file for first zeros, mzones and states, last 11 years

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("1st Zeros 2005-2015.RData", theFiles) == 1]
# theFiles <- theFiles[grepl(c("CA"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("CO"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("ID"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("MT"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("ND"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("NV"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("OR"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("SD"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("UT"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("WA"), theFiles) == 0]
# theFiles <- theFiles[grepl(c("WY"), theFiles) == 0]
theFiles <- theFiles[grepl(c("ZInf"), theFiles) == 0]
colVec <- brewer.pal(12,"Paired")

colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]

states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZones2 <- c("Zone 1","Zone 3","Zone 4","Zone 5","Zone 6","Zone 2 & 7")
textyZones <- c("Great Plains","Southern Great Basin","Snake River Plain","Northern Great Basin","Columbian Basin","Wyoming Basin & Colorado Plateau")

nStates <- length(states)
nMZones <- length(mZones)
textyStates <- c("California","Colorado","Idaho","Montana","North Dakota","Nevada","Oregon","South Dakota","Utah","Washington","Wyoming")

# readInAnalyticFiles(analDir)


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
#     units <- states
#     nUnits <- nStates
#     start <- 3
#     string <- 'State'
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
#         load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
      #       } else {
      #         load(paste0(outpDir,'/',))
      #       }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
#           dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
        allNLeks <- length(unique(dat1stZerosCore[[9]]$Lek_ID))
        file2 <- 'D'
        file3 <- 'Y05-15 - 1st Zeros - Core'
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
#           dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type  
        allNLeks <- length(unique(dat1stZerosNoco[[9]]$Lek_ID))
        file2 <- 'E'
        file3 <- 'Y05-15 - 1st Zeros - Non-Core'
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
#           dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type   
        allNLeks <- length(unique(dat1stZerosLeks[[9]]$Lek_ID))
        file2 <- 'F'
        file3 <- 'Y05-15 - 1st Zeros - All Leks'
      }
      
      if(h == 1 ){#| h == 2){
        dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
        
      }
      
      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else if(h == 2){
#         file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      } else {
        file <- 'Rangewide/Trend Plots'                                                             # get name of new folder 
      }
      
      if(h == 1 ){#| h == 2){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
        mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.a'   ,mean=thisOne[thisOne$X == 'mu.a'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a'   ,]$X95.)
        mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.b'   ,mean=thisOne[thisOne$X == 'mu.b'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b'   ,]$X95.)
        beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='beta'   ,mean=thisOne[thisOne$X == 'beta[1]',]$mean,X5.=thisOne[thisOne$X == 'beta[1]',]$X5.,X95.=thisOne[thisOne$X == 'beta[1]',]$X95.)
        sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sdnoise',mean=thisOne[thisOne$X == 'sdnoise',]$mean,X5.=thisOne[thisOne$X == 'sdnoise',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise',]$X95.)
        rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='rho'    ,mean=thisOne[thisOne$X == 'rho'    ,]$mean,X5.=thisOne[thisOne$X == 'rho'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho'    ,]$X95.)
        sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.a',mean=thisOne[thisOne$X == 'sigma.a',]$mean,X5.=thisOne[thisOne$X == 'sigma.a',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a',]$X95.)
        sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.b',mean=thisOne[thisOne$X == 'sigma.b',]$mean,X5.=thisOne[thisOne$X == 'sigma.b',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b',]$X95.)
      } else if (h == 3){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file3,'.csv'))
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
theSumms <- unique(theSumms)
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
write.csv(table1,paste0(manuDir,'/table1 - Y05-Y15 - 1st Zeros.csv'))









# build big file for all zeros, mzones and states, 11 years

theFiles <- list.files(outpDir)
theFiles <- theFiles[grepl("All Zeros 2005-2015.RData", theFiles) == 1]
colVec <- brewer.pal(12,"Paired")

colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]

states <- c("CA","CO","ID","MT","ND","NV","OR","SD","UT","WA","WY")
mZones <- c("MZone 1","MZone 3","MZone 4","MZone 5","MZone 6","MZone 8")
mZones2 <- c("Zone 1","Zone 3","Zone 4","Zone 5","Zone 6","Zone 2 & 7")
textyZones <- c("Great Plains","Southern Great Basin","Snake River Plain","Northern Great Basin","Columbian Basin","Wyoming Basin & Colorado Plateau")

nStates <- length(states)
nMZones <- length(mZones)
textyStates <- c("California","Colorado","Idaho","Montana","North Dakota","Nevada","Oregon","South Dakota","Utah","Washington","Wyoming")

# readInAnalyticFiles(analDir)


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
#     units <- states
#     nUnits <- nStates
#     start <- 3
#     string <- 'State'
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
#         load(paste0(outpDir,'/',unitFile[1]))                                                # get bayesian stuff -- state        
      }
      #       } else {
      #         load(paste0(outpDir,'/',))
      #       }
      if(j == 1){
        if(h == 1){
          dat <- dat1stZerosCore[[numUnit]]                                                  # get orig data
        } else {
#           dat <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Core'                                                                    # assign run type
        allNLeks <- length(unique(dat1stZerosCore[[9]]$Lek_ID))
        file2 <- 'D'
        file3 <- 'Y05-15 - All Zeros - Core'
      } else if(j == 2){
        if(h == 1){
          dat <- dat1stZerosNoco[[numUnit]]                                                  # get orig data   
        } else {
#           dat <- dat1stZerosNoco[[9]][dat1stZerosNoco[[9]]$State == theUnit,]                # get state data          
        }
        runType <- 'Non-Core'                                                                # assign run type  
        allNLeks <- length(unique(dat1stZerosNoco[[9]]$Lek_ID))
        file2 <- 'E'
        file3 <- 'Y05-15 - All Zeros - Non-Core'
      } else {
        if(h == 1){
          dat <- dat1stZerosLeks[[numUnit]]                                                  # get mzone data
          #dat <- datAllZerosLeks[[6]]
        } else { 
#           dat <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == theUnit,]                # get state data
        }
        runType <- 'All Leks'                                                                # assign run type   
        allNLeks <- length(unique(dat1stZerosLeks[[9]]$Lek_ID))
        file2 <- 'F'
        file3 <- 'Y05-15 - All Zeros - All Leks'
      }
      
      if(h == 1 ){#| h == 2){
        dat$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(dat$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))
        
      }
      
      # make a folder where needed, for files in tracDir, if it doesn't already exist.
      if(h == 1){
        file <- substr(unitFile[j],1,nchar(unitFile[j]) - 6)                                        # get name of new folder                
      } else if(h == 2){
#         file <- substr(unitFile[1],1,nchar(unitFile[1]) - 6)                                        # get name of new folder        
      } else {
        file <- 'Rangewide/Trend Plots'                                                             # get name of new folder 
      }
      
      if(h == 1 ){#| h == 2){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
        mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.a'   ,mean=thisOne[thisOne$X == 'mu.a'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a'   ,]$X95.)
        mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='mu.b'   ,mean=thisOne[thisOne$X == 'mu.b'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b'   ,]$X95.)
        beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='beta'   ,mean=thisOne[thisOne$X == 'beta[1]',]$mean,X5.=thisOne[thisOne$X == 'beta[1]',]$X5.,X95.=thisOne[thisOne$X == 'beta[1]',]$X95.)
        sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sdnoise',mean=thisOne[thisOne$X == 'sdnoise',]$mean,X5.=thisOne[thisOne$X == 'sdnoise',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise',]$X95.)
        rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='rho'    ,mean=thisOne[thisOne$X == 'rho'    ,]$mean,X5.=thisOne[thisOne$X == 'rho'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho'    ,]$X95.)
        sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.a',mean=thisOne[thisOne$X == 'sigma.a',]$mean,X5.=thisOne[thisOne$X == 'sigma.a',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a',]$X95.)
        sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=nrow(thisOne[substr(thisOne$X,1,1) == 'a',]),stat='sigma.b',mean=thisOne[thisOne$X == 'sigma.b',]$mean,X5.=thisOne[thisOne$X == 'sigma.b',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b',]$X95.)
      } else if (h == 3){
        thisOne <- read.csv(paste0(outpDir,'/',file,'/bayesSummary - ',file3,'.csv'))
        mu.a    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='mu.a0'   ,mean=thisOne[thisOne$X == 'mu.a0'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.a0'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.a0'   ,]$X95.)
        mu.b    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='mu.b0'   ,mean=thisOne[thisOne$X == 'mu.b0'   ,]$mean,X5.=thisOne[thisOne$X == 'mu.b0'   ,]$X5.,X95.=thisOne[thisOne$X == 'mu.b0'   ,]$X95.)
        beta    <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='beta0'   ,mean=thisOne[thisOne$X == 'beta0'   ,]$mean,X5.=thisOne[thisOne$X == 'beta0'   ,]$X5.,X95.=thisOne[thisOne$X == 'beta0'   ,]$X95.)
        sdnoise <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sdnoise0',mean=thisOne[thisOne$X == 'sdnoise0',]$mean,X5.=thisOne[thisOne$X == 'sdnoise0',]$X5.,X95.=thisOne[thisOne$X == 'sdnoise0',]$X95.)
        rho     <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='rho0'    ,mean=thisOne[thisOne$X == 'rho0'    ,]$mean,X5.=thisOne[thisOne$X == 'rho0'    ,]$X5.,X95.=thisOne[thisOne$X == 'rho0'    ,]$X95.)
        sigma.a <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sigma.a0',mean=thisOne[thisOne$X == 'sigma.a0',]$mean,X5.=thisOne[thisOne$X == 'sigma.a0',]$X5.,X95.=thisOne[thisOne$X == 'sigma.a0',]$X95.)
        sigma.b <- data.frame(focus=string,group=theUnit,cut=runType,Nleks=allNLeks,stat='sigma.b0',mean=thisOne[thisOne$X == 'sigma.b0',]$mean,X5.=thisOne[thisOne$X == 'sigma.b0',]$X5.,X95.=thisOne[thisOne$X == 'sigma.b0',]$X95.)
        
      }
      
      if(h != 2){
        summs <- rbind(mu.a,mu.b,beta,sdnoise,rho,sigma.a,sigma.b)
      }

      theSumms <- rbind(theSumms,summs)
      
    }
  }
}

theSumms <- theSumms[order(theSumms$focus,theSumms$stat,theSumms$group),]
theSumms <- unique(theSumms)
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
write.csv(table1,paste0(manuDir,'/table1 - Y05-Y15 - All Zeros.csv'))





























# make table 2 of b10s.
Rcore1st.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Core.csv'))
Rnoco1st.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Non-Core.csv'))
Rleks1st.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))

names(Rcore1st.CSV) <- c('X','core.mean','core.sd','coreX2.5','coreX5','coreX25','coreX50','coreX75','coreX95','coreX97.5')
names(Rnoco1st.CSV) <- c('X','noco.mean','noco.sd','nocoX2.5','nocoX5','nocoX25','nocoX50','nocoX75','nocoX95','nocoX97.5')
names(Rleks1st.CSV) <- c('X','leks.mean','leks.sd','leksX2.5','leksX5','leksX25','leksX50','leksX75','leksX95','leksX97.5')

Rcore1st.CSV2 <- Rcore1st.CSV[,c('X','core.mean','coreX5','coreX95')]
Rnoco1st.CSV2 <- Rnoco1st.CSV[,c('X','noco.mean','nocoX5','nocoX95')]
Rleks1st.CSV2 <- Rleks1st.CSV[,c('X','leks.mean','leksX5','leksX95')]

pre1stTable2   <- merge(Rleks1st.CSV2,Rcore1st.CSV2,by=c('X'))
pre1stTable2.B <- merge(pre1stTable2,Rnoco1st.CSV2,by=c('X'))
pre1stTable2.C <- pre1stTable2.B[1:41,]
table1st2 <- rbind(pre1stTable2.C[7:41,],pre1stTable2.C[1:6,] )
rownames(table1st2) <- NULL
write.csv(table1st2,paste0(manuDir,'/table2- 1st zeros.csv'))


RcoreAll.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Core.csv'))
RnocoAll.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Non-Core.csv'))
RleksAll.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))

names(RcoreAll.CSV) <- c('X','core.mean','core.sd','coreX2.5','coreX5','coreX25','coreX50','coreX75','coreX95','coreX97.5')
names(RnocoAll.CSV) <- c('X','noco.mean','noco.sd','nocoX2.5','nocoX5','nocoX25','nocoX50','nocoX75','nocoX95','nocoX97.5')
names(RleksAll.CSV) <- c('X','leks.mean','leks.sd','leksX2.5','leksX5','leksX25','leksX50','leksX75','leksX95','leksX97.5')

RcoreAll.CSV2 <- RcoreAll.CSV[,c('X','core.mean','coreX5','coreX95')]
RnocoAll.CSV2 <- RnocoAll.CSV[,c('X','noco.mean','nocoX5','nocoX95')]
RleksAll.CSV2 <- RleksAll.CSV[,c('X','leks.mean','leksX5','leksX95')]

preAllTable2   <- merge(RleksAll.CSV2,RcoreAll.CSV2,by=c('X'))
preAllTable2.B <- merge(preAllTable2,RnocoAll.CSV2,by=c('X'))
preAllTable2.C <- preAllTable2.B[1:41,]
tableAll2 <- rbind(preAllTable2.C[7:41,],preAllTable2.C[1:6,] )
rownames(tableAll2) <- NULL
write.csv(tableAll2,paste0(manuDir,'/table2 - All zeros.csv'))























# make plot, for each mzone, of core, non-core, all leks together - 1ST ZEROS

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
#   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
#   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
#   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
#   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
#   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
#   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  

  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
#   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
#   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
#   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
#   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
#   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
#   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
#   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
#   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
#   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
#   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
#   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
#   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
#   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
#   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
#   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
#   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
#   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
#   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)

  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - 1st Zeros - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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


# make plot, for each mzone, of core, non-core, all leks together - all ZEROS

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' All Zeros 1/bayesSummary - Model D ',units[i],' All Zeros 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' All Zeros 1/bayesSummary - Model E ',units[i],' All Zeros 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros 1/bayesSummary - Model F ',units[i],' All Zeros 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  #   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  #   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  #   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  #   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - All Zeros - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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
















# make plot, for each mzone, of core, non-core, all leks together - 1ST ZEROS - 11 years!!!!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model D ',units[i],' 1st Zeros 2005-2015.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model E ',units[i],' 1st Zeros 2005-2015.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model F ',units[i],' 1st Zeros 2005-2015.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  #   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  #   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  #   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  #   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - Y05-15 - 1st Zeros - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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


# make plot, for each mzone, of core, non-core, all leks together - all ZEROS - 11 years!!!!!!!!!!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' All Zeros 2005-2015/bayesSummary - Model D ',units[i],' All Zeros 2005-2015.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' All Zeros 2005-2015/bayesSummary - Model E ',units[i],' All Zeros 2005-2015.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros 2005-2015/bayesSummary - Model F ',units[i],' All Zeros 2005-2015.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  #   RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean  
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  #   Y.Rleks     <- exp(RleksMu.a)*(RleksBeta51^YearC)
  #   Y.Rleks.X5  <- exp(RleksMu.a.X5)*(RleksBeta51^YearC)
  #   Y.Rleks.X95 <- exp(RleksMu.a.X95)*(RleksBeta51^YearC)
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - Y05-15 - All Zeros - ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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



































# make plot, for each state, of all leks together -- 1st zeros!!!!!!!!!!!!!!!


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
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
#   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
#   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
#   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
#   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
#   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
#   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
#   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
#   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
#   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
#   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - 1st Zeros ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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







# make plot, for each state, of all leks together -- 1st zeros!!!!!!!!!!!!!!!


bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:11){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros State 1/bayesSummary - Model F ',units[i],' All Zeros State 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
  CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - All Zeros ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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




















































# make plot, for rangewide alone, that includes core, non-core, all leks -- 1st zeros!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:1){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - 1st Zeros - Rangewide.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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


# make plot, for rangewide alone, that includes core, non-core, all leks -- all zeros!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:1){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - All Zeros - Rangewide.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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












# make plot, for rangewide alone, that includes core, non-core, all leks -- 1st zeros!!!!   11 years!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:1){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - Y05-15 - 1st Zeros - Rangewide.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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


# make plot, for rangewide alone, that includes core, non-core, all leks -- all zeros!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:1){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - Core.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - Non-Core.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
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
  
  yMax <- 40
  
  png(filename=paste0(manuDir,'/Trend Plot of C, NC, AL - Y05-15 - All Zeros - Rangewide.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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

# png(filename=paste0(manuDir,'/B10 Rangewide Histograms - C, NC, AL.png'),width=12.5,height=9,units="in",res=600,pointsize=12)
# par(mfrow=c(2,2))
# min <- min(table2$core.mean)
# max <- max(table2$core.mean)
# hist(table2$core.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="Core",xlab="Estimates of Annual % Change",col="gray",las=1)
# abline(v=table2[41,]$core.mean,lwd=3)
# axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
# axis(2,at=c(0:10),labels=c(0:10),las=1)
# text(table2[41,]$core.mean-0.00005,9.9,"2005-2015",cex=0.8,pos=2)
# text(table2[41,]$core.mean-0.00005,9.5,"Annualized % Change",cex=0.8,pos=2)
# 
# min <- min(table2$noco.mean)
# max <- max(table2$noco.mean)
# hist(table2$noco.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="Periphery",xlab="Estimates of Annual % Change",col="gray",las=1)
# abline(v=table2[41,]$noco.mean,lwd=3)
# axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
# axis(2,at=c(0:10),labels=c(0:10),las=1)
# text(table2[41,]$noco.mean-0.0005,9.9,"2005-2015",cex=0.8,pos=2)
# text(table2[41,]$noco.mean-0.0005,9.5,"Annualized % Change",cex=0.8,pos=2)
# 
# min <- min(table2$leks.mean)
# max <- max(table2$leks.mean)
# hist(table2$leks.mean,breaks=seq(min,max,(max-min)/9),xaxt="n",yaxt="n",ylim=c(0,10),cex=0.8,main="All Leks Combined",xlab="Estimates of Annual % Change",col="gray",las=1)
# abline(v=table2[41,]$leks.mean,lwd=3)
# # ylab <- ifelse(seq(min,max,(max-min)/9) >= 1,seq(min,max,(max-min)/9)-1,seq(min,max,(max-min)/9)-1)
# axis(1,at=round(seq(min,max,(max-min)/9),4),labels=100*round(seq(min,max,(max-min)/9)-1,4))
# axis(2,at=c(0:10),labels=c(0:10),las=1)
# text(table2[41,]$leks.mean-0.0001,9.9,"2005-2015",cex=0.8,pos=2)
# text(table2[41,]$leks.mean-0.0001,9.5,"Annualized % Change",cex=0.8,pos=2)
# 
# plot.new()
# dev.off()
# 
# par(mfrow=c(1,1))



print(formatC(signif(round(seq(min,max,(max-min)/9),3))),digits=3)






# make the B10 plots -- old with all B10s on one graph
# Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
# Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
# Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
# 
# Year <- seq(1965,2015,1)
# YearC <- Year - 1964 - 26
# 
# RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
# RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
# RleksBeta51 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$mean
# 
# RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
# RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
# RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean
# 
# RYcore <- exp(RcoreMu.a)*RcoreBeta51^YearC
# RYnoco <- exp(RnocoMu.a)*RnocoBeta51^YearC
# RYleks <- exp(RleksMu.a)*RleksBeta51^YearC
# 
# yMin <- 0
# yMax <- 20
# 
# 
# # all leks
# png(filename=paste0(manuDir,'/B10 Plot of All Leks - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
# plot(Year,RYleks,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# 
# for(i in 1:41){
#   x <- Year[i:(i + 10)]                      # restrict to years of interest
#   B10 <- Rleks.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
#   mid <- median(x)                           # calculate temporal median of this B10 stretch
#   B10.mu.a <- RYleks[mid - 1964]             # this is the intercept to use for this 10-year stretch
#   B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
#    
#   par(new=TRUE)
#   plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# }
# 
# axis(side=1,labels=TRUE,seq(1965,2015,5))
# axis(side=2,labels=TRUE,seq(yMin,yMax,1))
# dev.off()
# 
# 
# # core
# png(filename=paste0(manuDir,'/B10 Plot of Core - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
# plot(Year,RYcore,type='l',col='darkgreen' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# 
# for(i in 1:41){
#   x <- Year[i:(i + 10)]                      # restrict to years of interest
#   B10 <- Rcore.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
#   mid <- median(x)                           # calculate temporal median of this B10 stretch
#   B10.mu.a <- RYcore[mid - 1964]             # this is the intercept to use for this 10-year stretch
#   B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
#   
#   par(new=TRUE)
#   plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# }
# 
# axis(side=1,labels=TRUE,seq(1965,2015,5))
# axis(side=2,labels=TRUE,seq(yMin,yMax,1))
# dev.off()
# 
# 
# # non-core
# png(filename=paste0(manuDir,'/B10 Plot of Non-Core - Rangewide.png'),width=8,height=6,units="in",res=600,pointsize=12)
# plot(Year,RYnoco,type='l',col='darkred' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# 
# for(i in 1:41){
#   x <- Year[i:(i + 10)]                      # restrict to years of interest
#   B10 <- Rnoco.CSV[100 - i,2]                # this is the slope term for this 10-year stretch
#   mid <- median(x)                           # calculate temporal median of this B10 stretch
#   B10.mu.a <- RYnoco[mid - 1964]             # this is the intercept to use for this 10-year stretch
#   B10.trend <- B10.mu.a*B10^seq(-5,5,1)      # estimate N sage grouse based on B10 parameters
#   
#   par(new=TRUE)
#   plot(x,B10.trend,axes=FALSE,type='l',frame.plot=TRUE,xlim=c(1965,2015),ylim=c(yMin,yMax),xlab='',ylab='',lwd=2)
# }
# 
# axis(side=1,labels=TRUE,seq(1965,2015,5))
# axis(side=2,labels=TRUE,seq(yMin,yMax,1))
# dev.off()













# make plot of all mzone-core, mzone-non-core, mzone-all leks, states-all leks















# make plot of all states - all leks together -- 1st zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - 1st Zeros - All States.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:11){
  if(i != 8){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
    
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


# make plot of all states - all leks together -- all zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - All Zeros - All States.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:11){
  if(i != 8){
    units <- paste0("MZone ",states)
    nUnits <- nStates
    #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
    #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
    bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros State 1/bayesSummary - Model F ',units[i],' All Zeros State 1.csv'))
    #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
    #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
    Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
    
    #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
    #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
    leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
    #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
    #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
    RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
    
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
    Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
    Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
    Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
    #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
    #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
    #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
    #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
    #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
    #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
    Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
    Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
    Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
    
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
    
    yMax <- 40
    
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

















# NEED TO RUN 11-YEAR MODELS FOR STATES IF WE WANT THESE.
# make plot of all states - all leks together -- 1st zeros --- 11 years!!!!
# 
# bsums.core.CSV <- vector("list",6)
# bsums.noco.CSV <- vector("list",6)
# bsums.leks.CSV <- vector("list",6)
# CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - Y05-15 - 1st Zeros - All States.png'),width=8,height=6,units="in",res=300,pointsize=12)
# 
# for(i in 1:11){
#   if(i != 8){
#     units <- paste0("MZone ",states)
#     nUnits <- nStates
#     #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
#     #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
#     bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' 1st Zeros 2005-2015/bayesSummary - Y05-15 - Model F ',units[i],' 1st Zeros 2005-2015.csv'))
#     #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
#     #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
#     Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - All Leks.csv'))
#     
#     #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
#     #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
#     leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
#     #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
#     #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
#     RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
#     
#     #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
#     #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
#     #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
#     #   
#     #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
#     #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
#     #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
#     
#     #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
#     #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
#     leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
#     #   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
#     #   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
#     RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
#     
#     #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
#     #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
#     leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
#     #   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
#     #   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
#     RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
#     
#     #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
#     #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
#     leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
#     #   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
#     #   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
#     RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
#     
#     
#     
#     Year <- data.frame(Year=seq(2005,2015,1))
#     YearC <- Year - 2004 - 6
#     
#     #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
#     #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
#     #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
#     #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
#     #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
#     #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
#     Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
#     Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
#     Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
#     #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
#     #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
#     #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
#     #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
#     #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
#     #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
#     Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
#     Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
#     Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
#     
#     plotYears <- data.frame(Year=Year,YearC=YearC, 
#                             #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
#                             #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
#                             Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
#                             #                           Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
#                             #                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
#                             Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
#     names(plotYears) <- c('Year','YearC',
#                           #                         'Y.core','Y.core.X5','Y.core.X95',
#                           #                         'Y.noco','Y.noco.X5','Y.noco.X95',
#                           'Y.leks','Y.leks.X5','Y.leks.X95',
#                           #                         'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
#                           #                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
#                           'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
#     
#     
#     x <- plotYears$Year
#     #   yA1 <- plotYears$Y.core
#     #   yA2 <- plotYears$Y.core.X5
#     #   yA3 <- plotYears$Y.core.X95
#     #   yB1 <- plotYears$Y.noco
#     #   yB2 <- plotYears$Y.noco.X5
#     #   yB3 <- plotYears$Y.noco.X95
#     yC1 <- plotYears$Y.leks
#     yC2 <- plotYears$Y.leks.X5
#     yC3 <- plotYears$Y.leks.X95
#     
#     #   yRA1 <- plotYears$Y.Rcore
#     #   yRA2 <- plotYears$Y.Rcore.X5
#     #   yRA3 <- plotYears$Y.Rcore.X95
#     #   yRB1 <- plotYears$Y.Rnoco
#     #   yRB2 <- plotYears$Y.Rnoco.X5
#     #   yRB3 <- plotYears$Y.Rnoco.X95
#     yRC1 <- plotYears$Y.Rleks
#     yRC2 <- plotYears$Y.Rleks.X5
#     yRC3 <- plotYears$Y.Rleks.X95
#     
#     yMax <- 40
#     
#     # state specific
#     #   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
#     #   par(new=TRUE)
#     #   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE)
#     #   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   
#     #   par(new=TRUE)
#     #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     if(i > 1){
#       par(new=TRUE)
#     }
#     if(i == 1){
#       plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0("All States - Combined Leks"))    
#     } else {
#       plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)     
#     }
#     if(i == 11){
#       par(new=TRUE)
#       plot(x,yRC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=5)
#     }
#     #   par(new=TRUE)
#     #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     # rangewide specific
#     #   par(new=TRUE)
#     #   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE)
#     #   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   
#     #   par(new=TRUE)
#     #   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     
#     #   par(new=TRUE)
#     #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   }
# }
# legend("topright",c(states[-8],"Rangewide"),lty=c(rep(1,10),1),lwd=c(rep(2,10),5),col=c(colVec[-8],'black'),bty="n",ncol=2)
# axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
# axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
# mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
# mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
# dev.off()
# 
# 
# # make plot of all states - all leks together -- all zeros -- 11 years!!!!!!!!
# 
# bsums.core.CSV <- vector("list",6)
# bsums.noco.CSV <- vector("list",6)
# bsums.leks.CSV <- vector("list",6)
# CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - Y05-15 - All Zeros - All States.png'),width=8,height=6,units="in",res=300,pointsize=12)
# 
# for(i in 1:11){
#   if(i != 8){
#     units <- paste0("MZone ",states)
#     nUnits <- nStates
#     #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
#     #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
#     bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros State 1/bayesSummary - Model F ',units[i],' All Zeros State 1.csv'))
#     #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
#     #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
#     Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
#     
#     #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
#     #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
#     leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
#     #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
#     #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
#     RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
#     
#     #   coreBeta51.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   nocoBeta51.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   leksBeta51.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X5.
#     #   RcoreBeta51.X5 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X5.
#     #   RnocoBeta51.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X5.
#     #   RleksBeta51.X5 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X5. 
#     #   
#     #   coreBeta51.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   nocoBeta51.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   leksBeta51.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$X95.
#     #   RcoreBeta51.X95 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$X95.
#     #   RnocoBeta51.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$X95.
#     #   RleksBeta51.X95 <- Rleks.CSV[Rleks.CSV$X == 'beta0',]$X95. 
#     
#     #   coreMu.a <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$mean
#     #   nocoMu.a <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$mean
#     leksMu.a <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$mean
#     #   RcoreMu.a <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$mean
#     #   RnocoMu.a <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$mean
#     RleksMu.a <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$mean 
#     
#     #   coreMu.a.X5 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X5.
#     #   nocoMu.a.X5 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X5.
#     leksMu.a.X5 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X5.
#     #   RcoreMu.a.X5 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X5.
#     #   RnocoMu.a.X5 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X5.
#     RleksMu.a.X5 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X5. 
#     
#     #   coreMu.a.X95 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.a',]$X95.
#     #   nocoMu.a.X95 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.a',]$X95.
#     leksMu.a.X95 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.a',]$X95.
#     #   RcoreMu.a.X95 <- Rcore.CSV[Rcore.CSV$X == 'mu.a0',]$X95.
#     #   RnocoMu.a.X95 <- Rnoco.CSV[Rnoco.CSV$X == 'mu.a0',]$X95.
#     RleksMu.a.X95 <- Rleks.CSV[Rleks.CSV$X == 'mu.a0',]$X95.
#     
#     
#     
#     Year <- data.frame(Year=seq(2005,2015,1))
#     YearC <- Year - 2004 - 6
#     
#     #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
#     #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
#     #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
#     #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
#     #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
#     #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
#     Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
#     Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
#     Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
#     #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
#     #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
#     #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
#     #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
#     #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
#     #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
#     Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
#     Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
#     Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
#     
#     plotYears <- data.frame(Year=Year,YearC=YearC, 
#                             #                           Y.core=Y.core , Y.core.X5=Y.core.X5 , Y.core.X95=Y.core.X95,
#                             #                           Y.noco=Y.noco , Y.noco.X5=Y.noco.X5 , Y.noco.X95=Y.noco.X95,
#                             Y.leks=Y.leks , Y.leks.X5=Y.leks.X5 , Y.leks.X95=Y.leks.X95,
#                             #                           Y.Rcore=Y.Rcore,Y.Rcore.X5=Y.Rcore.X5,Y.Rcore.X95=Y.Rcore.X95,
#                             #                           Y.Rnoco=Y.Rnoco,Y.Rnoco.X5=Y.Rnoco.X5,Y.Rnoco.X95=Y.Rnoco.X95,
#                             Y.Rleks=Y.Rleks,Y.Rleks.X5=Y.Rleks.X5,Y.Rleks.X95=Y.Rleks.X95)
#     names(plotYears) <- c('Year','YearC',
#                           #                         'Y.core','Y.core.X5','Y.core.X95',
#                           #                         'Y.noco','Y.noco.X5','Y.noco.X95',
#                           'Y.leks','Y.leks.X5','Y.leks.X95',
#                           #                         'Y.Rcore','Y.Rcore.X5','Y.Rcore.X95',
#                           #                         'Y.Rnoco','Y.Rnoco.X5','Y.Rnoco.X95',
#                           'Y.Rleks','Y.Rleks.X5','Y.Rleks.X95')
#     
#     
#     x <- plotYears$Year
#     #   yA1 <- plotYears$Y.core
#     #   yA2 <- plotYears$Y.core.X5
#     #   yA3 <- plotYears$Y.core.X95
#     #   yB1 <- plotYears$Y.noco
#     #   yB2 <- plotYears$Y.noco.X5
#     #   yB3 <- plotYears$Y.noco.X95
#     yC1 <- plotYears$Y.leks
#     yC2 <- plotYears$Y.leks.X5
#     yC3 <- plotYears$Y.leks.X95
#     
#     #   yRA1 <- plotYears$Y.Rcore
#     #   yRA2 <- plotYears$Y.Rcore.X5
#     #   yRA3 <- plotYears$Y.Rcore.X95
#     #   yRB1 <- plotYears$Y.Rnoco
#     #   yRB2 <- plotYears$Y.Rnoco.X5
#     #   yRB3 <- plotYears$Y.Rnoco.X95
#     yRC1 <- plotYears$Y.Rleks
#     yRC2 <- plotYears$Y.Rleks.X5
#     yRC3 <- plotYears$Y.Rleks.X95
#     
#     yMax <- 40
#     
#     # state specific
#     #   plot(x,yA1,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2015\n75% ",units[i]))
#     #   par(new=TRUE)
#     #   plot(x,yA2,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE)
#     #   plot(x,yA3,type='l',col='lightgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   
#     #   par(new=TRUE)
#     #   plot(x,yB1,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yB2,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yB3,type='l',col='pink'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     if(i > 1){
#       par(new=TRUE)
#     }
#     if(i == 1){
#       plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2,main=paste0("All States - Combined Leks"))    
#     } else {
#       plot(x,yC1,type='l',col=colVec[i] ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=2)     
#     }
#     if(i == 11){
#       par(new=TRUE)
#       plot(x,yRC1,type='l',col='black' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=1,lwd=5)
#     }
#     #   par(new=TRUE)
#     #   plot(x,yC2,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yC3,type='l',col='lightblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     # rangewide specific
#     #   par(new=TRUE)
#     #   plot(x,yRA1,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yRA2,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE)
#     #   plot(x,yRA3,type='l',col='darkgreen',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   
#     #   par(new=TRUE)
#     #   plot(x,yRB1,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     #   par(new=TRUE)
#     #   plot(x,yRB2,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yRB3,type='l',col='darkred'  ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     
#     
#     #   par(new=TRUE)
#     #   plot(x,yRC2,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     #   par(new=TRUE) 
#     #   plot(x,yRC3,type='l',col='darkblue' ,axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#   }
# }
# legend("topright",c(states[-8],"Rangewide"),lty=c(rep(1,10),1),lwd=c(rep(2,10),5),col=c(colVec[-8],'black'),bty="n",ncol=2)
# axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
# axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
# mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
# mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
# dev.off()















# make plot of all mzones - all leks together - 1ST ZEROS

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
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

# make plot of all mzones - all leks together - all ZEROS

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros 1/bayesSummary - Model F ',units[i],' All Zeros 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
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




# make plot of all mzones - all leks together - 1ST ZEROS - 11 years!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - Y05-15 - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model F ',units[i],' 1st Zeros 2005-2015.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
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

# make plot of all mzones - all leks together - all ZEROS - 11 years!!!!!

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of AL - Y05-15 - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros 2005-2015/bayesSummary - Model F ',units[i],' All Zeros 2005-2015.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 40
  
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



































# make plot of all mzones - core together -- 1st Zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of C - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
    bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
    Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Core.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
    coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
    RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
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
  
    Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
    Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
    Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
    Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
    Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
    Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
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
  
  yMax <- 40
  
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


# make plot of all mzones - core together -- all Zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of C - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' All Zeros 1/bayesSummary - Model D ',units[i],' All Zeros 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Core.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
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
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
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
  
  yMax <- 40
  
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


# make plot of all mzones - core together -- 1st Zeros - 11 years

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of C - Y05-15 - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model D ',units[i],' 1st Zeros 2005-2015.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Core.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
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
  
  yMax <- 40
  
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


# make plot of all mzones - core together -- all Zeros -- 11 years

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of C - Y05-15 - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' All Zeros 2005-2015/bayesSummary - Model D ',units[i],' All Zeros 2005-2015.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - Core.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  coreMuB <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'mu.b',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  RcoreMuB <- Rcore.CSV[Rcore.CSV$X == 'mu.b0',]$mean
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  Y.core      <- exp(coreMu.a)*((coreMuB+1)^YearC)
  Y.core.X5   <- exp(coreMu.a.X5)*((coreMuB+1)^YearC)
  Y.core.X95  <- exp(coreMu.a.X95)*((coreMuB+1)^YearC)
  #   Y.noco      <- exp(nocoMu.a)*(nocoBeta51^YearC)
  #   Y.noco.X5   <- exp(nocoMu.a.X5)*(nocoBeta51^YearC)
  #   Y.noco.X95  <- exp(nocoMu.a.X95)*(nocoBeta51^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  Y.Rcore     <- exp(RcoreMu.a)*((RcoreMuB+1)^YearC)
  Y.Rcore.X5  <- exp(RcoreMu.a.X5)*((RcoreMuB+1)^YearC)
  Y.Rcore.X95 <- exp(RcoreMu.a.X95)*((RcoreMuB+1)^YearC)
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
  
  yMax <- 40
  
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



































# make plot of all mzones - non-core together - 1st zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of NC - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - Non-Core.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
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
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
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
  
  yMax <- 40
  
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



# make plot of all mzones - non-core together - all zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of NC - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' All Zeros 1/bayesSummary - Model E ',units[i],' All Zeros 1.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - Non-Core.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
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
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
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
  
  yMax <- 40
  
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







# make plot of all mzones - non-core together - 1st zeros - 11 years

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of NC - Y05-15 - 1st Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' 1st Zeros 2005-2015/bayesSummary - Model E ',units[i],' 1st Zeros 2005-2015.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - 1st Zeros - Non-Core.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
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
  
  yMax <- 40
  
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



# make plot of all mzones - non-core together - all zeros

bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
CairoPNG(filename=paste0(manuDir,'/Trend Plot of NC - Y05-15 - All Zeros - All MZones.png'),width=8,height=6,units="in",res=300,pointsize=12)

for(i in 1:6){
  units <- mZones
  nUnits <- nMZones
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' All Zeros 2005-2015/bayesSummary - Model E ',units[i],' All Zeros 2005-2015.csv'))
  #   bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' Try 1/bayesSummary - Model F ',units[i],' Try 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Y05-15 - All Zeros - Non-Core.csv'))
  #   Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model F.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  nocoMuB <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'mu.b',]$mean
  #   leksBeta51 <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'beta[1]',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  RnocoMuB <- Rnoco.CSV[Rnoco.CSV$X == 'mu.b0',]$mean
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
  
  
  
  Year <- data.frame(Year=seq(2005,2015,1))
  YearC <- Year - 2004 - 6
  
  #   Y.core      <- exp(coreMu.a)*(coreBeta51^YearC)
  #   Y.core.X5   <- exp(coreMu.a.X5)*(coreBeta51^YearC)
  #   Y.core.X95  <- exp(coreMu.a.X95)*(coreBeta51^YearC)
  Y.noco      <- exp(nocoMu.a)*((nocoMuB+1)^YearC)
  Y.noco.X5   <- exp(nocoMu.a.X5)*((nocoMuB+1)^YearC)
  Y.noco.X95  <- exp(nocoMu.a.X95)*((nocoMuB+1)^YearC)
  #   Y.leks      <- exp(leksMu.a)*(leksBeta51^YearC)
  #   Y.leks.X5   <- exp(leksMu.a.X5)*(leksBeta51^YearC)
  #   Y.leks.X95  <- exp(leksMu.a.X95)*(leksBeta51^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  Y.Rnoco     <- exp(RnocoMu.a)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*((RnocoMuB+1)^YearC)
  Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*((RnocoMuB+1)^YearC)
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
  
  yMax <- 40
  
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



















































# make simple trends, with the 'jagged lines' of observed time series.

# over states!!!!

# make plot, for each state, of all leks together -- 1st zeros!!!!!!!!!!!!!!!


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
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - 1st Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 60
  
  df <- dat1stZerosLeks[[9]][dat1stZerosLeks[[9]]$State == states[i],]
  tsdf <- aggregate(data=df, Peak_Males ~ Year, mean)
  
  CairoPNG(filename=paste0(manuDir,'/Jagged Line Plots/Trend Plot of AL - 1st Zeros ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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
  
  par(new=TRUE)
  plot(tsdf$Year,tsdf$Peak_Males,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=3,lwd=0.5)
  legend("topright",c("Observed","Combined"),lty=c(3,1),lwd=c(0.5,2),col=c('black','black'),bty="n")
  axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
  axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
  mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
  mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
  dev.off()
  
}







# make plot, for each state, of all leks together -- all zeros!!!!!!!!!!!!!!!


bsums.core.CSV <- vector("list",6)
bsums.noco.CSV <- vector("list",6)
bsums.leks.CSV <- vector("list",6)
for(i in 1:11){
  units <- paste0("MZone ",states)
  nUnits <- nStates
  #   bsums.core.CSV[[i]] <- read.csv(paste0(outpDir,'/Model D ',units[i],' Try 1/bayesSummary - Model D ',units[i],' Try 1.csv'))
  #   bsums.noco.CSV[[i]] <- read.csv(paste0(outpDir,'/Model E ',units[i],' Try 1/bayesSummary - Model E ',units[i],' Try 1.csv'))
  bsums.leks.CSV[[i]] <- read.csv(paste0(outpDir,'/Model F ',units[i],' All Zeros State 1/bayesSummary - Model F ',units[i],' All Zeros State 1.csv'))
  #   Rcore.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model D.csv'))
  #   Rnoco.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - Model E.csv'))
  Rleks.CSV      <- read.csv(paste0(outpDir,'/Rangewide/Trend Plots/bayesSummary - All Zeros - All Leks.csv'))
  
  #   coreBeta51 <- bsums.core.CSV[[i]][bsums.core.CSV[[i]]$X == 'beta[1]',]$mean
  #   nocoBeta51 <- bsums.noco.CSV[[i]][bsums.noco.CSV[[i]]$X == 'beta[1]',]$mean
  leksMuB <- bsums.leks.CSV[[i]][bsums.leks.CSV[[i]]$X == 'mu.b',]$mean
  #   RcoreBeta51 <- Rcore.CSV[Rcore.CSV$X == 'beta0',]$mean
  #   RnocoBeta51 <- Rnoco.CSV[Rnoco.CSV$X == 'beta0',]$mean
  RleksMuB <- Rleks.CSV[Rleks.CSV$X == 'mu.b0',]$mean  
  
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
  Y.leks      <- exp(leksMu.a)*((leksMuB+1)^YearC)
  Y.leks.X5   <- exp(leksMu.a.X5)*((leksMuB+1)^YearC)
  Y.leks.X95  <- exp(leksMu.a.X95)*((leksMuB+1)^YearC)
  #   Y.Rcore     <- exp(RcoreMu.a)*(RcoreBeta51^YearC)
  #   Y.Rcore.X5  <- exp(RcoreMu.a.X5)*(RcoreBeta51^YearC)
  #   Y.Rcore.X95 <- exp(RcoreMu.a.X95)*(RcoreBeta51^YearC)
  #   Y.Rnoco     <- exp(RnocoMu.a)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X5  <- exp(RnocoMu.a.X5)*(RnocoBeta51^YearC)
  #   Y.Rnoco.X95 <- exp(RnocoMu.a.X95)*(RnocoBeta51^YearC)
  Y.Rleks     <- exp(RleksMu.a)*((RleksMuB+1)^YearC)
  Y.Rleks.X5  <- exp(RleksMu.a.X5)*((RleksMuB+1)^YearC)
  Y.Rleks.X95 <- exp(RleksMu.a.X95)*((RleksMuB+1)^YearC)
  
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
  
  yMax <- 60
  
  df <- datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == states[i],]
  tsdf <- aggregate(data=df, Peak_Males ~ Year, mean)
  
  CairoPNG(filename=paste0(manuDir,'/Jagged Line Plots/Trend Plot of AL - All Zeros ',units[i],'.png'),width=8,height=6,units="in",res=300,pointsize=12)
  
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
  par(new=TRUE)
  plot(tsdf$Year,tsdf$Peak_Males,type='l',col='black',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lty=3,lwd=0.5)
  legend("topright",c("Observed","Combined"),lty=c(3,1),lwd=c(0.5,2),col=c('black','black'),bty="n")
  axis(side=1,labels=TRUE,seq(1965,2015,5),cex=1.5)
  axis(side=2,labels=TRUE,seq(0,yMax,10),cex=1.5,las=1)
  mtext("Year", side=1, line=3, cex.lab=1.5,las=1, col="black")   # x-axis
  mtext("Peak Males / Lek Base Count", side=2, line=3, cex.lab=1.5, col="black")   # y-axis
  dev.off()
  
}
















