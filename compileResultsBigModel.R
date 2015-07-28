

modelLetter <- "P"

theFiles <- list.files(outpDir)[substr(list.files(outpDir),1,6) == 'Model ']


assign(paste0("ests",modelLetter),summarizeModels(modelLetter=modelLetter,theFiles))


# allResults <- rbind(estsC,estsD,estsF)

allResults <- estsP
rm(estsP)
names(allResults)[names(allResults) == 'Zone'] <- 'mZone'
allResults <- data.frame(allResults,ParmType=unlist(lapply(strsplit(as.character(droplevels(allResults$Parameter)),'[',fixed=TRUE), function(x) x[[1]][1])))

allResults <- allResults[order(allResults$mZone,allResults$modCode),]

pi0All <- allResults[allResults$Parameter == 'pi0',]
pi1All <- allResults[allResults$Parameter == 'pi1',]
beta00All <- allResults[allResults$Parameter == 'beta00',]
beta10All <- allResults[allResults$Parameter == 'beta10',]
gamma000All <- allResults[allResults$Parameter == 'gamma000',]
gamma100All <- allResults[allResults$Parameter == 'gamma100',]
taonoiseAll <- allResults[allResults$Parameter == 'taonoise',]

SumData <- read.csv(paste0(analDir,'/SumData.csv'))

yMax <- max(datList[[1]]$Peak_Males)

helper <- allResults[,c('Model','Zeros','Cut','mZone','modCode')]
helper <- unique(helper)
rownames(helper) <- NULL





doThese <- c('N')
nDoThese <- length(doThese)
table1 <- NULL
colVec <- brewer.pal(8,"Set1")

for(i in 9:9){
  
  #A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
  #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
  #C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
  
  #D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
  #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
  #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
  
  N <- readOGR(analDir,'Core1stZeros15p')@data        # N - read in 1st zeros, core data, all zones, 15% sample  
  
  for(j in 1:nDoThese){
    
    thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]
    
    if(doThese[j] %in% c('D','E','F')){
      thisOne <- get(doThese[j])
      #thisOne <- thisOne[thisOne$DupZero == 0,]  # dont think this is necessary when doing core
    } else {
      thisOne <- get(doThese[j])   
    }
    leks <- unique(thisOne$Lek_IDnum)
    mzones <- unique(thisOne$Mgmt_zone)
    thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
    thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID)) 
    thisOne$mZoneCls <- as.numeric(as.factor(thisOne$Mgmt_zone))
    nLeks <- length(leks)
    nmZones <- length(mzones)
  
    Model <- as.character(droplevels(thisHelper$Model))
    Zeros <- as.character(droplevels(thisHelper$Zeros))
    Cut   <- as.character(droplevels(thisHelper$Cut))
    mZone <- i
    
    # make individual lek trends
    for(k in 1:nLeks){    
      # do this for winbugs numbering (so on factor)
      xlek <- thisOne[thisOne$lekCls == k,]$Year - 1990
      ylek <- thisOne[thisOne$lekCls == k,]$Peak_Males
      the.lek  <- thisOne[thisOne$lekCls == k,]$Lek_ID[1]
      lek.pi0  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('pi0[',k,']'),]$mean
      lek.pi1  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('pi1[',k,']'),]$mean
      taunoise <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == 'taunoise',]$mean
      ypred1 <- exp(lek.pi0 + lek.pi1*xlek)
      ypred2 <- exp(lek.pi0 + lek.pi1*xlek + 0.5*taunoise)
      
      thisOne$origMZone <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(thisOne$Mgmt_zone)),"MZ ",fixed=TRUE))[c(FALSE,TRUE)]))
      origMZone <- thisOne[thisOne$lekCls == k,]$origMZone[1]
      
      png(filename=paste0(rsltDir,'/Lek Trends/',the.lek,' - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',origMZone,'.png'),width=10,height=6,units="in",res=500,pointsize=12)
      
      plot(xlek + 1990,ylek,col='gray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5,pch=19)
      par(new=TRUE)
      plot(xlek + 1990,ypred1,col=colVec[origMZone],axes=FALSE,lty=1,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Lek ",the.lek," in Management Zone ",origMZone,"\nYears 1965-2014, ",Model,', ',Zeros,', ',Cut))
      par(new=TRUE)
      plot(xlek + 1990,ypred2,col=colVec[origMZone],axes=FALSE,lty=3,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Lek ",the.lek," in Management Zone ",origMZone,"\nYears 1965-2014, ",Model,', ',Zeros,', ',Cut))
      axis(side=1,labels=TRUE,seq(1965,2015,5))
      axis(side=2,labels=TRUE,seq(0,yMax,10))
      
      dev.off()
    }
    
    # make individual mzone trends
    for(k in 1:nmZones){    
      # do this for winbugs numbering (so on factor)
      
      thisOnek <- thisOne[thisOne$mZoneCls == k,]
      thisOnek <- thisOnek[order(thisOnek$Year),]
      
      xlek <- unique(thisOnek$Year) - 1990
      obsMeans <- data.frame(MeanPMales=tapply(thisOnek$Peak_Males, list(thisOnek$Year), mean))
      ylek <- obsMeans[,1]
      #the.lek  <- thisOne[thisOne$lekCls == k,]$Lek_ID[1]
      mzone.beta0 <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('beta00[',k,']'),]$mean
      mzone.beta1 <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == paste0('beta10[',k,']'),]$mean
      taunoise <- allResults[allResults$mZone == i & allResults$modCode == doThese[j] & allResults$Parameter == 'taunoise',]$mean
      ypred1 <- exp(mzone.beta0 + mzone.beta1*xlek)
      ypred2 <- exp(mzone.beta0 + mzone.beta1*xlek + 0.5*taunoise)
      
      thisOne$origMZone <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(thisOne$Mgmt_zone)),"MZ ",fixed=TRUE))[c(FALSE,TRUE)]))
      origMZone <- thisOne[thisOne$mZoneCls == k,]$origMZone[1]
      
      png(filename=paste0(rsltDir,'/MZone Trends/',Model,' - ',Zeros,' - ',Cut,' - Zone ',origMZone,'.png'),width=10,height=6,units="in",res=500,pointsize=12)
      
      plot(xlek + 1990,ylek,col='gray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5,pch=19)
      par(new=TRUE)
      plot(xlek + 1990,ypred1,col=colVec[origMZone],axes=FALSE,lty=1,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Management Zone ",origMZone,"\nYears 1965-2014, ",Model,', ',Zeros,', ',Cut))
      par(new=TRUE)
      plot(xlek + 1990,ypred2,col=colVec[origMZone],axes=FALSE,lty=3,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,type='l',main=paste0("Temporal Trend of Management Zone ",origMZone,"\nYears 1965-2014, ",Model,', ',Zeros,', ',Cut))
      axis(side=1,labels=TRUE,seq(1965,2015,5))
      axis(side=2,labels=TRUE,seq(0,yMax,10))
      
      dev.off()
    }   
    
    
    
    # all trend
    xlek <- unique(thisOne$Year) - 1990
    obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
    ylek <- obsMeans[,1]
    
    plot(xlek,ylek)
    
    
    
    
    
    
    
    
    
    
   CairoPNG(filename=paste0(rsltDir,'/Trends + Yes Tao - No Leks - ',Model,' - ',Zeros,' - ',Cut,' - Zone ',mZone,'.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)

    # make gray individual lek trends
    for(k in 1:nLeks){
      
      # do this on k numbering
      x <- thisOne[thisOne$Lek_IDnum == leks[k],]$Year
      y <- thisOne[thisOne$Lek_IDnum == leks[k],]$Peak_Males
            
      yMax <- 100#max(thisOne$Peak_Males)
      
#       if(k == 1){
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
#       } else if(k > 1 & k < nLeks){
#         par(new=TRUE)
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=0.5)
#       } else if(k == nLeks){
#         par(new=TRUE)
#         plot(x,y,type='l',col='lightgray',axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='Year',ylab='Peak Males',lwd=0.5)        
#       }
    }
     
#     get mean trend info for this mzone and model
    mu.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$mean
    lo.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X5.
    hi.a    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.a',]$X95.
    mu.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$mean
    lo.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X5.
    hi.b    <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'mu.b',]$X95. 
    tao     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$mean
    lo.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X5.
    hi.tao  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'taunoise',]$X95.
    rho     <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$mean
    lo.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X5.
    hi.rho  <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'rho',]$X95.
    sigma.a <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.a',]$mean
    sigma.b <- allResults[allResults$mZone == i & allResults$modCode == doThese[j],][allResults[allResults$mZone == i & allResults$modCode == doThese[j],]$Parameter == 'sigma.b',]$mean

    thisRow <- data.frame(thisHelper,mu.a=mu.a,lo.a,hi.a,mu.b,lo.b,hi.b,tao,lo.tao,hi.tao,rho,lo.rho,hi.rho,sigma.a,sigma.b)
    table1 <- rbind(table1,thisRow)

    obsMeans <- data.frame(MeanPMales=tapply(thisOne$Peak_Males, list(thisOne$Year), mean))
    obsMeans$Year <- rownames(obsMeans)

    year <- seq(-25,25,1)
    est <- exp(mu.a + mu.b*year + 0.5*tao)
    se <- sqrt(sigma.a^2 + sigma.b^2 + 2*rho*sigma.a*sigma.b)
    lo <- exp(mu.a + mu.b*year + 0.5*tao) - 1.645*se 
    hi <- exp(mu.a + mu.b*year + 0.5*tao) + 1.645*se 
    x <- year + 1990
    
    #par(new=TRUE)
    plot(obsMeans$Year,obsMeans$MeanPMales,type='p',pch=19,col='darkgray',axes=FALSE,frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\n75% Core - Management Zone ",i))
    par(new=TRUE)
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
    par(new=TRUE)
    plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
    par(new=TRUE)
    plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)    
    axis(side=1,labels=TRUE,seq(1965,2015,5))
    axis(side=2,labels=TRUE,seq(0,yMax,10))
    legend("topleft",bty = "n",legend=c("Annual Mean","Temporal Trend & 90% Credible Intervals","Individual Lek Counts"),pch=c(19,NA,NA),lwd=c(NA,2,1),col=c('darkgray',colVec[i],"gray95"))
    dev.off()
  }
}

write.csv(table1,paste0(rsltDir,'/table1 + Yes Tao.csv'))

table1 <- read.csv(paste0(rsltDir,'/table1 + No Tao.csv'))


CairoPNG(filename=paste0(rsltDir,'/Trends + No Tao - All MZones - 1st Zeros - Core Leks.png'),width=10,height=6,units="in",res=500,quality=600,pointsize=12)

for(i in 1:8){
  thisRow <- table1[i,]
  year <- seq(-25,25,1)
  est <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao)
  se <- sqrt(thisRow$sigma.a^2 + thisRow$sigma.b^2 + 2*thisRow$rho*thisRow$sigma.a*thisRow$sigma.b)
  lo <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) - 1.645*se
  hi <- exp(thisRow$mu.a + thisRow$mu.b*year + 0.5*thisRow$tao) + 1.645*se
  x <- year + 1990
  if(i <= 7){
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2)
#     par(new=TRUE)
#     plot(x,lo,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)
#     par(new=TRUE)
#     plot(x,hi,type='l',col=colVec[i],axes=FALSE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=1)   
    par(new=TRUE)
  } else if(i == 8){
    plot(x,est,type='l',col=colVec[i],axes=FALSE, frame.plot=TRUE,xlim=c(1965,2015),ylim=c(0,yMax),xlab='',ylab='',lwd=2,main=paste0("Temporal Trend of Observed Peak Males, Years 1965-2014\nAll Management Zones"))
    axis(side=1,labels=TRUE,seq(1965,2015,5))
    axis(side=2,labels=TRUE,seq(0,yMax,10))
    legend("topright",bty="n",legend=c('I','II','III','IV','V','VI','VII','VIII'),col=colVec,lwd=rep(2,8))
  }

}
dev.off()





doThese <- c('A','C')
nDoThese <- length(doThese)
the.zeros <- NULL


# read in all the data
for(i in 9:9){
  
   A <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))@data        # A - read in all zeros, core data, ith mzone
  #B <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))@data    # B - read in all zeros, non-core data, ith mzone
   C <- readOGR(analDir,paste0('Zone ',i,' Both - All Zero'))@data           # C - read in all zeros, all data, ith mzone
  
  #D <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # D - read in 1st zeros, core data, ith mzone
  #E <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))@data    # E - read in 1st zeros, non-core data, ith mzone
  #F <- readOGR(analDir,paste0('Zone ',i,' Both - 1st Zero'))@data           # F - read in 1st zeros, all data, ith mzone
}
                                                                                                                                                                             # close out if
states <- as.character(droplevels(unique(A$State)))
nStates <- length(states)

for(i in 1:nStates){

  for(j in 1:nDoThese){
    
#   thisHelper <- helper[helper$modCode == doThese[j] & helper$mZone == i,]

    thisOne <- get(doThese[j])
    thisOne <- thisOne[thisOne$State == states[i],]
    
    thisOne$zeroPresent <- ifelse(thisOne$Peak_Males == 0,1,0)

    leks <- length(unique(thisOne$Lek_IDnum))
    leks.zero <- length(unique(thisOne[thisOne$zeroPresent == 1,]$Lek_IDnum))
    leks.zero.dupped <- length(unique(thisOne[thisOne$zeroPresent == 1 & thisOne$DupZero == 1,]$Lek_IDnum))
    
#     thisOne$yearCls = as.numeric(as.factor(thisOne$Year)) - 26
#     thisOne$lekCls <- as.numeric(as.factor(thisOne$Lek_ID))                     
#     nLeks <- length(leks)
#     
#     Model <- as.character(droplevels(thisHelper$Model))
#     Zeros <- as.character(droplevels(thisHelper$Zeros))
#     Cut   <- as.character(droplevels(thisHelper$Cut))
#     mZone <- i
    
    thisOne$Lek_ID <- as.character(droplevels(thisOne$Lek_ID))

    zeroInv <- data.frame(NDuppedZeros=tapply(thisOne$DupZero, thisOne$Lek_ID, FUN=sum))
    zeroInv$Lek <- rownames(zeroInv) 
    the.zero.mean <- data.frame(Model=doThese[j],NLeks=leks,NLeksZero=leks.zero,NLeksZeroDupped=leks.zero.dupped,AvgDuppedZeroMeanLek=mean(zeroInv[zeroInv$NDuppedZeros > 0,]$NDuppedZeros))
    the.zero.mean$State <- states[i]
    the.zero.mean$PercentLeksWithDup <- round(100*the.zero.mean$NLeksZeroDupped / the.zero.mean$NLeksZero,2)
    the.zeros <- rbind(the.zeros,the.zero.mean)
    write.csv(the.zeros[the.zeros$Model == 'C',],paste0(rsltDir,'/the.zeros.modelC.csv'))
  }
}




