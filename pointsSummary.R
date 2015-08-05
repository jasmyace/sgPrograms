pointsSummary <- function(){

  FleksAllZero <- FcoreAllZero <- FnocoAllZero <- Fleks1stZero <- Fcore1stZero <- Fnoco1stZero <- list("vector",9)
  SumFleksAllZero <- SumFcoreAllZero <- SumFnocoAllZero <- SumFleks1stZero <- SumFcore1stZero <- SumFnoco1stZero <- NULL
  
  for(i in 1:9){
    
    FleksAllZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both-75 - All Zero'))
    FcoreAllZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - All Zero'))
    FnocoAllZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - All Zero'))
    Fleks1stZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Both-75 - 1st Zero'))
    Fcore1stZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))
    Fnoco1stZero[[i]] <- readOGR(analDir,paste0('Zone ',i,' Non-Core-75 - 1st Zero'))
      
    nSumFleksAllZero <- data.frame(mZone=FleksAllZero[[i]]@data$mZoneNum[1],Zeros='All Zeros',Cut='All Leks',NData=nrow(FleksAllZero[[i]]@data),meanPeakMales=mean(FleksAllZero[[i]]@data$Peak_Males))
    nSumFcoreAllZero <- data.frame(mZone=FcoreAllZero[[i]]@data$mZoneNum[1],Zeros='All Zeros',Cut='Core Leks',NData=nrow(FcoreAllZero[[i]]@data),meanPeakMales=mean(FcoreAllZero[[i]]@data$Peak_Males))
    nSumFnocoAllZero <- data.frame(mZone=FnocoAllZero[[i]]@data$mZoneNum[1],Zeros='All Zeros',Cut='No-Core Leks',NData=nrow(FnocoAllZero[[i]]@data),meanPeakMales=mean(FnocoAllZero[[i]]@data$Peak_Males))
    nSumFleks1stZero <- data.frame(mZone=Fleks1stZero[[i]]@data$mZoneNum[1],Zeros='1st Zeros',Cut='All Leks',NData=nrow(Fleks1stZero[[i]]@data),meanPeakMales=mean(Fleks1stZero[[i]]@data$Peak_Males))
    nSumFcore1stZero <- data.frame(mZone=Fcore1stZero[[i]]@data$mZoneNum[1],Zeros='1st Zeros',Cut='Core Leks',NData=nrow(Fcore1stZero[[i]]@data),meanPeakMales=mean(Fcore1stZero[[i]]@data$Peak_Males))
    nSumFnoco1stZero <- data.frame(mZone=Fnoco1stZero[[i]]@data$mZoneNum[1],Zeros='1st Zeros',Cut='No-Core Leks',NData=nrow(Fnoco1stZero[[i]]@data),meanPeakMales=mean(Fnoco1stZero[[i]]@data$Peak_Males))
      
    SumFleksAllZero <- rbind(SumFleksAllZero,nSumFleksAllZero)
    SumFcoreAllZero <- rbind(SumFcoreAllZero,nSumFcoreAllZero)
    SumFnocoAllZero <- rbind(SumFnocoAllZero,nSumFnocoAllZero)
    SumFleks1stZero <- rbind(SumFleks1stZero,nSumFleks1stZero)
    SumFcore1stZero <- rbind(SumFcore1stZero,nSumFcore1stZero)
    SumFnoco1stZero <- rbind(SumFnoco1stZero,nSumFnoco1stZero)
  
    SumData <- rbind(SumFleksAllZero,SumFcoreAllZero,SumFnocoAllZero,SumFleks1stZero,SumFcore1stZero,SumFnoco1stZero)
    SumData <- SumData[order(SumData$mZone,SumData$Cut),]
  }
  
write.csv(SumData,paste0(analDir,"/SumData - 2015 Data - 20150731.csv"))

}