readInAnalyticFiles <- function(analDir){

  for(i in 1:7){
    
    if(i == 1){mZone <- 1}
    if(i == 2){mZone <- 3}
    if(i == 3){mZone <- 4}
    if(i == 4){mZone <- 5}
    if(i == 5){mZone <- 6}
    if(i == 6){mZone <- 8}
    if(i == 7){mZone <- 9}
    
    datAllZerosCore[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - All Zero'))@data        # read in all zeros, core data, ith mzone
    datAllZerosNoco[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - All Zero'))@data    # read in all zeros, non-core data, ith mzone
    datAllZerosLeks[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - All Zero'))@data        # read in all zeros, all data, ith mzone
    
    dat1stZerosCore[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
    dat1stZerosNoco[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Non-Core-75 - 1st Zero'))@data    # read in 1st zeros, non-core data, ith mzone
    dat1stZerosLeks[[mZone]] <<- readOGR(analDir,paste0('Zone ',mZone,' Both-75 - 1st Zero'))@data        # read in 1st zeros, all data, ith mzone
  }
#   datAllZerosCore
#   datAllZerosNoco
#   datAllZerosLeks
#   dat1stZerosCore
#   dat1stZerosNoco
#   dat1stZerosLeks
}