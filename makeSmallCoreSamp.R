makeSmallCoreSamp <- function(){

# combine regions 1, 6, 5

i <- 9
dat1stZerosCore <- vector("list",9) 
dat1stZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
dat1stZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosCore[[9]]$Mgmt_zone))))

mz1leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ I',]$Lek_ID)
mz5leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ V',]$Lek_ID)
mz6leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ VI',]$Lek_ID)

mz1sampLeks <- as.character(droplevels(sample(mz1leks,30)))
mz5sampLeks <- as.character(droplevels(sample(mz5leks,30)))
mz6sampLeks <- as.character(droplevels(sample(mz6leks,30)))

the1 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz1sampLeks,]
the5 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz5sampLeks,]
the6 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz6sampLeks,]

rm(dat1stZerosCore)

smallCoreSamp <- rbind(the1,the5,the6)

}
