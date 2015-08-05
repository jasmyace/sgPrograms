makeSmallCoreSamp <- function(){

# combine regions 1, 6, 5

i <- 9
dat1stZerosCore <- vector("list",9) 
dat1stZerosCore[[i]] <- readOGR(analDir,paste0('Zone ',i,' Core-75 - 1st Zero'))@data        # read in 1st zeros, core data, ith mzone
dat1stZerosCore[[9]]$Mgmt_zone <- as.factor(ifelse(dat1stZerosCore[[9]]$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(dat1stZerosCore[[9]]$Mgmt_zone))))

mz1leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ I',]$Lek_ID)
mz3leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ III',]$Lek_ID)
mz4leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ IV',]$Lek_ID)
mz5leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ V',]$Lek_ID)
mz6leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ VI',]$Lek_ID)
mz8leks <- unique(dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Mgmt_zone == 'MZ VIII',]$Lek_ID)

mz1sampLeks <- as.character(droplevels(sample(mz1leks,877/2)))
mz3sampLeks <- as.character(droplevels(sample(mz3leks,464/2)))
#mz4sampLeks <- as.character(droplevels(sample(mz4leks,30)))
mz5sampLeks <- as.character(droplevels(sample(mz5leks,358/2)))
#mz6sampLeks <- as.character(droplevels(sample(mz6leks,31)))
#mz8sampLeks <- as.character(droplevels(sample(mz8leks,30)))

the1 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz1sampLeks,]
the3 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz3sampLeks,]
#the4 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz4sampLeks,]
the5 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz5sampLeks,]
#the6 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz6sampLeks,]
#the8 <- dat1stZerosCore[[9]][dat1stZerosCore[[9]]$Lek_ID %in% mz8sampLeks,]

rm(dat1stZerosCore)

smallCoreSamp <- rbind(the1,the3,the5)#the4,the5,the6,the8)

}
