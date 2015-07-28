


# we want to bas sample 15% from each of the 6 management zones

getIt <- readOGR('//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Analysis Sets','Zone 9 Core-75 - 1st Zero')
getIt@data$Mgmt_zone <- as.factor(ifelse(getIt$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(getIt$Mgmt_zone))))
getIt@data$R_ID <- c(1:nrow(getIt@data))

oneLek <- getIt[!duplicated(getIt@data$Lek_ID),]

# now, grts sample from the oneLek spdf.  each lek has equal weight. 4,236 of them. 

#uniques <- unique(getIt@data[,c('Lek_ID','Mgmt_zone','Lat','Long')])



# run grts function to make shapefile Core1stZeros15p
# Core1stZeros15p has the LEKS we care about now.
# now we have to fetch their temporal data.

sample <- Core1stZeros15p@data[,c('Lek_ID','siteID')]

coreBAS15p <- getIt[getIt@data$Lek_ID %in% sample$Lek_ID,]


writeOGR(coreBAS15p,"//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Analysis Sets","Core1stZeros15p",overwrite_layer=TRUE,driver="ESRI Shapefile") 