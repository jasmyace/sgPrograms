
mZSpDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/Management Zones/SG_MgmtZones_ver2_20061018'
USSpDir <- '//LAR-FILE-SRV/Data/Jason/sage grouse/Data/Spatial/US States'

sg075Shp <- readOGR(polyDir,'MZone_9_the75')# get 75% core
sg100Shp <- readOGR(polyDir,'MZone_9_the100')# get 100% core
mzShp  <- readOGR(mZSpDir,'SG_MgmtZones_ver2_20061018')# get 75% core

usShp <- readOGR(USSpDir,"states")


mzShp <- spTransform(mzShp,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "))
usShp <- spTransform(usShp,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "))

usShp <- usShp[usShp@data$STATE_ABBR %in% c("WA","OR","CA","NV","UT","CO","WY","MT","ID","SD","ND"),]

mzCents <- gCentroid(mzShp,byid=TRUE)
mzNames <- c("2","4","3","7","1","5","6")

# core map
CairoPNG(filename=paste0(rsltDir,'/MZones + Core.png'),width=6,height=6,units="in",res=500,quality=600,pointsize=12)
par(mar=c(0,0,0,0))
plot(usShp,col="lightblue",border="white")
plot(mzShp,add=TRUE,lwd=2)
plot(sg100Shp,col="gray95",border="gray95",add=TRUE)
plot(sg075Shp,col="pink" ,border="pink",add=TRUE)
text(mzCents@coords[,1],mzCents@coords[,2],mzNames,cex=1.5)
legend("bottomright",c("Management Zone Boundary","US State","100% Contour","75% Contour 'Core'"),cex=0.8,col=c("black","lightblue","gray95","pink"),bty="n",pch=c(NA,15,15,15),lwd=c(2,NA,NA,NA))
dev.off()


# mzone map
CairoPNG(filename=paste0(rsltDir,'/MZones.png'),width=6,height=6,units="in",res=500,quality=600,pointsize=12)
par(mar=c(0,0,0,0))
plot(usShp,col="lightblue",border="white")
plot(mzShp,add=TRUE,lwd=2)
text(mzCents@coords[,1],mzCents@coords[,2],mzNames,cex=1.5)
legend("bottomright",c("Management Zone Boundary","US State"),cex=0.8,col=c("black","lightblue"),bty="n",pch=c(NA,15),lwd=c(2,NA))
dev.off()