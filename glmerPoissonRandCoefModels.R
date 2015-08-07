
require(lme4)

leks <- dat1stZerosCore[[9]]
leks$Mgmt_zone <- as.factor(ifelse(leks$Mgmt_zone %in% c('MZ II','MZ VII'),'MZ VIII',as.character(droplevels(leks$Mgmt_zone))))
leks$mZone_num <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(leks$Mgmt_zone))," ",fixed=TRUE))[c(FALSE,TRUE)]))

CYear <- leks$Year - rep(1964,nrow(leks)) - rep(26,nrow(leks))
test <- glmer(Peak_Males ~ 1 + CYear + Mgmt_zone + CYear*Mgmt_zone + (1 + CYear | Lek_ID),data=leks,family=poisson)

results <- vector("list",6)
for(i in 1:6){
  if(i == 1){mZone <- 1}
  if(i == 2){mZone <- 3}
  if(i == 3){mZone <- 4}
  if(i == 4){mZone <- 5}
  if(i == 5){mZone <- 6}
  if(i == 6){mZone <- 8}
  
  thisLeks <- leks[leks$mZone_num == mZone,]
  CYear <- thisLeks$Year - rep(1964,nrow(thisLeks)) - rep(26,nrow(thisLeks))
  results[[i]] <- glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=thisLeks,family=poisson,verbose=1)
}