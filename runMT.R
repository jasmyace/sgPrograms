runMT <- function(dat,progDir,BUGSDir,thisRun,z){
  
  #dat <- smallCoreSamp2


  
  # Add indices for model variables.
  dat$zoneCls = as.numeric(as.factor(droplevels(dat$Mgmt_zone)))                   # management zone
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
 
  est.m <- mean(data.frame(MeanPMales=log(tapply(dat[dat$Year == 1990,]$Peak_Males, list(dat[dat$Year == 1990,]$zoneCls), mean)))$MeanPMales)
  
  
  zoneCls <- sort(rep(unique(dat$zoneCls),length(sort(unique(dat$yearCls)))))
  yearCls <- rep(sort(unique(dat$yearCls)),length(unique(dat$zoneCls)))
  BackBone <- as.data.frame(cbind(zoneCls,yearCls))
  tempAgg <- aggregate(data=dat, lekCls ~ zoneCls + yearCls, function(x) length(unique(x)))
  preLekMat <- merge(BackBone,tempAgg,by=c('zoneCls','yearCls'),all.x=TRUE)
  preLekMat <- preLekMat[order(preLekMat$zoneCls,preLekMat$yearCls),]
  preLekMat[is.na(preLekMat)] <- 0
  
  lekMat <- matrix(NA,nrow=length(unique(dat$yearCls)),ncol=length(unique(dat$zoneCls)))
  mZoneClasses <- sort(unique(dat$zoneCls))
  for(i in 1:length(unique(zoneCls))){
    if(i == 1){
      getMZone <-mZoneClasses[1] 
      lekMat[,1] <- preLekMat[preLekMat$zoneCls == getMZone,]$lekCls
    } else {
      getMZone <-mZoneClasses[i] 
      lekMat[,i] <- preLekMat[preLekMat$zoneCls == getMZone,]$lekCls
    }
  }
  nLeksMat <- t(lekMat)
  nLeksVec <- colSums(nLeksMat)
  
#   # check.
#   thing <- rep(NA,51)
#   for(j in 1:51){
#     thing[j] <- sum(nLeksMat[1:3,j] * nLeksMat[1:3,j]) / nLeksVec[j]
#   }

  
  data <- list(pMales = dat$Peak_Males,                                            # count outcome
               lekCls = dat$lekCls,                                                # integers representing different leks
               yearCls = dat$yearCls,                                              # years from 1-whatever (maybe 0?)
               zoneCls = dat$zoneCls,                                              # zone class
               medYear = 26,                                                       # median year index (constant) [1990]
               nCounts = length(dat$Peak_Males),                                   # number of peak_male counts 
               nZones = length(unique(dat$Mgmt_zone)),                             # number of zones
               nYears = length(unique(dat$Year)),                                  # number of years
               nLeks = length(unique(dat$Lek_ID)),                                 # number of leks
               nLeksMat = matrix(nLeksMat,nrow=length(unique(dat$Mgmt_zone)),ncol=length(unique(dat$Year))),# zones x time of lek counts
               nLeksVec = matrix(nLeksVec,nrow=1,ncol=length(unique(dat$Year))))   # 1 x time of lek counts
               
  # the use of a f'n that pulls random numbers ensures each chain starts at a diff value
  inits <- function(){                                                             # initialize the gibbs sampler
    list(B=array(c(rnorm(data$nLeks,mean=3),rnorm(data$nLeks)),c(data$nLeks,2)),   # pull from multivariate normal 
         taunoise=rgamma(1,2),                                                     # pull from rand uni [0,1], positive
         mu.a=rnorm(1),                                                            # pull from normal of 1x1
         sigma.a=runif(1),                                                         # pull from rand uni [0,1], positive
         mu.b=rnorm(1),                                                            # pull from normal of 1x1
         sigma.b=runif(1),                                                         # pull from rand uni [0,1], positive
         m=rnorm(length(unique(dat$Mgmt_zone)),mean=est.m),
         n=rnorm(length(unique(dat$Mgmt_zone))),
         rho=runif(1,-0.2,0.2),                                                    # pull from rand uni [-0.2,0.2], positive
         noise=rep(0.01,length(dat$Peak_Males))
    )
  }
  
  parameters <- c(
                  #"a",
                  #"b",
                  "m",
                  "m.adj",
                  "n.adj",
                  "n",
                  "taunoise",
                  "sdnoise",
                  "mu.a",
                  "sigma.a",
                  "mu.b",
                  "sigma.b",
                  "rho",
                  "beta","beta0","N","N0",
                  "B10.05.15",
                  "B10.04.14",
                  "B10.03.13",
                  "B10.02.12",
                  "B10.01.11",
                  "B10.00.10",
                  "B10.99.09",
                  "B10.98.08",
                  "B10.97.07",
                  "B10.96.06",
                  "B10.95.05",
                  "B10.94.04",
                  "B10.93.03",
                  "B10.92.02",
                  "B10.91.01",
                  "B10.90.00",
                  "B10.89.99",
                  "B10.88.98",
                  "B10.87.97",
                  "B10.86.96",
                  "B10.85.95",
                  "B10.84.94",
                  "B10.83.93",
                  "B10.82.92",
                  "B10.81.91",
                  "B10.80.90",
                  "B10.79.89",
                  "B10.78.88",
                  "B10.77.87",
                  "B10.76.86",
                  "B10.75.85",
                  "B10.74.84",
                  "B10.73.83",
                  "B10.72.82",
                  "B10.71.81",
                  "B10.70.80",
                  "B10.69.79",
                  "B10.68.78",
                  "B10.67.77",
                  "B10.66.76",
                  "B10.65.75",                 
                  "BX10.05.15",
                  "BX10.04.14",
                  "BX10.03.13",
                  "BX10.02.12",
                  "BX10.01.11",
                  "BX10.00.10",
                  "BX10.99.09",
                  "BX10.98.08",
                  "BX10.97.07",
                  "BX10.96.06",
                  "BX10.95.05",
                  "BX10.94.04",
                  "BX10.93.03",
                  "BX10.92.02",
                  "BX10.91.01",
                  "BX10.90.00",
                  "BX10.89.99",
                  "BX10.88.98",
                  "BX10.87.97",
                  "BX10.86.96",
                  "BX10.85.95",
                  "BX10.84.94",
                  "BX10.83.93",
                  "BX10.82.92",
                  "BX10.81.91",
                  "BX10.80.90",
                  "BX10.79.89",
                  "BX10.78.88",
                  "BX10.77.87",
                  "BX10.76.86",
                  "BX10.75.85",
                  "BX10.74.84",
                  "BX10.73.83",
                  "BX10.72.82",
                  "BX10.71.81",
                  "BX10.70.80",
                  "BX10.69.79",
                  "BX10.68.78",
                  "BX10.67.77",
                  "BX10.66.76",
                  "BX10.65.75"
                  )                                             # the obj to a samp from post dist'n
  
  start <- Sys.time()
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/CorrRandSlopesIntsFixedMZoneBMatAdjM.R"),  # winbugs code
                #debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=1,                                                        # n chains >= 3
                n.burnin=76000,
                n.iter=80000,
                n.thin=1)                                                          # default burn-in tosses half
  end <- Sys.time()
  time <- end - start
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- "Hello."#glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=dat,family="poisson")
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,time,file=paste0(outpDir,"/Model T MZone ",z," ",thisRun,".RData"))
  
  list(data,inits,parameters,bayes,CYear,freqs,time)
  
}