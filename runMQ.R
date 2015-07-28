runMQ <- function(dat,progDir,BUGSDir,i){
  
  
  #   dat <- dat1stZerosCore[[9]]
  #   progDir <- progDir
  #   BUGSDir <- BUGSDir
  
  # -----------------------------------------------------------------------------------------------------------------------  
  
  
  # Add indices for model variables.
  dat$zoneCls = as.numeric(as.factor(droplevels(dat$Mgmt_zone)))                   # management zone
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  data <- list(pMales = dat$Peak_Males,                                            # count outcome
               lekCls = dat$lekCls,                                                # integers representing different leks
               yearCls = dat$yearCls,                                              # years from 1-whatever (maybe 0?)
               zoneCls = dat$zoneCls,                                              # zone class
               medYear = 26,                                                       # median year index (constant) [1990]
               nCounts = length(dat$Peak_Males),                                   # number of peak_male counts
               nZones = length(unique(dat$Mgmt_zone)),                             # number of zones
               nLeks = length(unique(dat$Lek_ID)))                                 # number of leks 
  
  # the use of a f'n that pulls random numbers ensures each chain starts at a diff value
  inits <- function(){                                                             # initialize the gibbs sampler
    list(taunoise=rgamma(1,2),                                                     # pull from rand uni [0,1], positive
         noise=rnorm(length(data$pMales),0.1),
         pi0=rnorm(length(unique(data$lekCls)),1),
         pi1=rnorm(length(unique(data$lekCls)),0),
         beta00=rnorm(length(unique(data$zoneCls)),2),
         beta10=rnorm(length(unique(data$zoneCls)),0),
         tauinv2a=rgamma(1,2),#rgamma(length(unique(data$lekCls)),2),
         tauinv2b=rgamma(1,2),#rgamma(length(unique(data$lekCls)),2),
         tauinv3a=rgamma(1,2),#rgamma(length(unique(data$lekCls)),2),
         tauinv3b=rgamma(1,2),#rgamma(length(unique(data$lekCls)),2),
         gamma000=rnorm(1,mean=3.5),
         gamma100=rnorm(1,mean=4.0)
    )
  }
  
  parameters <- c("taunoise",
                  "pi0","pi1","beta00","beta10","gamma000","gamma100",
                  "tauinv2a","tauinv2b","tauinv3a","tauinv3b")                     # the obj to a samp from post dist'n
  
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/TripleLevelRandomEffectMZoneCollapsed.R"),       # winbugs code
                #debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=3,                                                        # n chains >= 3
                n.burnin=50000,
                n.iter=60000,
                n.thin=1)                                                          # default burn-in tosses half
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- "Nothing for now."  
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,file=paste0(outpDir,"/Model Q MZone 9.RData"))
  
  list(data,inits,parameters,bayes,CYear,freqs)
  
}