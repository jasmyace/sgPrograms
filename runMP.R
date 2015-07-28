runMP <- function(dat,progDir,BUGSDir,i){
  
  
  #   dat <- Core1stZeros15p
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
    list(sigma.noise=rlnorm(1),                                                     # pull from rand uni [0,1], positive
         noise=rnorm(length(data$pMales),0.1),
         pi0=rnorm(length(unique(data$lekCls)),1),
         pi1=rnorm(length(unique(data$lekCls)),0),
         beta00=rnorm(length(unique(data$zoneCls)),2),
         beta10=rnorm(length(unique(data$zoneCls)),0),
         sigma.pi0=rlnorm(1),
         sigma.pi1=rlnorm(1),
         gamma000=rnorm(1,mean=3.5),
         gamma100=rnorm(1,mean=4.0)
    )
  }
  
  parameters <- c("sigma2.noise",
                  "pi0","pi1","beta00","beta10","gamma000","gamma100",
                  "sigma2.pi0","sigma2.pi1",
                  "pi0.adj","pi1.adj","beta00.adj","beta10.adj")                   # the obj to a samp from post dist'n
  
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/TripleLevelFixedEffectMZoneCollapsed.R"),       # winbugs code
                debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=3,                                                        # n chains >= 3
                n.burnin=36000,
                n.iter=40000,
                n.thin=1)                                                          # default burn-in tosses half
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- "Nothing for now."  
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,file=paste0(outpDir,"/Model P MZone 9.RData"))
  
  list(data,inits,parameters,bayes,CYear,freqs)
  
}