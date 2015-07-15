runMD <- function(dat,progDir,BUGSDir,i){
  
  # Add indices for model variables.
  dat$yearCls = as.numeric(as.factor(dat$Year))                                    # the year
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  data <- list(nCounts = length(dat$Peak_Males),                                   # number of peak_male counts
               pMales = dat$Peak_Males,                                            # count outcome
               medYear = 26,                                                       # median year index (constant) [1990]
               lekCls = dat$lekCls,                                                # integers representing different leks
               yearCls = dat$yearCls,                                              # years from 1-whatever (maybe 0?)
               nLeks = length(unique(dat$Lek_ID)))                                 # number of leks 
  
  # the use of a f'n that pulls random numbers ensures each chain starts at a diff value
  inits <- function(){                                                             # initialize the gibbs sampler
    list(B=array(rnorm(2*data$nLeks),c(data$nLeks,2)),                             # pull from multivariate normal 
         taunoise=rgamma(1,2),                                                     # pull from rand uni [0,1], positive
         mu.a=rnorm(1),                                                            # pull from normal of 1x1
         sigma.a=runif(1),                                                         # pull from rand uni [0,1], positive
         mu.b=rnorm(1),                                                            # pull from normal of 1x1
         sigma.b=runif(1),                                                         # pull from rand uni [0,1], positive
         rho=runif(1),                                                             # pull from rand uni [0,1], positive
         noise=rep(-0.1,length(data$pMales))
    )
  }
  
  parameters <- c("a","b",
                  "taunoise",
                  "mu.a","sigma.a",
                  "mu.b","sigma.b","rho")                                          # the obj to a samp from post dist'n
  
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/CorrRandSlopesInts.R"),                # winbugs code
                #debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=3,                                                        # n chains >= 3
                n.burnin=36000,
                n.iter=40000,
                n.thin=1)                                                          # default burn-in tosses half
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=dat,family="poisson")
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,file=paste0(outpDir,"/Model D MZone",i,".RData"))
  
  list(data,inits,parameters,bayes,CYear,freqs)

}