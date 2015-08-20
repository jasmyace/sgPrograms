runTestAn <- function(dat,progDir,BUGSDir,thisRun,z,burn,iter,thin){
  
  #dat <- dat1stZerosCore[[1]]
  
  # check if there are any years missing data -- affects B matrix calculations if missing
  if(dim(data.frame(table(dat$Year)))[1] != 51){
    look <- merge(data.frame(Year=seq(1965,2015,1)),data.frame(Year=unique(dat$Year),Here=rep(1,length(unique(dat$Year)))),by=c('Year'),all.x=TRUE)
    miss <- look[is.na(look$Here),]$Year                                           # extract the missing years
    miss.index <- miss - 1964                                                      # scale them to the factor index
    temp <- as.numeric(as.factor(c(miss,dat$Year)))                                # do the factor
    dat$yearCls = temp[-c(miss.index)]                                             # remove the fake years
  } else {                                                                       
    dat$yearCls = as.numeric(as.factor(dat$Year))                                  # the year    
  }
  
  # Add indices for model variables
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  data <- list(nCounts = length(dat$Peak_Males),                                   # number of peak_male counts
               nZones = 1,                                                         # number of management zones - set to 1?
               nLeks = length(unique(dat$Lek_ID)),                                 # number of leks 
               nYears = 51,                                                        # number of years
               pMales = dat$Peak_Males,                                            # count outcome
               medYear = 26,                                                       # median year index (constant) [1990]
               lekCls = dat$lekCls,                                                # integers representing different leks
               yearCls = dat$yearCls)                                              # years from 1-whatever (maybe 0?)
  
  # the use of a f'n that pulls random numbers ensures each chain starts at a diff value
  inits <- function(){                                                             # initialize the gibbs sampler
    list(B=array(rnorm(2*data$nLeks),c(data$nLeks,2)),                             # pull from multivariate normal 
         taunoise=rgamma(1,2),                                                     # pull from rand uni [0,1], positive
         mu.a=rnorm(1),                                                            # pull from normal of 1x1
         sigma.a=runif(1),                                                         # pull from rand uni [0,1], positive
         mu.b=rnorm(1),                                                            # pull from normal of 1x1
         sigma.b=runif(1),                                                         # pull from rand uni [0,1], positive
         rho=runif(1,-0.2,0.2),                                                    # pull from rand uni [-0.2,0.2]
         noise=runif(length(data$pMales),min=0,max=0.5)
    )
  }
  
  parameters <- c("a","b",
                  "taunoise","sdnoise",
                  "mu.a","sigma.a","mu.b","sigma.b","rho",
                  "N","beta",
                  "B10.05.15","B10.04.14","B10.03.13","B10.02.12","B10.01.11",
                  "B10.00.10","B10.99.09","B10.98.08","B10.97.07","B10.96.06",
                  "B10.95.05","B10.94.04","B10.93.03","B10.92.02","B10.91.01",
                  "B10.90.00","B10.89.99","B10.88.98","B10.87.97","B10.86.96",
                  "B10.85.95","B10.84.94","B10.83.93","B10.82.92","B10.81.91",
                  "B10.80.90","B10.79.89","B10.78.88","B10.77.87","B10.76.86",
                  "B10.75.85","B10.74.84","B10.73.83","B10.72.82","B10.71.81",
                  "B10.70.80","B10.69.79","B10.68.78","B10.67.77","B10.66.76",
                  "B10.65.75"
  )                                                                                # the obj to a samp from post dist'n
  
  start <- Sys.time()
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/WBugsAn.R"),                           # winbugs code
                #debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=1,                                                        # n chains >= 3
                n.burnin=burn,
                n.iter=iter,
                n.thin=thin)                                                       # default burn-in tosses half  
  end <- Sys.time()
  time <- end - start
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- "Hello." #glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=dat,family="poisson")
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,time,file=paste0(outpDir,"/Model D MZone ",z," ",thisRun,".RData"))
  list(data,inits,parameters,bayes,CYear,freqs,time)
  
}