runMD10 <- function(dat,progDir,BUGSDir,thisRun,z,j,burn,iter,thin){
  
  #dat <- ten1stZerosCore 
  
  # check if there are any years missing data -- affects B matrix calculations if missing
  if(dim(data.frame(table(dat$Year)))[1] != 11){
    look <- merge(data.frame(Year=seq(1964 + j,1964 + j + 10)),data.frame(Year=unique(dat$Year),Here=rep(1,length(unique(dat$Year)))),by=c('Year'),all.x=TRUE)
    miss <- look[is.na(look$Here),]$Year                                           # extract the missing years
    miss.index <- miss - 1963 - j                                                  # scale them to the factor index
    temp <- as.numeric(as.factor(c(miss,dat$Year)))                                # do the factor
    dat$yearCls = temp[temp != c(miss.index)]                                      # remove the fake years
  } else {                                                                       
    dat$yearCls = as.numeric(as.factor(dat$Year))                                  # the year    
  }
  
  # Add indices for model variables
  dat$lekCls = as.numeric(as.factor(droplevels(dat$Lek_ID)))                       # lek
  
  data <- list(nCounts = length(dat$Peak_Males),                                   # number of peak_male counts
               nZones = 1,                                                         # number of management zones - set to 1?
               nLeks = length(unique(dat$Lek_ID)),                                 # number of leks 
               nYears = 11,                                                        # number of years
               pMales = dat$Peak_Males,                                            # count outcome
               medYear = 6,                                                        # median year index (constant) [start year + 5]
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
                  "N","beta"
  )                                                                                # the obj to a samp from post dist'n
  
  start <- Sys.time()
  bayes <- bugs(data=data,                                                         # data to feed to winbugs
                inits=inits,                                                       # start chain with these
                parameters.to.save=parameters,                                     # monitor these
                model.file=paste0(progDir,"/CorrRandSlopesIntsBMat10.R"),          # winbugs code
               #debug=TRUE,
                bugs.directory=BUGSDir,
                n.chains=1,                                                        # n chains >= 3
                n.burnin=burn,
                n.iter=iter,
                n.thin=thin)                                                       # default burn-in tosses half  
  end <- Sys.time()
  time <- end - start
  
  CYear <- dat$yearCls - rep(26,nrow(dat))
  freqs <- seq(1964 + j,1964 + j + 10) #glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=dat,family="poisson")
  
  # output results
  save(data,inits,parameters,bayes,CYear,freqs,time,file=paste0(outpDir,"/Model D MZone ",z," ",thisRun,".RData"))
  list(data,inits,parameters,bayes,CYear,freqs,time)

}