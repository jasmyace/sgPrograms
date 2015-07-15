runAllZerosLeksModel <- function(shp){
  
    # Add indices for model variables.
    countDat[[i]]$zoneCls = as.numeric(as.factor(droplevels(countDat[[i]]$Mgmt_zone)))                   # management zone
    countDat[[i]]$yearCls = as.numeric(as.factor(countDat[[i]]$Year))                                    # the year
    countDat[[i]]$lekCls = as.numeric(as.factor(droplevels(countDat[[i]]$Lek_ID)))                       # lek
    
    data[[i]] <- list(nCounts = length(countDat[[i]]$Peak_Males),                                        # number of peak_male counts
                      pMales = countDat[[i]]$Peak_Males,                                                 # count outcome
                      medYear = 26,                                                                      # median year index (constant) [1990]
                      lekCls = countDat[[i]]$lekCls,                                                     # integers representing different leks
                      yearCls = countDat[[i]]$yearCls,                                                   # years from 1-whatever (maybe 0?)
                      #                     nYears = length(unique(countDat[[i]]$yrCls)),                                      # number of years
                      nLeks = length(unique(countDat[[i]]$Lek_ID)))                                      # number of leks 
    
    # the use of a f'n that pulls random numbers ensures each chain starts at a diff value
    inits[[i]] <- function(){                                                                            # initialize the gibbs sampler
      list(B=array(rnorm(2*data[[i]]$nLeks),c(data[[i]]$nLeks,2)),                                       # pull from multivariate normal 
           taunoise=rgamma(1,2),                                                                         # pull from rand uni [0,1], positive
           mu.a=rnorm(1),                                                                                # pull from normal of 1x1
           sigma.a=runif(1),                                                                             # pull from rand uni [0,1], positive
           mu.b=rnorm(1),                                                                                # pull from normal of 1x1
           sigma.b=runif(1),                                                                             # pull from rand uni [0,1], positive
           rho=runif(1),                                                                                 # pull from rand uni [0,1], positive
           noise=rep(-0.1,length(countDat[[i]]$Peak_Males))
      )
    }
    
    parameters[[i]] <- c("a","b",
                         "taunoise",
                         "mu.a","sigma.a",
                         "mu.b","sigma.b","rho")                                                         # the obj to a samp from post dist'n
    
    bayes[[i]] <- bugs(data=data[[i]],                                                                   # data to feed to winbugs
                       inits=inits[[i]],                                                                 # start chain with these
                       parameters.to.save=parameters[[i]],                                               # monitor these
                       model.file=paste0(progDir,"/CorrRandSlopesInts.R"),                               # winbugs code
                       #debug=TRUE,
                       bugs.directory="C:/Users/jmitchell/WinBUGS_14/winbugs14/WinBUGS14",
                       n.chains=3,                                                                       # n chains >= 3
                       n.burnin=36000,
                       n.iter=40000,
                       n.thin=1)                                                                        # default burn-in tosses half
    
    
    CYear <- countDat[[i]]$yearCls - rep(26,nrow(countDat[[i]]))
    freqs[[i]] <- glmer(Peak_Males ~ 1 + CYear + (1 + CYear | Lek_ID),data=countDat[[i]],family="poisson")
    
    # output results
    sBayes <- bayes[[i]]
    sFreqs <- freqs[[i]]
    save(sFreqs,sBayes,file=paste0(outpDir,"/MZone",i,".RData"))
    rm(sBayes,sFreqs)
    
  }
}