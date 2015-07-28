# t :: nCounts = total number of observations
# u :: nLeks = total number of leks
# v :: nZones = total number of mzones

{
  # level-1 model
  for(t in 1:nCounts){   
    pMales[t] ~ dpois(lambda[t])                         # count (number of peak males) distribution
    log(lambda[t]) <-                                    # expected number of peaked males 
      pi0[lekCls[t]] +                                   # lek intercept
      pi1[lekCls[t]]*(yearCls[t] - medYear) +            # lek slope
      noise[t]                                           # overdispersion / noise
    noise[t] ~ dnorm(0.0,taunoise)
  }

  # level-2 model
  for(u in 1:nLeks){
    pi0[u] ~ dnorm(beta00[zoneCls[u]],tauinv2a[zoneCls[u]])  # random effect for lek intercept
    pi1[u] ~ dnorm(beta10[zoneCls[u]],tauinv2b[zoneCls[u]])  # random effect for lek slope    
  }
  
  # level-3 model
  for(v in 1:nZones){
    beta00[v] ~ dnorm(gamma000,tauinv3a)                # random effect for mzone intercept
    beta10[v] ~ dnorm(gamma100,tauinv3b)                # random effect for mzone slope
  }
  
  # prior specification for fixed effects
  gamma000 ~ dnorm(0,1)
  gamma100 ~ dnorm(0,1)
  
  # prior specification for the level-1 precision parameter
  taunoise ~ dgamma(0.1,0.1)                           # variation of residual error distribution
  sdnoise <- 1 / pow(taunoise, 0.5)                      # sd
  
  # prior specification for the level-2 precision parameter
  for(v in 1:nZones){
    tauinv2a[v] ~ dgamma(0.1,0.1)   
    tauinv2b[v] ~ dgamma(0.1,0.1)   
    tauvar2a[v] <- 1/tauinv2a[v]
    tauvar2b[v] <- 1/tauinv2b[v]
  }
  
  # prior specification for the level-3 precision parameters
  tauinv3a ~ dgamma(0.1,0.1)   
  tauinv3b ~ dgamma(0.1,0.1)   
  tauvar3a <- 1/tauinv3a
  tauvar3b <- 1/tauinv3b
}