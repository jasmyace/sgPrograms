# t :: nCounts = total number of observations
# u :: nLeks = total number of leks
# v :: nZones = total number of mzones

{
  # combined model
  for(t in 1:nCounts){   
    pMales[t] ~ dpois(lambda[t])                          # count (number of peak males) distribution
    log(lambda[t]) <-                                     # expected number of peaked males 
      gamma000 +
      beta00[zoneCls[t]] + 
      pi0[lekCls[t]] +                                    # lek intercept
      gamma100*(yearCls[t] - medYear) +
      beta10[zoneCls[t]]*(yearCls[t] - medYear) + 
      pi1[lekCls[t]]*(yearCls[t] - medYear) +             # lek slope
      noise[t]                                            # overdispersion / noise
    noise[t] ~ dnorm(0.0,tau.noise)
  }
  
  # leks
  for(u in 1:nLeks){
    #pi0[u] ~ dnorm(beta00[zoneCls[u]],tauinv2a)          # random effect for lek intercept
    #pi1[u] ~ dnorm(beta10[zoneCls[u]],tauinv2b)          # random effect for lek slope   
    pi0[u] ~ dnorm(0,tau.pi0)                            # random effect for lek intercept
    pi1[u] ~ dnorm(0,tau.pi1)                            # random effect for lek slope   
    pi0.adj[u] <- pi0[u] - mean(pi0[])
    pi1.adj[u] <- pi1[u] - mean(pi1[])
  }
  
  # management zones
  for(v in 1:nZones){    
    #beta00[v] ~ dnorm(gamma000,0.001)                    # fixed effect for mzone intercept
    #beta10[v] ~ dnorm(gamma100,0.001)                    # fixed effect for mzone slope    
    beta00[v] ~ dnorm(0,0.001)                            # fixed effect for mzone intercept
    beta10[v] ~ dnorm(0,0.001)                            # fixed effect for mzone slope   
    beta00.adj[v] <- beta00[v] - mean(beta00[])
    beta10.adj[v] <- beta10[v] - mean(beta10[])
  }
  
  # prior specification for fixed effects
  gamma000 ~ dnorm(0,0.001)
  gamma100 ~ dnorm(0,0.001)
  
  # prior specification for the level-1 precision parameter
  tau.noise <- pow(sigma.noise, -2)                       # variation of residual error distribution
  sigma.noise ~ dunif(0,100)                              # sd
  sigma2.noise <- pow(sigma.noise, 2)                     # variance for additive lognormal effect
  
  # prior specification for the level-2 precision parameter
  tau.pi0 <- pow(sigma.pi0,-2)   
  tau.pi1 <- pow(sigma.pi1,-2)   
  sigma.pi0 ~ dunif(0,100)
  sigma.pi1 ~ dunif(0,100)
  sigma2.pi0 <- pow(sigma.pi0, 2)
  sigma2.pi1 <- pow(sigma.pi1, 2)
  
}