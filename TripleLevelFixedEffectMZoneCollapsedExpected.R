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
    pi0[u] ~ dnorm(beta00[zoneCls[u]],tau.pi0)          # random effect for lek intercept
    pi1[u] ~ dnorm(beta10[zoneCls[u]],tau.pi1)          # random effect for lek slope   
    e.pi0[u] <- pi0[u] - beta00[zoneCls[u]] 
    e.pi1[u] <- pi1[u] - beta10[zoneCls[u]] 
    f.pi0[u] <- pi0[u] - e.beta00[zoneCls[u]] 
    f.pi1[u] <- pi1[u] - e.beta10[zoneCls[u]] 
  }
  
  # management zones
  for(v in 1:nZones){    
    beta00[v] ~ dnorm(gamma000,0.001)                    # fixed effect for mzone intercept
    beta10[v] ~ dnorm(gamma100,0.001)                    # fixed effect for mzone slope      
    e.beta00[v] <- beta00[v] - gamma000
    e.beta10[v] <- beta10[v] - gamma100
  }
  
  # prior specification for fixed effects
  gamma000 ~ dnorm(0,0.001)
  gamma100 ~ dnorm(0,0.001)
  e.gamma000 <- gamma000 + mean(beta00[])
  e.gamma100 <- gamma100 + mean(beta10[])
  f.gamma000 <- gamma000 + mean(beta00[]) + mean(pi0[])
  f.gamma100 <- gamma100 + mean(beta10[]) + mean(pi1[])
  
  # prior specification for the level-1 precision parameter
  tau.noise <- pow(sigma.noise, -2)                       # variation of residual error distribution
  sigma.noise ~ dunif(0,5)                                # sd
  sigma2.noise <- pow(sigma.noise, 2)                     # variance for additive lognormal effect
  
  # prior specification for the level-2 precision parameter
  tau.pi0 <- pow(sigma.pi0,-2)   
  tau.pi1 <- pow(sigma.pi1,-2)   
  sigma.pi0 ~ dunif(0,5)
  sigma.pi1 ~ dunif(0,5)
  sigma2.pi0 <- pow(sigma.pi0, 2)
  sigma2.pi1 <- pow(sigma.pi1, 2)
  
  # make the B-matrices
  
  ##### Calculate the scaling factor for the BBS data.  NOTE INDEXING FOR TIME:
  ##### Only covers years of West survey from years 40 to 44 (2006-2010) #####
  for (i in 1:nZones){
    for (j in 1:nYears) {
      n[i,j] <- exp(gamma000 +
                    beta00[zoneCls[i]] +                                
                    gamma100*(yearCls[j] - medYear) +
                    beta10[zoneCls[i]]*(yearCls[j] - medYear) +
                    0.5*sigma2.noise)
    }
  }
  for (i in 1:nZones) {
    beta[i] <- pow(n[i,nYears] / n[i,1], 1/(nYears-1))
  }
  for (j in 1:nYears) {
    n_0[j] <- sum(n[1:nZones,j])
  }
  beta0 <- pow(n_0[nYears] / n_0[1], 1/(nYears-1))

}