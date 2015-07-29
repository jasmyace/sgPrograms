
# random slopes & intercepts model, with fixed management zone effect and slope.
# always specify dist'ns and formulas separately.

{
  for(k in 1:nCounts) {                                  # specify model for each obs in the df
    pMales[k] ~ dpois(lambda[k])                         # count (number of peak males) distribution
    log(lambda[k]) <-                                    # expected number of peaked males
      
      m[zoneCls[k]] +                                    # fixed zone effect
      n[zoneCls[k]]*(yearCls[k] - medYear) +             # fixed zone slope
      
      a[lekCls[k]] +                                     # random lek intercept
      b[lekCls[k]]*(yearCls[k] - medYear) +              # random lek slope
      noise[k]                                           # residual error         
    noise[k] ~ dnorm(0.0, taunoise)                      # residual error distribution -- non-informative prior?
  }
  taunoise ~ dgamma(0.001,0.001)                         # variation of residual error distribution
  sdnoise <- 1 / pow(taunoise, 0.5)                      # sd
  
  for(z in 1:nZones){
    m[z] ~ dnorm(0,.0001)                              # for the fixed-effect zone intercepts
    n[z] ~ dnorm(0,.0001)                              # for the fixed-effect zone slopes
  }
  
  for(j in 1:nLeks){
    a[j] <- B[j,1]                                     # put random intercepts in 1st col of B matrix
    b[j] <- B[j,2]                                     # put random slopes in 2nd col of B matrix
    B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])             # sample from mulitvariate normal
    B.hat[j,1] <- mu.a                                 # make a matrix of overall intercept
    B.hat[j,2] <- mu.b                                 # make a matrix of overall slope
  }
  
  mu.a ~ dnorm(0,.0001)                                # like an overall intercept -- uninformative prior
  mu.b ~ dnorm(0,.0001)                                # like an overall slope -- uninformative prior
  Tau.B[1:2,1:2] <- inverse(Sigma.B[,])                # take inverse of covariance matrix to get precision matrix
  Sigma.B[1,1]<-pow(sigma.a,2)                         # square estimate of sd(int) and put in [1,1] spot in Sigma.B matrix
  sigma.a ~ dunif(0,100)                               # positive uninformative prior for estimate of sigma.a
  Sigma.B[2,2] <- pow(sigma.b,2)                       # square estimate of sd(slope) and put in [2,2] spot in Sigma.B matrix
  sigma.b ~ dunif(0,100)                               # positive uninformative prior for estimate of sigma.b
  
  #correlation
  Sigma.B[1,2] <- rho*sigma.a*sigma.b                  # calculate covariance of intercepts and slopes and place in [1,2]
  Sigma.B[2,1] <- Sigma.B[1,2]                         # and in [2,1] of covariance matrix of Sigma.B
  rho ~ dunif(-1,1)                                    # uninformative prior on correlation
   
  # make the B-matrices
  for (i in 1:nZones){
    for (j in 1:nYears) {
      N[i,j] <- exp(m[i] +                                    # fixed zone effect
                    n[i]*(yearCls[j] - medYear) +             # fixed zone slope                        
                    mu.a +                                    # grand intercept
                    mu.b*(yearCls[j] - medYear) +             # grand slope
                    0.5*sdnoise*sdnoise)
    }
  }
  for (i in 1:nZones) {
    beta[i] <- pow(N[i,nYears] / N[i,1], 1/(nYears-1))
  }
  for (j in 1:nYears) {
    N0[j] <- sum(N[1:nZones,j])
  }
  beta0 <- pow(N0[nYears] / N0[1], 1/(nYears-1)) 
}