
  # random slopes & intercepts model
  # always specify dist'ns and formulas separately.
  
{
  for(k in 1:nCounts) {                                  # specify model for each obs in the df
    pMales[k] ~ dpois(lambda[k])                         # count (number of peak males) distribution
    log(lambda[k]) <-                                    # expected number of peaked males
    
      a[lekCls[k]] +                                     # random lek intercept
      b[lekCls[k]]*(yearCls[k] - medYear) +              # random lek slope
      noise[k]                                           # residual error         
    noise[k] ~ dnorm(0.0, taunoise)                      # residual error distribution -- non-informative prior?
  }
  taunoise ~ dgamma(0.001,0.001)                         # variation of residual error distribution
  sdnoise <- 1 / pow(taunoise, 0.5)                      # sd

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
}


#####################################################
##### SUMMARY STATISTICS  
#####################################################

# make the B-vectors
for (j in 1:nYears) {
  n[j] <- exp(mu.a + 
              mu.b*(yearCls[j] - medYear) +                             
              0.5*sd.noise*sd.noise)
}

# B-vector: compare first year to all other years
for (i in 1:nZones) {
  beta[i] <- pow(n[i,nYears] / n[i,1], 1/(nYears-1))
}

# build up the N_t and \bar{B}
for (j in 1:nYears) {
  n_0[j] <- sum(n[1:nZones,j])
}
beta0 <- pow(n_0[nYears] / n_0[1], 1/(nYears-1))
