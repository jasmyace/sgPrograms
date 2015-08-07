
  # random slopes & intercepts model
  # always specify dist'ns and formulas separately.
  
{
  for(k in 1:nCounts) {                                # specify model for each obs in the df
    #w[k] ~ dbern(psi)
    #pMales[k] ~ dpois(eff.lambda[k])                   # count (number of peak males) distribution
    pMales[k] ~ dpois(lambda[k])
    #eff.lambda[k] <- w[k]*lambda[k]
    log(lambda[k]) <-                                  # expected number of peaked males
      
      a[lekCls[k]] +                                   # random lek intercept
      b[lekCls[k]]*(yearCls[k] - medYear) +            # random lek slope
      noise[k]                                         # residual error         
    noise[k] ~ dnorm(0.0, taunoise)                    # residual error distribution -- non-informative prior?
  }
  taunoise ~ dgamma(0.001,0.001)                       # variation of residual error distribution
  sdnoise <- 1 / pow(taunoise, 0.5)                    # sd

#   psi ~ dunif(0,1)
#   R.lpsi <- logit(1 - psi)  
  
  for(j in 1:nLeks){                                   # loop over leks
    a[j]~dnorm(a.hat[j],tau.a)
    b[j]~dnorm(b.hat[j],tau.b)
    a.hat[j]<-mu.a
    b.hat[j]<-mu.b
#     a[j] <- B[j,1]                                     # put random intercepts in 1st col of B matrix
#     b[j] <- B[j,2]                                     # put random slopes in 2nd col of B matrix
#     B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])             # sample from mulitvariate normal
#     B.hat[j,1] <- mu.a                                 # make a vector of overall intercept
#     B.hat[j,2] <- mu.b                                 # make a vector of overall slope
  }

  mu.a ~ dnorm(0,.0001)                                # like an overall intercept -- uninformative prior
  mu.b ~ dnorm(0,.0001)                                # like an overall slope -- uninformative prior

  tau.a<-pow(sigma.a,-2)
  sigma.a~dunif(0,100)
  tau.b<-pow(sigma.b,-2)
  sigma.b~dunif(0,100)

#   Tau.B[1:2,1:2] <- inverse(Sigma.B[,])                # take inverse of covariance matrix to get precision matrix
#   Sigma.B[1,1]<-pow(sigma.a,2)                         # square estimate of sd(int) and put in [1,1] spot in Sigma.B matrix
#   sigma.a ~ dunif(0,100)                               # positive uninformative prior for estimate of sigma.a
#   Sigma.B[2,2] <- pow(sigma.b,2)                       # square estimate of sd(slope) and put in [2,2] spot in Sigma.B matrix
#   sigma.b ~ dunif(0,100)                               # positive uninformative prior for estimate of sigma.b
  
  #correlation
#   Sigma.B[1,2] <- rho*sigma.a*sigma.b                  # calculate covariance of intercepts and slopes and place in [1,2]
#   Sigma.B[2,1] <- Sigma.B[1,2]                         # and in [2,1] of covariance matrix of Sigma.B
#   rho ~ dunif(-1,1)                                    # uninformative prior on correlation

  # make the B-matrices
  for (d in 1:nZones){
    for (e in 1:nYears) {
#       N[i,j] <- exp(mu.a +                             # mean intercept for mzone 
#                     mu.b*(yearCls[j] - medYear) +      # mean slope for mzone
#                     0.5*sdnoise*sdnoise)               # log-normal adjustment
       N[d,e] <- exp(mu.a +                        # mean intercept for mzone 
                      mu.b*(yearCls[e] - medYear) +    # mean slope for mzone
                      0.5*sdnoise*sdnoise)            # log-normal adjustment
#        O[d,e] <- psi*exp(mu.a +                        # mean intercept for mzone 
#                          mu.b*(yearCls[e] - medYear) +    # mean slope for mzone
#                          0.5*sdnoise*sdnoise)            # log-normal adjustment
#        P[d,e] <- (1-psi)*exp(mu.a +                        # mean intercept for mzone 
#                              mu.b*(yearCls[e] - medYear) +    # mean slope for mzone
#                              0.5*sdnoise*sdnoise)            # log-normal adjustment
#        
    }
  }

  # all-year estimates.
  for (d in 1:nZones) {
    beta[d] <- pow(N[d,nYears] / N[d,1], 1/(nYears-1))
  }
  
#   # weight the expected Ns, for all mzones and time periods -- build up to wgt average
#   for(i in 1:nZones){
#     for(j in 1:nYears){
#       preN0[i,j] <- nLeksMat[i,j] * N[i,j]
#     }
#   }
  
#   # sum over each time period, and divide by total number of leks -- gets wgt average
#   for(j in 1:nYears){
#     N0[j] <- sum(preN0[1:nZones,j]) / nLeksVec[1,j]
#   }
  
#   # calculate grand weighted-mean slope
#   beta0 <- pow(N0[nYears] / N0[1], 1/(nYears-1)) 
  
  # 10-year Bs
  for (d in 1:nZones) {
    B10.05.15[d] <- pow(N[d,51] / N[d,41], 1/(51-41))
    B10.04.14[d] <- pow(N[d,50] / N[d,40], 1/(50-40))
    B10.03.13[d] <- pow(N[d,49] / N[d,39], 1/(49-39))
    B10.02.12[d] <- pow(N[d,48] / N[d,38], 1/(48-38))
    B10.01.11[d] <- pow(N[d,47] / N[d,37], 1/(47-37))
    B10.00.10[d] <- pow(N[d,46] / N[d,36], 1/(46-36))
    B10.99.09[d] <- pow(N[d,45] / N[d,35], 1/(45-35))
    B10.98.08[d] <- pow(N[d,44] / N[d,34], 1/(44-34))
    B10.97.07[d] <- pow(N[d,43] / N[d,33], 1/(43-33))
    B10.96.06[d] <- pow(N[d,42] / N[d,32], 1/(42-32))
    B10.95.05[d] <- pow(N[d,41] / N[d,31], 1/(41-31))
    B10.94.04[d] <- pow(N[d,40] / N[d,30], 1/(40-30))
    B10.93.03[d] <- pow(N[d,39] / N[d,29], 1/(39-29))
    B10.92.02[d] <- pow(N[d,38] / N[d,28], 1/(38-28))
    B10.91.01[d] <- pow(N[d,37] / N[d,27], 1/(37-27))
    B10.90.00[d] <- pow(N[d,36] / N[d,26], 1/(36-26))
    B10.89.99[d] <- pow(N[d,35] / N[d,25], 1/(35-25))
    B10.88.98[d] <- pow(N[d,34] / N[d,24], 1/(34-24))
    B10.87.97[d] <- pow(N[d,33] / N[d,23], 1/(33-23))
    B10.86.96[d] <- pow(N[d,32] / N[d,22], 1/(32-22))
    B10.85.95[d] <- pow(N[d,31] / N[d,21], 1/(31-21))
    B10.84.94[d] <- pow(N[d,30] / N[d,20], 1/(30-20))
    B10.83.93[d] <- pow(N[d,29] / N[d,19], 1/(29-19))
    B10.82.92[d] <- pow(N[d,28] / N[d,18], 1/(28-18))
    B10.81.91[d] <- pow(N[d,27] / N[d,17], 1/(27-17))
    B10.80.90[d] <- pow(N[d,26] / N[d,16], 1/(26-16))
    B10.79.89[d] <- pow(N[d,25] / N[d,15], 1/(25-15))
    B10.78.88[d] <- pow(N[d,24] / N[d,14], 1/(24-14))
    B10.77.87[d] <- pow(N[d,23] / N[d,13], 1/(23-13))
    B10.76.86[d] <- pow(N[d,22] / N[d,12], 1/(22-12))
    B10.75.85[d] <- pow(N[d,21] / N[d,11], 1/(21-11))
    B10.74.84[d] <- pow(N[d,20] / N[d,10], 1/(20-10))
    B10.73.83[d] <- pow(N[d,19] / N[d,9], 1/(19-9))
    B10.72.82[d] <- pow(N[d,18] / N[d,8], 1/(18-8))
    B10.71.81[d] <- pow(N[d,17] / N[d,7], 1/(17-7))
    B10.70.80[d] <- pow(N[d,16] / N[d,6], 1/(16-6))
    B10.69.79[d] <- pow(N[d,15] / N[d,5], 1/(15-5))
    B10.68.78[d] <- pow(N[d,14] / N[d,4], 1/(14-4))
    B10.67.77[d] <- pow(N[d,13] / N[d,3], 1/(13-3))
    B10.66.76[d] <- pow(N[d,12] / N[d,2], 1/(12-2))
    B10.65.75[d] <- pow(N[d,11] / N[d,1], 1/(11-1))
  }
}


















