
# random slopes & intercepts model, with fixed management zone effect and slope.
# always specify dist'ns and formulas separately.

{
  for(k in 1:nCounts) {                                  # specify model for each obs in the df
    pMales[k] ~ dpois(lambda[k])                         # count (number of peak males) distribution
    log(lambda[k]) <-                                    # expected number of peaked males
      
      m.adj[zoneCls[k]] +                                # fixed zone effect
      n.adj[zoneCls[k]]*(yearCls[k] - medYear) +         # fixed zone slope
      
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
    m.adj[z] <- m[z] - mu.a
    n.adj[z] <- n[z] - mu.b
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
      N[i,j] <- exp(m.adj[i] +                         # fixed zone effect
                    n.adj[i]*(yearCls[j] - medYear) +  # fixed zone slope                        
                    mu.a +                             # grand intercept
                    mu.b*(yearCls[j] - medYear) +      # grand slope
                    0.5*sdnoise*sdnoise) 
    }
  }
  
  # all-year estimates.
  for (i in 1:nZones) {
    beta[i] <- pow(N[i,nYears] / N[i,1], 1/(nYears-1))
  }
  
  # weight the expected Ns, for all mzones and time periods -- build up to wgt average
  for(i in 1:nZones){
    for(j in 1:nYears){
      preN0[i,j] <- nLeksMat[i,j] * N[i,j]
    }
  }
  
  # sum over each time period, and divide by total number of leks -- gets wgt average
  for(j in 1:nYears){
    N0[j] <- sum(preN0[1:nZones,j]) / nLeksVec[1,j]
  }

  # calculate grand weighted-mean slope
  beta0 <- pow(N0[nYears] / N0[1], 1/(nYears-1)) 
   
  # 10-year Bs
  for (i in 1:nZones) {
    B10.05.15[i] <- pow(N[i,51] / N[i,41], 1/(51-41))
    B10.04.14[i] <- pow(N[i,50] / N[i,40], 1/(50-40))
    B10.03.13[i] <- pow(N[i,49] / N[i,39], 1/(49-39))
    B10.02.12[i] <- pow(N[i,48] / N[i,38], 1/(48-38))
    B10.01.11[i] <- pow(N[i,47] / N[i,37], 1/(47-37))
    B10.00.10[i] <- pow(N[i,46] / N[i,36], 1/(46-36))
    B10.99.09[i] <- pow(N[i,45] / N[i,35], 1/(45-35))
    B10.98.08[i] <- pow(N[i,44] / N[i,34], 1/(44-34))
    B10.97.07[i] <- pow(N[i,43] / N[i,33], 1/(43-33))
    B10.96.06[i] <- pow(N[i,42] / N[i,32], 1/(42-32))
    B10.95.05[i] <- pow(N[i,41] / N[i,31], 1/(41-31))
    B10.94.04[i] <- pow(N[i,40] / N[i,30], 1/(40-30))
    B10.93.03[i] <- pow(N[i,39] / N[i,29], 1/(39-29))
    B10.92.02[i] <- pow(N[i,38] / N[i,28], 1/(38-28))
    B10.91.01[i] <- pow(N[i,37] / N[i,27], 1/(37-27))
    B10.90.00[i] <- pow(N[i,36] / N[i,26], 1/(36-26))
    B10.89.99[i] <- pow(N[i,35] / N[i,25], 1/(35-25))
    B10.88.98[i] <- pow(N[i,34] / N[i,24], 1/(34-24))
    B10.87.97[i] <- pow(N[i,33] / N[i,23], 1/(33-23))
    B10.86.96[i] <- pow(N[i,32] / N[i,22], 1/(32-22))
    B10.85.95[i] <- pow(N[i,31] / N[i,21], 1/(31-21))
    B10.84.94[i] <- pow(N[i,30] / N[i,20], 1/(30-20))
    B10.83.93[i] <- pow(N[i,29] / N[i,19], 1/(29-19))
    B10.82.92[i] <- pow(N[i,28] / N[i,18], 1/(28-18))
    B10.81.91[i] <- pow(N[i,27] / N[i,17], 1/(27-17))
    B10.80.90[i] <- pow(N[i,26] / N[i,16], 1/(26-16))
    B10.79.89[i] <- pow(N[i,25] / N[i,15], 1/(25-15))
    B10.78.88[i] <- pow(N[i,24] / N[i,14], 1/(24-14))
    B10.77.87[i] <- pow(N[i,23] / N[i,13], 1/(23-13))
    B10.76.86[i] <- pow(N[i,22] / N[i,12], 1/(22-12))
    B10.75.85[i] <- pow(N[i,21] / N[i,11], 1/(21-11))
    B10.74.84[i] <- pow(N[i,20] / N[i,10], 1/(20-10))
    B10.73.83[i] <- pow(N[i,19] / N[i,9], 1/(19-9))
    B10.72.82[i] <- pow(N[i,18] / N[i,8], 1/(18-8))
    B10.71.81[i] <- pow(N[i,17] / N[i,7], 1/(17-7))
    B10.70.80[i] <- pow(N[i,16] / N[i,6], 1/(16-6))
    B10.69.79[i] <- pow(N[i,15] / N[i,5], 1/(15-5))
    B10.68.78[i] <- pow(N[i,14] / N[i,4], 1/(14-4))
    B10.67.77[i] <- pow(N[i,13] / N[i,3], 1/(13-3))
    B10.66.76[i] <- pow(N[i,12] / N[i,2], 1/(12-2))
    B10.65.75[i] <- pow(N[i,11] / N[i,1], 1/(11-1))
  }
  BX10.05.15 <- pow(N0[51] / N0[41], 1/(51-41))
  BX10.04.14 <- pow(N0[50] / N0[40], 1/(50-40))
  BX10.03.13 <- pow(N0[49] / N0[39], 1/(49-39))
  BX10.02.12 <- pow(N0[48] / N0[38], 1/(48-38))
  BX10.01.11 <- pow(N0[47] / N0[37], 1/(47-37))
  BX10.00.10 <- pow(N0[46] / N0[36], 1/(46-36))
  BX10.99.09 <- pow(N0[45] / N0[35], 1/(45-35))
  BX10.98.08 <- pow(N0[44] / N0[34], 1/(44-34))
  BX10.97.07 <- pow(N0[43] / N0[33], 1/(43-33))
  BX10.96.06 <- pow(N0[42] / N0[32], 1/(42-32))
  BX10.95.05 <- pow(N0[41] / N0[31], 1/(41-31))
  BX10.94.04 <- pow(N0[40] / N0[30], 1/(40-30))
  BX10.93.03 <- pow(N0[39] / N0[29], 1/(39-29))
  BX10.92.02 <- pow(N0[38] / N0[28], 1/(38-28))
  BX10.91.01 <- pow(N0[37] / N0[27], 1/(37-27))
  BX10.90.00 <- pow(N0[36] / N0[26], 1/(36-26))
  BX10.89.99 <- pow(N0[35] / N0[25], 1/(35-25))
  BX10.88.98 <- pow(N0[34] / N0[24], 1/(34-24))
  BX10.87.97 <- pow(N0[33] / N0[23], 1/(33-23))
  BX10.86.96 <- pow(N0[32] / N0[22], 1/(32-22))
  BX10.85.95 <- pow(N0[31] / N0[21], 1/(31-21))
  BX10.84.94 <- pow(N0[30] / N0[20], 1/(30-20))
  BX10.83.93 <- pow(N0[29] / N0[19], 1/(29-19))
  BX10.82.92 <- pow(N0[28] / N0[18], 1/(28-18))
  BX10.81.91 <- pow(N0[27] / N0[17], 1/(27-17))
  BX10.80.90 <- pow(N0[26] / N0[16], 1/(26-16))
  BX10.79.89 <- pow(N0[25] / N0[15], 1/(25-15))
  BX10.78.88 <- pow(N0[24] / N0[14], 1/(24-14))
  BX10.77.87 <- pow(N0[23] / N0[13], 1/(23-13))
  BX10.76.86 <- pow(N0[22] / N0[12], 1/(22-12))
  BX10.75.85 <- pow(N0[21] / N0[11], 1/(21-11))
  BX10.74.84 <- pow(N0[20] / N0[10], 1/(20-10))
  BX10.73.83 <- pow(N0[19] / N0[9], 1/(19-9))
  BX10.72.82 <- pow(N0[18] / N0[8], 1/(18-8))
  BX10.71.81 <- pow(N0[17] / N0[7], 1/(17-7))
  BX10.70.80 <- pow(N0[16] / N0[6], 1/(16-6))
  BX10.69.79 <- pow(N0[15] / N0[5], 1/(15-5))
  BX10.68.78 <- pow(N0[14] / N0[4], 1/(14-4))
  BX10.67.77 <- pow(N0[13] / N0[3], 1/(13-3))
  BX10.66.76 <- pow(N0[12] / N0[2], 1/(12-2))
  BX10.65.75 <- pow(N0[11] / N0[1], 1/(11-1))
}