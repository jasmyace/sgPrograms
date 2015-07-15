
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

#-----------------------------------------------------------------------------

# for(i in 1:nLeks) {
#   for(t in 1:nYears) {
# 
#     n[i,t] <- exp(a[i] + b[i] * (t - medYear) + 0.5*sdnoise*sdnoise)
# 
#     ##### EXPAND WEST DENSITY TO POPULATION SIZE FOR STRATA #####
#     N_W[i,t] <- nDen[i,t]*areaweight[i] # (l) abundance by BCR and year
#     # ##### DERIVE POPULATION INDICIES FROM THE BBS SURVEY #####
#     # n[i,t] <- nonzeroweight[i]*exp(strata[i]+beta[i]*(t-fixedyear)+yreff[t,i]+ 0.5*sdnoise*sdnoise+0.5*sdobs*sdobs)
#     # N[i,t] <- areaweight[i]*n[i,t]/totareaweight
#   }
#   # compadj[i]<-sum(nDen[i,((nyears-nyears_W)+1):nyears])/sum(n[i,((nyears-nyears_W)+1):nyears])
#   # (l) adjustment for BBS data
#   # #---------------------------------------------------------------------------------------------------------#
#   # ##### Calculate the adjusted time series for the BBS data.  Prior to year 40, all we had was BBS data, so the composite estimate is the adjusted BBS data. #####
#   # for( t in 1 : (nyears-nyears_W) ) {
#   # n[i,t] <- nonzeroweight[i]*exp(strata[i]+beta[i]*(t-fixedyear)+yreff[t,i]+ 0.5*sdnoise*sdnoise+0.5*sdobs*sdobs)
#   # N[i,t] <- areaweight[i]*n[i,t]/totareaweight
#   # nadj[i,t] <- n[i,t]*compadj[i]
#   # Ncomp[i,t] <- areaweight[i]*nadj[i,t]
#   # }
#   # #---------------------------------------------------------------------------------------------------------#
#   # ##### Calculate the adjusted BBS estimates for the overlap years (years 40-45), and then calculate the composite estimate as the mean of the adjusted BBS and the West data #####
#   # for( t in ((nyears-nyears_W)+1) : nyears ) {
#   # nadj[i,t] <- n[i,t]*compadj[i]
#   # Ncomp[i,t] <- areaweight[i]*((nadj[i,t]+nDen[i,t])/2)
#   # }
# }
# # compadjave <- sum(compadj[1:nOVstrata])/4
# 
# # ##########################################################
# # ##### SUMMARY STATISTICS  FOR NON OVERLAP BCR ######
# # ##########################################################
# # for( i in (nOVstrata+1) : nstrata) {
# # for( t in 1 : nyears ) {
# # n[i,t] <- nonzeroweight[i]*exp(strata[i]+beta[i]*(t-fixedyear)+yreff[t,i]+ 0.5*sdnoise*sdnoise+0.5*sdobs*sdobs)
# # N[i,t] <- areaweight[i]*n[i,t]/totareaweight
# # nadj[i,t] <- n[i,t]*compadjave
# # Ncomp[i,t] <- areaweight[i]*nadj[i,t]
# # }
# # }
# 
# 
# # ################################################################################
# # ######################### ESTIMATE TREND STATISTICS ##########################
# # ################################################################################
# # ##### STRATA-SPECIFIC OVERALL TRENDS #####
# # for( i in 1 : nstrata ) {
# # B[i] <- pow(Ncomp[i,nyears]/Ncomp[i,2],1/(nyears-2))
# # }
# # ##### TOTAL ANNUAL INDICES #####
# # for( t in 1 : nyears ) {
# # CompIndex[t] <- sum(Ncomp[1:12,t])
# # BBSComInd[t] <- sum(N[1:12,t])
# # }
# # ##### TOTAL ANNUAL INDICES OVERLAP STRATA #####
# # for( t in ((nyears-nyears_W)+1) : nyears ) {
# # CompIndOV[t] <- sum(Ncomp[1:4,t])
# # }
# # ##### OVERALL TREND FROM OVERLAP BCR?S AND ALL BCR'S #####
# # BbarOV <- pow(CompIndOV[nyears]/CompIndOV[((nyears-nyears_W)+1)],1/(nyears_W-1))
# # Bbar <- pow(CompIndex[nyears]/CompIndex[2],1/(nyears-2))
# 
# # (l) To obtain an overall trend (by BCR) for all eagles:
# # (l) The following is based on equation at bottom of page 1441 in Millsap et al. 2013:
# for (i in 1:nOVstrata) {
#   beta[i] <- pow(N_W[i,nyears_W] / N_W[i,1], 1/(nyears_W-1))
# }
# # (l) for study area
# for (i in 1:nyears_W) {
#   N_W0[i] <- sum(N_W[1:nOVstrata,i])
# }
# beta0 <- pow(N_W0[nyears_W] / N_W0[1], 1/(nyears_W-1))
# }