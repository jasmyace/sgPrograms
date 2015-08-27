
require(MASS)
require(lattice)



# MANAGEMENT ZONE 6

# make simulated data with correlated random slopes and random intercepts.  31 LEKS
n.groups <- 31
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.45
intercept.sd <- 0.9
slope.mean <- -0.28
slope.sd <- 0.05
rho <- -0.1
sdnoise <- 0.62
intercept.slope.covariance <- intercept.sd*slope.sd*rho


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

# intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
# slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
# all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)

lin.pred <- Xmat[,] %*% all.effects + noise.effects

C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


library("lattice")
xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simPoisRandCoefCorrData31Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)






# make simulated data with correlated random slopes and random intercepts.  5 LEKS
n.groups <- 5
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.45
intercept.sd <- 0.9
slope.mean <- -0.28
slope.sd <- 0.05
rho <- -0.1
sdnoise <- 0.62
intercept.slope.covariance <- intercept.sd*slope.sd*rho


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

# intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
# slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
# all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)

lin.pred <- Xmat[,] %*% all.effects + noise.effects

C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


library("lattice")
xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simPoisRandCoefCorrData5Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)







# make simulated data with UNcorrelated random intercepts and slopes.  5 LEKS

n.groups <- 5
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.45
intercept.sd <- 0.9
slope.mean <- -0.28
slope.sd <- 0.05
intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
all.effects <- c(intercept.effects,slope.effects)


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)
lin.pred <- Xmat[,] %*% all.effects + noise.effects
C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


library("lattice")
xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simPoisRandCoefUncorrData5Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)








# make simulated data with UNcorrelated random intercepts and slopes.  31 LEKS

n.groups <-31
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.45
intercept.sd <- 0.9
slope.mean <- -0.28
slope.sd <- 0.05
intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
all.effects <- c(intercept.effects,slope.effects)


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)
lin.pred <- Xmat[,] %*% all.effects + noise.effects
C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


library("lattice")
xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simPoisRandCoefUncorrData31Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)




save(simPoisRandCoefCorrData31Leks ,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefCorrData31Leks.RData")
save(simPoisRandCoefCorrData5Leks  ,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefCorrData5Leks.RData")
save(simPoisRandCoefUncorrData31Leks,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefUncorrData31Leks.RData")
save(simPoisRandCoefUncorrData5Leks,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefUncorrData5Leks.RData")




























# MANAGEMENT ZONE 1

# make simulated data with correlated random slopes and random intercepts.  877 LEKS
n.groups <- 877
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.72
intercept.sd <- 0.80
slope.mean <- -0.034
slope.sd <- 0.035
rho <- -0.22
sdnoise <- 0.79
intercept.slope.covariance <- intercept.sd*slope.sd*rho


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

# intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
# slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
# all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)

lin.pred <- Xmat[,] %*% all.effects + noise.effects

C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


# library("lattice")
# xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simMZ1PoisRandCoefCorrData877Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)




# 
# 
# # make simulated data with correlated random slopes and random intercepts.  5 LEKS
# n.groups <- 5
# n.years <- 51
# n <- n.groups * n.years
# pop <- gl(n=n.groups,k=n.years)
# 
# original.year <- rep(1:n.years, n.groups)
# year <- (original.year-1)/51
# 
# Xmat <- model.matrix(~pop*year - 1 - year)
# 
# intercept.mean <- 2.72
# intercept.sd <- 0.80
# slope.mean <- -0.034
# slope.sd <- 0.035
# rho <- -0.22
# sdnoise <- 0.79
# intercept.slope.covariance <- intercept.sd*slope.sd*rho
# 
# 
# mu.vector <- c(intercept.mean,slope.mean)
# var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
# effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
# effects
# apply(effects,2,mean)
# var(effects)
# intercept.effects <- effects[,1]
# slope.effects <- effects[,2]
# all.effects <- c(intercept.effects,slope.effects)
# 
# # intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
# # slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
# # all.effects <- c(intercept.effects,slope.effects)
# 
# noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)
# 
# lin.pred <- Xmat[,] %*% all.effects + noise.effects
# 
# C <- rpois(n=n,lambda=exp(lin.pred))
# hist(C[C<100],col="grey")
# 
# pMales <- C
# lek_ID <- pop
# Year <- original.year
# 
# 
# library("lattice")
# xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")
# 
# simMZ1PoisRandCoefCorrData5Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)
# 
# 
# 
# 
# 
# 
# 
# # make simulated data with UNcorrelated random intercepts and slopes.  5 LEKS
# 
# n.groups <- 5
# n.years <- 51
# n <- n.groups * n.years
# pop <- gl(n=n.groups,k=n.years)
# 
# original.year <- rep(1:n.years, n.groups)
# year <- (original.year-1)/51
# 
# Xmat <- model.matrix(~pop*year - 1 - year)
# 
# intercept.mean <- 2.72
# intercept.sd <- 0.80
# slope.mean <- -0.034
# slope.sd <- 0.035
# sdnoise <- 0.79
# intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
# slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
# all.effects <- c(intercept.effects,slope.effects)
# 
# 
# mu.vector <- c(intercept.mean,slope.mean)
# var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
# effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
# effects
# apply(effects,2,mean)
# var(effects)
# intercept.effects <- effects[,1]
# slope.effects <- effects[,2]
# all.effects <- c(intercept.effects,slope.effects)
# 
# noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)
# lin.pred <- Xmat[,] %*% all.effects + noise.effects
# C <- rpois(n=n,lambda=exp(lin.pred))
# hist(C[C<100],col="grey")
# 
# pMales <- C
# lek_ID <- pop
# Year <- original.year
# 
# 
# library("lattice")
# xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")
# 
# simMZ1PoisRandCoefUncorrData5Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)
# 
# 
# 
# 




# make simulated data with UNcorrelated random intercepts and slopes.  877 LEKS

n.groups <-877
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)

intercept.mean <- 2.72
intercept.sd <- 0.80
slope.mean <- -0.034
slope.sd <- 0.035
sdnoise <- 0.79
intercept.effects <- rnorm(n=n.groups,mean=intercept.mean,sd=intercept.sd)
slope.effects <- rnorm(n=n.groups,mean=slope.mean,sd=slope.sd)
all.effects <- c(intercept.effects,slope.effects)


mu.vector <- c(intercept.mean,slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2,intercept.slope.covariance,intercept.slope.covariance, slope.sd^2),2,2)
effects <- mvrnorm(n=n.groups,mu=mu.vector,Sigma=var.cova.matrix)
effects
apply(effects,2,mean)
var(effects)
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
all.effects <- c(intercept.effects,slope.effects)

noise.effects <- rnorm(n=length(pop),mean=0,sd=sdnoise)
lin.pred <- Xmat[,] %*% all.effects + noise.effects
C <- rpois(n=n,lambda=exp(lin.pred))
hist(C[C<100],col="grey")

pMales <- C
lek_ID <- pop
Year <- original.year


# library("lattice")
# xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")

simMZ1PoisRandCoefUncorrData877Leks <- data.frame(Peak_Males=C,Lek_ID=pop,Year=original.year)




save(simMZ1PoisRandCoefCorrData877Leks ,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simMZ1PoisRandCoefCorrData877Leks.RData")
# save(simMZ1PoisRandCoefCorrData5Leks  ,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simMZ1PoisRandCoefCorrData5Leks.RData")
save(simMZ1PoisRandCoefUncorrData877Leks,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simMZ1PoisRandCoefUncorrData877Leks.RData")
# save(simMZ1PoisRandCoefUncorrData5Leks,file="//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simMZ1PoisRandCoefUncorrData5Leks.RData")
