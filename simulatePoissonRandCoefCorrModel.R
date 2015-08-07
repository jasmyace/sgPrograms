
require(MASS)
require(lattice)

n.groups <- 8
n.years <- 51
n <- n.groups * n.years
pop <- gl(n=n.groups,k=n.years)

original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/51

Xmat <- model.matrix(~pop*year - 1 - year)
print(Xmat[1:91,], dig=2)

intercept.mean <- 2.5
intercept.sd <- 0.84
slope.mean <- -1
slope.sd <- 0.25
rho <- -0.4
sdnoise <- 0.7
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

fakeData <- data.frame(Peak_Males=C,lekCls=pop,yearCls=original.year)
save(fakeData,'fakeData.')


library("lattice")
xyplot(C ~ original.year | pop, ylab="Peak Male counts",xlab="Year")