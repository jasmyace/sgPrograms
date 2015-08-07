make90pCredInt <- function(bayes){
  
  blist <- bayes$sims.list
  bsums <- bayes$summary
  
#   # correcting duplicate multiple mzones for states -- duplicate info
#   dups <- substr(rownames(bsums),nchar(rownames(bsums))-2,nchar(rownames(bsums))) %in% c('[2]','[3]','[4]') & !(substr(rownames(bsums),1,1) %in% c('a','b')) 
#   bsums <- bsums[!dups,]
#   
#   # correcting duplicate multiple mzones for states -- duplicate info
#   dups <- substr(rownames(bsums),1,3) %in% c('N[2','N[3','N[4') | substr(rownames(bsums),1,6) %in% c('beta[2','beta[3','beta[4')
#   bsums <- bsums[!dups,]
  
  get5.95 <- function(metric){
    if(is.matrix(blist[[metric]])){         # matrix
      the.5.95 <- data.frame(X5.=apply(blist[[metric]],2,function(x) quantile(x,0.05)),X95.=apply(blist[[metric]],2,function(x) quantile(x,0.95)))
    } else if(is.array(blist[[metric]])){   # array
      the.5.95 <- data.frame(X5.=apply(blist[[metric]][,1,],2,function(x) quantile(x,0.05)),X95.=apply(blist[[metric]][,1,],2,function(x) quantile(x,0.95)))
    } else {                                # vector
      the.5.95 <- data.frame(X5.=quantile(blist[[metric]],0.05),X95.=quantile(blist[[metric]],0.95))      
    }
    the.5.95
  }
  
  ans5.95 <- NULL
  for(x in 1:length(bayes$sims.list)){
    metric <- attributes(bayes$sims.list[x])$names[1]
    ans <- get5.95(metric)
    ans5.95 <- rbind(ans5.95,ans)
  }
  
  newAns <- cbind(bsums,ans5.95)
  newAns <- newAns[,c('mean','sd','2.5%','X5.','25%','50%','75%','X95.','97.5%')]
  names(newAns)[names(newAns) == 'X5.'] <- '5%'
  names(newAns)[names(newAns) == 'X95.'] <- '95%'
  
  newAns
}