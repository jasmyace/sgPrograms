summarizeModels <- function(modelLetter){
  
  # load model 
  these <- grep(paste0("Model ",modelLetter), theFiles , ignore.case=FALSE, fixed=TRUE)
  zoneString <- unlist(lapply(strsplit(theFiles,".",fixed=TRUE),function(x) strsplit(x,".",fixed=TRUE)[[1]][1]))
  zones <- as.numeric(substr(zoneString,nchar(zoneString),nchar(zoneString)))
  loadThese <- paste0(outpDir,"/",theFiles[these])
  assign(paste0("nModel",modelLetter),length(loadThese))
  
  # compile data
  assign(paste0("ests",modelLetter),NULL)
  for(i in 1:get(paste0("nModel",modelLetter))){
    load(loadThese[i])   
    rm(data, inits, freqs, CYear, parameters)
    
    if(modelLetter == 'C'){
      df <- data.frame(modCode=modelLetter,Model="Ind MZone",Zeros='All Zeros',Cut='All Leks')
    } else if(modelLetter == 'D'){
      df <- data.frame(modCode=modelLetter,Model="Ind MZone",Zeros='1st Zeros',Cut='Core Leks')      
    } else if(modelLetter == 'F'){
      df <- data.frame(modCode=modelLetter,Model="Ind MZone",Zeros='1st Zeros',Cut='All Leks')      
    }
    assign(paste0("ests",modelLetter,i),data.frame(df,Zone=zones[these][i],Parameter=rownames(bayes$summary),bayes$summary))
    
    
    get5.95 <- function(metric){
      if(is.matrix(bayes$sims.list[[metric]])){
        the.5.95 <- data.frame(X5.=apply(bayes$sims.list[[metric]],2,function(x) quantile(x,0.05)),X95.=apply(bayes$sims.list[[metric]],2,function(x) quantile(x,0.95)))
      } else {
        the.5.95 <- data.frame(X5.=quantile(bayes$sims.list[[metric]],0.05),X95.=quantile(bayes$sims.list[[metric]],0.95))
      }
      the.5.95
    }
    
    the.a        <- data.frame(df,Zone=zones[these][i],Parameter=paste0('a[',rownames(get5.95('a')),']'),get5.95('a'))
    the.b        <- data.frame(df,Zone=zones[these][i],Parameter=paste0('b[',rownames(get5.95('b')),']'),get5.95('b'))
    the.taunoise <- data.frame(df,Zone=zones[these][i],Parameter='taunoise',get5.95('taunoise'))
    the.mu.a     <- data.frame(df,Zone=zones[these][i],Parameter='mu.a',get5.95('mu.a'))
    the.mu.b     <- data.frame(df,Zone=zones[these][i],Parameter='mu.b',get5.95('mu.b'))
    the.sigma.a  <- data.frame(df,Zone=zones[these][i],Parameter='sigma.a',get5.95('sigma.a'))
    the.sigma.b  <- data.frame(df,Zone=zones[these][i],Parameter='sigma.b',get5.95('sigma.b'))
    the.rho      <- data.frame(df,Zone=zones[these][i],Parameter='rho',get5.95('rho'))
    the.deviance <- data.frame(df,Zone=zones[these][i],Parameter='deviance',get5.95('deviance'))

    buildN <- data.frame(sorter=c(1:nrow(bayes$summary)))
    build0 <- cbind(rbind(the.a,the.b,the.taunoise,the.mu.a,the.mu.b,the.sigma.a,the.sigma.b,the.rho,the.deviance),buildN)
    build0$modCode <- NULL
    build1 <- merge(get(paste0("ests",modelLetter,i)),build0,by=c('Model','Zeros','Cut','Zone','Parameter'),all.x=TRUE)
    build2 <- build1[,c('Model','Zeros','Cut','Zone','Parameter','modCode','mean','sd','X2.5.','X5.','X25.','X50.','X75.','X95.','X97.5.','Rhat','n.eff','sorter')]
    build3 <- build2[order(build2$sorter),]
    
    assign(paste0("ests",modelLetter),rbind(get(paste0("ests",modelLetter)),build3))
    rm(the.a,the.b,the.taunoise,the.mu.a,the.mu.b,the.sigma.a,the.sigma.b,the.rho,the.deviance,buildN,build0,build1,build2,build3,bayes)
  }
  
  get(paste0("ests",modelLetter))
  
}