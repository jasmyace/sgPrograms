makeSummary <- function(paramList,theZones){
  
  b <- length(paramList)        # get number of models to summarize
  
  summary <- NULL
  
  for(i in 1:b){
    if(i == 1){
      summary <- data.frame(tail(paramList[[1]]$summary,7),model=theZones[1])
    } else {
      preSummary <- data.frame(tail(paramList[[i]]$summary,7),model=theZones[i])
      summary <- rbind(summary,preSummary)
    }
  }
  
  summary$param <- ifelse(is.na(as.numeric(substr(rownames(summary),nchar(rownames(summary)),nchar(rownames(summary))))),substr(rownames(summary),1,nchar(rownames(summary))),substr(rownames(summary),1,nchar(rownames(summary)) - 1))
  summary$mZoneNum <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(summary$model))," ",fixed=TRUE))[seq(2,length(unlist(strsplit(as.character(droplevels(summary$model))," ",fixed=TRUE))),2)]))
  summary <- summary[order(summary$param,summary$mZoneNum),]
  
  summary$exp.mean <- ifelse(summary$param %in% c('mu.a','mu.b'),round(exp(summary$mean),4),'N/A')
  rownames(summary) <- NULL
  summary <- summary[,c('param','mZoneNum','model','mean','exp.mean','sd','X2.5.','X25.','X50.','X75.','X97.5.','Rhat','n.eff')]
  
  summary
}