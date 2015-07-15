

setwd <- outpDir
# read in the results

theFiles <- list.files(outpDir)

# load model C
these <- grep("Model C", theFiles , ignore.case=FALSE, fixed=TRUE)
zoneString <- unlist(lapply(strsplit(theFiles,".",fixed=TRUE),function(x) strsplit(x,".",fixed=TRUE)[[1]][1]))
zones <- as.numeric(substr(zoneString,nchar(zoneString),nchar(zoneString)))
loadThese <- paste0(outpDir,"/",theFiles[these])
nModelC <- length(loadThese)

estsC <- NULL
for(i in 1:nModelC){
  load(loadThese[i])                
  assign(paste0("estsC",i),data.frame(modCode="C",model="Ind MZone, All Zeros",Zone=zones[these][i],Parameter=rownames(bayes$summary),bayes$summary))
  estsC <- rbind(estsC,get(paste0("estsC",i)))
  rm(data, inits, freqs, bayes, CYear, parameters)
}



# load model F
these <- grep("Model F", theFiles , ignore.case=FALSE, fixed=TRUE)
zoneString <- unlist(lapply(strsplit(theFiles,".",fixed=TRUE),function(x) strsplit(x,".",fixed=TRUE)[[1]][1]))
zones <- as.numeric(substr(zoneString,nchar(zoneString),nchar(zoneString)))
loadThese <- paste0(outpDir,"/",theFiles[these])
nModelF <- length(loadThese)

estsF <- NULL
for(i in 1:nModelF){
  load(loadThese[i])                
  assign(paste0("estsF",i),data.frame(modCode="F",model="Ind MZone, 1st Zeros",Zone=zones[these][i],Parameter=rownames(bayes$summary),bayes$summary))
  estsF <- rbind(estsF,get(paste0("estsF",i)))
  rm(data, inits, freqs, bayes, CYear, parameters)
}