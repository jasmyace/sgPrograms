prepareForAnalysis <- function(df){
   
  # combine mzones 2 and 7 into new mzone 8
  zone2or7 <- df[df$Mgmt_zone %in% c('MZ II','MZ VII'),]
  zone2or7$Mgmt_zone <- "MZ VIII"
  df <- rbind(df,zone2or7)
  
  # combine mzones 1, 2, 3, 4, 5, 6, 7, into new mzone 9
  zoneall <- df[df$Mgmt_zone %in% c('MZ I','MZ II','MZ III','MZ IV','MZ V','MZ VI','MZ VII'),]
  zoneall$Mgmt_zone <- "MZ IX"
  df <- rbind(df,zoneall)  
  
  df$mZoneNum <- as.numeric(as.roman(unlist(strsplit(as.character(droplevels(df$Mgmt_zone))," ",fixed=TRUE))[seq(2,length(unlist(strsplit(as.character(droplevels(df$Mgmt_zone))," ",fixed=TRUE))),2)]))
  
  
  # identify strings of zeros
  df$Peak_Males1 <- shift(df$Peak_Males,-1)
  df$Peak_Males2 <- shift(df$Peak_Males,-2)
  df$Lek_IDnum <- as.numeric(df$Lek_ID)
  df$Lek_IDnum1 <- shift(df$Lek_IDnum,-1)
  df$DupZero <- ifelse( ( df$Peak_Males == 0 & df$Peak_Males1 == 0 ) & ( df$Lek_IDnum == df$Lek_IDnum1),1,0)                
  
  # identify strings of 1s (< 2 males for two or more years)
  df$P_M12 <- ifelse(df$Peak_Males < 2 & df$Peak_Males > 0,1,0)
  df$P_M12.1 <- shift(df$P_M12,-1)
  df$P_M12.N1 <- shift(df$P_M12,1)
  # last                                   # first                               # middle
  df$LT2 <- ifelse( df$P_M12 == 1 & ( (df$P_M12.1 == 1 & df$P_M12.N1 == 0) | ( df$P_M12.1 == 0 & df$P_M12.N1 == 1 ) | ( df$P_M12.1 == 1 & df$P_M12.N1 == 1 )),1,0)
  
  df$Peak_Males1 <- df$Peak_Males2 <- df$P_M12 <- df$P_M12.1 <- df$P_M12.N1 <- NULL
  
                                                                            # pre-combo of mzones 2 and 7
  length(unique(df$Lek_ID))                                                 #  6,594
  dim(df)[1]                                                                # 92,031
  df.noZero.0 <- df[df$DupZero == 0,]                                       # 75,726
  
  df.1 <- df[df$Year <= 1989 & df$DupZero == 0 ,]            
  length(unique(df.1$Lek_ID))                                               #  3,180
  dim(df.1)[1]                                                              # 18,211
  
  df.2 <- df[df$Year >= 1990 & df$DupZero == 0 ,]            
  length(unique(df.2$Lek_ID))                                               #  6,349
  dim(df.2)[1]                                                              # 57,515
  
  length(unique(df[df$Year == 1965 & df$Peak_Males != 0,]$Lek_ID))          #    267
  length(unique(df[df$Year == 2014 & df$Peak_Males != 0,]$Lek_ID))          #  3,034
  
  df.weird <- df[df$Peak_Males != 0,]
  
  df.list <- list(df,df.1,df.2,df.weird)
}