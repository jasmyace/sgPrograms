
      
    readInAnalyticFiles(analDir)
    
    colVec <- brewer.pal(12,"Paired")
    colVec <- colVec[c(1,3,5,7,9,2,4,6,8,10,12)]
    

    load("//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefCorrData5Leks.RData")
    load("//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefUncorrData31Leks.RData")
    load("//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefUncorrData5Leks.RData")
    
    
    
    load("//LAR-FILE-SRV/Data/Jason/sage grouse/Output/simPoisRandCoefCorrData31Leks.RData")    
    

    
    library("lattice")
    xyplot(simPoisRandCoefCorrData31Leks$Peak_Males ~ simPoisRandCoefCorrData31Leks$Year | simPoisRandCoefCorrData31Leks$Lek_ID, ylab="Peak Male counts",xlab="Year")
    
    
  
    
    
      
compileResults("Model D MZone 6 Test All-An-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test All-Bn-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test All-Cn-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test All-Dn-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test All-En-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test All-Fn-4k-76k-1" ,datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}

compileResults("Model D MZone 6 Test 1st-An-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test 1st-Bn-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test 1st-Cn-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test 1st-Dn-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test 1st-En-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test 1st-Fn-4k-76k-1" ,dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 

compileResults("Model D MZone 6 Test All-An-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test All-Bn-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test All-Cn-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test All-Dn-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test All-En-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test All-Fn-4k-304k-4",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}

compileResults("Model D MZone 6 Test 1st-An-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test 1st-Bn-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test 1st-Cn-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test 1st-Dn-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test 1st-En-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test 1st-Fn-4k-304k-4",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 

compileResults("Model D MZone 6 Test All-An-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test All-Bn-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test All-Cn-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test All-Dn-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test All-En-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test All-Fn-4k-316k-1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}

compileResults("Model D MZone 6 Test 1st-An-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test 1st-Bn-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test 1st-Cn-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test 1st-Dn-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test 1st-En-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test 1st-Fn-4k-316k-1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
    
    
    

    
    
    
    

compileResults("Model D MZone 6 Test simCorr31Lek-An-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr31Lek-Bn-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simCorr31Lek-Cn-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simCorr31Lek-Dn-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr31Lek-En-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr31Lek-Fn-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simCorr5Lek-An-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr5Lek-Bn-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test simCorr5Lek-Cn-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simCorr5Lek-Dn-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr5Lek-En-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr5Lek-Fn-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
compileResults("Model D MZone 6 Test simUncorr31Lek-An-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr31Lek-Bn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simUncorr31Lek-Cn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simUncorr31Lek-Dn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr31Lek-En-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr31Lek-Fn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simUncorr5Lek-An-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr5Lek-Bn-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test simUncorr5Lek-Cn-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simUncorr5Lek-Dn-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr5Lek-En-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr5Lek-Fn-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
    
    
    
    
    

    
compileResults("Model D MZone 6 Test simCorr31Lek-An-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr31Lek-Bn-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simCorr31Lek-Cn-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simCorr31Lek-Dn-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr31Lek-En-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr31Lek-Fn-4k-304k-4",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simCorr5Lek-An-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr5Lek-Bn-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test simCorr5Lek-Cn-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simCorr5Lek-Dn-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr5Lek-En-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr5Lek-Fn-4k-304k-4",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
compileResults("Model D MZone 6 Test simUncorr31Lek-An-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr31Lek-Bn-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simUncorr31Lek-Cn-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simUncorr31Lek-Dn-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr31Lek-En-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr31Lek-Fn-4k-304k-4",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simUncorr5Lek-An-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr5Lek-Bn-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test simUncorr5Lek-Cn-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simUncorr5Lek-Dn-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr5Lek-En-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr5Lek-Fn-4k-304k-4",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
    
    

    
    
    
    
    
compileResults("Model D MZone 6 Test simCorr31Lek-An-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr31Lek-Bn-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simCorr31Lek-Cn-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simCorr31Lek-Dn-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr31Lek-En-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr31Lek-Fn-4k-316k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6',"Core")  # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simCorr5Lek-An-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simCorr5Lek-Bn-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")  # Try - fix mu_a
compileResults("Model D MZone 6 Test simCorr5Lek-Cn-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simCorr5Lek-Dn-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simCorr5Lek-En-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simCorr5Lek-Fn-4k-316k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
compileResults("Model D MZone 6 Test simUncorr31Lek-An-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr31Lek-Bn-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_a
compileResults("Model D MZone 6 Test simUncorr31Lek-Cn-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix mu_b
compileResults("Model D MZone 6 Test simUncorr31Lek-Dn-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr31Lek-En-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr31Lek-Fn-4k-316k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6',"Core")   # all zeros - fix 0.5*sigma^2_{\epsilon}
    
compileResults("Model D MZone 6 Test simUncorr5Lek-An-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop 0.5*sigma^2_{\epsilon}
compileResults("Model D MZone 6 Test simUncorr5Lek-Bn-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_a
compileResults("Model D MZone 6 Test simUncorr5Lek-Cn-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix mu_b
compileResults("Model D MZone 6 Test simUncorr5Lek-Dn-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a, mu_b correlation
compileResults("Model D MZone 6 Test simUncorr5Lek-En-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - drop mu_a altogether (regression through the origin)
compileResults("Model D MZone 6 Test simUncorr5Lek-Fn-4k-316k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6',"Core")   # Try - fix 0.5*sigma^2_{\epsilon} 
    
    
    
    
compileResults("Model D MZone 6 Test simCorr31Lek-Gn-4k-76k-1",simPoisRandCoefCorrData31Leks,'Management Zone','mZone 6','Core')
compileResults("Model D MZone 6 Test simCorr5Lek-Gn-4k-76k-1",simPoisRandCoefCorrData5Leks,'Management Zone','mZone 6','Core')    
compileResults("Model D MZone 6 Test simUncorr31Lek-Gn-4k-76k-1",simPoisRandCoefUncorrData31Leks,'Management Zone','mZone 6','Core')    
compileResults("Model D MZone 6 Test simUncorr5Lek-Gn-4k-76k-1",simPoisRandCoefUncorrData5Leks,'Management Zone','mZone 6','Core')    
    
compileResults("Model D MZone 1 Test simCorr877Lek-Gn-4k-76k-1",simMZ1PoisRandCoefCorrData877Leks,'Management Zone','mZone 1','Core')    
compileResults("Model D MZone 6 Try 1",dat1stZerosCore[[6]],'Management Zone','mZone 6','Core')
compileResults("Model D MZone 1 Try 1",dat1stZerosCore[[1]],'Management Zone','mZone 1','Core')
    
    


# 1st Zeros analysis 
    
compileResults("Model D MZone 1 Try 1",dat1stZerosCore[[1]],'Management Zone','mZone 1',"Core")
compileResults("Model E MZone 1 Try 1",dat1stZerosNoco[[1]],'Management Zone','mZone 1',"Periphery")
compileResults("Model F MZone 1 Try 1",dat1stZerosLeks[[1]],'Management Zone','mZone 1',"All Leks")
    
compileResults("Model D MZone 3 Try 1",dat1stZerosCore[[3]],'Management Zone','mZone 3',"Core")
compileResults("Model E MZone 3 Try 1",dat1stZerosNoco[[3]],'Management Zone','mZone 3',"Periphery")
compileResults("Model F MZone 3 Try 1",dat1stZerosLeks[[3]],'Management Zone','mZone 3',"All Leks")
    
compileResults("Model D MZone 4 Try 1",dat1stZerosCore[[4]],'Management Zone','mZone 4',"Core")
compileResults("Model E MZone 4 Try 1",dat1stZerosNoco[[4]],'Management Zone','mZone 4',"Periphery")
compileResults("Model F MZone 4 Try 1",dat1stZerosLeks[[4]],'Management Zone','mZone 4',"All Leks")
    
compileResults("Model D MZone 5 Try 1",dat1stZerosCore[[5]],'Management Zone','mZone 5',"Core")
compileResults("Model E MZone 5 Try 1",dat1stZerosNoco[[5]],'Management Zone','mZone 5',"Periphery")
compileResults("Model F MZone 5 Try 1",dat1stZerosLeks[[5]],'Management Zone','mZone 5',"All Leks")
    
compileResults("Model D MZone 6 Try 1",dat1stZerosCore[[6]],'Management Zone','mZone 6',"Core")
compileResults("Model E MZone 6 Try 1",dat1stZerosNoco[[6]],'Management Zone','mZone 6',"Periphery")
compileResults("Model F MZone 6 Try 1",dat1stZerosLeks[[6]],'Management Zone','mZone 6',"All Leks")
    
compileResults("Model D MZone 8 Try 1",dat1stZerosCore[[8]],'Management Zone','mZone 8',"Core")
compileResults("Model E MZone 8 Try 1",dat1stZerosNoco[[8]],'Management Zone','mZone 8',"Periphery")
compileResults("Model F MZone 8 Try 1",dat1stZerosLeks[[8]],'Management Zone','mZone 8',"All Leks")   
    
    
    
    
    
# all zeros analysis 
    
compileResults("Model D MZone 1 All Zeros 1",datAllZerosCore[[1]],'Management Zone','mZone 1',"Core")
compileResults("Model E MZone 1 All Zeros 1",datAllZerosNoco[[1]],'Management Zone','mZone 1',"Periphery")
compileResults("Model F MZone 1 All Zeros 1",datAllZerosLeks[[1]],'Management Zone','mZone 1',"All Leks")
    
compileResults("Model D MZone 3 All Zeros 1",datAllZerosCore[[3]],'Management Zone','mZone 3',"Core")
compileResults("Model E MZone 3 All Zeros 1",datAllZerosNoco[[3]],'Management Zone','mZone 3',"Periphery")
compileResults("Model F MZone 3 All Zeros 1",datAllZerosLeks[[3]],'Management Zone','mZone 3',"All Leks")

compileResults("Model D MZone 4 All Zeros 1",datAllZerosCore[[4]],'Management Zone','mZone 4',"Core")
compileResults("Model E MZone 4 All Zeros 1",datAllZerosNoco[[4]],'Management Zone','mZone 4',"Periphery")
compileResults("Model F MZone 4 All Zeros 1",datAllZerosLeks[[4]],'Management Zone','mZone 4',"All Leks")

compileResults("Model D MZone 5 All Zeros 1",datAllZerosCore[[5]],'Management Zone','mZone 5',"Core")
compileResults("Model E MZone 5 All Zeros 1",datAllZerosNoco[[5]],'Management Zone','mZone 5',"Periphery")
compileResults("Model F MZone 5 All Zeros 1",datAllZerosLeks[[5]],'Management Zone','mZone 5',"All Leks")

compileResults("Model D MZone 6 All Zeros 1",datAllZerosCore[[6]],'Management Zone','mZone 6',"Core")
compileResults("Model E MZone 6 All Zeros 1",datAllZerosNoco[[6]],'Management Zone','mZone 6',"Periphery")
compileResults("Model F MZone 6 All Zeros 1",datAllZerosLeks[[6]],'Management Zone','mZone 6',"All Leks")
    
compileResults("Model D MZone 8 All Zeros 1",datAllZerosCore[[8]],'Management Zone','mZone 8',"Core")
compileResults("Model E MZone 8 All Zeros 1",datAllZerosNoco[[8]],'Management Zone','mZone 8',"Periphery")
compileResults("Model F MZone 8 All Zeros 1",datAllZerosLeks[[8]],'Management Zone','mZone 8',"All Leks")
    

compileResults("Model F MZone MT All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'MT',],'State','MT',"All Leks")    
compileResults("Model F MZone CA All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'CA',],'State','CA',"All Leks")  
compileResults("Model F MZone ID All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'ID',],'State','ID',"All Leks")  
compileResults("Model F MZone WA All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'WA',],'State','WA',"All Leks")  
compileResults("Model F MZone WY All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'WY',],'State','WY',"All Leks")  
compileResults("Model F MZone UT All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'UT',],'State','UT',"All Leks")  
compileResults("Model F MZone NV All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'NV',],'State','NV',"All Leks")  
compileResults("Model F MZone OR All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'OR',],'State','OR',"All Leks")  
compileResults("Model F MZone CO All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'CO',],'State','CO',"All Leks")  
compileResults("Model F MZone SD All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'SD',],'State','SD',"All Leks")  
compileResults("Model F MZone ND All Zeros State 1",datAllZerosLeks[[9]][datAllZerosLeks[[9]]$State == 'ND',],'State','ND',"All Leks")  
    
    
    
    
    
    # 10 years analysis
compileResults("Model D MZone 1 All Zeros 2005-2015",datAllZerosCore[[1]][datAllZerosCore[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','Core')
compileResults("Model E MZone 1 All Zeros 2005-2015",datAllZerosNoco[[1]][datAllZerosNoco[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','Periphery')
compileResults("Model F MZone 1 All Zeros 2005-2015",datAllZerosLeks[[1]][datAllZerosLeks[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','All Leks')
compileResults("Model D MZone 1 1st Zeros 2005-2015",dat1stZerosCore[[1]][dat1stZerosCore[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','Core')
compileResults("Model E MZone 1 1st Zeros 2005-2015",dat1stZerosNoco[[1]][dat1stZerosNoco[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','Periphery')
compileResults("Model F MZone 1 1st Zeros 2005-2015",dat1stZerosLeks[[1]][dat1stZerosLeks[[1]]$Year %in% seq(2005,2015),],'Management Zone','mZone 1','All Leks')    

compileResults("Model D MZone 3 All Zeros 2005-2015",datAllZerosCore[[3]][datAllZerosCore[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','Core')
compileResults("Model E MZone 3 All Zeros 2005-2015",datAllZerosNoco[[3]][datAllZerosNoco[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','Periphery')
compileResults("Model F MZone 3 All Zeros 2005-2015",datAllZerosLeks[[3]][datAllZerosLeks[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','All Leks')
compileResults("Model D MZone 3 1st Zeros 2005-2015",dat1stZerosCore[[3]][dat1stZerosCore[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','Core')
compileResults("Model E MZone 3 1st Zeros 2005-2015",dat1stZerosNoco[[3]][dat1stZerosNoco[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','Periphery')
compileResults("Model F MZone 3 1st Zeros 2005-2015",dat1stZerosLeks[[3]][dat1stZerosLeks[[3]]$Year %in% seq(2005,2015),],'Management Zone','mZone 3','All Leks')    

compileResults("Model D MZone 4 All Zeros 2005-2015",datAllZerosCore[[4]][datAllZerosCore[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','Core')
compileResults("Model E MZone 4 All Zeros 2005-2015",datAllZerosNoco[[4]][datAllZerosNoco[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','Periphery')
compileResults("Model F MZone 4 All Zeros 2005-2015",datAllZerosLeks[[4]][datAllZerosLeks[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','All Leks')
compileResults("Model D MZone 4 1st Zeros 2005-2015",dat1stZerosCore[[4]][dat1stZerosCore[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','Core')
compileResults("Model E MZone 4 1st Zeros 2005-2015",dat1stZerosNoco[[4]][dat1stZerosNoco[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','Periphery')
compileResults("Model F MZone 4 1st Zeros 2005-2015",dat1stZerosLeks[[4]][dat1stZerosLeks[[4]]$Year %in% seq(2005,2015),],'Management Zone','mZone 4','All Leks')    

compileResults("Model D MZone 5 All Zeros 2005-2015",datAllZerosCore[[5]][datAllZerosCore[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','Core')
compileResults("Model E MZone 5 All Zeros 2005-2015",datAllZerosNoco[[5]][datAllZerosNoco[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','Periphery')
compileResults("Model F MZone 5 All Zeros 2005-2015",datAllZerosLeks[[5]][datAllZerosLeks[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','All Leks')
compileResults("Model D MZone 5 1st Zeros 2005-2015",dat1stZerosCore[[5]][dat1stZerosCore[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','Core')
compileResults("Model E MZone 5 1st Zeros 2005-2015",dat1stZerosNoco[[5]][dat1stZerosNoco[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','Periphery')
compileResults("Model F MZone 5 1st Zeros 2005-2015",dat1stZerosLeks[[5]][dat1stZerosLeks[[5]]$Year %in% seq(2005,2015),],'Management Zone','mZone 5','All Leks')    

compileResults("Model D MZone 6 All Zeros 2005-2015",datAllZerosCore[[6]][datAllZerosCore[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','Core')
compileResults("Model E MZone 6 All Zeros 2005-2015",datAllZerosNoco[[6]][datAllZerosNoco[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','Periphery')
compileResults("Model F MZone 6 All Zeros 2005-2015",datAllZerosLeks[[6]][datAllZerosLeks[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','All Leks')
compileResults("Model D MZone 6 1st Zeros 2005-2015",dat1stZerosCore[[6]][dat1stZerosCore[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','Core')
compileResults("Model E MZone 6 1st Zeros 2005-2015",dat1stZerosNoco[[6]][dat1stZerosNoco[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','Periphery')
compileResults("Model F MZone 6 1st Zeros 2005-2015",dat1stZerosLeks[[6]][dat1stZerosLeks[[6]]$Year %in% seq(2005,2015),],'Management Zone','mZone 6','All Leks')    
    
compileResults("Model D MZone 8 All Zeros 2005-2015",datAllZerosCore[[8]][datAllZerosCore[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','Core')
compileResults("Model E MZone 8 All Zeros 2005-2015",datAllZerosNoco[[8]][datAllZerosNoco[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','Periphery')
compileResults("Model F MZone 8 All Zeros 2005-2015",datAllZerosLeks[[8]][datAllZerosLeks[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','All Leks')
compileResults("Model D MZone 8 1st Zeros 2005-2015",dat1stZerosCore[[8]][dat1stZerosCore[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','Core')
compileResults("Model E MZone 8 1st Zeros 2005-2015",dat1stZerosNoco[[8]][dat1stZerosNoco[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','Periphery')
compileResults("Model F MZone 8 1st Zeros 2005-2015",dat1stZerosLeks[[8]][dat1stZerosLeks[[8]]$Year %in% seq(2005,2015),],'Management Zone','mZone 8','All Leks')    

    
    
    
    
    
    
    
makeHistogramPlots(datAllZerosCore[[6]],theUnit='Management Zone mZone 6',tracDir=paste0(outpDir,'/','Model D MZone 6 Try 1','/Zeros Plots'),file='Model D MZone 6 Try 1')  
makeHistogramPlots(datAllZerosCore[[1]],theUnit='Management Zone mZone 1',tracDir=paste0(outpDir,'/','Model D MZone 1 Try 1','/Zeros Plots'),file='Model D MZone 1 Try 1')  
    
    
    
    
    compileResults <- function(file,dat,string,theUnit,runType){
      
#       file <- "Model E MZone 3 All Zeros 2005-2015"
#       dat <- datAllZerosNoco[[3]][datAllZerosNoco[[3]]$Year %in% seq(2005,2015),]
#       string <- 'mZone 3'
#       theUnit <- "mZone 3"
#       runType <- 'Periphery'
    
      load(paste0(outpDir,'/',file,".RData"))  
      
      ifelse(!dir.exists(file.path(outpDir,file)), dir.create(file.path(outpDir,file)), FALSE)      # make new folder
      file.copy(paste0(outpDir,"/",file,".RData"),paste0(outpDir,"/",file))                         # copy bayes output to new folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trace Plots'))), FALSE)        # make new trace plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Posterior Plots'))), FALSE)# make new posterior plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Trend Plots'))), FALSE)        # make new mzone trend plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Zeros Plots'))), FALSE)        # make new mzone trend plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Lek Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Lek Plots'))), FALSE)            # make new mzone lek plots folder
      ifelse(!dir.exists(file.path(paste0(outpDir,'/',file,'/Random Lek Plots'))), dir.create(file.path(paste0(outpDir,'/',file,'/Random Lek Plots'))), FALSE)            # make new mzone lek plots folder
      
      # make trace, posterior, histogram, lek-trend, random-lek plots
      nParms <- dim(bayes$sims.array)[3]
      parmList <- dimnames(bayes$sims.array)[[3]]
      #makeTracePlots(nParms,parmList,paste0(outpDir,'/',file,'/Trace Plots'),file,bayes)    
      #makePosteriorPlots(nParms,parmList,paste0(outpDir,'/',file,'/Posterior Plots'),file,bayes)
      #makeHistogramPlots(dat,paste0(" - ",string," ",theUnit),paste0(outpDir,'/',file,'/Zeros Plots'),file)
      #makeLekTrendPlots(dat,paste0(outpDir,'/',file,'/Lek Plots'),bayes)

      bsums90 <- make90pCredInt(bayes)
      makeRandomLekPlots(dat,paste0(outpDir,'/',file,'/Random Lek Plots'),file,bayes,bsums90)     
      
      # make bayes summary file of estimates
      #write.csv(bsums90,paste0(outpDir,'/',file,'/bayesSummary - ',file,'.csv'))
      
      # make trend plots
      #makeTrendPlots(dat,runType,paste0(" - ",string," ",theUnit),1,paste0(outpDir,'/',file,'/Trend Plots'),file,bayes)
      
      rm(nParms,parmList,bayes)
    }  
    
    