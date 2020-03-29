
mars <- function(stat, geno, subsize = 50, causalCount = 2, NCP = 5.7, gamma = 0.01, 
                 simNum = 1000, fast = TRUE, threshold = 5e-6, setseed = 1, cores = 1){
  
  start <- Sys.time()
  
  stopifnot(nrow(stat) == nrow(geno))
  
  cat("[",format(Sys.time()),"]"," - calculating MARS alt LRT\n",sep = "")
  alt <- computeLRT(stat, geno, subsize, causalCount, NCP, gamma)
  
  cat("[",format(Sys.time()),"]"," - make null data\n",sep = "")
  nulldat <- nullmars(geno, simNum = simNum, setseed = setseed, fast = fast)
  
  cat("[",format(Sys.time()),"]"," - calculating MARS null LRT\n",sep = "")
  if(.Platform$OS.type == "windows") {
    null <- data.table::rbindlist(pbapply::pblapply(1:simNum, function(i){
      computeLRT(nulldat[[2]][,i, drop=F], geno)
    }))
  }else{
    if(cores == 1){
      null <- data.table::rbindlist(pbapply::pblapply(1:simNum, function(i){
        computeLRT(nulldat[[2]][,i, drop=F], geno)
      }))
    }else{
      null <- data.table::rbindlist(pbmcapply::pbmclapply(1:simNum, function(i){
        computeLRT(nulldat[[2]][,i, drop=F], geno)
      },mc.cores = cores))    
    }
  }
  
  nullLRT = as.matrix(null[,1])
  nullUNI = as.matrix(null[,2])
  altLRT = alt[[1]]
  altUNI = alt[[2]]
  
  if(fast){
    w = nulldat[[1]]
    ### pthreshold ###
    sum = 0;
    wsum =  0;
    for(i in c(1:simNum)){
      if(nullUNI[i]<threshold) wsum=wsum+w[i]
      sum=sum+w[i]
    }
    pvalue_threshold=wsum/sum;
    ### UNI p pvalue ###
    sum = 0;
    wsum =  0;
    for(i in c(1:simNum)){
      if(nullUNI[i]<altUNI) wsum=wsum+w[i]
      sum=sum+w[i]
    }
    UNI_pvalue = wsum/sum
    ### UNI LRT value ###	
    sum = 0;
    wsum =  0;
    for(i in c(1:simNum)){
      if(nullLRT[i]>altLRT) wsum=wsum+w[i]
      sum=sum+w[i]
    }
    LRT_pvalue = wsum/sum
    res <- data.frame("LRT_pvalue" = LRT_pvalue, "univariate_pvalue" = UNI_pvalue, "threshold_pvalue" = pvalue_threshold, "threshold_UNI" = threshold,
                      "significance_LRT" = (pvalue_threshold > LRT_pvalue), "significance_UNI" = (pvalue_threshold > UNI_pvalue))    
  }else{
    LRT_pvalue <- length(which(abs(nullLRT)>abs(altLRT)))/simNum
    UNI_pvalue <- length(which(nullUNI<altUNI))/simNum
    quantile <- length(which(nullUNI<threshold))
    LRTthreshold <- sort(nullLRT)[simNum-quantile+1]
    res <- data.frame("LRT_score" = altLRT, "univariate_pvalue" = altUNI, "threshold_LRT" = LRTthreshold, "threshold_UNI" = threshold,
                      "significance_LRT" = (altLRT > LRTthreshold), "significance_UNI" = (altUNI < threshold))
  }
  timdiff <- unlist(strsplit(format(Sys.time() - start)," "))
  
  cat("Total analysis time:", round(as.numeric(timdiff[1]),2),timdiff[2]," \n",sep = " ")
  final_res <- list(alt = alt, results = res)
  return(final_res)
}
