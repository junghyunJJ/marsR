
colmat2listR <- function(dat){
  res <- lapply(1:ncol(dat),function(i){
    as.matrix(dat[,i])
  })
  return(res)
}

mat2listR_marsR <- function(mat){
  list_length <- ncol(mat)
  out_list <- vector("list", list_length)
  for(i in 1:list_length) out_list[[i]] <- mat[,i]
  out_list
}

nCrR_marsR <- function(n, r, v = 1:n) {
  if (r == 0) 
    warning("please set r option")
  else if (r == 1) 
    matrix(v, n, 1)
  else if (r == n) 
    matrix(v, 1, n)
  else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n - 1, r, v[-1]))
}

combnR_marsR <- function(n, k){
  raw_tot_causalIndex <- list()
  for( i in 1:k) raw_tot_causalIndex <- append(raw_tot_causalIndex, mat2listR_marsR(t(nCrR_marsR(n,i))))
  tot_causalIndex <- lapply(raw_tot_causalIndex,function(i){i-1})
  return(tot_causalIndex)
}


pbpmcapply_win <- function(dat, cl, X, FUN, ...) {
  # suppressMessages(library(foreach))
  # suppressMessages(library(parallel))
  
  TIME <- Sys.time()
  cl <- parallel::makeCluster(cl)
  parallel::clusterExport(cl, dat)
  pb <- txtProgressBar(max=length(X), style = 3)
  on.exit(close(pb))
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  res <- foreach(i=X, .combine='rbind', .options.snow=opts) %dopar% {FUN(i, ...)}
  timediff <- Sys.time() - TIME
  cat("\n")
  print(timediff)
  stopCluster(cl)
  return(res)
}
