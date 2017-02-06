mclapply.win <- function(X, FUN, mc.cores, ...) {
  cl <- makeCluster(mc.cores)
  res <- parLapply(cl, X, FUN, ...)
  stopCluster(cl)
  return(res)
}
