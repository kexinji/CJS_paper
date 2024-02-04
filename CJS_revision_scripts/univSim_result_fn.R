getResult <- function(file.path, true) {
  files <- list.files(path = file.path, pattern = "*.txt", full.names = T)
  myfiles <- do.call(cbind, lapply(files, function(x) read.table(x, header = T, 
                                                                 colClasses = c(rep("NULL", 1), rep("double", 1)), 
                                                                 sep = ",")))
  myfiles.new <- myfiles[, which(is.na(myfiles[37,]) == F)]
  
  res <- rowMeans(myfiles.new[1:94, ])
  
  theta <- res[1:7]
  betaHat <- res[8]
  fHat <- res[9:36]
  varBeta <- res[37]
  varF <- res[38:65]
  biasBeta <- res[66]
  biasF <- res[67:94]
  
  # LATEX table
  Parameters <- c("betaHat", "tau", "sigma.b", "sigma", "rho", "a0", "a1", "a2")
  estim <- c(betaHat, theta)
  relative_bias <- round(abs(true-estim)/true, 4)
  se <- sapply(1:8, function(i) return(sd(myfiles.new[i,])))
  empirical_se <- round(c(se[8], se[1:7]), 4)
  table <- data.frame(Parameters, true, estim, relative_bias, empirical_se)
  
  # param <- c("betaHat", "sdBeta", "biasBeta")
  # table <- data.frame(param, c(betaHat, varBeta, biasBeta))
  
  return(table)
}
