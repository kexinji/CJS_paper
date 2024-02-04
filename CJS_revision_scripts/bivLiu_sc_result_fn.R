#============================================================================================================
getTable <- function(result.path, result.name, s, b) {
  # Extract model parameters from the txt file and put them into a data.frame
  #
  # Args: 
  #   result.path: the directory containing the result file
  #   result.name: the name of the result file
  # 	s: 	number of patameters 
  # 	b: 	number of fixed effect parameters
  #
  # Returns:
  # 	Output the data.frame 
  
  liuResults <- read.table(paste0(result.path, result.name), header = T, sep = ",")
  res <- round(liuResults[,2], 4)

  
  Parameters <- c("betaHat1", "betaHat2", "betaHat3", "betaHat4", "betaHat5", "betaHat6",
                  "tau1", "tau2", "sigma.b1", "sigma.b2", "sigma.b3", "sigma1", "sigma2",  
                  "rho1", "a10", "a11", "a12", "rho2", "a20", "a21", "a22", "rho3")
  
  estim <- c(res[(s+1):(s+b)], res[1:s])
  
  
  sdTheta <- res[which(liuResults[,1] == "sdTheta1"):which(liuResults[,1] == "sdTheta16")]
  sdBeta <- res[which(liuResults[,1] == "sdBeta1"):which(liuResults[,1] == "sdBeta6")]
  se <- c(sdBeta, sdTheta)
  
  sep <- noquote(rep("&", s+b))
  lastSep <- noquote(rep("\\", s+b))
  table <- data.frame(Parameters, sep, estim, sep, se, lastSep)
  
  return(table)
}

#============================================================================================================


