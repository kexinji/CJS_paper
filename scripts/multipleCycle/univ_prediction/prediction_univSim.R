source("univFns_2c.R")

# n.sim is the number of simulations
# m is the number of subjects in a simulation
prediction_univSim <- function(n.sim, m) {
  # Store pse and mu from simulation
  pse_sim <- matrix(0, 7, n.sim) # record the average pse values for each simulation
  mu_sim <- list() # record the smallest, median, 2nd largest and largest mu's for each simulation 
  pseSort <- matrix(0, 4, n.sim) # record the smallest, median, 2nd largest and largest pse for each simulation 
  
  for(k in 1:n.sim) {
    # print simulation progress
    cat(paste0("Simulation: ", k, "\n"))
    
    t <- (1:56 %% 28)/10
    ij <- 29:56
    ij_1 <- 28:55
    
    theta <- c(rep(1, 3), 0.1, -3, -5, -2)
    alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))
    
    # For each simulation, generate new dataset
    set.seed(k)
    data <- newData(theta, beta=-0.5, m, 28, process="NOU", cycle=2)
    res <- results(data, data$Y, data$age, random = 1, process = "NOU", data$day, data$id, tol = 0.1, cap = 10)
    
    betaHat <- res$betaHat
    fHat <- res$fHat
    
    # Store pse for each subj to a column
    pse_sub <- matrix(0, 7, m) 
    
    # Initializations, each column is a subject
    mu <- matrix(0, 28, m) 
    mu_pc <- matrix(0, 28, m)
    mu_pcv <- matrix(0, 28, m)
    mu_pre <- matrix(0, 28, m)
    
    mu25 <- matrix(0, 28, m)
    mu50 <- matrix(0, 28, m)
    mu75 <- matrix(0, 28, m)
    
    resp <- matrix(0, 28, m)
    
    for (i in 1:m){

      # indexes for each subj
      inc <- 56*(i-1)
      ind <- ij+inc
      ind_1 <- ij_1 + inc
      
      #----------------------------------------------------------------
      # Average the true simulated data for plot
      resp[,i] <- data[ind,3]
      
      #----------------------------------------------------------------
      # point estimates conditional on prior true simulated data, more accurate
      mu[,i] <- data[ind,4]*betaHat + fHat + alpha*(data[ind_1,3] - data[ind_1,4]*betaHat - c(fHat[28], fHat[1:27]))
      pse_sub[1, i] <- mean((mu[,i] - data[ind, 3])^2) # conditional on simulated previous true value
      
      #----------------------------------------------------------------
      # point estimates based on previous cycle, y_1,...,y_28
      mu_pc[,i] <- data[ind,4]*betaHat + fHat + alpha*(c(data[(28 + inc), 3],data[(1:27+inc),3]) - data[ind_1,4]*betaHat - c(fHat[28], fHat[1:27]))
      mu_pcv[,i] <- data[ind,4]*betaHat + fHat + alpha*(data[(1:28 + inc),3] - data[ind_1,4]*betaHat - c(fHat[28], fHat[1:27]))
      pse_sub[2,i] <- mean((mu_pc[,i] - data[ind, 3])^2)
      pse_sub[3,i] <- mean((mu_pcv[,i] - data[ind, 3])^2)
      
      #----------------------------------------------------------------
      # point estimates conditional on prior predictive obs'n
      mu_pre[1,i] <- data[(29+inc),4]*betaHat + fHat[1] + alpha[1]*(data[(28+inc), 3] - data[(28+inc),4]*betaHat - fHat[28])
      for (j in 2:28) {
        mu_pre[j,i] <- data[(28+j+inc),4]*betaHat + fHat[j] + alpha[j]*(mu_pre[(j-1),i] - data[(27+j+inc),4]*betaHat - fHat[j-1])
      }
      pse_sub[4,i] <- mean((mu_pre[,i] - data[ind, 3])^2)
      
      #----------------------------------------------------------------
      # Weighted version
      
      mu25[,i] <- 0.25*mu_pre[,i] + 0.75*mu_pc[,i]
      mu50[,i] <- 0.50*mu_pre[,i] + 0.50*mu_pc[,i]
      mu75[,i] <- 0.75*mu_pre[,i] + 0.25*mu_pc[,i]
      
      pse_sub[5,i] <- mean((mu25[,i] - data[ind, 3])^2)
      pse_sub[6,i] <- mean((mu50[,i] - data[ind, 3])^2)
      pse_sub[7,i] <- mean((mu75[,i] - data[ind, 3])^2)
    }
    
    # Record the average pse over all m subjects for simulation k 
    pse_sim[,k] <- rowMeans(pse_sub)

    # Index identifying subj's with the smallest, second smallest, median and largest pse values
    pseIn <- c(which.min(pse_sub[2,]), match(sort(pse_sub[2,], T)[m-1], pse_sub[2,]), match(n.median(pse_sub[2,]), pse_sub[2,]), which.max(pse_sub[2,]))
    mu_sim[[k]] <- cbind(resp[,pseIn], mu[,pseIn], mu_pc[,pseIn], mu_pcv[,pseIn], mu_pre[,pseIn], mu25[,pseIn], mu50[,pseIn], mu75[,pseIn])
    
    # Record the smallest, second smallest, median and largest pse values for kth simulation
    pseSort[,k] <- pse_sub[2,][pseIn]
  }
  
  # Sort the smallest, second smallest, medianand largest pse values over all simulations
  pseSortF <- c(min(pseSort[,1]), min(pseSort[2,]), n.median(pseSort[3,]), max(pseSort[4,]))
  
  # Index of simulation that have the smallest, second smallest, median and largest pse values
  pseFI <- sapply(1:4, function(v) {return(match(pseSortF[v], pseSort[v,]))})
  
  # Find the corresponding mu's according to the pse index
  muRank <- lapply(1:4, function(v) return(mu_sim[[pseFI[v]]][,seq(v, 32, by = 4)]))
  
  # smallest <- mu_sim[[pseFI[1]]][,seq(1, 32, by = 4)]
  # median <- mu_sim[[pseFI[2]]][,seq(2, 32, by = 4)]
  # sedLarg <- mu_sim[[pseFI[3]]][,seq(3, 32, by = 4)]
  # larg <- mu_sim[[pseFI[4]]][,seq(4, 32, by = 4)]
  
  # newlist <- list("pse" = pse, "smallest" = muRank[[1]], "sedSm" = muRank[[2]], "median" = muRank[[3]], "larg" = muRank[[4]])
  # return(newlist)
  
  pse <- c(rowMeans(pse_sim), NA)
  predRes <- rbind(pse, do.call(rbind, muRank))
  
  write.table(predRes, file='pred_univSim.txt', row.names=FALSE, col.names=FALSE) 
}

# return the first of middle elements when the vector is even number
n.median = function(x) {
  if (length(x) %% 2 != 0) {
    x <- median(x)
  } else if (length(x) %% 2 == 0) {
    a <- sort(x)[length(x)/2]
    x <- a[1]
  }
  return(x)
}
