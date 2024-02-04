source("univFns_2c.R")

theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU
m <- 1000 # number of subjects/simulations

process <- "NOU"
# process <- "OU"
cycle <- 2
set.seed(1)
beta <- -0.5
data <- newData(theta, beta, m, 28, process, cycle)


#======================================================================================
# Prediction of 2nd cycle, j = 29,...,56
#======================================================================================
t <- (1:56 %% 28)/10
ij <- 29:56
ij_1 <- 28:55

alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))

# Initializations, each column is a simulation
mu <- matrix(0, 28, m) 
mu_pc <- matrix(0, 28, m)
mu_pcv <- matrix(0, 28, m)
mu_pre <- matrix(0, 28, m)

mu25 <- matrix(0, 28, m)
mu50 <- matrix(0, 28, m)
mu75 <- matrix(0, 28, m)

resp <- matrix(0, 28, m)

pse <- matrix(0, m)
pse_pc_loop <- matrix(0, m)
pse_pcv_loop <- matrix(0, m)
pse_pre_loop <- matrix(0, m)

pse25_loop <- matrix(0, m)
pse50_loop <- matrix(0, m)
pse75_loop <- matrix(0, m)


for (i in 1:m) {
  inc <- 56*(i-1)
  ind <- ij+inc
  ind_1 <- ij_1 + inc

  #--------------------------------------------------------------------------
  # point estimates conditional on prior true simulated data, more accurate
  mu[,i] <- data[ind,4]*beta + sine(2.8, t[ij]) + alpha*(data[ind_1,3] - data[ind_1,4]*beta - sine(2.8, t[ij_1]))
  pse[i] <- mean((mu[,i] - data[ind, 3])^2) # conditional on simulated previous true value
  
  #--------------------------------------------------------------------------
  # point estimates based on previous cycle, y_1,...,y_28
  mu_pc[,i] <- data[ind,4]*beta + sine(2.8, t[ij]) + alpha*(c(data[(28 + inc), 3],data[(1:27+inc),3]) - data[ind_1,4]*beta - sine(2.8, t[ij_1]))
  mu_pcv[,i] <- data[ind,4]*beta + sine(2.8, t[ij]) + alpha*(data[(1:28 + inc),3] - data[ind_1,4]*beta - sine(2.8, t[ij_1]))
  pse_pc_loop[i] <- mean((mu_pc[,i] - data[ind, 3])^2)
  pse_pcv_loop[i] <- mean((mu_pcv[,i] - data[ind, 3])^2)
  
  #--------------------------------------------------------------------------
  # point estimates conditional on prior predictive obs'n
  mu_pre[1,i] <- data[(29+inc),4]*beta + sine(2.8, t[29]) + alpha[1]*(data[(28+inc), 3] - data[(28+inc),4]*beta - sine(2.8, t[28]))
  for (j in 2:28) {
    mu_pre[j,i] <- data[(28+j+inc),4]*beta + sine(2.8, t[28+j]) + alpha[j]*(mu_pre[(j-1),i] - data[(27+j+inc),4]*beta - sine(2.8, t[27+j]))
  }
  pse_pre_loop[i] <- mean((mu_pre[,i] - data[ind, 3])^2)
  
  #--------------------------------------------------------------------------
  # Weighted version
  
  mu25[,i] <- 0.25*mu_pre[,i] + 0.75*mu_pc[,i]
  mu50[,i] <- 0.50*mu_pre[,i] + 0.50*mu_pc[,i]
  mu75[,i] <- 0.75*mu_pre[,i] + 0.25*mu_pc[,i]
  
  pse25_loop[i] <- mean((mu25[,i] - data[ind, 3])^2)
  pse50_loop[i] <- mean((mu50[,i] - data[ind, 3])^2)
  pse75_loop[i] <- mean((mu75[,i] - data[ind, 3])^2)
  
  #--------------------------------------------------------------------------
  # Average the true simulated data for plot
  resp[,i] <- data[ind,3]
}

pseTru <- mean(pse) # conditional on simulated previous true value
pse_pc <- mean(pse_pc_loop)
pse_pcv <- mean(pse_pcv_loop)
pse_pre <- mean(pse_pre_loop)

pse25 <- mean(pse25_loop)
pse50 <- mean(pse50_loop)
pse75 <- mean(pse75_loop)

pseTru
pse_pc
pse_pcv
pse_pre

pse25
pse50
pse75

matplot(1:28, cbind(rowMeans(resp), rowMeans(mu),  rowMeans(mu_pc), rowMeans(mu_pcv), rowMeans(mu_pre)), type="l", lty = c(1:5, 6, 6), col=c("black", "orange", "blue", "red", "green"))

#--------------------------------------------------------------------------
# Variance
expo <- theta[5] + (theta[6]*(s(28/3, t[ij_1]) + s(28/3, t[ij])) + theta[7]*(s(56/3, t[ij_1]) + s(56/3, t[ij])))/2
variance <- (1 + alpha^2)*1 + (1 + alpha^2)*1 + exp(expo)*(1 - alpha^2)

# Confidence Intervals
sd <- sqrt(variance)
CI_lower <- rowMeans(mu) - 1.96*sd
CI_upper <- rowMeans(mu) + 1.96*sd

matplot(1:28, cbind(rowMeans(resp), rowMeans(mu),  rowMeans(mu_pc), rowMeans(mu_pcv), rowMeans(mu_pre), CI_lower, CI_upper), type="l", lty = c(1:5, 6, 6), col=c("black", "orange", "blue", "red", "green", "pink", "pink"))
