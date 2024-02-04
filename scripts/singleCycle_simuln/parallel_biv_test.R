#======================================================================================
# TESTING
#======================================================================================

source("bivFns_sim.R")
# theta <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)
# theta <- c(rep(1, 3), -0.5, rep(1, 3), -log(0.0963), rep(0, 3), -log(0.8), rep(0, 3))

# OU
# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.9), sqrt(-log(0.9)*2)), 2))

# NOU
# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.4), -5, 1.5, -0.1), 2))
# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.2, -0.44, 0.3, -0.2, 0.15, -1.6, 0.3, -0.1)

# Wiener
theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.1, 0.1)

t <- 1:28
exponent1 <- log(2*theta[8]) + theta[9] + theta[10]*t + theta[11]*(t^2)
exponent2 <- log(2*theta[12]) + theta[13] + theta[14]*t + theta[15]*(t^2)
theta13 <- rep(exp(exponent1/2), 30)
theta23 <- rep(exp(exponent2/2), 30)
u1 <- sapply(1:840, function(v) return(rsOU(1, theta=c(0, theta[8], theta13[v]))))
u2 <- sapply(1:840, function(v) return(rsOU(1, theta=c(0, theta[12], theta23[v]))))
range(u1)
range(u2)

process <- "Wiener"
# process <- "NOU"
# process <- "OU"
set.seed(16)
data <- newData(theta, beta = c(1, 0.75), 30, 28, process)

# NORMALLY DISTRIBUTED.
hist(data$Y1)
hist(data$Y2)

# NOU range
range(data$Y1)
# [1] 12.79824 50.04450
range(data$Y2)
# [1] 25.28612 70.94049

# OU range
# > range(data$Y1)
# [1] 11.75678 49.89749
# > range(data$Y2)
# [1] 24.95699 71.55520

# Wiener range, sigma = 0.5
# > range(data$Y1)
# [1]  9.353117 51.007229
# > range(data$Y2)
# [1] 10.80251 40.60409



time <- data$day
id <- data$id
# dim <- 15
# dim <- 11
dim <- 9

helper <- helperFn(data, time = data$day, id=data$id, fixed = data$age, response=cbind(data$Y1, data$Y2))
niVec <- helper$niVec
K <- helper$K
N1 <- helper$N1
N2 <- helper$N2
r <- helper$r
B1dou <- helper$B1dou
B2dou <- helper$B2dou

X <- helper$X
Y <- helper$Y

# V1 <- matrixVi(theta, ni = niVec[1], ti = data$day[which(data$id==1)], process)

res <- est(theta, X, Y, N1, N2, K, niVec, r, time = data$day, id=data$id, process)
res$lik
# [1,] -930.8878
plot(res$f1Hat)
plot(res$f2Hat)

# PV1NOU <- PVi(theta, ni = niVec[1], ti = data$day[which(data$id==1)], process)

# thetaNew <- findPar(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, process, dim)
# fs <- fisher.scoring(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, tol=0.1, cap=10, process, dim)
# res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "Wiener", time = data$day, id = data$id, tol = 0.1, cap = 10)


#======================================================================================
# PARALLEL COMPUTING FOR BIV. SIMULATION
#======================================================================================
eachRes <- function(i) {
  source("bivFns_lik.R")
  
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), -2.1563, 3.6710, -0.145), 2)) # NOU
  # > range(u1)
  # [1] -65234.08 121079.08
  
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), -2, 2, -0.1), 2)) # NOU
  # > range(u1)
  # [1] -59.05207  71.27276
  
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), -2, 1.5, -0.1), 2)) # NOU
  # > range(u1)
  # [1] -11.64795  10.96343
  
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), 0, 0, 0), 2)) # NOU
  # > range(u1)
  # [1] -2.205472  2.053775	
  
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.2, -0.44, 0.3, -0.2, 0.15, -1.6, 0.3, -0.1)
  # > range(u1)
  # [1] -1.595918  2.096956
  # > range(u2)
  # [1] -1.3781716  0.8454772
  
  #------------------------------------------------------------------------------------
  # Parameter values used in simulating datasets.
  
  # OU
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.9), sqrt(-log(0.9)*2)), 2))
  
  # NOU 
  # theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.4), -5, 1.5, -0.1), 2)) 
  
  # Wiener
  theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(0.1, 2))
  
  set.seed(i)
  data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "Wiener")
  
  #------------------------------------------------------------------------------------
  # Obtain moadel parameter estimates using the specified parameters
  res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "Wiener", time = data$day, id = data$id, tol = 0.001, cap = 50)
  
  write.csv(unlist(res), paste('bivW100_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(100)

parLapply(cl, 1:100, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)



#----------------------------------------------------------------------------------------------------------
# Use for loop instead of parallel computing
# source("bivFns.R")
# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), -5, 1.5, -0.1), 2)) # NOU
# for (i in 1:100) {
# print(paste('Iteration', i))
# set.seed(i)
# data <- newData(theta, beta = c(0.03, 0.07), 30, 28, process = "NOU")

# res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "NOU", time = data$day, id = data$id, tol = 0.001, cap = 200)

# write.csv(unlist(res), paste('NPbivNOU_t001c200_' , i, '.txt', sep=''))
# }
