# PARALLEL COMPUTING FOR BIV. SIMULATION

eachRes <- function(i) {
	source("bivFns.R")
	
	set.seed(i)
	theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.0963), -2.1563, 3.6710, -0.145), 2)) # NOU
	data <- newData(theta, beta = c(0.03, 0.07), 30, 28, process = "NOU")
	
	res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "NOU", time = data$day, id = data$id, tol = 0.001, cap = 200)
	
	write.csv(unlist(res), paste('bivNOU64_t001c200_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(64)

parLapply(cl, 1:64, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

# result in error message - need to figure out why
# > source("parallel_biv.R")
# Error in checkForRemoteErrors(val) : 
  # 52 nodes produced errors; first error: infinite or missing values in 'x'


#======================================================================================
# TESTING
#======================================================================================

# source("bivFns.R")
# theta <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)
# theta <- c(rep(1, 3), -0.5, rep(1, 3), -log(0.8), rep(0, 3), -log(0.8), rep(0, 3))

# process <- "NOU"
# data <- newData(theta, beta = c(0.03, 0.07), 30, 28, process="NOU")
# time <- data$day
# id <- data$id
# dim <- 15

# helper <- helperFn(data, time = data$day, id=data$id, fixed = data$age, response=cbind(data$Y1, data$Y2))
# niVec <- helper$niVec
# K <- helper$K
# N1 <- helper$N1
# N2 <- helper$N2
# r <- helper$r
# B1dou <- helper$B1dou
# B2dou <- helper$B2dou

# X <- helper$X
# Y <- helper$Y

# V1 <- matrixVi(theta, ni = niVec[1], ti = data$day[which(data$id==1)], process)


# res <- est(theta, X, Y, N1, N2, K, niVec, r, time = data$day, id=data$id, process="NOU")
# res$betaHat
# plot(res$f1Hat)
# plot(res$f2Hat)

# PV1NOU <- PVi(theta, ni = niVec[1], ti = data$day[which(data$id==1)], process)

# thetaNew <- findPar(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, process="NOU", dim)
# fs <- fisher.scoring(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, tol=0.9, cap=1, process="NOU", dim)