# PARALLEL COMPUTING CODE FOR SIMULATION

eachRes <- function(i) {
	
	source("biv_fns.R")
	
	# parameter length of 11
	# theta0 = c(tau1, tau2, phi1, phi2, phi3, sigma1, sigma2, theta12, theta13, theta22, theta23)
	param <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)
	
	set.seed(i)
	data <- newData(param, beta11 = 0.03, beta21 = 0.07, m=30, tp=28)
	Y <- data$Y
	X <- data$X
	N1 <- data$N1
	N2 <- data$N2
	K <- matrixK(h=1, tp=28)
		
	# Estimates of beta, f, theta, bias and Covariance from one simulation
	res <- results(param, Y, tol = rep(0.1, 11), cap = 10, m=30, tp=28, K, X, N1, N2, niVec=rep(28, 30))
	
	write.csv(unlist(res), paste('t01_c100_gen_' , i, '.txt', sep=''))	
}

library(parallel)
cl <- makeCluster(2)

parLapply(cl, 1:2, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)


# > detectCores()
# [1] 4

# 2 CLUSTERS, parLapply(cl, 1:2, eachRes), OFFICE COMPUTER
# > system.time(source("parallel.r"))
   # user  system elapsed 
  # 0.134   0.103 556.735 
# > 556.735/60
# [1] 9.278917

# cpu139 100 CLUSTERS, parLappy(cl, 1:100, eachRes)
# > system.time(source("parallel.R"))
   # user  system elapsed 
  # 0.167   0.253 496.821 