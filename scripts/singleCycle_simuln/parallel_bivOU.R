eachRes <- function(i) {
	source("bivFns_lik.R")
	
	# # OU
	# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.9), sqrt(-log(0.9)*2)), 2))
	
	# OU testing for equality with NOU
	theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.25), sqrt(-2*log(0.25)*exp(-1.02))), 2))
	
  set.seed(i)
	data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "OU")
	
	# NOU 
	# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.3), -5, 1.5, -0.1), 2)) 
	# set.seed(i)
	# data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "NOU")
	
	#------------------------------------------------------------------------------------
	# Fit OU on NOU datasets
	res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "OU", time = data$day, id = data$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('bivOUeq_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(64)

parLapply(cl, 1:300, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

