eachRes <- function(i) {
	source("bivFns_lik.R")

  # Wiener
	theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(0.1, 2))
	
	set.seed(i)
	data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "Wiener")
	
	#------------------------------------------------------------------------------------
	# Obtain moadel parameter estimates using the specified parameters
	res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "Wiener", time = data$day, id = data$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('bivW300.05_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(100)

parLapply(cl, 1:300, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)



