# PARALLEL COMPUTING FOR UNIV. SIMULATION

eachRes <- function(i) {
  source("univSim_fn.R")
  
	# theta0 <- c(tau, sigma.b, sigma, theta2, theta3)
	# theta <- c(rep(1, 3), -log(0.0963), sqrt(-log(0.0963)*2))
  theta <- c(rep(1, 3), 0.2, -0.44, 0.3, -0.2)
  # theta <- c(rep(1, 3), 0.15, -1.60, 0.3, -0.1)

	set.seed(i)
	simdata <- newData(theta, 1, 30, 28, process = "NOU")
	
	res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "NOU", simdata$day, simdata$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('paper_univSim_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(100)

parLapply(cl, 1:500, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

#----------------------------------------------------------------------------------------------------------
# LOOP, NOT parallel
# setwd("C:/Users/s7486232/Downloads")
# source("univSim_fn.R")
# 
# theta <- c(rep(1, 3), 0.2, -0.44, 0.3, -0.2)
# for (i in 1:500) {
# 	print(paste('Simulation', i))
# 	set.seed(i)
# 	simdata <- newData(theta, 1, 30, 28, process = "NOU")
# 	
# 	res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "NOU", simdata$day, simdata$id, tol = 0.005, cap = 50)
# 
# 	write.csv(unlist(res), paste('paper_univSim_' , i, '.txt', sep=''))
# }

