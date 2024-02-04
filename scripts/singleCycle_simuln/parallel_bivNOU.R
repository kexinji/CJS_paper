eachRes <- function(i) {
	source("bivFns_sim.R")
	# NOU
  # old initializations
	# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(0.1, -5, 1.5, -0.1), 2)) 
  
  # initializations from liu analysis
  theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.2, -0.44, 0.3, -0.2, 0.15, -1.6, 0.3, -0.1)
	
	# NOU testing for equality with OU
	# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.25), -1.02, 0, 0), 2))
	
	set.seed(i)
	data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "NOU")
	
	#------------------------------------------------------------------------------------
	# Obtain moadel parameter estimates using the specified parameters
	res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "NOU", time = data$day, id = data$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('bivNOUliu_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(100)

parLapply(cl, 1:600, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

#------------------------------------------------------------------------------------
# for loop

# source("bivFns_lik.R")
# for (i in c(55, 92:95, 158:160, 182:185, 190, 218:220, 377:380, 438:440)) {
	# print(i)
	# # NOU 
	# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.4), -5, 1.5, -0.1), 2)) 
	# set.seed(440)
	# data <- newData(theta, beta = c(1, 0.75), 30, 28, process = "NOU")
	
	# #------------------------------------------------------------------------------------
	# # Obtain moadel parameter estimates using the specified parameters
	# res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "NOU", time = data$day, id = data$id, tol = 0.001, cap = 50)
	# write.csv(unlist(res), paste('bivNOU500_440.txt', sep=''))

	# write.csv(unlist(res), paste('bivNOU500_' , i, '.txt', sep=''))
# }

