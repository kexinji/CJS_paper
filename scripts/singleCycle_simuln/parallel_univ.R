# PARALLEL COMPUTING FOR UNIV. SIMULATION

eachRes <- function(i) {
	source("univFns_lik.R")
	
	# OU, theta0 <- c(tau, sigma.b, sigma, theta2, theta3)
	# theta <- c(rep(1, 3), -log(0.0963), sqrt(-log(0.0963)*2))
	
	# NOU
	# theta <- c(rep(1, 3), -log(0.0963), -5, 1.5, -0.1)	

	# Wiener
	theta <- c(rep(1, 3), 0.7)
	
	process <- "Wiener"

	set.seed(i)
	simdata <- newData(theta, -0.5, 30, 28, process)
	
	res <- results(simdata, simdata$Y, simdata$age, random = 1, process, simdata$day, simdata$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('univW_07_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(64)

parLapply(cl, 1:64, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

#----------------------------------------------------------------------------------------------------------
# LOOP, NOT parallel
# source("univFns_lik.R")
# theta <- c(rep(1, 3), -log(0.8), sqrt(-log(0.8)*2))
# for (i in 1:100) {
	# print(paste('Iteration', i))
	# set.seed(i)
	# simdata <- newData(theta, -0.5, 30, 28, process = "OU")
	
	# res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "OU", simdata$day, simdata$id, tol = 0.01, cap = 100)
	
	# write.csv(unlist(res), paste('univ_OUt01c100_' , i, '.txt', sep=''))

# res <- results(data, data$Y, data$age, random = 1, process, data$day, data$id, tol = 0.001, cap = 50)
# write.csv(unlist(res), paste('univNP_NOUt001c50_400.txt'))
# }



#======================================================================================
# TESTING
#======================================================================================
# # source("univFns_lik.R")
# theta <- c(rep(1, 3), -log(0.0963), -5, 1.5, -0.1) # NOU
# # # theta <- c(rep(1, 3), -log(0.0963), sqrt(-log(0.0963)*2)) # OU - Code works!
# # theta <- c(rep(1, 3), 0.7) # Wiener 

# # process <- "Wiener"
# # # process <- "NOU"
# # # process <- "OU"
# # set.seed(1)
# # data <- newData(theta, -0.5, 30, 28, process)
# # time <- data$day
# # id <- data$id
# # # dim <- 7
# # # dim <- 5
# # dim <- 4

# helper <- helperFn(data, time, id)
# X <- as.matrix(data$age)
# Y <- as.matrix(data$Y)
# N <- helper$N
# K <- helper$K
# r <- helper$r
# Bdou <- helper$Bdou
# niVec <- helper$niVec

# # res <- est(theta, X, Y, N, K, niVec, r, time, id, process)
# # theta1 <- findPar(theta, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)


# # fs <- fisher.scoring(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol = 0.001, cap = 10, process, dim)
# # res <- results(data, data$Y, data$age, random = 1, process, data$day, data$id, tol = 0.001, cap = 50)
