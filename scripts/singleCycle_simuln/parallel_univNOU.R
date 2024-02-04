# PARALLEL COMPUTING FOR UNIV. SIMULATION

eachRes <- function(i) {
	source("univFns_lik.R")
	
	# theta0 <- c(tau, sigma.b, sigma, theta2, theta3)
	# theta <- c(rep(1, 3), -log(0.0963), sqrt(-log(0.0963)*2))
	theta <- c(rep(1, 3), -log(0.0963), -5, 1.5, -0.1)	

# TESTING INITIATION VALUES FOR THETA, NOU 	
# t <- 1:28
# exponent <- log(2*theta[4]) + theta[5] + theta[6]*t + theta[7]*(t^2)
# t3 <- exp(exponent/2)
# sapply(1:28, function(v) return(rsOU(1, theta=c(0, theta[4], t3[v]))))

	set.seed(i)
	simdata <- newData(theta, -0.5, 30, 28, process = "NOU")
	
	res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "NOU", simdata$day, simdata$id, tol = 0.001, cap = 50)
	
	write.csv(unlist(res), paste('univNOU_t001c50_' , i, '.txt', sep=''))
}

library(parallel)
cl <- makeCluster(100)

parLapply(cl, 1:100, eachRes)

# Once we are done we need to close the cluster so that resources such as memory are returned to the operating system.
stopCluster(cl)

#----------------------------------------------------------------------------------------------------------
# LOOP, NOT parallel
# theta <- c(rep(1, 3), -log(0.8), sqrt(-log(0.8)*2))
# for (i in 1:100) {
	# print(paste('Iteration', i))
	# set.seed(i)
	# simdata <- newData(theta, -0.5, 30, 28, process = "OU")
	
	# res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "OU", simdata$day, simdata$id, tol = 0.01, cap = 100)
	
	# write.csv(unlist(res), paste('univ_OUt01c100_' , i, '.txt', sep=''))
# }



#======================================================================================
# TESTING
#======================================================================================
source("univFns_lik.R")
# theta <- c(rep(1, 3), -log(0.0963), -5, 1.5, -0.1) # NOU
# theta <- c(rep(1, 3), -log(0.0963), sqrt(-log(0.0963)*2)) # OU - Code works!

theta <- c(0.001553081, 2.947410649, 1.747088038, 1.802835685, 8.770661877, -8.655713148, 2.269095144)

process <- "NOU"
# process <- "OU"
set.seed(2)
data <- newData(theta, -0.5, 30, 28, process)
time <- data$day
id <- data$id
dim <- 7
# dim <- 5

helper <- helperFn(data, time, id)
X <- as.matrix(data$age)
Y <- as.matrix(data$Y)
N <- helper$N
K <- helper$K
r <- helper$r
Bdou <- helper$Bdou
niVec <- helper$niVec

res <- est(theta, X, Y, N, K, niVec, r, time, id, process)

thetanew <- findPar(theta, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)

thetaupdated <- findPar(thetanew, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)


fs <- fisher.scoring(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol = 0.001, cap = 50, process, dim)

res <- results(data, data$Y, data$age, random = 1, process = "NOU", data$day, data$id, tol = 0.001, cap = 50)