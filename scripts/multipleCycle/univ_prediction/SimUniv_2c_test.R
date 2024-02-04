source("univFns_2c.R")

theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU

process <- "NOU"
# process <- "OU"
cycle <- 2
set.seed(2)
beta <- -0.5
data <- newData(theta, beta, 30, 28, process, cycle)

#--------------------------------------------------------------------------
# Check simulated data, whether it's cyclic, normal etc.

# # NOU range
# range(data$Y)
# 
# plot(data[1:56,3], type = "l")
# 
# # NORMALLY DISTRIBUTED.
# hist(data$Y)

#--------------------------------------------------------------------------
# prediction of 2nd cycle, i = 1, j = 29,...,56
t <- (1:56 %% 28)/10
ij <- 29:56
ij_1 <- 28:55

alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))
# point estimates based on simulated data, more accurate
mu <- 24*beta + sine(2.8, t[ij]) + alpha*(data[ij_1,3] - 24*beta - sine(2.8, t[ij_1]))

# point estimates based on previous predictive obs'n
mu2 <- rep(0, 28)
mu2[1] <- 24*beta + sine(2.8, t[29]) + alpha[1]*(data[28,3] - 24*beta - sine(2.8, t[28]))
for (i in 2:28) {
  mu2[i] <- 24*beta + sine(2.8, t[28+i]) + alpha[i]*(mu2[i-1] - 24*beta - sine(2.8, t[28+i-1]))
}

# point estimates based on previous cycle, y_1,...,y_28
mu_pc <- 24*beta + sine(2.8, t[ij]) + alpha*(c(data[28, 3],data[1:27,3]) - 24*beta - sine(2.8, t[ij_1]))
mu_pc1 <- 24*beta + sine(2.8, t[ij]) + alpha*(data[1:28,3] - 24*beta - sine(2.8, t[ij_1]))  # vector form


# Variance
expo <- theta[5] + (theta[6]*(s(28/3, t[ij_1]) + s(28/3, t[ij])) + theta[7]*(s(56/3, t[ij_1]) + s(56/3, t[ij])))/2
variance <- (1 + alpha^2)*1 + (1 + alpha^2)*1 + exp(expo)*(1 - alpha^2)

# Confidence Intervals
sd <- sqrt(variance)
CI_lower <- mu - 1.96*sd
CI_upper <- mu + 1.96*sd

matplot(1:28, cbind(data[ij, 3], mu_pc, mu_pc1, mu2, CI_lower, CI_upper), type="l", lty = c(1:4, 5, 5), col=c("black", "red", "blue", "orange", "green", "green"))

# PSE
pse_pc <- mean((mu_pc - data[ij, 3])^2) # conditional on previous cycle
pse_pc1 <- mean((mu_pc1 - data[ij, 3])^2) # conditional on previous cycle, vector form
pse_p <- mean((mu2 - data[ij, 3])^2) # conditional on predictive values
pse <- mean((mu - data[ij, 3])^2) # conditional on simulated previous true value

# fitted values from proposed model
res <- results(data, data$Y, data$age, random = 1, process = "NOU", data$day, data$id, tol = 0.1, cap = 3)
yHat <- res$betaHat*24 + res$fHat + res$U[29:56] + res$b[1]

yHat_1 <- res$betaHat*24 + res$fHat + res$U[28:55] + res$b[1]
mu3 <- 24*beta + sine(2.8, t[ij]) + alpha*(yHat_1 - 24*beta - sine(2.8, t[ij_1]))

#--------------------------------------------------------------------------
# Obtain relevant parameter values for est(), findPar(), and  without FS

# time <- data$day
id <- data$id
dim <- 7

# helperFn() requires time = data$day
helper <- helperFn(data, time=data$day, id)
niVec <- helper$niVec
K <- helper$K
N <- helper$N
r <- helper$r
Bdou <- helper$Bdou
tprime <- helper$tprime

X <- as.matrix(data$age)
Y <- as.matrix(data$Y)

# est(), findPar(), fisher.scoring() requires time = tprime
time <- tprime

#--------------------------------------------------------------------------
# Obtain betaHat and fHat without FS

res <- est(theta, X, Y, N, K, niVec, r, time, id, process)

# Check whether the likelihood function is inf or -inf
res$lik

# Test whether fHat is continuous
plot(rep(res$fHat, 2), type = "l")

plot(res$fHat, type = "l")

theta1 <- findPar(theta, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)

fs <- fisher.scoring(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol = 0.1, cap = 2, process, dim)
#--------------------------------------------------------------------------
# Test run results()

response <- data$Y
fixed <- data$age 
random <- 1
time <- data$day
id <- data$id 
tol <- 0.1 
cap <- 2

res <- results(data, response, fixed, random, process, time, id, tol, cap)

res <- results(data, data$Y, data$age, random = 1, process = "NOU", data$day, data$id, tol = 0.1, cap = 3)