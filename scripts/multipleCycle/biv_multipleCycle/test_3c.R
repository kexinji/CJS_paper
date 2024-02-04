dataComp <- read.table("data_3cycle.txt", header = T)

dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL
#------------------------------------------------------------------------------
# m = 10, randomly select 30 subjects with 3 concecutive cycles
id <- dataComp[,1]
m <- 10

set.seed(42)
s <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data)

# Centre the covariates, standday and age
# data[,2] <- data[,2]/10

medianAge <- median(data[,5][!duplicated(data[,c(1, 5)])])
data[,5] <- (data[,5] - medianAge)/100

#------------------------------------------------------------------------------
# parameters in results()

response <- cbind(log(data$adjpdg2), log(data$adje1c2))
fixed <- cbind(data$age, data$underWeight, data$overWeight)
random <- cbind(1, 1)
process <- "NOU"
time <- data$standday
id <- data$womanid
tol <- 0.1
cap <- 10

dim <- 15

#--------------------------------------------------------------------------	
source("bivFns_3c.R")
# Obtain values from helper functions
helper <- helperFn(data, time=data[,2], id, fixed, response)
niVec <- helper$niVec
K <- helper$K
N1 <- helper$N1
N2 <- helper$N2
r <- helper$r
B1dou <- helper$B1dou
B2dou <- helper$B2dou
tprime <- helper$tprime

X <- helper$X
Y <- helper$Y

# theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(0.1, -5, 1.5, -0.1), 2)) # res$lik = inf
# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.2, -0.44, 0.3, -0.2, 0.15, -1.6, 0.3, -0.1)  # res$lik = inf

# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.4, -1, 0.1, -0.1, 0.4, -1, 0.1, -0.1)  # res$lik = inf
# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.4, -1, 0.1, -0.1, 0.1, -5, 1.5, -0.1) # res$lik = inf
# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.4, -1, 0.1, -0.1, 0.15, -1.6, 0.3, -0.1) # res$lik = inf

# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.4, -1, 0.1, -0.1, 0.15, -0.6, 0.3, -0.1) # res$lik = inf

# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.1, -5, 1.5, -0.1, 0.2, -0.44, 0.3, -0.2)  # res$lik = inf
# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.1, -5, 1.5, -0.1, 0.15, -1.6, 0.3, -0.1)  # res$lik = inf

# theta <- c(rep(1, 3), -0.5, rep(1, 3), 0.4, -1, 0.1, -0.1, 0.1, -2, 1, -0.1) # res$lik = inf

theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(0.07, -1.32, -4.8, 0.6), 2)) 

ti <- time[which(id==uniqId[10])]
ni <- niVec[10]
log(det(matrixVi(theta, ni, ti, process)$Vi))

G1 <- outer(ti, ti, nouVar, theta[8], theta[9], theta[10], theta[11])	
G2 <- outer(ti, ti, nouVar, theta[12], theta[13], theta[14], theta[15])	


res <- est(theta, X, Y, N1, N2, K, niVec, r, time=tprime, id, process)
res$lik

theta1 <- theta
lik1 <- res$lik

theta0 <- theta1
lik0 <- lik1

theta1 <- findPar(theta0, X, Y, N1, N2, K, niVec, r, time=tprime, id, B1dou, B2dou, process, dim)
fs <- fisher.scoring(theta, X, Y, N1, N2, K, niVec, r, time=tprime, id, B1dou, B2dou, tol=0.1, cap=10, process, dim)
res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age, data$underWeight, data$overWeight), random = cbind(1, 1), process = "NOU", time = data$standday, id = data$womanid, tol = 0.1, cap = 10)


#==========================================================================================================
library(sde) # rsOU()
library(MASS) # mvrnorm()
library(psych) # tr()
library(magic) # adiag()
library(Rcpp) # TO C++
sourceCpp("matrixmulti_1.cpp") # NEED TO INSTALL Package 'RcppEigen'
library(e1071) # rwiener()

# install.packages("sde")
# install.packages("MASS")
# install.packages("psych")
# install.packages("magic")
# install.packages("Rcpp")
# install.packages("e1071")
# install.packages("RcppEigen")

#--------------------------------------------------------------------------	
# niVec
uniqId <- unique(id)
m <- length(uniqId)
niVec <- sapply(1:m, function(v) return(nrow(data[which(id == uniqId[v]),])))
n <- length(id)

#--------------------------------------------------------------------------
# r, different than before
tprime <- time %% 28 # mod(time, 28)
uniquet <- sort(unique(tprime))
r <- length(uniquet) 

#--------------------------------------------------------------------------
# N1 & N2
# colIndex <- sapply(time, function(v) return(which(v == uniquet))) # slower
colIndex <- sapply(tprime, function(v) return(match(v, uniquet)))
N <- matrix(0, n, r)
N[cbind(1:n, colIndex)] <- 1
N1 <- do.call(rbind, lapply(1:n, function(v) return(rbind(N[v,], rep(0, r)))))
N2 <- do.call(rbind, lapply(1:n, function(v) return(rbind(rep(0, r), N[v,]))))

#--------------------------------------------------------------------------
# K
h <- c(uniquet[2:r], 28) - uniquet	

# testing
# r <- 5
# h <- c(0.5, rep(0.1,3), 0.5)

# Q 
Q <- matrix(0, r, r)

Q[1,1] <- -1/h[1] - 1/h[r]
diag(Q) <- c(Q[1,1], (-1/h[1:(r-1)] - 1/h[2:r]))

Q[1,r] <- 1/h[r]
Q[r,1] <- Q[1,r]

diag(Q[-r, -1]) <- 1/h[1:(r-1)]
diag(Q[-1, -r]) <- 1/h[1:(r-1)]

# R 
R <- matrix(0, r, r)

R[1,1] <- (h[1] + h[r])/3
diag(R) <- c(R[1,1], (h[1:(r-1)] + h[2:r])/3)

R[1,r] <- h[r]/6
R[r,1] <- R[1,r]

diag(R[-r, -1]) <- h[1:(r-1)]/6
diag(R[-1, -r]) <- h[1:(r-1)]/6

K <- eigenMapMatMult3(Q, ginv(R), t(Q))	

#--------------------------------------------------------------------------
# B
B <- Q %*% ginv(crossprod(Q)) %*% t(chol(R))

# LTranspose <- chol(K, pivot=T, LDL = T) # pivot = T to handle positive-semi-definite
# L <- t(LTranspose)
# B <- L %*% ginv(crossprod(L))

# B1* B2*, TO BE USED IN findPar, SCORE
B1dou <- tcrossprod(N1 %*% B)
B2dou <- tcrossprod(N2 %*% B)