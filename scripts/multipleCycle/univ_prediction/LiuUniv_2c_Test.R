#======================================================================================
# Import data with with 4 concecutive cycles, but only use 2 cycles
#======================================================================================

dataComp <- read.table("data_2cycle_old.txt", header = T)

# dataComp <- subset(dataComp, select = c(womanid, standday, adjpdg2, age, BMI))

# create two indicator variables
dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# m = 10, randomly select 10 subjects from the dataset
id <- dataComp[,1]
m <- 10

set.seed(42)
samp <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == samp[v]),])))
data <- do.call(rbind, data)

# Centre the covariates, standday and age
# data[,2] <- (data[,2] %% 28)/10

medianAge <- median(data[,5][!duplicated(data[,c(1, 5)])])
data[,5] <- (data[,5] - medianAge)/100

# medianAge <- median(data[,6][!duplicated(data[,c(1, 6)])])
# data[,6] <- (data[,6] - medianAge)/100

niVec <- sapply(1:m, function(v) return(nrow(data[which(id == samp[v]),])))

#======================================================================================
# Test
#======================================================================================
source("univFns_2c.R")

# parameters in results()

response <- log(data$adjpdg2)
fixed <- cbind(data$age, data$underWeight, data$overWeight)
random <- 1
time <- data$standday
id <- data$womanid
tol <- 0.1
cap <- 10

process <- "NOU"
theta <- c(rep(1, 3), 0.4, -1, 0.1, -0.1)	
dim <- 7

#--------------------------------------------------------------------------	
# Obtain values from helper functions
helper <- helperFn(data, time, id)
niVec <- helper$niVec
K <- helper$K
N <- helper$N
r <- helper$r
Bdou <- helper$Bdou

X <- as.matrix(fixed)
Y <- as.matrix(response)


res <- est(theta, X, Y, N, K, niVec, r, time, id, process)


theta1 <- theta
lik1 <- res$lik

theta0 <- theta1
lik0 <- lik1

theta1 <- findPar(theta, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)

fs <- fisher.scoring(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol, cap, process, dim)

system.time(res <- results(data, response = log(data$adjpdg2), fixed = cbind(data$age, data$underWeight, data$overWeight), random = 1, process = "NOU", time = data$standday, id = data$womanid, tol=0.1, cap=10))

#======================================================================================
# prediction for subject 1, 55 obs'ns, seed = 36
#======================================================================================
theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU

t <- (data[1:55, 2] %% 28)/10
ij <- 31:55
ij_1 <- 30:54

alpha31 <- exp(abs(t[31] - t[30])*log(theta[4]))
# point estimates based on simulated data, more accurate
mu31 <- res$betaHat[1]*(-0.005) + res$fHat[1] + alpha31*(data[30,3] - res$betaHat[1]*(-0.005) - res$fHat[28])

diff31 <- data[31,3] - mu31

alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))

# Find the index where the 2nd cycle has the same time point as the 1st cycle
index <- sapply(ij_1, function(v) return(match(t[v], t[1:30])))

mu <- res$betaHat[1]*(-0.005) + res$fHat[c((t*10)[31:54],28)] + alpha*(data[index,3] - res$betaHat[1]*(-0.005) - res$fHat[c(28, (t*10)[31:54])])


matplot(1:25, cbind(data[ij, 3], mu), type="l", lty = 1:2, col=c("black", "red"))
pse <- mean((mu - data[ij, 3])^2)

splinefun(data[1:30,2], data[1:30,3])

x <- data[1:30, 2]
x <- 1:30
y <- data[1:30, 3]

plot(x, y)
lines(spline(x, y))
splinefun(x, y)

#======================================================================================
# prediction for subject 1, 55 obs'ns, seed = 42
#======================================================================================
theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU

t <- (data[1:55, 2] %% 28)/10
ij <- 31:55
ij_1 <- 30:54

alpha31 <- exp(abs(t[31] - t[30])*log(theta[4]))
# point estimates based on simulated data, more accurate
mu31 <- res$betaHat[1]*(-0.005) + res$fHat[1] + alpha31*(data[30,3] - res$betaHat[1]*(-0.005) - res$fHat[28])

diff31 <- data[31,3] - mu31

mun <- matrix(0, 25)
alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))

mun[1] <- mu31
mun[2] <- res$betaHat[1]*(-0.005) + res$fHat[2] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[3] <- res$betaHat[1]*(-0.005) + res$fHat[3] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[4] <- res$betaHat[1]*(-0.005) + res$fHat[4] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[5] <- res$betaHat[1]*(-0.005) + res$fHat[5] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[6] <- res$betaHat[1]*(-0.005) + res$fHat[6] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[7] <- res$betaHat[1]*(-0.005) + res$fHat[7] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[8] <- res$betaHat[1]*(-0.005) + res$fHat[8] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[9] <- res$betaHat[1]*(-0.005) + res$fHat[9] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[10] <- res$betaHat[1]*(-0.005) + res$fHat[10] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[11] <- res$betaHat[1]*(-0.005) + res$fHat[11] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[12] <- res$betaHat[1]*(-0.005) + res$fHat[12] + alpha[2]*(data[1,3] - res$betaHat[1]*(-0.005) - res$fHat[1])
mun[14] <- res$betaHat[1]*(-0.005) + res$fHat[15] + alpha[14]*(mean(data[14,3], data[15,3]) - res$betaHat[1]*(-0.005) - res$fHat[14])





alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))

# Find the index where the 2nd cycle has the same time point as the 1st cycle
tempp <- c((t[31:54] - 0.1), 2.7)
temp <- c(seq(0, 1.2, by = 0.1), 1.4, seq(1.6, 2.2, by = 0.1), seq(2.4, 2.7, by = 0.1))
index <- sapply(1:25, function(v) return(match(temp[v], t[1:30])))

index[4] <- 3
index[7] <-6
index[8] <- 7
index[16] <- 18
index[18] <- 20

mu <- res$betaHat[1]*(-0.005) + res$fHat[c((t*10)[31:54],28)] + alpha*(data[index,3] - res$betaHat[1]*(-0.005) - res$fHat[c(28, (t*10)[31:54])])


matplot(1:25, cbind(data[ij, 3], mu), type="l", lty = 1:2, col=c("black", "red"))
pse <- mean((mu - data[ij, 3])^2)

#======================================================================================
# results
#======================================================================================
result <- read.table("LiuUniv_2c.txt", header = T, sep = ",")

p <- 7 # number of smoothing and variance parmaeters
r <- 28 # number of distinct time points
b <- 3 # number of fixed effect parameters, only age

fhat <- result[(p+b+1) : (p+b+r),2]
plot(fhat, type = "l")

sdF <- result[(p+2*b+r+1):(p+2*b+2*r), 2]
lowerF <- fhat - 1.96*sdF
upperF <- fhat + 1.96*sdF

#------------------------------------------------------------------------------
# Organize the dataset according to observed timepoints
tp <- sort(unique(data[,2]))
toltp <- length(tp)
logPro <- lapply(1:toltp, function(v) return(log(data[which(data$standday == tp[v]), 3])))
day <- lapply(1:toltp, function(v) return(rep(tp[v], length(logPro[[v]]))))

pro <- cbind(unlist(day), unlist(logPro))

#------------------------------------------------------------------------------
# plots
pdf("LiuUniv_2c_f.pdf")
plot(pro, pch = 39, ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n")
par(new = T)
plot(tp, rep(fhat, 2), type = "l",ylim = c(-3, 4), lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(lowerF, 2), type = "l", col = "red", ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(upperF, 2), type = "l", col = "green", ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n", yaxt="n")
title(xlab="Days in standardized menstrual cycle", ylab="Log progesterone")
axis(1, at = seq(0, 56, by = 7), labels = c(seq(0, 28, by = 7), seq(7, 28, by = 7)))
dev.off()
