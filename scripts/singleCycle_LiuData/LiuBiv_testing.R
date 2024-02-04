#======================================================================================
# TESTING
#======================================================================================
source("bivFns_lik.R")

#------------------------------------------------------------------------------

# NOU
theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(0.1, -2, 1.5, -0.1), 2)) 

# OU
theta <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.9), sqrt(-log(0.9)*2)), 2))

#------------------------------------------------------------------------------
# Test whether OU == NOU

# NOU
thetaNOU <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(0.5, -0.3266, 0, 0), 2)) 

# OU 
thetaOU <- c(rep(1, 3), -0.5, rep(1, 3), rep(c(-log(0.5), sqrt(-2*log(0.5)*exp(-0.3266))), 2))

#------------------------------------------------------------------------------
# test whether Vou == Vnou
Vou <- matrixVi(thetaOU, 5, 1:5, process = "OU") 

Vnou <- matrixVi(thetaNOU, 5, 1:5, process = "NOU") 

Vou$Vi
Vnou$Vi

Vou$Gammai[1,]
Vnou$Gammai[1,]
#------------------------------------------------------------------------------
helper <- helperFn(data, time = data$standday, id = data$womanid, fixed = cbind(data$age, data$BMI), response = cbind(data$adjpdg2, data$adje1c2))
niVec <- helper$niVec
K <- helper$K
N1 <- helper$N1
N2 <- helper$N2
r <- helper$r
B1dou <- helper$B1dou
B2dou <- helper$B2dou

X <- helper$X
Y <- helper$Y

time <- data$standday
id <- data$womanid
dim <- 15

# V1 <- matrixVi(theta, ni = niVec[1], ti = data$standday[which(data$womanid== uniqId[1])], process)

process <- "NOU"
process <- "OU"

res <- est(theta, X, Y, N1, N2, K, niVec, r, time, id, process)
res$lik

thetaNew <- findPar(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, process, dim)
fs <- fisher.scoring(theta, X, Y, N1, N2, K, niVec, r, time, id, B1dou, B2dou, tol=0.1, cap=10, process, dim)
res <- results(data, response = cbind(data$Y1, data$Y2), fixed = data$age, random = cbind(1, 1), process = "Wiener", time = data$day, id = data$id, tol = 0.1, cap = 10)

matW <- matrixW(param, niVec)

# tol = 0.05 cap = 50 - NOT converge
# > system.time(source("LiuBiv_userFrdly.R"))
     # user    system   elapsed 
# 29158.227   477.292 40204.045 


# OU 

#======================================================================================
# Input data
#======================================================================================

# m = 307 is the total number of women in the dataset (dataset w/ complete age BMI)

#------------------------------------------------------------------------------
# want <- read.table("ageBMI_want.txt", header = T)
# data <- read.table("ageBMI.txt", header = T)

# INPUT the FIRST 34 SUBJ. TO REDUCE DIMENSION, r is still 56.
# m = 34
# want <- read.table("ageBMI_want.txt", header = T)[1:34,]
# data <- read.table("ageBMI.txt", header = T)[1:992,]

# m = 50
# want <- read.table("ageBMI_want.txt", header = T)[1:50,]
# data <- read.table("ageBMI.txt", header = T)[1:1463,]

# m = 100
# want <- read.table("ageBMI_want.txt", header = T)[1:100,]
# data <- read.table("ageBMI.txt", header = T)[1:2810,]

#------------------------------------------------------------------------------
# create two indicator variables

dataComp <- read.table("ageBMI.txt", header = T)

dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# m = 50, randomly select 50 subjects with max 28 in standday
id <- dataComp[,1]
m <- max(id)
niMax <- sapply(1:m, function(v) return(max(dataComp[which(id == v), 2])))
niMax <- cbind(matrix(1:m, m), as.matrix(niMax))

# subject id with 28 standday
niMax28 <- niMax[which(niMax[,2] == 28)]

set.seed(100)
s <- sample(niMax28, 50)

data <- lapply(1:50, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data)

# Centre the covariates standday age. 
data[,2] <- (data[,2] - 14)/10
data[,5] <- (data[,5] - 33)/100




#------------------------------------------------------------------------------
# Run 

system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age, data$underWeight, data$overWeight), random = cbind(1, 1), process = "OU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))

write.csv(unlist(res), "LiuBiv50_OUeq_cat.txt")


# Notes: 135 NOU (result's out), 133 OU (awaiting for result)
