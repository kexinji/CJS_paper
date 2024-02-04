#======================================================================================
# Import data with with 4 concecutive cycles, but only use 2 cycles
#======================================================================================

# dataComp <- read.table("data_2cycle_old.txt", header = T)

dataComp <- read.table("data_3cycle.txt", header = T)

# dataComp <- subset(dataComp, select = c(womanid, standday, adjpdg2, age, BMI))

# create two indicator variables
dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# randomly select m subjects from the dataset
m <- 50

set.seed(42)
samp <- sample(1:max(dataComp[,1]), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(dataComp[,1] == samp[v]),])))
data <- do.call(rbind, data)

# Centre the covariates, standday and age
# data[,2] <- (data[,2] %% 28)/10

niVec <- sapply(1:m, function(v) return(nrow(data[which(dataComp[,1] == samp[v]),])))

# niVec <- sapply(1:94, function(v) return(nrow(dataComp[which(dataComp[,1] == v),]))) # for 94 women

age_index <- sapply(1:(m-1), function(v) return(1+sum(niVec[1:v])))

# age range of 94 subj is 23-44 yo, median = 36
# median(dataComp[c(1, age_index), 6]) 

medianAge <- median(data[c(1, age_index), 6]) 
data$age <- (data$age - medianAge)/100

# medianAge <- median(data[,6][!duplicated(data[,c(1, 6)])])
# data[,6] <- (data[,6] - medianAge)/100

source("univFns_2c.R")
system.time(res <- results(data, response = log(data$adjpdg2), fixed = cbind(data$age, data$underWeight, data$overWeight), random = 1, process = "NOU", time = data$standday, id = data$womanid, tol=0.1, cap=10))

# length of 3rd cycle 
n3i <- sapply(1:m, function(v) return(nrow(data[which(data[,1] == samp[v] & data[,2] == 3),])))
#======================================================================================
# prediction for the 3rd cycle
#======================================================================================
theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU

# Initializations
pse_sub <- matrix(0, 7, m) # Store PSE, column for each subj, row for naive, pre, pc PSE

mu_true <- matrix(NA, max(n3i), m) 
mu_naive <- matrix(NA, max(n3i), m) 
mu_pc <- matrix(NA, max(n3i), m) 
mu_pc1 <- matrix(NA, max(n3i), m) 
mu_pre <- matrix(NA, max(n3i), m) 

mu25 <- matrix(NA, max(n3i), m)
mu50 <- matrix(NA, max(n3i), m)
mu75 <- matrix(NA, max(n3i), m)

CI_lower <- matrix(NA, max(n3i), m)
CI_upper <- matrix(NA, max(n3i), m)
 
standay <- matrix(NA, max(n3i), m)

for (i in 1:m) {
  temp <- data[which(data$womanid == samp[i]),]
  temp$index <- 1:dim(temp)[1]
  
  t <- (temp$standday %% 28)/10
  
  # indexes of time points
  ij <- temp$index[which(temp$CYCLE == 3)]
  ij_1 <- ij - 1
  ij_pc <- temp$index[which(temp$CYCLE == 2)]
  ij_pc1 <- temp$index[which(temp$CYCLE == 1)]
  cl <- 1:length(ij)
  
  alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))
  # Variance
  expo <- theta[5] + (theta[6]*(s(28/3, t[ij_1]) + s(28/3, t[ij])) + theta[7]*(s(56/3, t[ij_1]) + s(56/3, t[ij])))/2
  variance <- (1 + alpha^2)*1 + (1 + alpha^2)*1 + exp(expo)*(1 - alpha^2)
  
  
  tAge <- unique(temp$age)
  uW <- unique(temp$underWeight)
  oW <- unique(temp$overWeight)
  
  f_ij <- res$fHat[temp$standday[ij]]
  f_ij1 <- res$fHat[temp$standday[ij_1]]
  
  # Record response, standday for plotting
  standay[cl, i] <- temp$standday[ij]
  mu_true[cl, i] <- log(data[ij, 4])
  #--------------------------------------------------------------------------
  # prediction est conditional on true prior data
  mu_naive[cl, i] <- res$betaHat[1]*tAge + res$betaHat[2]*uW + res$betaHat[3]*oW + f_ij + alpha*(log(data[ij_1,4]) - res$betaHat[1]*tAge - res$betaHat[2]*uW - res$betaHat[3]*oW - f_ij1)
  pse_sub[1,i] <- mean((mu_naive[cl, i] - log(data[ij, 4]))^2)
  
  #--------------------
  # Variance
  expo <- theta[5] + (theta[6]*(s(28/3, t[ij_1]) + s(28/3, t[ij])) + theta[7]*(s(56/3, t[ij_1]) + s(56/3, t[ij])))/2
  variance <- (1 + alpha^2)*1 + (1 + alpha^2)*1 + exp(expo)*(1 - alpha^2)
  
  # Confidence Intervals
  sd <- sqrt(variance)
  CI_lower[cl, i] <- mu_naive[cl, i] - 1.96*sd
  CI_upper[cl, i] <- mu_naive[cl, i] + 1.96*sd
  
  #--------------------------------------------------------------------------
  # prediction results conditional (mostly) on cycle 2
  
  if (all(t[ij] %in% t[ij_pc])) { # Check if time points in C3 is a subset of C2
    
    # data from previous cycle, where t[ij][k] = t[ij_pc][k] for all k
    data_pc <- log(temp[match(t[ij_1][-1], t[ij_pc]),4]) 
    
  } else {
    
    # linear interpolation from data in Cycle 2, produce NA for end points
    int <- approxfun(temp$standday[which(temp$CYCLE == 2)], temp[ij_pc,4]) 
    
    # cubic spline interpolation
    # int_sp <- splinefun(t[ij_pc], temp[ij_pc,4]) # cubic spline interpolation
    
    data_pc <- log(int(temp$standday[which(temp$CYCLE == 3)][1:(length(ij)-1)]))
    
    if (sum( is.na( data_pc ) ) > 0) {
      # check the position of NA
      na.pos <- which(is.na(data_pc))

      # obtain interpolation for NA produced from the first cycle, hopefully
      # Still need to come up with ways to interpolations if neither of the past two cycle data works
      int.na <- approxfun(temp$standday[which(temp$CYCLE == 1)], temp[ij_pc1,4])

      # cubic spline interpolation
      # int_sp <- splinefun(t[ij_pc], temp[ij_pc,4]) # cubic spline interpolation

      data_na <- log(int.na(temp$standday[which(temp$CYCLE == 3)][na.pos]))
      
      # merge the interpolation for NA with the rest of Non-NA interpolation produced from 2nd cycle.
      data_pc[na.pos] <- data_na
    }
  }
  
  mu_pc[cl, i] <- res$betaHat[1]*tAge + res$betaHat[2]*uW + res$betaHat[3]*oW + f_ij + alpha*(c(log(data[ij_1[1],4]),data_pc) - res$betaHat[1]*tAge - res$betaHat[2]*uW - res$betaHat[3]*oW - f_ij1)
  
  pse_sub[2,i] <- mean((mu_pc[cl, i] - log(data[ij, 4]))^2)
  
  
  #--------------------------------------------------------------------------
  # prediction results conditional (mostly) on cycle 1
  
  # Obtain observations from C1, at time points from C3
  if (all(t[ij] %in% t[ij_pc1])) { # Check if time points in C3 is a subset of C1
    
    # data from previous cycle, where t[ij][k] = t[ij_pc][k] for all k
    data_pc1 <- log(temp[match(t[ij_1][-1], t[ij_pc1]),4]) 
  } else {
    # linear interpolation from data in Cycle 1, produce NA for end points
    int1 <- approxfun(temp$standday[which(temp$CYCLE == 1)], temp[ij_pc1,4]) 
    
    data_pc1 <- log(int1(temp$standday[which(temp$CYCLE == 3)][1:(length(ij)-1)]))
    
    if (sum( is.na( data_pc1 ) ) > 0) {
      # check the position of NA
      na.pos1 <- which(is.na(data_pc1))

      # obtain interpolation for NA produced from cycle 2, hopefully
      # Still need to come up with ways to interpolations if neither of the past two cycle data works
      int.na1 <- approxfun(temp$standday[which(temp$CYCLE == 2)], temp[ij_pc,4])

      # cubic spline interpolation
      # int_sp <- splinefun(t[ij_pc], temp[ij_pc,4]) # cubic spline interpolation

      data_na1 <- log(int.na1(temp$standday[which(temp$CYCLE == 3)][na.pos]))
      
      # merge the interpolation for NA with the rest of Non-NA interpolation produced from cycle 1.
      data_pc1[na.pos] <- data_na1
    }
  }
  
  
  mu_pc1[cl, i] <- res$betaHat[1]*tAge + res$betaHat[2]*uW + res$betaHat[3]*oW + f_ij + alpha*(c(log(data[ij_1[1],4]),data_pc1) - res$betaHat[1]*tAge - res$betaHat[2]*uW - res$betaHat[3]*oW - f_ij1)
  
  pse_sub[3,i] <- mean((mu_pc1[cl, i] - log(data[ij, 4]))^2)
  #--------------------------------------------------------------------------
  # prediction est conditional on predictive prior data
  mu_pre[1, i] <- res$betaHat[1]*tAge + res$betaHat[2]*uW + res$betaHat[3]*oW + f_ij[1] + alpha[1]*(log(data[ij_1[1],4]) - res$betaHat[1]*tAge - res$betaHat[2]*uW - res$betaHat[3]*oW - f_ij1[1])
  for (j in 2:length(ij)) {
    mu_pre[j, i] <- res$betaHat[1]*tAge + res$betaHat[2]*uW + res$betaHat[3]*oW + f_ij[j] + alpha[j]*(mu_pre[j-1, i] - res$betaHat[1]*tAge - res$betaHat[2]*uW - res$betaHat[3]*oW - f_ij1[j])
  }

  pse_sub[4,i] <- mean((mu_pre[cl, i] - log(data[ij, 4]))^2)
  
  #----------------------------------------------------------------
  # Weighted version of past two cycles
  
  mu25[cl, i] <- 0.25*mu_pre[cl, i] + 0.75*mu_pc[cl, i]
  mu50[cl, i] <- 0.50*mu_pre[cl, i] + 0.50*mu_pc[cl, i]
  mu75[cl, i] <- 0.75*mu_pre[cl, i] + 0.25*mu_pc[cl, i]
  
  pse_sub[5,i] <- mean((mu25[cl, i] - log(data[ij, 4]))^2)
  pse_sub[6,i] <- mean((mu50[cl, i] - log(data[ij, 4]))^2)
  pse_sub[7,i] <- mean((mu75[cl, i] - log(data[ij, 4]))^2)
}
# Record the average pse over all m subjects for the dataset
pse <- rowMeans(pse_sub)
round(pse, 4)
pse_sub

#--------------------------------------------------------------------------
# Index identifying subj's with the smallest, second smallest, median and largest PSE values
sortpse <- sort(pse_sub[2,])
pseIn <- c(which.min(pse_sub[2,]), match(sortpse[5], pse_sub[2,]), match(sortpse[25], pse_sub[2,]), match(sortpse[44], pse_sub[2,]))

# Index identifying subj's with the smallest, 10%, median and 90% PSE values
# pseIn <- c(which.min(pse_sub[2,]), match(quantile(pse_sub[2,], .1), pse_sub[2,]), match(quantile(pse_sub[2,], .5), pse_sub[2,]), quantile(pse_sub[2,], .9))


mu <- cbind(mu_true[,pseIn], mu_naive[,pseIn], mu_pc[,pseIn], mu_pc1[,pseIn], mu_pre[,pseIn], mu25[,pseIn], mu50[,pseIn], mu75[,pseIn], CI_lower[,pseIn], CI_upper[,pseIn])
sday <- standay[,pseIn]


# index of non-na standdaysay
si <- sapply(1:4, function(v) return(1:sum(!is.na(sday[,v]))))

# plot smallest, second smallest, median and largest PSE values
dev.new(width=8, height=2)
par(mfrow=c(2,2))

matplot(sday[si[[1]], 1], mu[si[[1]], seq(1, 40, by = 4)], type="l", lty = 1:10, col=c("black", "pink", "blue", "red", "green", "orange", "yellow", "purple", "gray", "gray"), ylim = c(-5, 5), xlab = "Days in the 3rd standardized cycle", ylab = "Log progesterone", xaxt="n", main = "Smallest PSE ")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(sday[si[[2]], 2], mu[si[[2]], seq(2, 40, by = 4)], type="l", lty = 1:10, col=c("black", "pink", "blue", "red", "green", "orange", "yellow", "purple", "gray", "gray"), ylim = c(-5, 5), xlab = "Days in the 3rd standardized cycle", ylab = "Log progesterone", xaxt="n", main = "10% Smallest PSE ")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(sday[si[[3]], 3], mu[si[[3]], seq(3, 40, by = 4)], type="l", lty = 1:10, col=c("black", "pink", "blue", "red", "green", "orange", "yellow", "purple", "gray", "gray"), ylim = c(-5, 5), xlab = "Days in the 3rd standardized cycle", ylab = "Log progesterone", xaxt="n", main = "Median PSE ")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(sday[si[[4]], 4], mu[si[[4]], seq(4, 40, by = 4)], type="l", lty = 1:10, col=c("black", "pink", "blue", "red", "green", "orange", "yellow", "purple", "gray", "gray"), ylim = c(-5, 5), xlab = "Days in the 3rd standardized cycle", ylab = "Log progesterone", xaxt="n", main = "90th PSE ")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

#--------------------------------------------------------------------------

# return the first of middle elements when the vector is even number
n.median = function(x) {
  if (length(x) %% 2 != 0) {
    x <- median(x)
  } else if (length(x) %% 2 == 0) {
    a <- sort(x)[length(x)/2]
    x <- a[1]
  }
  return(x)
}


#======================================================================================
# prediction for subject id 77, 55 obs'ns, seed = 42
#======================================================================================
theta <- c(rep(1, 3), 0.1, -3, -5, -2) # NOU

t <- (data[1:55, 3] %% 28)/10
ij <- 31:55
ij_1 <- 30:54

# # point estimates based on true prior data, most accurate
# # prediction for immediate next point, id 77
# alpha31 <- exp(abs(t[31] - t[30])*log(theta[4]))
# mu31 <- res$betaHat[1]*(-0.005) + res$fHat[1] + alpha31*(data[30,4] - res$betaHat[1]*(-0.005) - res$fHat[28])
# diff31 <- data[31,4] - mu31
# 
# # prediction for 2nd cycle, for subject id 77
# 
alpha77 <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))
mu77 <-  res$betaHat[1]*(-0.005) + res$fHat[c(t[31:54]*10, 28)] + alpha77*(data[ij_1,4] - res$betaHat[1]*(-0.005) - res$fHat[c(28,t[31:54]*10)])
pse <- mean((mu77 - data[ij, 4])^2)
# matplot(1:25, cbind(data[ij, 4], mu77), type="l", lty = 1:2, col=c("black", "red"))