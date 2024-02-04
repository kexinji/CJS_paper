# Bivaraite simulations, one cycle, known variances & smoothing parameters

#==============================================================================
# Initialization
#==============================================================================

library(sde) # rsOU
library(MASS) # mvrnorm

n = 30 # number of subjects
tp = 28 # 28 time points for one period

#------------------------------------------------------------------------------
# Variance initialization 

# theta0 = c(lambda1, lambda2, phi1, phi2, phi3, sigma1, sigma2, theta12, theta13, theta22, theta23)

par = c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)

D = matrix(c(par[1], par[2], par[2], par[3]), 2, 2) # variance for bi
var_eps = c(par[4], par[5])
sigma = diag(var_eps) # variance for epsilon
# theta12 = 2 # par[6]
# theta13 = 3 # par[7]
# theta22 = 2 # par[8]
# theta23 = 5 # par[9]


# # Smoothing parameters initialization
# lambda1 = 0.15 # par[10]
# lambda2 = 2 # par[11]

#------------------------------------------------------------------------------
# periodic function with period length equal to 30 days
t = seq(1, 28) # 
f1 = 5*sin((2*pi/tp)*t)
f2 = 3*cos((2*pi/tp)*t)
#plot(t, f_1, type = "l")
plot(t, f2, type = "l")

#------------------------------------------------------------------------------
# X
set.seed(34)
age = sample(20:44, size = n, replace = T)

X = matrix(0, 2*tp*n, 2)
for (i in 1:n) {
	# the first column of 1st subject of Xi is age[1], 0, age[2], 0...
	X[((2*tp)*(i-1)+1):(2*tp*i), 1] = rep(c(age[i], 0), tp) 	
	
	# the second column of 1st subject of Xi is 0, age[1], 0, age[2]...
	X[((2*tp)*(i-1)+1):(2*tp*i), 2] = rep(c(0, age[i]), tp) 
}


#------------------------------------------------------------------------------
# Initializaiton of beta1
beta11 = 0.03
beta21 = 0.07


#==============================================================================
# Matrices that does not depend on Y
#==============================================================================

# A function that return a block matrix, code found online.
blockMatrixDiagonal <- function(...) {
	matrixList <- list(...)
	if(is.list(matrixList[[1]])) matrixList <- matrixList[[1]]
	
	dimensions <- sapply(matrixList, FUN = function(x) dim(x)[1])
	finalDimension <- sum(dimensions)
  	finalMatrix<-matrix(0, nrow=finalDimension, ncol=finalDimension)
	index <- 1
	for (k in 1:length(dimensions)) {
		finalMatrix[index:(index+dimensions[k]-1), index:(index + dimensions[k] -1)] <- matrixList[[k]]
		index <- index + dimensions[k]
	}
	finalMatrix
}



# the variance 
Zi = matrix(0, 2*tp, 2) # initialize matrix Zi
Zi[ ,1] = rep(c(1, 0), tp) # the first column of Zi is 1, 0, 1, 0...
Zi[ ,2] = rep(c(0, 1), tp) # the second column of Zi is 0, 1, 0, 1...
Z  = do.call("rbind", replicate(n, Zi, simplify = F))


sigma_i = diag(rep(var_eps, tp))

gamma_i = matrix(0, 2*tp, 2*tp)
for (i in 1:tp) {
	gamma_i[2*i - 1, ] = rep(c(9/4, 0), tp) 
	gamma_i[2*i, ] = rep(c(0, 25/4), tp)
	}

Vi = Zi %*% D %*% t(Zi) + sigma_i + gamma_i
Wi = solve(Vi)
V = do.call("blockMatrixDiagonal", replicate(n, Vi, simplify = F))
W = do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))





# Matrix K, does not depend on theta

# Define matrix Q, double checked.
Q = matrix(0, tp, tp-2)
for(j in 1:(tp-2)) {
	Q[j, j] = 1/2
	Q[j + 1, j] = -1
	Q[j + 2, j] = 1/2
}

# Define matrix R 
R = matrix(0, tp-2, tp-2)
for (i in 1:(tp-2)) {
	R[i, i] = (1/3)*4
}
for (i in 1:(tp-3)) {
	R[i, i + 1] = (1/6)*2
	R[i + 1, i] = (1/6)*2
}

K = Q %*% ginv(R) %*% t(Q)




# Incidence matrices N1 & N2

	
Ni <- diag(28)
A1i <- matrix(0, 2*tp, tp) # initialize matrix A1i
A2i <- matrix(0, 2*tp, tp) # initialize matrix A2i
for (j in 1:tp){
	A1i[2*j-1, j] <- 1
	A2i[2*j, j] <- 1 
}
		
N1i <- A1i %*% Ni
N2i <- A2i %*% Ni

N1  = do.call("rbind", replicate(n, N1i, simplify = F))
N2  = do.call("rbind", replicate(n, N2i, simplify = F))


	




# NTemp = diag(1, tp, tp) 
# z = matrix(0, tp, tp) # matrix of zeros 
# N1i = matrix(0, 2*tp, tp) 
# for(i in 1:tp) {
	# N1i[(2*i-1):(2*i),] = rbind(NTemp[i, ], z[i, ])
# }
# N2i = matrix(0, 2*tp, tp) 
# for(i in 1:tp) {
	# N2i[(2*i-1):(2*i),] = rbind(z[i, ], NTemp[i, ])
# }
# N1  = do.call("rbind", replicate(n, N1i, simplify = F))
# N2  = do.call("rbind", replicate(n, N2i, simplify = F))



# INVERSE of MATRIX C

lambda1 = 0.15
lambda2 = 2

# The 1st row of matrix C
C11 = t(X) %*% W %*% X
C12 = t(X) %*% W %*% N1
C13 = t(X) %*% W %*% N2
CRow1 = cbind(C11, C12, C13)
# The 2nd row of matrix C
C21 = t(N1) %*% W %*% X
C22 = t(N1) %*% W %*% N1 + lambda1 * K
C23 = t(N1) %*% W %*% N2
CRow2 = cbind(C21, C22, C23)
# The 2nd row of matrix C
C31 = t(N2) %*% W %*% X
C32 = t(N2) %*% W %*% N1
C33 = t(N2) %*% W %*% N2 + lambda2 * K
CRow3 = cbind(C31, C32, C33)
# Inverse of coefficient matrix
C = rbind(CRow1, CRow2, CRow3)
invC = ginv(C)




#==============================================================================
# Simulation
#==============================================================================

sim = 100 # number of simulation

# Store simulation results, each column corresponds to a simulation result. 
betaHat = matrix(0, 2, sim)
f1Hat = matrix(0, tp, sim)
f2Hat = matrix(0, tp, sim)
bHat = matrix(0, 2, sim)

# Store the final simulation estimates into matrices, average of beta and fhat
beta.sim = matrix(0, 2, 1)
f1.sim = matrix(0, tp, 1)
f2.sim = matrix(0, tp, 1)

for (j in 1:sim) {
	
	# random intercepts;
	b = mvrnorm(n, c(0, 0), D) 

	# measurement error - epsilon
	eps = mvrnorm(n*tp, c(0, 0), sigma)
	eps1 = matrix(eps[,1], n, tp)
	eps2 = matrix(eps[,2], n, tp)
	
	# bivariate Gaussian field
	tempU_1 = rsOU(n*tp, theta=c(0, par[6], par[7])) # theta12 = 2, theta13 = 3
	tempU_2 = rsOU(n*tp, theta=c(0, par[8], par[9])) # theta22 = 2, theta23 = 5
	u1 = matrix(tempU_1, n, tp)
	u2 = matrix(tempU_2, n, tp)
	
	# bivariate Longitudinal data Y
	tempResponse = matrix(0, 2*n, tp)
	tempR = matrix(0, tp, 2) # for each subject
	bivY = matrix(0, n*tp, 2)
	Y <- matrix(0, 2*n*tp, 1)
	for (i in 1:n) {
		# every two rows are one subject
		tempResponse[(2*i-1),] = beta11*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
		tempResponse[2*i,] = beta21*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
		tempR <- t(tempResponse[(2*i-1):(2*i),])
		bivY[((i-1)*tp+1):(i*tp), ] <- tempR
	}
	
	for (i in 1:(n*tp)) {
		Y[2*i - 1] = bivY[i,1]
		Y[2*i] = bivY[i,2]
	} 
	
	
	# ESTIMATES FOR beta, f1 AND f2
	temp = rbind(t(X) %*% W %*% Y, t(N1) %*% W %*% Y, t(N2) %*% W %*% Y)
	est = invC %*% temp
	betaHat[,j] = est[1:2] 
	f1Hat[,j] = matrix(est[3:30])
	f2Hat[,j] = matrix(est[31:58])
	
	
	res <- Y - X %*% betaHat[,j] - N1 %*% f1Hat[,j] - N2 %*% f2Hat[,j]
	b <- D %*% t(Z) %*% W %*% res


}


#==============================================================================
# Simulation ESTIMATES, beta hat, f hat, PLOTS
#==============================================================================

# beta Hat
beta.sim[1] = mean(betaHat[1,])
beta.sim[2] = mean(betaHat[2,])


# f Hat
for (k in 1:tp) {
	f1.sim[k, ] = mean(f1Hat[k,])
	f2.sim[k, ] = mean(f2Hat[k,])
}

# plots of f Hat and Beta Hat
par(mfrow = c(2, 2))
plot(f1.sim, type = "l", ylab = "f1 Hat", xlab = "days")
plot(f2.sim, type = "l", ylab = "f2 Hat", xlab = "days")
hist(betaHat[1,], xlab = expression(beta_1), main = "Histogram of Beta_1")
hist(betaHat[2,], xlab = expression(beta_2), main = "Histogram of Beta_2")


# # plot of f hat for all simulations 
# for (m in 1:sim) {
	# plot(f1Hat[,m], type = "l", ylim = c(-10, 10), ylab = "f1Hat")
	# par(new = T)
# }

# for (m in 1:sim) {
	# plot(f2Hat[,m], type = "l", ylim = c(-10, 10), ylab = "f2Hat")
	# par(new = T)
# }


#==============================================================================
# BIASES and COVARIANCES
#==============================================================================

# Weight matrices
W1inv <- ginv(C22)
W2inv <- ginv(C33)

a1 <- W %*% N1
W1 <- W - a1 %*% W1inv %*% t(a1)
a2 <- W %*% N2
W2 <- W - a2 %*% W2inv %*% t(a2)


invX <- t(N2) %*% W1 %*% N2 + lambda2 * K
WxInv <- ginv(invX)
invF1 <- t(X) %*% W2 %*% X
Wf1Inv <- ginv(invF1)
invF2 <- t(X) %*% W1 %*% X
Wf2Inv <- ginv(invF2)

Wx <- W1 - W1 %*% N2 %*% WxInv %*% t(N2) %*% W1
Wf1 <- W2 - W2 %*% X %*% Wf1Inv %*% t(X) %*% W2
Wf2 <- W1 - W1 %*% X %*% Wf2Inv %*% t(X) %*% W1


# BIASES & COVARIANCES - beta, f1 AND f2

invB <- t(X) %*% Wx %*% X
betaInv <- ginv(invB) 
hatB <-  betaInv %*% t(X) %*% Wx
biasBeta<- hatB %*% (N1 %*% f1 + N2 %*% f2)
covBeta <- hatB %*% V %*% t(hatB)

invF1 <- t(N1) %*% Wf1 %*% N1 + lambda1 * K
f1Inv <- solve(invF1)
hatF1 <- f1Inv %*% t(N1) %*% Wf1
covF1 <- hatF1 %*% V %*% t(hatF1)

invF2 <- t(N2) %*% Wf2 %*% N2 + lambda2 * K
f2Inv <- solve(invF1)
hatF2 <- f2Inv %*% t(N2) %*% Wf2
covF2 <- hatF2 %*% V %*% t(hatF2)








# findPar = function(sim, theta) {
	
	# for (j in 1:sim) {
		
		# # random intercepts;
		# b = mvrnorm(n, c(0, 0), D) 

		# # measurement error - epsilon
		# eps = mvrnorm(n*tp, c(0, 0), sigma)
		# eps1 = matrix(eps[,1], n, tp)
		# eps2 = matrix(eps[,2], n, tp)
		
		# # bivariate Gaussian field
		# tempU_1 = rsOU(n*tp, theta=c(0, par[6], par[7])) # theta12 = 2, theta13 = 3
		# tempU_2 = rsOU(n*tp, theta=c(0, par[8], par[9])) # theta22 = 2, theta23 = 5
		# u1 = matrix(tempU_1, n, tp)
		# u2 = matrix(tempU_2, n, tp)
		
		# # bivariate Longitudinal data Y
		# tempResponse = matrix(rep(0, 2*n*tp), 2*n, tp)
		# for (i in 1:n) {
			# # every two rows are one subject
			# tempResponse[(2*i-1),] = beta11*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
			# tempResponse[2*i,] = beta21*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
		# }
		# Y = as.vector(t(tempResponse))
		# Y = as.matrix(round(Y))
 
		
		# # ESTIMATES FOR beta, f1 AND f2
		# temp = rbind(t(X) %*% W %*% Y, t(N1) %*% W %*% Y, t(N2) %*% W %*% Y)
		# est = invC %*% temp
		# betaHat[,j] = est[1:2] 
		# f1Hat[,j] = matrix(est[3:30])
		# f2Hat[,j] = matrix(est[31:58])
		
		
		# res <- Y - X %*% betaHat[,j] - N1 %*% f1Hat[,j] - N2 %*% f2Hat[,j]
		# b <- D %*% t(Z) %*% W %*% res
	
	
	# }
	
	# # beta Hat
	# beta.sim[1] = mean(betaHat[1,])
	# beta.sim[2] = mean(betaHat[2,])
	
	
	# # f Hat
	# for (k in 1:tp) {
		# f1.sim[k, ] = mean(f1Hat[k,])
		# f2.sim[k, ] = mean(f2Hat[k,])
	# }
	
	# # histogram for beta Hat
	# par(mfrow = c(2, 1))
	# hist(betaHat[1,], xlab = expression(beta_11), main = "Histogram of Beta_1")
	# hist(betaHat[2,], xlab = expression(beta_21), main = "Histogram of Beta_2")
	
	# # plot of f hat for all simulations 
	# for (m in 1:sim) {
		# plot(f1Hat[,m], type = "l", ylim = c(-10, 10), ylab = "f1Hat")
		# par(new = T)
	# }
	
	# for (m in 1:sim) {
		# plot(f2Hat[,m], type = "l", ylim = c(-10, 10), ylab = "f2Hat")
		# par(new = T)
	# }
	
	# # plot of f hat
	# plot(f1.sim, type = "l")
	# par(new = T)
	# plot(f2.sim, type = "l")
		
	# return(list(beta.sim, f1.sim, f2.sim, betaHat, f1Hat, f2Hat))
	
# }








