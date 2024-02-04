# Zhang et al 2000


#==============================================================================
# Initialization
#==============================================================================
n = 30 # number of subjects
#tp = 30 # 15 time points for one period; 30 time points for two periods
tp = 15 # one period first


# Initialization of variance and smoothing parameters
sigma.b = 1.5 # variance for bi
sigma = 0.8 # variance for epsilon
library(sde) 
# xi0 = 0.01
# xi1 = 0.02
# xi2 = 0.03
# sigma.ou = exp(xi0 + xi1*t + xi2*t^2)
theta2 = 2
theta3 = 0.5
lambda = 22

# Initialization of regression coefficients
# beta0 = 27.05 # initiation of beta0
beta1 = -0.5 # initiation of beta1




#==============================================================================
# Non-parametric function & other covariates
#==============================================================================

# Periodic function with period length equal to 30 days
# use the cyclic spline, alternatively; and see what I can do..
t = seq(1, 30, by = 2) # 60 for two periods, 30 for one period.
f = 5*sin((2*pi/30)*t)
# plot(t, f, type = "l")


# # Other covariates
# id = rep(1:n, each = tp)
# day = rep(t, n)


#==============================================================================
# Smoothing matrix K - nonnegative definite, Green & Silverman (1994).
# K does not depend on simulation.
#==============================================================================

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


# Check whether R is strictly positive definite

# library(matrixcalc)
# is.positive.definite(R) # True. 


# Define matrix K. Use ginv() instead of solve() for inverse matrix.
K = Q %*% ginv(R) %*% t(Q)

#==============================================================================
# Incidence matrix N
# N does not depend on simulation.
#==============================================================================

Ni = diag(1, tp, tp) 

# N = rbind(Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni, Ni)
# The above is the naive way. The below is simpler way.
# Stack Ni for n subjects.
N = do.call("rbind", replicate(n, Ni, simplify = F))

#==============================================================================
# Covariance matrix V
# Depend on parameter initialization; does not depend on simulation. 
#==============================================================================

Zi = as.matrix(rep(1, tp))

Gammai = matrix(rep(0, tp*tp), tp, tp)
for (i in 1:15) {
	for (j in 1:15) {
		Gammai[i,j] = (theta3^2 / (2*theta2)) * exp(-theta2*abs(i-j))
	}
}


Vi = Zi %*% sigma.b %*% t(Zi) + Gammai + sigma * diag(tp) 

Wi = solve(Vi)

# from internet http://www.r-bloggers.com/block-diagonal-matrices-in-r/
# builds a block matrix whose diagonals are the square matrices provided.
# m1=matrix(runif(10*10),nrow=10,ncol=10)
# m2=matrix(runif(5*5),nrow=5,ncol=5)
# blockMatrix<-blockMatrixDiagonal(m1,m2,m2,m1)
# or
# blockMatrix<-blockMatrixDiagonal(list(m1,m2,m2,m1))
# C.Ladroue
 
blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
 
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
    }
    finalMatrix
  }


# W = blockMatrixDiagonal(Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi, Wi)

# Simpler code.

W = do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))



#==============================================================================
# Covariate matrix X - one covariate, one intercept.
#==============================================================================

# Covariate - age
age = sample(35:60, size = n, replace = T)

		
# Covariate matrix
X = as.matrix(rep(age, each = tp))



#==============================================================================
# The MPLE estimates of beta and f 
#==============================================================================


#==============================================================================
# beta hat
# Depend on both simulation, X and Y, and lambda.
#==============================================================================


beta.est = function(lambda) {
	invWx = solve(t(N) %*% W %*% N + lambda * K)
	Wx = W - W %*% N %*% invWx %*% t(N) %*% W
	betaHat = ginv(t(X) %*% Wx %*% X) %*% t(X) %*% Wx %*% Y
	return(betaHat)
}

#==============================================================================
# f hat 
# Depend on both simulation, X and Y, and lambda.
#==============================================================================



f.est = function(lambda) {
	invWf = solve(t(X) %*% W %*% X) 
	Wf = W - W %*% X %*% invWf %*% t(X) %*% W 
	fHat = ginv(t(N) %*% Wf %*% N + lambda * K) %*% t(N) %*% Wf %*% Y
	return(fHat)
}



 

#==============================================================================
# Simulation
#==============================================================================

# Number of simulation
sim = 100 

# Final simulation results, averge of all simulation rounds. 
beta = matrix(0, 2, sim)
fhat = matrix(0, 15, sim)

# Store simulation estimates into matrices
beta.sim = matrix(0, 2, 1)
f.sim = matrix(0, 15, 1)

# Simulations
simulation = function(n, sim, sigma.b, sigma, theta2, theta3, lambda) {
	for (i in 1:sim) {
		
		# independent random intercepts; sigma.b > sigma; random time effect corresponding to different subject;
		bi = rnorm(n, 0, sigma.b) 
		
		# measurement error - epsilon
		eps.vec = rnorm(n*tp, 0, sigma)
		epsilon = matrix(eps.vec, n, tp)
		
		# Gaussian process
		tempU = rsOU(n*tp, theta=c(0, theta2, theta3))
		u = matrix(tempU, n, tp)
		
		# Longitudinal data Y
		tempResponse = matrix(rep(0, n*tp), n, tp)
		for(i in 1:n) {
			# each row is one subject
			tempResponse[i,] = beta1*age[i] + f + bi[i] + u[i,] + epsilon[i,]
		}
		Y = as.vector(t(tempResponse))
		Y = as.matrix(round(Y))
	
		
		# Beta and f estimates from one simulation
 		beta[,i] = beta.est(lambda)
 		fhat[,i] = f.est(lambda)
	}
	
	# Average of simulation results. 
	beta.sim[1, ] = mean(beta[1, ])
	beta.sim[2, ] = mean(beta[2, ])
	
	for (i in 1:tp) {
		f.sim[i, ] = mean(fhat[i,])
	}
	
	return(list(beta.sim, f.sim))
}



# Note: I do not need to run simulation now, i.e. simulate sim number of Y, thus betaHat and fHat; then average.


# # plot of the response - verify that the response is indeed cyclic
# for (i in 1:n) {
	# plot(t, tempResponse[i,], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	# par(new = TRUE)
# }




#==============================================================================
# Estimation of lambda and variance components 
#==============================================================================

#==============================================================================
# Matrix C
#==============================================================================



# Inverse of Coefficient matrix 
# Depend on simulation X

invC = function(lambda, X) {
	# The 1st row of block matrix C
	C11 = t(X) %*% W %*% X
	C12 = t(X) %*% W %*% N
	CRow1 = cbind(C11, C12)

	# The 2nd row of block matrix C
 	C21 = t(N) %*% W %*% X
	C22 = t(N) %*% W %*% N + lambda * K
	CRow2 = cbind(C21, C22)

	# Inverse of coefficient matrix
	C = rbind(CRow1, CRow2)
	invC = ginv(C)
	return(invC)
}




# # Check whether matrix [X, NT] is of full rank.
# # This is a necessary and sufficient condition for C to be positive definite.
# oneT = rep(1, tp)
# T = cbind(1, t)
# Xstar = cbind(X, N %*% T)

# # Rank of Xstar
# rankMatrix(Xstar)

# Alternatively, can use qr(Xstar)$rank 
# rank = 3; full rank = 4



#==============================================================================
# Matrix P*
#==============================================================================

# Matrix chi in P*
chi = cbind(X, N)

# Matrix P*
PstarFunction = function (lambda) {
	invC = invC(lambda, X)
	Pstar = W - W %*% chi %*% invC %*% t(chi) %*% W
	return(Pstar)
}



#==============================================================================
# Matrix B*
# Does not depend on simulation or lambda
#==============================================================================

# Cholesky decomposition of K to obtain matrix L
# Matrix L is tp * (tp - 2) rull-rank matrix

#> chol(K)
#Error in chol.default(K) : 
#  the leading minor of order 15 is not positive definite

# Use option pivot = T for sparse matrix
LTranspose = chol(K, pivot = T, LDL = T) # upper triangular matrix
L = t(LTranspose) # lower triangular matrix


# Matrix B
B = L %*% ginv((t(L) %*% L))

# Matrix Bstar
Bstar = N %*% B





#==============================================================================
# Score functions for tau and theta
#==============================================================================


# Partial derivative of V with respect to the variance of random intercept b
parV_b = do.call("blockMatrixDiagonal", replicate(n, Zi %*% t(Zi), simplify = F))

# Partial derivative of V with respect to the variance of the measurement error epsilon
parV_eps = do.call("blockMatrixDiagonal", replicate(n, diag(tp), simplify = F))

# Partial derivative of V with respect to theta2 in the OU process
parV_theta2i = matrix(rep(0, tp*tp), tp, tp)
for (i in 1:15) {
	for (j in 1:15) {
		parV_theta2i[i,j] = -(theta3^2 / (2*theta2^2)) * exp(-theta2*abs(i-j)) - abs(i-j)*(theta3^2 / (2*theta2)) * exp(-theta2*abs(i-j)) 
	}
}
parV_theta2 = do.call("blockMatrixDiagonal", replicate(n, parV_theta2i, simplify = F))

# Partial derivative of V with respect to theta3 in the OU process
parV_theta3i = matrix(rep(0, tp*tp), tp, tp)
for (i in 1:15) {
	for (j in 1:15) {
		parV_theta3i[i,j] = (theta3 / theta2) * exp(-theta2*abs(i-j))
	}
}
parV_theta3 = do.call("blockMatrixDiagonal", replicate(n, parV_theta3i, simplify = F))



# To use tr() for trace of a function
library(psych)

# Store score functions to score
score = matrix(0, 5, 1)


# Score function
scoreFunction = function(lambda) {
	Pstar = PstarFunction(lambda)
	
	# score function for tau
	score_tau = -(1/2) * tr(Pstar %*% Bstar %*% t(Bstar)) + (1/2)*t(Y - X %*% betaHat - N %*% fHat) %*% W %*% Bstar %*% t(Bstar) %*% W %*% (Y - X %*% betaHat - N %*% fHat)
	score[1, 1] = score_tau
	
	# score function for the variance of random intercept b
	score_b = -(1/2) * tr(Pstar %*% parV_b) + (1/2)*t(Y - X %*% betaHat - N %*% fHat) %*% W %*% parV_b %*% W %*% (Y - X %*% betaHat - N %*% fHat)
	score[2, 1] = score_b
	
	# score function for the variance of the measurement error epsilon
	score_eps = -(1/2) * tr(Pstar %*% parV_eps) + (1/2)*t(Y - X %*% betaHat - N %*% fHat) %*% W %*% parV_eps %*% W %*% (Y - X %*% betaHat - N %*% fHat)
	score[3, 1] = score_eps

	# score function for theta2
	score_theta2 = -(1/2) * tr(Pstar %*% parV_theta2) + (1/2)*t(Y - X %*% betaHat - N %*% fHat) %*% W %*% parV_theta2 %*% W %*% (Y - X %*% betaHat - N %*% fHat)
	score[4, 1] = score_theta2
	
	# score function for theta3
	score_theta3 = -(1/2) * tr(Pstar %*% parV_theta3) + (1/2)*t(Y - X %*% betaHat - N %*% fHat) %*% W %*% parV_theta3 %*% W %*% (Y - X %*% betaHat - N %*% fHat)
	score[5, 1] = score_theta3

	return(score)
}



#==============================================================================
# Information matrices for tau and theta
#==============================================================================

# Store all entries of Fisher Info matrix into fisherInfo
fisherInfo = matrix(0, 5, 5)

# Information matrix 
infoFunction = function(lambda) {
	Pstar = PstarFunction(lambda)
	
	# Informatin matrix of tau 
	info_tau = (1/2) * tr(Pstar %*% Bstar %*% t(Bstar) %*% Pstar %*% Bstar %*% t(Bstar))
	fisherInfo[1, 1] = info_tau
	
	# Information matrix of tau and theta
	info_taub = (1/2) * tr(Pstar %*% Bstar %*% t(Bstar) %*% Pstar %*% parV_b)
	fisherInfo[1, 2] = fisherInfo[2, 1] = info_taub
	
	info_tauEps = (1/2) * tr(Pstar %*% Bstar %*% t(Bstar) %*% Pstar %*% parV_eps)
	fisherInfo[1, 3] = fisherInfo[3, 1] = info_tauEps
	
	info_tauTheta2 = (1/2) * tr(Pstar %*% Bstar %*% t(Bstar) %*% Pstar %*% parV_theta2)
	fisherInfo[1, 4] = fisherInfo[4, 1] = info_tauTheta2
		
	info_tauTheta3 = (1/2) * tr(Pstar %*% Bstar %*% t(Bstar) %*% Pstar %*% parV_theta3)
	fisherInfo[1, 5] = fisherInfo[5, 1] = info_tauTheta3
	
	
	# Information matrices of var_b and theta
	info_bb = (1/2) * tr(Pstar %*% parV_b %*% Pstar %*% parV_b)
	fisherInfo[2, 2] = info_bb

	info_bEps = (1/2) * tr(Pstar %*% parV_b %*% Pstar %*% parV_eps)
	fisherInfo[2, 3] = info_bEps
	
	info_bTheta2 = (1/2) * tr(Pstar %*% parV_b %*% Pstar %*% parV_theta2)
	fisherInfo[2, 4] = info_bTheta2
	
	info_bTheta3 = (1/2) * tr(Pstar %*% parV_b %*% Pstar %*% parV_theta3)
	fisherInfo[2, 5] = info_bTheta3
	
	
	# Information matrices of var_eps and theta
	info_epsb = (1/2) * tr(Pstar %*% parV_eps %*% Pstar %*% parV_b)
	fisherInfo[3, 2] = info_epsb

	info_epsEps = (1/2) * tr(Pstar %*% parV_eps %*% Pstar %*% parV_eps)
	fisherInfo[3, 3] = info_epsEps
	
	info_epsTheta2 = (1/2) * tr(Pstar %*% parV_eps %*% Pstar %*% parV_theta2)
	fisherInfo[3, 4] = info_epsTheta2
	
	info_epsTheta3 = (1/2) * tr(Pstar %*% parV_eps %*% Pstar %*% parV_theta3)
	fisherInfo[3, 5] = info_epsTheta3
	
	 
	# Information matrices of theta2 and theta
	info_theta2b = (1/2) * tr(Pstar %*% parV_theta2 %*% Pstar %*% parV_b)
	fisherInfo[4, 2] = info_theta2b

	info_theta2Eps = (1/2) * tr(Pstar %*% parV_theta2 %*% Pstar %*% parV_eps)
	fisherInfo[4, 3] = info_theta2Eps


	info_theta2Theta2 = (1/2) * tr(Pstar %*% parV_theta2 %*% Pstar %*% parV_theta2)
	fisherInfo[4, 4] = info_theta2Theta2


	info_theta2Theta3 = (1/2) * tr(Pstar %*% parV_theta2 %*% Pstar %*% parV_theta3)
	fisherInfo[4, 5] = info_theta2Theta3
	
	
	 
	# Information matrices of theta3 and theta
	info_theta3b = (1/2) * tr(Pstar %*% parV_theta3 %*% Pstar %*% parV_b)
	fisherInfo[5, 2] = info_theta3b

	info_theta3Eps = (1/2) * tr(Pstar %*% parV_theta3 %*% Pstar %*% parV_eps)
	fisherInfo[5, 3] = info_theta3Eps

	info_theta3Theta2 = (1/2) * tr(Pstar %*% parV_theta3 %*% Pstar %*% parV_theta2)
	fisherInfo[5, 4] = info_theta3Theta2
	
	info_theta3Theta3 = (1/2) * tr(Pstar %*% parV_theta3 %*% Pstar %*% parV_theta3)
	fisherInfo[5, 5] = info_theta3Theta3

	return(fisherInfo)
}






#==============================================================================
# Fisher scoring algorithm for tau and theta
#==============================================================================

iternum = 1
tau = 1/lambda

lambda0 = matrix(tau, sigma.b, sigma, theta2, theta3)

sigma.b = 1.5 # variance for bi
sigma = 0.8 # variance for epsilon
library(sde) 
# xi0 = 0.01
# xi1 = 0.02
# xi2 = 0.03
# sigma.ou = exp(xi0 + xi1*t + xi2*t^2)
theta2 = 2
theta3 = 0.5
lambda = 22

# Fisher scoring iterations for smoothing parameter lambda
fisher.scoring = function(lambda0, eps = rep(0.01, 5)) {
	lambda = lambda0
	diff = 1
	while(diff > eps) {
		
		iternum = iternum + 1
		
		lambda.old = lambda 
		betaHat = beta.est(1/lambda[1, 1])
		fHat = f.est(1/lambda[1, 1])
		Pstar = PstarFunction(1/lambda[1, 1])
		
		score = scoreFunction(1/lambda[1, 1]) 	# score function
		info = infoFunction(1/lambda[1, 1])		# information matrix 
		
		lambda = lambda + solve(info) %*% score
		# lambda
		diff = abs(lambda - lambda.old)
	}
	list(lambda)
}






