# Zhang et al 2000

#===============================================================================
# PART I - DATA GENERATION
#===============================================================================

n = 30 # number of subjects
tp = 15 # one period first, 30 time points for two periods

#------------------------------------------------------------------------------
# Parameter initializaiton

# Variance initializations
sigma.b = 1.5 # variance for bi
sigma = 0.8 # variance for epsilon
theta2 = 2
theta3 = 0.5
lambda = 0.1389352  # based on Zhang et al(1998), lambda = 1/tau

tau = 1/lambda
theta = matrix(c(tau, sigma.b, sigma, theta2, theta3), 5, 1)


# Initializaiton of beta1
beta1 = -0.5
#------------------------------------------------------------------------------


# Nonparametric function f, period = 30 days 
t = seq(1, 30, by = 2) # 60 for two periods, 30 for one period.
f = 5*sin((2*pi/30)*t)

# Covariate matrix - age
age = sample(35:60, size = n, replace = T)
X = as.matrix(rep(age, each = tp))
Zi = as.matrix(rep(1, tp)) # random effect covariates

# Random intercepts; sigma.b > sigma; random time effect corresponding to different subject;
bi = rnorm(n, 0, sigma.b) 

		
# Measurement error - epsilon
eps.vec = rnorm(n*tp, 0, sigma)
epsilon = matrix(eps.vec, n, tp)
		
# Gaussian process
library(sde) 
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




#===============================================================================
# PART II - PARAMETER ESTIMATION, BETA & f 
#===============================================================================

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



#------------------------------------------------------------------------------
# THREE FOUNDAMENTAL MATRICES, W, K and N, TO OBTAIN C

# Covariance matrix V
# Depend on parameter initialization; does not depend on simulation. 

W = function(theta) {
	
	# Variance matrix for OU process
	Gammai = matrix(rep(0, tp*tp), tp, tp)
	for (i in 1:15) {
		for (j in 1:15) {
			Gammai[i,j] = (theta[5]^2 / (2*theta[4])) * exp(-theta[4]*abs(i-j))
		}
	}
	
	Vi = Zi %*% theta[2] %*% t(Zi) + Gammai + theta[3] * diag(tp)
	Wi = solve(Vi)
	invV = do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))
	return(invV)
}




# Matrix K 

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

# Matrix N, incidence matrix
Ni = diag(1, tp, tp) 
N = do.call("rbind", replicate(n, Ni, simplify = F))


# Inverse of Matrix C
invC = function(theta) {
	# The 1st row of block matrix C
	C11 = t(X) %*% W(theta) %*% X
	C12 = t(X) %*% W(theta) %*% N
	CRow1 = cbind(C11, C12)

	# The 2nd row of block matrix C
 	C21 = t(N) %*% W(theta) %*% X
	C22 = t(N) %*% W(theta) %*% N + (1/theta[1]) * K
	CRow2 = cbind(C21, C22)

	# Inverse of coefficient matrix
	C = rbind(CRow1, CRow2)
	Cinv = ginv(C)
	return(Cinv)
}


#------------------------------------------------------------------------------
# Obtain Estimates for beta and f

est = function(theta) {
	temp = rbind(t(X) %*% (W(theta)) %*% Y, t(N) %*% (W(theta)) %*% Y)
	estimates = invC(theta)%*% temp
	newlist = list("betaHat" = estimates[1], "fHat" = matrix(estimates[2:16]))
	return(newlist)
}





#===============================================================================
# PART III - VARIANCE ESTIMATION
#===============================================================================


#------------------------------------------------------------------------------
# TWO MORE MATRICES P* AND B*

# P*
chi = cbind(X, N)
Pstar = function (theta) {
	P = W(theta) - W(theta) %*% chi %*% invC(theta) %*% t(chi) %*% W(theta)
	return(P)
}


# B*

# Matrix B
LTranspose = chol(K, pivot = T, LDL = T) # upper triangular matrix
L = t(LTranspose) # lower triangular matrix
B = L %*% ginv((t(L) %*% L))

# Bstar
Bstar = N %*% B


#------------------------------------------------------------------------------
# PARTIAL DERIVATIVES

par = function(theta) {
	
	# Partial derivative of V with respect to the variance of random intercept b
	parV_b = do.call("blockMatrixDiagonal", replicate(n, Zi %*% t(Zi), simplify = F))

	# Partial derivative of V with respect to the variance of the measurement error epsilon
	parV_eps = do.call("blockMatrixDiagonal", replicate(n, diag(tp), simplify = F))
	
	# Partial derivative of V wrt theta2 in the OU process
	parV_theta2i = matrix(rep(0, tp*tp), tp, tp)
	for (i in 1:15) {
		for (j in 1:15) {
			parV_theta2i[i,j] = -(theta[5]^2 / (2*theta[4]^2)) * exp(-theta[4]*abs(i-j)) - abs(i-j)*(theta[5]^2 / (2*theta[4])) * exp(-theta[4]*abs(i-j)) 
		}
	}
	parV_theta2 = do.call("blockMatrixDiagonal", replicate(n, parV_theta2i, simplify = F))

	
	# Partial derivative of V wrt theta3 in the OU process
	parV_theta3i = matrix(rep(0, tp*tp), tp, tp)
	for (i in 1:15) {
		for (j in 1:15) {
			parV_theta3i[i,j] = (theta[5] / theta[4]) * exp(-theta[4]*abs(i-j))
		}
	}
	parV_theta3 = do.call("blockMatrixDiagonal", replicate(n, parV_theta3i, simplify = F))
	
	parVec = cbind(parV_b, parV_eps, parV_theta2, parV_theta3)
	return(parVec)
}


#------------------------------------------------------------------------------
# SCORE AND FISHER INFO


# Score
score = matrix(0, 5, 1)

# py = W * (Y - x*betaHat - N * fHat); tpy 
p = function(theta){
	betaHat = est(theta)$betaHat
	fHat = est(theta)$fHat
	return(list("py" = W(theta) %*% (Y - X %*% betaHat - N %*% fHat), "tpy" = t(Y - X %*% betaHat - N %*% fHat) %*% W(theta)))
}


library(psych) # to use function tr()

# Store some repeated appearing values
Pst = Pstar(theta)
parV = par(theta)
py = p(theta)$py
tpy = p(theta)$tpy
Bdou = Bstar %*% t(Bstar) 


# score function for tau
score[1] = tpy %*% Bdou %*% py - tr(Pst %*% Bdou)

for (i in 2:5) {
	score[i] = tpy %*% parV[, (450*(i-2)+1): (450*(i-1))] %*% py - tr(Pst %*% parV[, (450*(i-2)+1): (450*(i-1))])
}


# Fisher Info 
fInfo = matrix(0, 5, 5)

# I_tauTau
fInfo[1, 1] = tr(Pst %*% Bdou %*% Pst %*% Bdou)

# I_tauTheta
for (i in 2:5) {
	fInfo[1, i] = tr(Pst %*% Bdou %*% Pst %*% parV[, (450*(i-2)+1): (450*(i-1))])
}

# I_thetaTheta
for (i in 2:5) {
	for (j in 2:5) {
		fInfo[i, j] = tr(Pst %*% parV[, (450*(i-2)+1): (450*(i-1))] %*% Pst %*% parV[, (450*(j-2)+1): (450*(j-1))])
	}
}

score = score/2
fInfo = fInfo/2
Finv = ginv(fInfo)


#------------------------------------------------------------------------------
# FISHER SCORING ITERATION

iternum = 1

theta = theta + Finv %*% score

# For each iteration, need to update P*, C, betaHat, fHat, W, i.e. py and tpy, partial derivatives w.r.t. OU process
