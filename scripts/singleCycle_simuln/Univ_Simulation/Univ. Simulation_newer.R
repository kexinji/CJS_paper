# Zhang et al 2000

#===============================================================================
# PART I - DATA GENERATION
#===============================================================================

n = 30 # number of subjects
tp = 15 # one period first, 30 time points for two periods

#------------------------------------------------------------------------------
# Parameter initializaiton

# Variance initializations
sigma.b = 2 # variance for bi, initially 1.5
sigma =  1 # variance for epsilon, initially 2
theta2 = 2
theta3 = 0.5
lambda = 0.1389352  # based on Zhang et al(1998), lambda = 1/tau

tau = 1/lambda


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
# PART II - PARAMETER & VARIANCE ESTIMATION
#===============================================================================

#------------------------------------------------------------------------------
# A function that gives block matrix, found online.

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
# Matrices that does not depend on theta

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

# Matrix N, incidence matrix
Ni = diag(1, tp, tp) 
N = do.call("rbind", replicate(n, Ni, simplify = F))


# B*

# Matrix B
LTranspose = chol(K, pivot = T, LDL = T) # upper triangular matrix
L = t(LTranspose) # lower triangular matrix
B = L %*% ginv((t(L) %*% L))

# Bstar
Bstar = N %*% B



#------------------------------------------------------------------------------
# THE MAIN FUNCTION
findPar = function(theta) {
	
	# VARIANCE FOR OU PROCESS
	Gammai = matrix(rep(0, tp*tp), tp, tp)
	for (i in 1:15) {
		for (j in 1:15) {
			Gammai[i,j] = (theta[5]^2 / (2*theta[4])) * exp(-theta[4]*abs(i-j))
		}
	}
	
	Vi = Zi %*% theta[2] %*% t(Zi) + Gammai + theta[3] * diag(tp)
	Wi = solve(Vi)
	W = do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))
	
	# INVERSE of MATRIX C
	# The 1st row of block matrix C
	C11 = t(X) %*% W %*% X
	C12 = t(X) %*% W %*% N
	CRow1 = cbind(C11, C12)
	# The 2nd row of block matrix C
 	C21 = t(N) %*% W %*% X
	C22 = t(N) %*% W %*% N + (1/theta[1]) * K
	CRow2 = cbind(C21, C22)
	# Inverse of coefficient matrix
	C = rbind(CRow1, CRow2)
	invC = ginv(C)
	
	# ESTIMATES FOR beta AND f
	temp = rbind(t(X) %*% W %*% Y, t(N) %*% W %*% Y)
	estimates = invC%*% temp
	betaHat = estimates[1] 
	fHat = matrix(estimates[2:16])
	
	# P*
	chi = cbind(X, N)	
	Pst = W - W %*% chi %*% invC %*% t(chi) %*% W

	# PARTIAL DERIVATIVES
	# Partial derivative of V wrt the variance of b
	parV_b = do.call("blockMatrixDiagonal", replicate(n, Zi %*% t(Zi), simplify = F))
	# Partial derivative of V wrt to the variance of epsilon
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
	parV = cbind(parV_b, parV_eps, parV_theta2, parV_theta3)
	
	# Values needed in score, py = W * (Y - x*betaHat - N * fHat); tpy 
	py = W %*% (Y - X %*% betaHat - N %*% fHat)
	tpy = t(Y - X %*% betaHat - N %*% fHat) %*% 	W
	
	# SCORE; FISHER INFO
	score = matrix(0, 5, 1) 
	fInfo = matrix(0, 5, 5)

	library(psych) # to use function tr()
	Bdou = Bstar %*% t(Bstar) 

	score[1] = tpy %*% Bdou %*% py - tr(Pst %*% Bdou)
	for (i in 2:5) {
		score[i] = tpy %*% parV[, (450*(i-2)+1): (450*(i-1))] %*% py - tr(Pst %*% parV[, (450*(i-2)+1): (450*(i-1))])
	}
	
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

	thetaNew = theta + Finv %*% score	# Parameter estimate after one iteration
	
	
	newlist = list("thetaNew" = theta + Finv %*% score, "betaHat" = estimates[1], "fHat" = matrix(estimates[2:16]))
	return(newlist)
}

#------------------------------------------------------------------------------
# FISHER SCORING ITERATION

iternum = 0
diff = matrix(rep(1, 5), 1, 5)
tol = 0.005 # take 52 iterations for tol = 0.01

theta = matrix(c(tau, sigma.b, sigma, theta2, theta3), 5, 1)

theta1 = theta # initialization

fisher.scoring = function(theta1) {

while (diff[1] > tol | diff[2] > tol | diff[3] > tol | diff[4] > tol | diff[5] > tol) {
	iternum = iternum + 1
	
    theta0 = theta1
	theta1 = findPar(theta0)$thetaNew
	
	for (i in 1:5){
		while (theta1[i]<0) {
			theta1[i] = (theta1[i] + theta0[i])/2
		}
	}
	
	for (i in 1:5) {
		diff[i] = abs(theta1[i] - theta0[i])
	}
	
	print(iternum)
	print (diff)
	print(theta1)
}
return(theta1)
}

