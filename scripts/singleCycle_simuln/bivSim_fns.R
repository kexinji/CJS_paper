# THIS R FILE IS BASED ON bivSim_unknownVar.r, USED FOR PARALLEL COMPUTING.

#==============================================================================
# PART I - INITIALIZATION, AND SET UP 
#==============================================================================

library(sde) # rsOU()
library(MASS) # mvrnorm()
library(psych) # tr()
library(Rcpp) # TO C++
sourceCpp("matrixmulti_1.cpp") # NEED TO INSTALL Package 'RcppEigen'

n <- 30 # number of subjects
tp <- 28 # 28 time points for one period

m = 2*n*tp

#------------------------------------------------------------------------------
# X
set.seed(34)
age <- sample(20:44, size = n, replace = T)

X = matrix(0, 2*tp*n, 2)
for (i in 1:n) {
	# the first column of 1st subject of Xi is age[1], 0, age[2], 0...
	X[((2*tp)*(i-1)+1):(2*tp*i), 1] <- rep(c(age[i], 0), tp) 	
	
	# the second column of 1st subject of Xi is 0, age[1], 0, age[2]...
	X[((2*tp)*(i-1)+1):(2*tp*i), 2] <- rep(c(0, age[i]), tp) 
}

#------------------------------------------------------------------------------
# Initializaiton of fixed effects, beta1
beta11 <- 0.03
beta21 <- 0.07

#------------------------------------------------------------------------------
# INCIDENCE MATRICES N1 & N2
Ni <- diag(tp)
A1i <- matrix(0, 2*tp, tp) # initialize matrix A1i
A2i <- matrix(0, 2*tp, tp) # initialize matrix A2i
for (j in 1:tp){
	A1i[2*j-1, j] <- 1
	A2i[2*j, j] <- 1 
}
		
N1i <- A1i %*% Ni
N2i <- A2i %*% Ni

N1 <- do.call("rbind", replicate(n, N1i, simplify = F))
N2 <- do.call("rbind", replicate(n, N2i, simplify = F))

#------------------------------------------------------------------------------
# Periodic function with period length equal to 28 days
t <- seq(1, 28)  
f1 <- 5*sin((2*pi/tp)*t)
f2 <- 3*cos((2*pi/tp)*t)
# plot(t, f1, type = "l")
# plot(t, f2, type = "l")


#------------------------------------------------------------------------------
# Z
Zi <- matrix(0, 2*tp, 2) # initialize matrix Zi
Zi[ ,1] <- rep(c(1, 0), tp) # the first column of Zi is 1, 0, 1, 0...
Zi[ ,2] <- rep(c(0, 1), tp) # the second column of Zi is 0, 1, 0, 1...
Z <- do.call("rbind", replicate(n, Zi, simplify = F))


#===============================================================================
# PART II - RELEVANT FUNCTIONS
#===============================================================================

#------------------------------------------------------------------------------
# Matrices that does not depend on Y

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


#------------------------------------------------------------------------------
# Matrix K, does not depend on theta; assuming h is equally spaced

matrixK <- function(h, r) {
	# Define matrix Q, double checked.
	Q <- matrix(0, r, r-2)
	for(j in 1:(r-2)) {
		Q[j, j] <- 1/h
		Q[j + 1, j] <- -2/h
		Q[j + 2, j] <- 1/h
	}
	
	# Define matrix R 
	R = matrix(0, r-2, r-2)
	for (i in 1:(tp-2)) {
		R[i, i] <- (2/3)*h
	}
	for (i in 1:(r-3)) {
		R[i, i + 1] <- (1/6)*h
		R[i + 1, i] <- (1/6)*h
	}
	
	K <- eigenMapMatMult3(Q, ginv(R), t(Q))
	
	
	#------------------------------------------------------------------------------
	# MATRICES B1* B2*
	
	# Matrix B
	LTranspose <- chol(K, pivot = T, LDL = T) # upper triangular matrix
	L <- t(LTranspose) # lower triangular matrix
	B <- L %*% ginv(crossprod(L))
	
	# Bstar
	B1star <- N1 %*% B
	B2star <- N2 %*% B
	
	# TO BE USED IN findPar
	B1dou <- tcrossprod(B1star)
	B2dou <- tcrossprod(B2star)


	newlist <- list("K" = K, "B1dou" = B1dou, "B2dou" = B2dou)
	return(newlist)
}


#------------------------------------------------------------------------------
# VARIANCE FUNCTION, OU PROCESS
ouVar <- function(theta2, theta3, x) {
	ouV <- (theta3^2/(2*theta2)) * exp(-theta2*abs(x))
	return(ouV)
}

# PARTIAL DERIVATIVES OF VARIANCE FUNCTION WRT THETA2, OU PROCESS
ouPv1 <- function(theta2, theta3, x) {
	ouP1 <- (-1/theta2 - abs(x)) * ouVar(theta2, theta3, x)
	return(ouP1)
}


# PARTIAL DERIVATIVES OF VARIANCE FUNCTION WRT THETA3, OU PROCESS
ouPv2 <- function(theta2, theta3, x) {
	ouP2 <- (theta3/theta2) * exp(-theta2*abs(x))
	return(ouP2)
}


#------------------------------------------------------------------------------
# A FUNCTION FOR DATA GENERATION
# GENERATE A SET OF NEW DATA FOR EACH SIMULATION.

newData <- function(theta) {
	D <- matrix(c(theta[3], theta[4], theta[4], theta[5]), 2, 2) # variance for bi
	# random intercepts;
	b <- mvrnorm(n, c(0, 0), D) 
	
	
	# measurement error - epsilon
	var_eps <- c(theta[6], theta[7])
	sigma <- diag(var_eps) # variance for epsilon
	
	eps <- mvrnorm(n*tp, c(0, 0), sigma)
	eps1 <- matrix(eps[,1], n, tp)
	eps2 <- matrix(eps[,2], n, tp)
	
	
	# bivariate Gaussian field
	tempU_1 <- rsOU(n*tp, theta=c(0, theta[8], theta[9])) 
	tempU_2 <- rsOU(n*tp, theta=c(0, theta[10], theta[11]))
	u1 <- matrix(tempU_1, n, tp)
	u2 <- matrix(tempU_2, n, tp)
	
	# bivariate Longitudinal data Y
	tempResponse <- matrix(0, 2*n, tp)
	tempR <- matrix(0, tp, 2) # for each subject
	bivY <- matrix(0, n*tp, 2)
	Y <- matrix(0, 2*n*tp, 1)
	for (i in 1:n) {
		# every two rows are one subject
		tempResponse[(2*i-1),] <- beta11*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
		tempResponse[2*i,] <- beta21*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
		tempR <- t(tempResponse[(2*i-1):(2*i),])
		bivY[((i-1)*tp+1):(i*tp), ] <- tempR
	}
	
	for (i in 1:(n*tp)) {
		Y[2*i - 1] <- bivY[i,1]
		Y[2*i] <- bivY[i,2]
	} 

	return(Y)	
}


#------------------------------------------------------------------------------
# A FUNCTION TO COMPUTE ESTIMATES FOR BETA and F, GIVEN theta & Y.

est = function(theta, Y, h, r) {
	
	K <- matrixK(h, r)$K
	
	#--------------------------------------------------------------------------
	# VARIANCE for bi
	D <- matrix(c(theta[3], theta[4], theta[4], theta[5]), 2, 2) 
	
	# VARIANCE OF MEASUREMENT ERROR
	sigma_i <- diag(rep(c(theta[6], theta[7]), tp))
	
	# VARIANCE FOR OU PROCESS
	# USE VARIANCE FUNCTION, ASSIGN ONLY THE UPPER TRIANGLE, TAKE ACCOUNT OF SYMMETRY OF GAMMAi
	
	Gammai <- matrix(0, 2*tp, 2*tp)
	for (i in 1:tp) {
		for (j in i:tp) {
			Gammai[2*i-1, 2*j-1] <- ouVar(theta[8], theta[9], (i-j))
			Gammai[2*i, 2*j] <- ouVar(theta[10], theta[11], (i-j))
			Gammai[2*j-1, 2*i-1] <- t(Gammai[2*i-1, 2*j-1])
			Gammai[2*j, 2*i] <- t(Gammai[2*i, 2*j])
		}
	}
		
	Vi <- Zi %*% D %*% t(Zi) + sigma_i + Gammai
	V <- do.call("blockMatrixDiagonal", replicate(n, Vi, simplify = F))

	Wi <- solve(Vi)
	W <- do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))
	
		
	#--------------------------------------------------------------------------
	# INVERSE of MATRIX C
	
	# HELP IMPROVE EFFICIENCY OF THE CODE
	WX <- W %*% X
	WN1 <- eigenMapMatMult(W, N1)
	WN2 <- eigenMapMatMult(W, N2)
	
	# The 1st row of matrix C	
	C11 <- crossprod(X, WX)
	C12 <- crossprod(X, WN1)
	C13 <- crossprod(X, WN2)
	CRow1 <- cbind(C11, C12, C13)
	# The 2nd row of matrix C
	C21 <- crossprod(N1, WX)
	C22 <- crossprod(N1, WN1) + (1/theta[1]) * K  # lambda1 = 1/theta[1]
	C23 <- crossprod(N1, WN2)
	CRow2 <- cbind(C21, C22, C23)
	# The 2nd row of matrix C
	C31 <- crossprod(N2, WX)
	C32 <- crossprod(N2, WN1)
	C33 <- crossprod(N2, WN2) + (1/theta[2]) * K  # lambda2 = 1/theta[2]
	CRow3 <- cbind(C31, C32, C33)
	
	# Inverse of coefficient matrix
	C <- rbind(CRow1, CRow2, CRow3)
	invC <- ginv(C)
	
	
	#--------------------------------------------------------------------------
	# ESTIMATES FOR beta, f1 AND f2, TO BE USED IN SCORE IN FISHER-SCORING
	
	WY <- W %*% Y
	temp <- rbind(crossprod(X, WY), crossprod(N1, WY), crossprod(N2, WY))
	
	res <- invC %*% temp
	
	#--------------------------------------------------------------------------
	# RETURN RESULTS
	estimates <- list("betaHat" = res[1:2], "f1Hat" = matrix(res[3:30]), "f2Hat" = matrix(res[31:58]), "V" = V, "W" = W, "WN1" = WN1, "WN2" = WN2, "C22" = C22, "C33" = C33, "invC" = invC)
	return(estimates)
}



#===============================================================================
# PART III - FISHER SCORING
#===============================================================================

#------------------------------------------------------------------------------
# A FUNCTION OF ONE ITERATION OF FISHER-SCORING - INPUT AN OLD THETA, RETURN A NEW THETA

findPar = function(theta, Y, h, r){
	
	res <- est(theta, Y, h, r)
	betaHat <- res$betaHat
	f1Hat <- res$f1Hat
	f2Hat <- res$f2Hat
	W <- res$W
	invC <- res$invC
	
	resK <- matrixK(h, r)
	B1dou <- resK$B1dou
	B2dou <- resK$B2dou
	
	
	#--------------------------------------------------------------------------
	# P*
	chi <- cbind(X, N1, N2)	
	temp <- eigenMapMatMult(W, chi)
	Pst <- W - eigenMapMatMult3(temp, invC, t(temp))

	#--------------------------------------------------------------------------
	# Values needed in score, py = W * (Y - x*betaHat - N1 * f1Hat - N2 * f2Hat); tpy 
	py <- W %*% (Y - X %*% betaHat - N1 %*% f1Hat - N2 %*% f2Hat)
	tpy <- t(py)
	
	#--------------------------------------------------------------------------
	# PARTIAL DERIVATIVES, TO BE USED IN SCORE AND FISHER INFO


	# PARTIAL DERIVATIVES WRT b, 3:5 
	pV_phi1 <- do.call("blockMatrixDiagonal", replicate(n, Zi %*% matrix(c(1, 0, 0, 0), 2, 2) %*% t(Zi), simplify = F))
	pV_phi2 <- do.call("blockMatrixDiagonal", replicate(n, Zi %*% matrix(c(0, 1, 1, 0), 2, 2) %*% t(Zi), simplify = F))
	pV_phi3 <- do.call("blockMatrixDiagonal", replicate(n, Zi %*% matrix(c(0, 0, 0, 1), 2, 2) %*% t(Zi), simplify = F))
	
	
	# PARTIAL DERIVATIVES WRT epsilon, 6:7
	pV_sigma1 <- do.call("blockMatrixDiagonal", replicate(n, diag(rep(c(1, 0), tp)), simplify = F))
	pV_sigma2 <- do.call("blockMatrixDiagonal", replicate(n, diag(rep(c(0, 1), tp)), simplify = F))
	
	
	# PARTIAL DERIVATIVES WRT theta12, theta13, theta22, theta23, 8:11
	parV_theta12i <- matrix(0, 2*tp, 2*tp)
	parV_theta13i <- matrix(0, 2*tp, 2*tp)
	parV_theta22i <- matrix(0, 2*tp, 2*tp)
	parV_theta23i <- matrix(0, 2*tp, 2*tp)

	for (i in 1:tp) {
		for (j in i:tp) {
			parV_theta12i[2*i-1, 2*j-1] <- ouPv1(theta[8], theta[9], (i-j))			
			parV_theta12i[2*j-1, 2*i-1] <- t(parV_theta12i[2*i-1, 2*j-1])
			
			parV_theta13i[i,j] <- ouPv2(theta[8], theta[9], (i-j))
			parV_theta13i[j,i] <- t(parV_theta13i[i,j])
			
			parV_theta22i[2*i-1, 2*j-1] <- ouPv1(theta[10], theta[11], (i-j))
			parV_theta22i[2*j-1, 2*i-1] <- t(parV_theta22i[2*i-1, 2*j-1])
			
			parV_theta23i[i,j] <- ouPv1(theta[10], theta[11], (i-j))
			parV_theta23i[j,i] <- t(parV_theta23i[i,j])
		}
	}
	pV_theta12 <- do.call("blockMatrixDiagonal", replicate(n, parV_theta12i, simplify = F))
	pV_theta13 <- do.call("blockMatrixDiagonal", replicate(n, parV_theta13i, simplify = F))	
	pV_theta22 <- do.call("blockMatrixDiagonal", replicate(n, parV_theta22i, simplify = F))
	pV_theta23 <- do.call("blockMatrixDiagonal", replicate(n, parV_theta23i, simplify = F))
	

	#--------------------------------------------------------------------------	
	# SCORE
	score <- matrix(0, 11, 1) 
	
	# TEMP 1:11, NAIVE VERSION
	temp1 <- eigenMapMatMult(Pst, B1dou)
	temp2 <- eigenMapMatMult(Pst, B2dou)
	
	temp3 <- eigenMapMatMult(Pst, pV_phi1)
	temp4 <- eigenMapMatMult(Pst, pV_phi2)
	temp5 <- eigenMapMatMult(Pst, pV_phi3)
	
	temp6 <- eigenMapMatMult(Pst, pV_sigma1)
	temp7 <- eigenMapMatMult(Pst, pV_sigma2)
	
	temp8 <- eigenMapMatMult(Pst, pV_theta12)
	temp9 <- eigenMapMatMult(Pst, pV_theta13)
	temp10 <- eigenMapMatMult(Pst, pV_theta22)
	temp11 <- eigenMapMatMult(Pst, pV_theta23)
	
	
	# SCORE
	score[1] <- eigenMapMatMult3(tpy, B1dou, py) - tr(temp1)
	score[2] <- eigenMapMatMult3(tpy, B2dou, py) - tr(temp2)
	
	score[3] <- eigenMapMatMult3(tpy, pV_phi1, py) - tr(temp3)
	score[4] <- eigenMapMatMult3(tpy, pV_phi2, py) - tr(temp4)
	score[5] <- eigenMapMatMult3(tpy, pV_phi3, py) - tr(temp5)
	
	score[6] <- eigenMapMatMult3(tpy, pV_sigma1, py) - tr(temp6)
	score[7] <- eigenMapMatMult3(tpy, pV_sigma2, py) - tr(temp7)
	
	score[8] <- eigenMapMatMult3(tpy, pV_theta12, py) - tr(temp8)
	score[9] <- eigenMapMatMult3(tpy, pV_theta13, py) - tr(temp9)
	score[10] <- eigenMapMatMult3(tpy, pV_theta22, py) - tr(temp10)
	score[11] <- eigenMapMatMult3(tpy, pV_theta23, py) - tr(temp11)

	
	#--------------------------------------------------------------------------
	# FISHER INFO
	fInfo = matrix(0, 11, 11)
	
	# UPPER LEFT 4 ELEMENTS
	fInfo[1, 1] <- tr(eigenMapMatMult(temp1, temp1))	
	fInfo[1, 2] <- tr(eigenMapMatMult(temp1, temp2))	
	fInfo[2, 1] <- t(fInfo[1, 2])	
	fInfo[2, 2] <- tr(eigenMapMatMult(temp2, temp2))	
	
	# I_tau1Theta 
	fInfo[1, 3] <- tr(eigenMapMatMult(temp1, temp3))
	fInfo[1, 4] <- tr(eigenMapMatMult(temp1, temp4))
	fInfo[1, 5] <- tr(eigenMapMatMult(temp1, temp5))
	fInfo[1, 6] <- tr(eigenMapMatMult(temp1, temp6))
	fInfo[1, 7] <- tr(eigenMapMatMult(temp1, temp7))
	fInfo[1, 8] <- tr(eigenMapMatMult(temp1, temp8))
	fInfo[1, 9] <- tr(eigenMapMatMult(temp1, temp9))
	fInfo[1, 10] <- tr(eigenMapMatMult(temp1, temp10))
	fInfo[1, 11] <- tr(eigenMapMatMult(temp1, temp11))
	fInfo[3:11, 1] <- t(fInfo[1, 3:11]) # SYMMETRIC PART
	
	# I_tau2Theta
	fInfo[2, 3] <- tr(eigenMapMatMult(temp2, temp3))
	fInfo[2, 4] <- tr(eigenMapMatMult(temp2, temp4))
	fInfo[2, 5] <- tr(eigenMapMatMult(temp2, temp5))
	fInfo[2, 6] <- tr(eigenMapMatMult(temp2, temp6))
	fInfo[2, 7] <- tr(eigenMapMatMult(temp2, temp7))
	fInfo[2, 8] <- tr(eigenMapMatMult(temp2, temp8))
	fInfo[2, 9] <- tr(eigenMapMatMult(temp2, temp9))
	fInfo[2, 10] <- tr(eigenMapMatMult(temp2, temp10))
	fInfo[2, 11] <- tr(eigenMapMatMult(temp2, temp11))
	fInfo[3:11, 2] <- t(fInfo[2, 3:11]) # SYMMETRIC PART
	
	fInfo[3, 3] <- tr(eigenMapMatMult(temp3, temp3))
	fInfo[3, 4] <- tr(eigenMapMatMult(temp3, temp4))
	fInfo[3, 5] <- tr(eigenMapMatMult(temp3, temp5))
	fInfo[3, 6] <- tr(eigenMapMatMult(temp3, temp6))
	fInfo[3, 7] <- tr(eigenMapMatMult(temp3, temp7))
	fInfo[3, 8] <- tr(eigenMapMatMult(temp3, temp8))
	fInfo[3, 9] <- tr(eigenMapMatMult(temp3, temp9))
	fInfo[3, 10] <- tr(eigenMapMatMult(temp3, temp10))
	fInfo[3, 11] <- tr(eigenMapMatMult(temp3, temp11))
	
	fInfo[4, 3] <- t(fInfo[3, 4])
	fInfo[4, 4] <- tr(eigenMapMatMult(temp4, temp4))
	fInfo[4, 5] <- tr(eigenMapMatMult(temp4, temp5))
	fInfo[4, 6] <- tr(eigenMapMatMult(temp4, temp6))
	fInfo[4, 7] <- tr(eigenMapMatMult(temp4, temp7))
	fInfo[4, 8] <- tr(eigenMapMatMult(temp4, temp8))
	fInfo[4, 9] <- tr(eigenMapMatMult(temp4, temp9))
	fInfo[4, 10] <- tr(eigenMapMatMult(temp4, temp10))
	fInfo[4, 11] <- tr(eigenMapMatMult(temp4, temp11))
	
	fInfo[5, 3:4] <- t(fInfo[3:4, 5])
	fInfo[5, 5] <- tr(eigenMapMatMult(temp5, temp5))
	fInfo[5, 6] <- tr(eigenMapMatMult(temp5, temp6))
	fInfo[5, 7] <- tr(eigenMapMatMult(temp5, temp7))
	fInfo[5, 8] <- tr(eigenMapMatMult(temp5, temp8))
	fInfo[5, 9] <- tr(eigenMapMatMult(temp5, temp9))
	fInfo[5, 10] <- tr(eigenMapMatMult(temp5, temp10))
	fInfo[5, 11] <- tr(eigenMapMatMult(temp5, temp11))
	
	fInfo[6, 3:5] <- t(fInfo[3:5, 6])
	fInfo[6, 6] <- tr(eigenMapMatMult(temp6, temp6))
	fInfo[6, 7] <- tr(eigenMapMatMult(temp6, temp7))
	fInfo[6, 8] <- tr(eigenMapMatMult(temp6, temp8))
	fInfo[6, 9] <- tr(eigenMapMatMult(temp6, temp9))
	fInfo[6, 10] <- tr(eigenMapMatMult(temp6, temp10))
	fInfo[6, 11] <- tr(eigenMapMatMult(temp6, temp11))

	fInfo[7, 3:6] <- t(fInfo[3:6, 7])
	fInfo[7, 7] <- tr(eigenMapMatMult(temp7, temp7))
	fInfo[7, 8] <- tr(eigenMapMatMult(temp7, temp8))
	fInfo[7, 9] <- tr(eigenMapMatMult(temp7, temp9))
	fInfo[7, 10] <- tr(eigenMapMatMult(temp7, temp10))
	fInfo[7, 11] <- tr(eigenMapMatMult(temp7, temp11))
	
	fInfo[8, 3:7] <- t(fInfo[3:7, 8])
	fInfo[8, 8] <- tr(eigenMapMatMult(temp8, temp8))
	fInfo[8, 9] <- tr(eigenMapMatMult(temp8, temp9))
	fInfo[8, 10] <- tr(eigenMapMatMult(temp8, temp10))
	fInfo[8, 11] <- tr(eigenMapMatMult(temp8, temp11))
	
	fInfo[9, 3:8] <- t(fInfo[3:8, 9])
	fInfo[9, 9] <- tr(eigenMapMatMult(temp9, temp9))
	fInfo[9, 10] <- tr(eigenMapMatMult(temp9, temp10))
	fInfo[9, 11] <- tr(eigenMapMatMult(temp9, temp11))
	
	fInfo[10, 3:9] <- t(fInfo[3:9, 10])
	fInfo[10, 10] <- tr(eigenMapMatMult(temp10, temp10))
	fInfo[10, 11] <- tr(eigenMapMatMult(temp10, temp11))

	fInfo[11, 3:10] <- t(fInfo[3:10, 11])
	fInfo[11, 11] <- tr(eigenMapMatMult(temp11, temp11))
	
	#--------------------------------------------------------------------------
	score <- score/2
	fInfo <- fInfo/2
	Finv <- ginv(fInfo)

	thetaNew <- theta + Finv %*% score	# Parameter estimate after one iteration
		
	return(thetaNew)
}

 
#------------------------------------------------------------------------------
# FISHER-SCORING ALGORITHM

iternum <- 0
# initialization of difference of parameters
diff <- matrix(1, 1, 11) # a vector of 1 of length 11


# tol is a vector of tolerance levels for different parameters.
# cap is the maximum runs of iteration, after which the iteration will stop.
fisher.scoring <- function(par0, Y, tol, cap) {
	
	theta1 <- par0
	
	while (diff[1] > tol[1] | diff[2] > tol[2] | diff[3] > tol[3] | diff[4] > tol[4] | diff[5] > tol[5] |
			diff[6] > tol[6] | diff[7] > tol[7] | diff[8] > tol[8] | diff[9] > tol[9] 
			| diff[10] > tol[10] | diff[11] > tol[11]) {
		
		iternum <- iternum + 1
		
	    theta0 <- theta1
		theta1 <- findPar(theta0, Y)
		
		for (i in 1:11){
			while (theta1[i]<0) {
				theta1[i] <- (theta1[i] + theta0[i])/2
			}
		}
		
		for (i in 1:11) {
			diff[i] <- abs(theta1[i] - theta0[i])
		}
		
		# Abort this iteration if iternum > cap, and set the parameters to be zero
		if (iternum > cap) {
			theta1 <- matrix(0, 11)
			return (theta1)
		}
	
		print(iternum)
		print (diff)
		print(theta1)
	}
	
		
	estNew <- est(theta1, Y) # PARAMETER ESTIMATES BASED ON VARIANCE theta ESTIMATE

	newlist <- list("theta" = theta1, "betaHat" = estNew$betaHat, "f1Hat" = estNew$f1Hat, "f2Hat" = estNew$f2Hat)
	return(newlist)
}


#======================================================================================
# PART IV - COVARIANCES for beta, f1 & f2
#======================================================================================

# INPUT A VECTOR OF PARAMETERS theta, RETURN COVARIANCES for BETA & f.
# theta COMES FROM fisher.scoring()

# system.time(theta <- fisher.scoring(theta0, Y, rep(0.1, 5), cap = 10))

biasCov <- function(theta, h, r){
	
	res <- est(theta, Y) 
	V <- res$V
	W <- res$W
	WN1 <- res$WN1	
	WN2 <- res$WN2
	
	tX <- t(X)
	tN1 <- t(N1)
	tN2 <- t(N2)
	
	K <- matrixK(h, r)$K
	K1 <- (1/theta[1]) * K
	K2 <- (1/theta[2]) * K

	# WEIGHT MATRICES, W1 & W2	

	W1 <- W - eigenMapMatMult3(WN1, ginv(res$C22), t(WN1))
	W2 <- W - eigenMapMatMult3(WN2, ginv(res$C33), t(WN2))

	# WEIGHT MATRICES, Wx, Wf1 & Wf2 	
	
	W1N2 <- eigenMapMatMult(W1, N2)
	W2X <- eigenMapMatMult(W2, X)
	W1X <- eigenMapMatMult(W1, X)
	
	WxInv <- ginv(eigenMapMatMult3(tN2, W1, N2) + K2)
	Wf1Inv <- ginv(eigenMapMatMult3(tX, W2, X))
	Wf2Inv <- ginv(eigenMapMatMult3(tX, W1, X))
	
	Wx <- W1 - eigenMapMatMult3(W1N2, WxInv, t(W1N2))
	Wf1 <- W2 - eigenMapMatMult3(W2X, Wf1Inv, t(W2X))
	Wf2 <- W1 - eigenMapMatMult3(W1X, Wf2Inv, t(W1X))
	
	# COVARIANCES - beta
	
	betaInv <- ginv(eigenMapMatMult3(tX, Wx, X))
	hatB <- eigenMapMatMult3(betaInv, tX, Wx)
	
	N1f1 <- eigenMapMatMult(N1, f1)
	N2f2 <- eigenMapMatMult(N2, f2)
	Nf <- N1f1 + N2f2
	
	biasBeta <- eigenMapMatMult(hatB, Nf)
	covBeta <- eigenMapMatMult3(hatB, V, t(hatB))
	
	
	# COVARIANCES - f1

	tN1Wf1 <- eigenMapMatMult(tN1, Wf1)
	f1Inv <- ginv(eigenMapMatMult(tN1Wf1, N1) + K1) 
	hatF1 <- eigenMapMatMult(f1Inv, tN1Wf1)
	
	temp1 <- eigenMapMatMult(tN1Wf1, N2f2) - eigenMapMatMult(K1, f1)
	
	biasF1 <- eigenMapMatMult(f1Inv, temp1)
	covF1 <- eigenMapMatMult3(hatF1, V, t(hatF1))
	
	
	# COVARIANCES - f2
 
	tN2Wf2 <- eigenMapMatMult(tN2, Wf2)
	f2Inv <- ginv(eigenMapMatMult(tN2Wf2, N2) + K2)
	hatF2 <- eigenMapMatMult(f2Inv, tN2Wf2)
	
	temp2 <- eigenMapMatMult(tN2Wf2, N1f1) - eigenMapMatMult(K2, f2)
	
	biasF2 <- eigenMapMatMult(f2Inv, temp2)
	covF2 <- eigenMapMatMult3(hatF2, V, t(hatF2))
	
		
	newlist <- list("biasBeta" = biasBeta, "covBeta" = covBeta, "biasF1" = biasF1, "biasF2" = biasF2, "covF1" = covF1, "covF2" = covF2)
	return(newlist)
}







