# THIS R FILE IS BASED ON bivSim_fns.r, USED FOR BOTH SIMULATION AND REAL DATASET ANALYSIS.

library(sde) # rsOU()
library(MASS) # mvrnorm()
library(psych) # tr()
library(Rcpp) # TO C++
sourceCpp("matrixmulti_1.cpp") # NEED TO INSTALL Package 'RcppEigen'

#==============================================================================
# PART I - FUNCTIONS USED IN SIMULATION
#==============================================================================

#------------------------------------------------------------------------------
# A FUNCTION that RETURNS X, Y, N1 & N2 in A BIVARIATE LONGITUDINAL DATASET in a SIMULATION
# m is the number of subjects
# tp is the number of distinct points

newData <- function(theta, beta11, beta21, m, tp) {
	
	tp2 <- 2*tp
	n <- m*tp
	
	#------------------------------------------------------------------------------
	# X
	
	X <- matrix(0, tp2*m, 2)
	for (i in 1:m) {
		# set.seed(34)
		age <- sample(20:44, size = m, replace = T)
		
		index <- (tp2*(i-1)+1):(tp2*i)
		X[index, 1] <- rep(c(age[i], 0), tp) # the 1st column of Xi for subj. i is age[i], 0, age[i], 0...
		X[index, 2] <- rep(c(0, age[i]), tp) # the 2nd column of Xi for subj. i is 0, age[i], 0, age[i]...
	}	
	#------------------------------------------------------------------------------
	# N1 & N2 
	
	Ni <- diag(tp)
	A1i <- matrix(0, tp2, tp) # initialize matrix A1i
	A2i <- matrix(0, tp2, tp) # initialize matrix A2i	
	for (j in 1:tp){		
		A1i[2*j-1, j] <- 1
		A2i[2*j, j] <- 1 
	}
	N1  <- do.call("rbind", replicate(m, A1i %*% Ni, simplify = F))
	N2  <- do.call("rbind", replicate(m, A2i %*% Ni, simplify = F))

	#------------------------------------------------------------------------------
	# BIVARIATE CYCLIC RESPONSE Y 
		
	t <- seq(1, 28)
	f1 <- sine(tp, t)
	f2 <- cosine(tp, t)
	
	# random intercepts
	D <- matrix(c(theta[3], theta[4], theta[4], theta[5]), 2, 2) # variance for bi
	b <- mvrnorm(m, c(0, 0), D) 
	
	# measurement error - epsilon
	var_eps <- c(theta[6], theta[7])
	sigma <- diag(var_eps) # variance for epsilon
	
	eps <- mvrnorm(n, c(0, 0), sigma)
	eps1 <- matrix(eps[,1], m, tp)
	eps2 <- matrix(eps[,2], m, tp)
	
	# bivariate Gaussian field
	tempU_1 <- rsOU(n, theta=c(0, theta[8], theta[9])) 
	tempU_2 <- rsOU(n, theta=c(0, theta[10], theta[11]))
	u1 <- matrix(tempU_1, m, tp)
	u2 <- matrix(tempU_2, m, tp)
	
	# bivariate Longitudinal data Y
	tempResponse <- matrix(0, 2*m, tp)
	tempR <- matrix(0, tp, 2) # for each subject
	bivY <- matrix(0, n, 2)
	Y <- matrix(0, 2*n, 1)
	for (i in 1:m) {
		# every two rows are one subject
		tempResponse[(2*i-1),] <- beta11*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
		tempResponse[2*i,] <- beta21*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
		tempR <- t(tempResponse[(2*i-1):(2*i),])
		bivY[((i-1)*tp+1):(i*tp), ] <- tempR
	}
	for (i in 1:n) {
		Y[2*i - 1] <- bivY[i,1]
		Y[2*i] <- bivY[i,2]
	} 
	return(list("X" = X, "Y" = Y, "N1" = N1, "N2" = N2))	
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS A PERIODIC FUNCTION

sine <- function(tp, x) {
	f1 <- 5*sin((2*pi/tp)*x)
	return(f1)
}

cosine <- function(tp, x) {
	f2 <-  3*cos((2*pi/tp)*x)
	return(f2)
}
#------------------------------------------------------------------------------
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
#===============================================================================
# PART II - LIU data
#===============================================================================
# INPUT ONLY FIRST 34 SUBJ. TO REDUCE DIMENSION, r is still 56.
# want <- read.table("ageBMI_want.txt", header = T)
# data <- read.table("ageBMI.txt", header = T)

# m = 34
# want <- read.table("ageBMI_want.txt", header = T)[1:34,]
# data <- read.table("ageBMI.txt", header = T)[1:992,]

# m = 50
# want <- read.table("ageBMI_want.txt", header = T)[1:50,]
# data <- read.table("ageBMI.txt", header = T)[1:1463,]

# m = 100
want <- read.table("ageBMI_want.txt", header = T)[1:100,]
data <- read.table("ageBMI.txt", header = T)[1:2810,]
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS BIVARIATE RESPONSE Y, X, N1 & N2 from LIU DATASET
liuData <- function(m, niVec){
	
	n <- sum(niVec)
	s <- 2*n
	
	#------------------------------------------------------------------------------
	# Y
	Y <- matrix(0, s)
	for (i in 1:n) {
		i2 <- 2*i
		Y[i2-1] = log(data$adjpdg2[i])
		Y[i2] = log(data$adje1c2[i])
	}
	#------------------------------------------------------------------------------
	# X, N1 & N2 
	age <- want$age
	BMI <- want$BMI
	day <- data$standday
	uniqueDay <- sort(unique(day))
	r <- length(uniqueDay)
	
	# Initialization for X, N1 & N2
	X <- matrix(0, s, 4) # initialization
	N1 <- matrix(0, s, r)
	N2 <- matrix(0, s, r)
	for (i in 1:m) {
		
		ni = max[i]
		index <- ((2*sum(max[0:(i-1)])) + 1):(2*sum(max[1:i])) 
		
		# X
		X[index, 1] = rep(c(age[i], 0), ni) # the 1st column of Xi is age[1], 0, age[1], 0...
		X[index, 2] = rep(c(BMI[i], 0), ni) # the 2nd column of Xi is BMI[1], 0, BMI[1], 0...
		X[index, 3] = rep(c(0, age[i]), ni) # the 3rd column of Xi is 0, age[1], 0, age[1]...
		X[index, 4] = rep(c(0, BMI[i]), ni) # the 4th column of Xi is 0, BMI[1], 0, BMI[1]...
		
		# N1 & N2 
		dayIndex <- (sum(max[0:(i-1)]) + 1):sum(max[1:i]) 
		dayi <- day[dayIndex]
		
		Ni <- matrix(0, ni, r)
		for (j in 1:ni){
			for (l in 1:r){
				if(dayi[j] == uniqueDay[l]) {
					Ni[j, l] <- 1
				}
			}
		}
		A1i <- matrix(0, 2*ni, ni) # initialize matrix A1i
		A2i <- matrix(0, 2*ni, ni) # initialize matrix A2i
		for (j in 1:ni){
			A1i[2*j-1, j] <- 1
			A2i[2*j, j] <- 1 
		}	
		N1[index,] <- A1i %*% Ni
		N2[index,] <- A2i %*% Ni	
	}
	return(list("X" = X, "Y" = Y, "N1" = N1, "N2" = N2))
}
#===============================================================================
# PART III - COMMON FUNCTIONS
#===============================================================================

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS Zi

matrixZi <- function(ni) {
	
	Zi <- matrix(0, 2*ni, 2) # initialize matrix Zi
	Zi[ ,1] <- rep(c(1, 0), ni) # the first column of Zi is 1, 0, 1, 0...
	Zi[ ,2] <- rep(c(0, 1), ni) # the second column of Zi is 0, 1, 0, 1...
	return(Zi)
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS Vi 
# ni is the number of observations for subject i

matrixVi <- function(theta, ni) {
	
	Zi <- matrixZi(ni) 

	# VARIANCE for bi
	D <- matrix(c(theta[3], theta[4], theta[4], theta[5]), 2, 2) 
	
	# VARIANCE OF MEASUREMENT ERROR
	sigma_i <- diag(rep(c(theta[6], theta[7]), ni))
	
	# VARIANCE FOR OU PROCESS
	# USE VARIANCE FUNCTION, ASSIGN ONLY THE UPPER TRIANGLE, TAKE ACCOUNT OF SYMMETRY OF GAMMAi
	# Gamma_i
	Gammai <- matrix(0, 2*ni, 2*ni)
	for (i in 1:ni) {
		for (j in i:ni) {
			i2 <- 2*i
			j2 <- 2*j
			ij <- i - j
			Gammai[i2-1, j2-1] <- ouVar(theta[8], theta[9], ij)
			Gammai[i2, j2] <- ouVar(theta[10], theta[11], ij)
			Gammai[j2-1, i2-1] <- t(Gammai[i2-1, j2-1])
			Gammai[j2, i2] <- t(Gammai[i2, j2])
		}
	}		
	Vi <- Zi %*% D %*% t(Zi) + sigma_i + Gammai	
	return("Vi" = Vi)
}
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS W & V in SIMULATION and LIU dataset
# niVec is the vector of ni for all subjects 

matrixW <- function(theta, m, tp, niVec) {
	
	# SIMULATION 
	if (identical(niVec, rep(tp, m))){
		Vi <- matrixVi(theta, tp)
		Wi <- solve(Vi)
		V <- do.call("blockMatrixDiagonal", replicate(m, Vi, simplify = F))
		W <- do.call("blockMatrixDiagonal", replicate(m, Wi, simplify = F))
		return(list("V" = V, "W" = W))
	} else { # LIU DATASET
		s <- 2*sum(niVec)
		W <- diag(0, s) # initialization of W
		V <- diag(0, s) # initialization of V
		for (i in 1:m) {
			ni <- niVec[i]
			Vi <- matrixVi(theta, ni)
			# index for W, nrow = ncol
			index <- ((2*sum(max[0:(i-1)])) + 1):(2*sum(max[1:i])) 
			V[index, index] <- Vi
			W[index, index] <- solve(Vi)
		}
		return(list("V" = V, "W" = W))
	}
}
#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS PVi 
PVi <- function(theta, ni) {
	
	Zi <- matrixZi(ni)
	tZi <- t(Zi)

	pV_phi1i <- Zi %*% matrix(c(1, 0, 0, 0), 2, 2) %*% tZi
	pV_phi2i <- Zi %*% matrix(c(0, 1, 1, 0), 2, 2) %*% tZi
	pV_phi3i <- Zi %*% matrix(c(0, 0, 0, 1), 2, 2) %*% tZi

	pV_sigma1i <- diag(rep(c(1, 0), ni))
	pV_sigma2i <- diag(rep(c(0, 1), ni))
	
	# PARTIAL DERIVATIVES WRT theta12, theta13, theta22, theta23, 8:11
	ni2 <- 2*ni
	parV_theta12i <- matrix(0, ni2, ni2)
	parV_theta13i <- matrix(0, ni2, ni2)
	parV_theta22i <- matrix(0, ni2, ni2)
	parV_theta23i <- matrix(0, ni2, ni2)
	for (i in 1:ni) {
		for (j in i:ni) {
			
			i2 <- 2*i - 1
			j2 <- 2*j - 1
			ij <- i - j
			
			parV_theta12i[i2, j2] <- ouPv1(theta[8], theta[9], ij)			
			parV_theta12i[j2, i2] <- t(parV_theta12i[i2, j2])
			
			parV_theta13i[i,j] <- ouPv2(theta[8], theta[9], ij)
			parV_theta13i[j,i] <- t(parV_theta13i[i,j])
			
			parV_theta22i[i2, j2] <- ouPv1(theta[10], theta[11], ij)
			parV_theta22i[j2, i2] <- t(parV_theta22i[i2, j2])
			
			parV_theta23i[i,j] <- ouPv1(theta[10], theta[11], ij)
			parV_theta23i[j,i] <- t(parV_theta23i[i,j])
		}
	}
	pvlist <- list("pV_phi1i" = pV_phi1i, "pV_phi2i" = pV_phi2i, "pV_phi3i" = pV_phi3i, "pV_sigma1i" = pV_sigma1i, "pV_sigma2i" = pV_sigma2i, "pV_theta12i" = parV_theta12i, "pV_theta13i" = parV_theta13i, "pV_theta22i" = parV_theta22i, "pV_theta23i" = parV_theta23i)
	return(pvlist)
}

#------------------------------------------------------------------------------
# A FUNCTION that RETURNS PARTIAL DERIVATIVES wrt THETA

PVar <- function(theta, m, tp, niVec) {
	
	# SIMULATION 
	if (identical(niVec, rep(tp, m))){
		
		parVi <- PVi(theta, tp)

		# PARTIAL DERIVATIVES WRT b, 3:5 
		pV_phi1 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_phi1i, simplify = F))
		pV_phi2 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_phi2i, simplify = F))
		pV_phi3 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_phi3i, simplify = F))
		
		# PARTIAL DERIVATIVES WRT epsilon, 6:7
		pV_sigma1 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_sigma1i, simplify = F))
		pV_sigma2 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_sigma2i, simplify = F))
	
		# PARTIAL DERIVATIVES WRT theta12, theta13, theta22, theta23, 8:11
		pV_theta12 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_theta12i, simplify = F))
		pV_theta13 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_theta13i, simplify = F))	
		pV_theta22 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_theta22i, simplify = F))
		pV_theta23 <- do.call("blockMatrixDiagonal", replicate(m, parVi$pV_theta23i, simplify = F))
		
		pvlist <- list("pV_phi1" = pV_phi1, "pV_phi2" = pV_phi2, "pV_phi3" = pV_phi3, "pV_sigma1" = pV_sigma1, "pV_sigma2" = pV_sigma2, "pV_theta12" = pV_theta12, "pV_theta13" = pV_theta13, "pV_theta22" = pV_theta22, "pV_theta23" = pV_theta23)
		return(pvlist)
		
	} else { # LIU DATASET		
		s <- 2*sum(niVec)
		for (i in 1:m) {					
			pV_phi1 <- matrix(0, s, s)
			pV_phi2 <- matrix(0, s, s)
			pV_phi3 <- matrix(0, s, s)
			pV_sigma1 <- matrix(0, s, s)
			pV_sigma2 <- matrix(0, s, s)
			pV_theta12 <- matrix(0, s, s)
			pV_theta13 <- matrix(0, s, s)
			pV_theta22 <- matrix(0, s, s)
			pV_theta23 <- matrix(0, s, s)
			
			ni <- niVec[i]
			parVi <- PVi(theta, ni)
			index <- ((2*sum(max[0:(i-1)])) + 1):(2*sum(max[1:i]))
			
			pV_phi1[index, index] <- parVi$pV_phi1i
			pV_phi2[index, index] <- parVi$pV_phi2i
			pV_phi3[index, index] <- parVi$pV_phi3i
			pV_sigma1[index, index] <- parVi$pV_sigma1i
			pV_sigma2[index, index] <- parVi$pV_sigma2i
			pV_theta12[index, index] <- parVi$pV_theta12i
			pV_theta13[index, index] <- parVi$pV_theta13i
			pV_theta22[index, index] <- parVi$pV_theta22i
			pV_theta23[index, index] <- parVi$pV_theta23i			
		}
		pvlist <- list("pV_phi1" = pV_phi1, "pV_phi2" = pV_phi2, "pV_phi3" = pV_phi3, "pV_sigma1" = pV_sigma1, "pV_sigma2" = pV_sigma2, "pV_theta12" = pV_theta12, "pV_theta13" = pV_theta13, "pV_theta22" = pV_theta22, "pV_theta23" = pV_theta23)
		return(pvlist)
	}
}
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS K
# tp is the number of dictinct time points in one period; 
# h is the difference between two adjacent time points, i.e. h = t_(i+1) - t_i.

matrixK <- function(h, tp) {
	
	# Define matrix Q, double checked.
	Q <- matrix(0, tp, tp-2)
	for(j in 1:(tp-2)) {
		Q[j, j] <- 1/h
		Q[j + 1, j] <- -2/h
		Q[j + 2, j] <- 1/h
	}	
	# Define matrix R 
	R = matrix(0, tp-2, tp-2)
	for (i in 1:(tp-2)) {
		R[i, i] <- (2/3)*h
	}
	for (i in 1:(tp-3)) {
		R[i, i + 1] <- (1/6)*h
		R[i + 1, i] <- (1/6)*h
	}
	K <- eigenMapMatMult3(Q, ginv(R), t(Q))
	return(K)
}
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS VARIANCE FUNCTION for OU PROCESS
ouVar <- function(theta2, theta3, x) {
	ouV <- (theta3^2/(2*theta2)) * exp(-theta2*abs(x))
	return(ouV)
}

# A FUNCTION that RETURNS PARTIAL DERIVATIVES OF VARIANCE FUNCTION WRT THETA2, OU PROCESS
ouPv1 <- function(theta2, theta3, x) {
	ouP1 <- (-1/theta2 - abs(x)) * ouVar(theta2, theta3, x)
	return(ouP1)
}


# A FUNCTION that RETURNS PARTIAL DERIVATIVES OF VARIANCE FUNCTION WRT THETA3, OU PROCESS
ouPv2 <- function(theta2, theta3, x) {
	ouP2 <- (theta3/theta2) * exp(-theta2*abs(x))
	return(ouP2)
}
#======================================================================================
# PART IV - MAIN FUNCTIONS; FISHER-SCORING & RESULTS
#======================================================================================
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS ESTIMATES FOR BETA and F.

est = function(theta, Y, m, tp, K, X, N1, N2, niVec) {
	
	matW <- matrixW(theta, m, tp, niVec)
	W <- matW$W
	V <- matW$V
	
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
	
	WY <- W %*% Y
	temp <- rbind(crossprod(X, WY), crossprod(N1, WY), crossprod(N2, WY))
	
	# ESTIMATES FOR beta, f1 AND f2, TO BE USED IN SCORE IN FISHER-SCORING
	res <- invC %*% temp
	
	# RETURN RESULTS
	if (identical(niVec, rep(tp, m))){
		# SIMULATION
		estimates <- list("betaHat" = res[1:2], "f1Hat" = matrix(res[3:30]), "f2Hat" = matrix(res[31:58]), "WN1" = WN1, "WN2" = WN2, "C22" = C22, "C33" = C33, "invC" = invC, "V" = V, "W" = W)
		return(estimates)
	}
	# LIU 
	estimates <- list("betaHat" = res[1:4], "f1Hat" = matrix(res[5:60]), "f2Hat" = matrix(res[61:116]), "WN1" = WN1, "WN2" = WN2, "C22" = C22, "C33" = C33, "invC" = invC, "V" = V, "W" = W)
	return(estimates)
}
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS A NEW THETA, after ONE ITERATION OF FISHER-SCORING

findPar = function(theta, Y, m, tp, K, X, N1, N2, niVec){
	
	res <- est(theta, Y, m, tp, K, X, N1, N2, niVec)
	betaHat <- res$betaHat
	f1Hat <- res$f1Hat
	f2Hat <- res$f2Hat
	invC <- res$invC
	W <- res$W
		
	#--------------------------------------------------------------------------
	# P*
	chi <- cbind(X, N1, N2)	
	temp <- eigenMapMatMult(W, chi)
	Pst <- W - eigenMapMatMult3(temp, invC, t(temp))

	#--------------------------------------------------------------------------
	# Values needed in score, py = W * (Y - x*betaHat - N1 * f1Hat - N2 * f2Hat); tpy 
	py <- W %*% (Y - X %*% betaHat - N1 %*% f1Hat - N2 %*% f2Hat)
	tpy <- t(py)
		
	#------------------------------------------------------------------------------
	# MATRICES B1* B2*
	
	# Matrix B
	LTranspose <- chol(K, pivot = T, LDL = T) # upper triangular matrix
	L <- t(LTranspose) # lower triangular matrix
	B <- L %*% ginv(crossprod(L))
		
	# TO BE USED IN SCORE
	B1dou <- tcrossprod(N1 %*% B)
	B2dou <- tcrossprod(N2 %*% B)
		
	#--------------------------------------------------------------------------	
	# SCORE
	
	pv <- PVar(theta, m, tp, niVec)

	# TEMP 1:11, NAIVE VERSION
	temp1 <- eigenMapMatMult(Pst, B1dou)
	temp2 <- eigenMapMatMult(Pst, B2dou)
	temp3 <- eigenMapMatMult(Pst, pv$pV_phi1)
	temp4 <- eigenMapMatMult(Pst, pv$pV_phi2)
	temp5 <- eigenMapMatMult(Pst, pv$pV_phi3)
	temp6 <- eigenMapMatMult(Pst, pv$pV_sigma1)
	temp7 <- eigenMapMatMult(Pst, pv$pV_sigma2)
	temp8 <- eigenMapMatMult(Pst, pv$pV_theta12)
	temp9 <- eigenMapMatMult(Pst, pv$pV_theta13)
	temp10 <- eigenMapMatMult(Pst, pv$pV_theta22)
	temp11 <- eigenMapMatMult(Pst, pv$pV_theta23)
	
	# SCORE
	
	score <- matrix(0, 11, 1) 
	score[1] <- eigenMapMatMult3(tpy, B1dou, py) - tr(temp1)
	score[2] <- eigenMapMatMult3(tpy, B2dou, py) - tr(temp2)
	score[3] <- eigenMapMatMult3(tpy, pv$pV_phi1, py) - tr(temp3)
	score[4] <- eigenMapMatMult3(tpy, pv$pV_phi2, py) - tr(temp4)
	score[5] <- eigenMapMatMult3(tpy, pv$pV_phi3, py) - tr(temp5)
	score[6] <- eigenMapMatMult3(tpy, pv$pV_sigma1, py) - tr(temp6)
	score[7] <- eigenMapMatMult3(tpy, pv$pV_sigma2, py) - tr(temp7)
	score[8] <- eigenMapMatMult3(tpy, pv$pV_theta12, py) - tr(temp8)
	score[9] <- eigenMapMatMult3(tpy, pv$pV_theta13, py) - tr(temp9)
	score[10] <- eigenMapMatMult3(tpy, pv$pV_theta22, py) - tr(temp10)
	score[11] <- eigenMapMatMult3(tpy, pv$pV_theta23, py) - tr(temp11)

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
	
	# I_thetaTheta
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
# A FUNCTION that RETURNS THETA, using FISHER-SCORING ALGORITHM
# tol is a vector of tolerance levels for different parameters.
# cap is the maximum runs of iteration, after which the iteration will stop.

fisher.scoring <- function(par0, Y, tol, cap, m, tp, K, X, N1, N2, niVec) {
	
	theta1 <- par0
	iternum <- 0
	diff <- matrix(1, 1, 11) # a vector of 1 of length 11
	
	while (diff[1] > tol[1] | diff[2] > tol[2] | diff[3] > tol[3] | diff[4] > tol[4] | diff[5] > tol[5] |
			diff[6] > tol[6] | diff[7] > tol[7] | diff[8] > tol[8] | diff[9] > tol[9] 
			| diff[10] > tol[10] | diff[11] > tol[11]) {
		
		iternum <- iternum + 1
		
	    theta0 <- theta1
		theta1 <- findPar(theta0, Y, m, tp, K, X, N1, N2, niVec)
		
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
	return(theta1)
}
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS BIAS & COVARIANCES for BETA & f.

# theta COMES FROM fisher.scoring()
# system.time(theta <- fisher.scoring(theta0, Y, rep(0.1, 5), cap = 10))

results <- function(theta, Y, tol, cap, m, tp, K, X, N1, N2, niVec){
	
	newTheta <- fisher.scoring(theta, Y, tol, cap, m, tp, K, X, N1, N2, niVec)
	
	if (identical(newTheta, matrix(0, 11))) {
		return (NA)
	} else {
		res <- est(newTheta, Y, m, tp, K, X, N1, N2, niVec)
		WN1 <- res$WN1	
		WN2 <- res$WN2
		V <- res$V
		W <- res$W
		
		tX <- t(X)
		tN1 <- t(N1)
		tN2 <- t(N2)
		
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
		
		# COVARIANCE - beta	
		betaInv <- ginv(eigenMapMatMult3(tX, Wx, X))
		hatB <- eigenMapMatMult3(betaInv, tX, Wx)
		covBeta <- eigenMapMatMult3(hatB, V, t(hatB))
		
		# COVARIANCE - f1
		tN1Wf1 <- eigenMapMatMult(tN1, Wf1)
		f1Inv <- ginv(eigenMapMatMult(tN1Wf1, N1) + K1) 
		hatF1 <- eigenMapMatMult(f1Inv, tN1Wf1)
		covF1 <- eigenMapMatMult3(hatF1, V, t(hatF1))
		
		# COVARIANCE - f2
		tN2Wf2 <- eigenMapMatMult(tN2, Wf2)
		f2Inv <- ginv(eigenMapMatMult(tN2Wf2, N2) + K2)
		hatF2 <- eigenMapMatMult(f2Inv, tN2Wf2)
		covF2 <- eigenMapMatMult3(hatF2, V, t(hatF2))
		
		# FOR SIMULATION, THERE IS BIAS FOR beta, f1 & f2
		if (identical(niVec, rep(tp, m))) {
			
			# BIAS - beta
			t <- seq(1, 28)
			f1 <- sine(tp, t)
			f2 <- cosine(tp, t)
			N1f1 <- eigenMapMatMult(N1, f1)
			N2f2 <- eigenMapMatMult(N2, f2)
			Nf <- N1f1 + N2f2
			biasBeta <- eigenMapMatMult(hatB, Nf)
			
			# BIAS - f1
			temp1 <- eigenMapMatMult(tN1Wf1, N2f2) - eigenMapMatMult(K1, f1)
			biasF1 <- eigenMapMatMult(f1Inv, temp1)
			
			# BIAS - f2
			temp2 <- eigenMapMatMult(tN2Wf2, N1f1) - eigenMapMatMult(K2, f2)
			biasF2 <- eigenMapMatMult(f2Inv, temp2)
			
			newlist <- list("theta" = newTheta, "betaHat" = res$betaHat, "f1Hat" = res$f1Hat, "f2Hat" = res$f2Hat, "covBeta" = diag(covBeta), "covF1" = diag(covF1), "covF2" = diag(covF2), "biasBeta" = biasBeta, "biasF1" = biasF1, "biasF2" = biasF2)
			return(newlist)
		}
		# newlist <- list("theta" = newTheta, "betaHat" = res$betaHat, "f1Hat" = res$f1Hat, "f2Hat" = res$f2Hat, "covBeta" = covBeta, "covF1" = covF1, "covF2" = covF2)
		newlist <- list("theta" = newTheta, "betaHat" = res$betaHat, "f1Hat" = res$f1Hat, "f2Hat" = res$f2Hat, "covBeta" = diag(covBeta), "covF1" = diag(covF1), "covF2" = diag(covF2))
		return(newlist)
	}
}
