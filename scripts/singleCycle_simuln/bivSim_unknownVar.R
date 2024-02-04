# Bivaraite simulations, one cycle, unknown variances & smoothing parameters

#==============================================================================
# PART I - INITIALIZATION, AND SET UP 
#==============================================================================

library(sde) # rsOU()
library(MASS) # mvrnorm()
library(psych) # tr()

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
Ni <- diag(28)
A1i <- matrix(0, 2*tp, tp) # initialize matrix A1i
A2i <- matrix(0, 2*tp, tp) # initialize matrix A2i
for (j in 1:tp){
	A1i[2*j-1, j] <- 1
	A2i[2*j, j] <- 1 
}
		
N1i <- A1i %*% Ni
N2i <- A2i %*% Ni

N1  <- do.call("rbind", replicate(n, N1i, simplify = F))
N2  <- do.call("rbind", replicate(n, N2i, simplify = F))



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
# Matrix K, does not depend on theta

# Define matrix Q, double checked.
Q <- matrix(0, tp, tp-2)
for(j in 1:(tp-2)) {
	Q[j, j] <- 1/2
	Q[j + 1, j] <- -1
	Q[j + 2, j] <- 1/2
}

# Define matrix R 
R = matrix(0, tp-2, tp-2)
for (i in 1:(tp-2)) {
	R[i, i] <- (1/3)*4
}
for (i in 1:(tp-3)) {
	R[i, i + 1] <- (1/6)*2
	R[i + 1, i] <- (1/6)*2
}

# K <- Q %*% ginv(R) %*% t(Q)
K <- eigenMapMatMult3(Q, ginv(R), t(Q))


# # > microbenchmark(Q %*% ginv(R) %*% t(Q), eigenMapMatMult3(Q, ginv(R), t(Q)))
# Unit: microseconds
                               # expr     min      lq     mean   median      uq      max neval
             # Q %*% ginv(R) %*% t(Q) 362.309 367.408 430.5811 396.4275 411.042 1893.638   100
 # eigenMapMatMult3(Q, ginv(R), t(Q)) 347.985 355.422 380.6851 384.8550 394.874  436.923   100
 
#------------------------------------------------------------------------------
# MATRICES B1* B2*

# Matrix B
LTranspose <- chol(K, pivot = T, LDL = T) # upper triangular matrix
L <- t(LTranspose) # lower triangular matrix
B <- L %*% ginv((t(L) %*% L))

# Bstar
B1star <- N1 %*% B
B2star <- N2 %*% B


# TO BE USED IN findPar
# B1dou <- B1star %*% t(B1star) 
# B2dou <- B2star %*% t(B2star) 

# COMPUTATION TIME IS HALVED
B1douNew <- tcrossprod(B1star)
B2douNew <- tcrossprod(B2star)


# > system.time(B1dou <- B1star %*% t(B1star))
   # user  system elapsed 
  # 0.105   0.016   0.121 
# > system.time(B1douNew <- tcrossprod(B1star))
   # user  system elapsed 
  # 0.048   0.002   0.050 
# > identical(B1dou, B1douNew)
# [1] TRUE


#===============================================================================
# PART II - FISHER SCORING
#===============================================================================

#------------------------------------------------------------------------------
# Variance initialization 

# parameter length of 11
# theta0 = c(tau1, tau2, phi1, phi2, phi3, sigma1, sigma2, theta12, theta13, theta22, theta23)

param <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)


#------------------------------------------------------------------------------
# ONE SIMULATION OF Y, TO TRY OUT FISHER SCORING


D <- matrix(c(param[3], param[4], param[4], param[5]), 2, 2) # variance for bi
# random intercepts;
b <- mvrnorm(n, c(0, 0), D) 


# measurement error - epsilon
var_eps <- c(param[6], param[7])
sigma <- diag(var_eps) # variance for epsilon

eps <- mvrnorm(n*tp, c(0, 0), sigma)
eps1 <- matrix(eps[,1], n, tp)
eps2 <- matrix(eps[,2], n, tp)


# bivariate Gaussian field
tempU_1 <- rsOU(n*tp, theta=c(0, param[8], param[9])) 
tempU_2 <- rsOU(n*tp, theta=c(0, param[10], param[11]))
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
# ONE ITERATION OF FISHER-SCORING, INPUT AN OLD THETA, GIVES A NEW THETA
# GIVE ESTIMATES FOR BETA, F, AND THETA
findPar = function(theta, Y){
	
	#--------------------------------------------------------------------------
	# VARIANCE for bi
	D <- matrix(c(theta[3], theta[4], theta[4], theta[5]), 2, 2) 
	
	# VARIANCE OF MEASUREMENT ERROR
	sigma_i <- diag(rep(c(theta[6], theta[7]), tp))
	

	# VARIANCE FOR OU PROCESS
	
	# NAIVE WAY, ASSIGNING VALUES TO EACH ENTRY
	# Gammai <- matrix(0, 2*tp, 2*tp)
	# for (i in 1:tp) {
		# for (j in 1:tp) {
			# Gammai[2*i-1, 2*j-1] <- (theta[9]^2 / (2*theta[8])) * exp(-theta[8]*abs(i-j))
			# Gammai[2*i, 2*j] <- (theta[11]^2 / (2*theta[10])) * exp(-theta[10]*abs(i-j))
		# }
	# }
	
	
	
# # 	# ASSIGN ONLY THE UPPER TRIANGLE, TAKE ACCOUNT OF SYMMETRY OF GAMMAi
	# Gammai <- matrix(0, 2*tp, 2*tp)
	# for (i in 1:tp) {
		# for (j in i:tp) {
			# Gammai[2*i-1, 2*j-1] <- (theta[9]^2 / (2*theta[8])) * exp(-theta[8]*abs(i-j))
			# Gammai[2*i, 2*j] <- (theta[11]^2 / (2*theta[10])) * exp(-theta[10]*abs(i-j))
			# Gammai[2*j-1, 2*i-1] <- t(Gammai[2*i-1, 2*j-1])
			# Gammai[2*j, 2*i] <- t(Gammai[2*i, 2*j])
		# }
	# }
	
	
		
	# USE VARIANCE FUNCTION, ASSIGN ONLY THE UPPER TRIANGLE, TAKE ACCOUNT OF SYMMETRY OF GAMMAi
	Gammai <- matrix(0, 2*tp, 2*tp)
	for (i in 1:tp) {
		for (j in i:tp) {
			Gammai[2*i-1, 2*j-1] <- ouVar(param[8], param[9], (i-j))
			Gammai[2*i, 2*j] <- ouVar(param[10], param[11], (i-j))
			Gammai[2*j-1, 2*i-1] <- t(Gammai[2*i-1, 2*j-1])
			Gammai[2*j, 2*i] <- t(Gammai[2*i, 2*j])
		}
	}
	
	
	
	
	Vi <- Zi %*% D %*% t(Zi) + sigma_i + Gammai
	Wi <- solve(Vi)
	W <- do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))
	

		
	
	#--------------------------------------------------------------------------
	# INVERSE of MATRIX C
	
	# HELP IMPROVE EFFICIENCY OF THE CODE
	temp4 <- W %*% X
	temp5 <- W %*% N1
	temp6 <- W %*% N2
	
	# temp4N <- crossprod(t(W), X)
	# temp5 <- crossprod(t(W), N1)
	# temp5 <- crossprod(t(W), N2)

# USUAL MULTIPLICATION IS FASTER
# # > system.time(temp4 <- W %*% X)
   # user  system elapsed 
  # 0.017   0.000   0.017 
# > system.time(temp4N <- crossprod(t(W), X))
   # user  system elapsed 
  # 0.043   0.002   0.045
  
# > system.time(eigenMapMatMult(W, X))
   # user  system elapsed 
  # 0.007   0.000   0.007 
# > system.time(W %*% X)
   # user  system elapsed 
  # 0.007   0.000   0.006 
# > system.time(eigenMapMatMult(W, N1))
   # user  system elapsed 
  # 0.022   0.001   0.021 
# > system.time(W %*% N1)
   # user  system elapsed 
  # 0.068   0.000   0.067 
# > system.time(eigenMapMatMult(W, N2))
   # user  system elapsed 
  # 0.018   0.001   0.019 
# > system.time(W %*% N2)
   # user  system elapsed 
  # 0.065   0.000   0.065 

	
	
	# # The 1st row of matrix C
	# C11 <- t(X) %*% temp4
	# C12 <- t(X) %*% temp5
	# C13 <- t(X) %*% temp6
	# CRow1 <- cbind(C11, C12, C13)
	# # The 2nd row of matrix C
	# C21 <- t(N1) %*% temp4
	# C22 <- t(N1) %*% temp5 + (1/theta[1]) * K  # lambda1 = 1/theta[1]
	# C23 <- t(N1) %*% temp6
	# CRow2 <- cbind(C21, C22, C23)
	# # The 2nd row of matrix C
	# C31 <- t(N2) %*% temp4
	# C32 <- t(N2) %*% temp5
	# C33 <- t(N2) %*% temp6 + (1/theta[2]) * K  # lambda2 = 1/theta[2]
	# CRow3 <- cbind(C31, C32, C33)
	
	
# crossprod is faster		
# > microbenchmark(t(X) %*% temp4, crossprod(X, temp4))
# Unit: microseconds
                # expr   min      lq     mean median      uq     max neval
      # t(X) %*% temp4 48.75 50.1795 54.42733 51.387 53.2400 129.721   100
 # crossprod(X, temp4) 34.84 35.7295 40.75947 36.457 40.8325 132.107   100
# > 
# > microbenchmark(t(N1) %*% temp4, crossprod(N1, temp4))
# Unit: microseconds
                 # expr     min       lq     mean   median       uq      max neval
      # t(N1) %*% temp4 321.354 609.7745 664.1395 641.0890 681.8870 1486.638   100
 # crossprod(N1, temp4) 142.786 144.4185 173.1373 152.3155 180.4615  574.017   100

	
	
	# The 1st row of matrix C	
	C11 <- crossprod(X, temp4)
	C12 <- crossprod(X, temp5)
	C13 <- crossprod(X, temp6)
	CRow1 <- cbind(C11, C12, C13)
	# The 2nd row of matrix C
	C21 <- crossprod(N1, temp4)
	C22 <- crossprod(N1, temp5) + (1/theta[1]) * K  # lambda1 = 1/theta[1]
	C23 <- crossprod(N1, temp6)
	CRow2 <- cbind(C21, C22, C23)
	# The 2nd row of matrix C
	C31 <- crossprod(N2, temp4)
	C32 <- crossprod(N2, temp5)
	C33 <- crossprod(N2, temp6) + (1/theta[2]) * K  # lambda2 = 1/theta[2]
	CRow3 <- cbind(C31, C32, C33)
	
	# Inverse of coefficient matrix
	C <- rbind(CRow1, CRow2, CRow3)
	invC <- ginv(C)
	
	
	
	


	#--------------------------------------------------------------------------
	# ESTIMATES FOR beta, f1 AND f2, TO BE USED IN SCORE IN FISHER-SCORING
	
	# HELP IMPROVE EFFICIENCY OF THE CODE
	temp7 <- W %*% Y
	
	# SLOWER
	# temp <- rbind(t(X) %*% temp7, t(N1) %*% temp7, t(N2) %*% temp7)
	
	# FASTER VERSION
	temp <- rbind(crossprod(X, temp7), crossprod(N1, temp7), crossprod(N2, temp7))
	est <- invC %*% temp
	betaHat <- est[1:2] 
	f1Hat <- matrix(est[3:30])
	f2Hat <- matrix(est[31:58])
	
# COMPUTATION TIME COMPARISON	
# > microbenchmark(rbind(t(X) %*% temp7, t(N1) %*% temp7, t(N2) %*% temp7), rbind(crossprod(X, temp7), crossprod(N1, temp7), crossprod(N2, temp7)))
# Unit: microseconds
                                                                        # expr     min       lq      mean   median        uq      max neval
                     # rbind(t(X) %*% temp7, t(N1) %*% temp7, t(N2) %*% temp7) 649.157 1196.302 1220.7773 1231.100 1259.9935 1717.917   100
 # rbind(crossprod(X, temp7), crossprod(N1, temp7), crossprod(N2,      temp7)) 196.915  199.130  229.4711  222.196  252.5065  587.163   100
 
	
	#--------------------------------------------------------------------------
	# P*
	chi <- cbind(X, N1, N2)	
	# HELP IMPROVE EFFICIENCY OF CODE, BY INTRODUCING TEMP0
	temp0 <- W %*% chi
	
# > system.time(W %*% chi)
   # user  system elapsed 
  # 0.159   0.001   0.159 
# > system.time(eigenMapMatMult(W, chi))
   # user  system elapsed 
  # 0.035   0.000   0.035 
	
	Pst <- W - temp0 %*% invC %*% t(temp0)
	
# > system.time( W - temp0 %*% invC %*% t(temp0))
   # user  system elapsed 
  # 0.120   0.003   0.121 
# > system.time(W - eigenMapMatMult3(temp1, invC, t(temp1)))
   # user  system elapsed 
  # 0.055   0.020   0.074
	
	# Pst <- W - temp0 %*% tcrossprod(invC, temp0)
	
	# # NAIVE Pst	
	# oldPst <- W - W %*% chi %*% invC %*% t(chi) %*% W

# COMPUTATION TIME COMPARISON - tcrossprod is not faster here
# > microbenchmark(W - temp0 %*% invC %*% t(temp0), W - temp0 %*% tcrossprod(invC, temp0))
# Unit: milliseconds
                                  # expr      min       lq     mean   median       uq      max neval
       # W - temp0 %*% invC %*% t(temp0) 236.0391 237.9314 244.7776 239.0673 250.8327 299.4106   100
 # W - temp0 %*% tcrossprod(invC, temp0) 236.2239 238.3861 250.0274 239.7389 250.4865 438.0305   100

	#--------------------------------------------------------------------------
	# Values needed in score, py = W * (Y - x*betaHat - N1 * f1Hat - N2 * f2Hat); tpy 
	py <- W %*% (Y - X %*% betaHat - N1 %*% f1Hat - N2 %*% f2Hat)
	tpy <- t(py)
	# oldtpy <- t(Y - X %*% betaHat - N1 %*% f1Hat - N2 %*% f2Hat) %*% W
	
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
	
	
	
	# PARTIAL DERIVATIVES 3:11
	pV <- cbind(pV_phi1, pV_phi2, pV_phi3, pV_sigma1, pV_sigma2, pV_theta12, pV_theta13, pV_theta22, pV_theta23)




	#--------------------------------------------------------------------------	
	# SCORE
	score <- matrix(0, 11, 1) 
	
	temp1 <- Pst %*% B1dou
	
# > system.time(Pst %*% B1dou)
   # user  system elapsed 
  # 4.487   0.048   4.489 
# > system.time(eigenMapMatMult(Pst, B1dou))
   # user  system elapsed 
  # 0.843   0.020   0.857 

	temp2 <- Pst %*% B2dou
	score[1] <- tpy %*% B1dou %*% py - tr(temp1)
	score[2] <- tpy %*% B2dou %*% py - tr(temp2)


	
	# ## TEMP 3:11, NAIVE VERSION
	# temp3 <- Pst %*% pV_phi1
	# temp4 <- Pst %*% pV_phi2
	# temp5 <- Pst %*% pV_phi3
	# temp6 <- Pst %*% pV_sigma1
	# temp7 <- Pst %*% pV_sigma2
	# temp8 <- Pst %*% pV_theta12
	# temp9 <- Pst %*% pV_theta13
	# temp10 <- Pst %*% pV_theta22
	# temp11 <- Pst %*% pV_theta23
	
	# # SCORE 3:11, NAIVE VERSION, W/ TEMP 3:11
	# score[3] <- tpy %*% pV_phi1 %*% py - tr(temp3)
	# score[4] <- tpy %*% pV_phi2 %*% py - tr(temp4)
	# score[5] <- tpy %*% pV_phi3 %*% py - tr(temp5)
	
	# score[6] <- tpy %*% pV_sigma1 %*% py - tr(temp6)
	# score[7] <- tpy %*% pV_sigma2 %*% py - tr(temp7)
	
	# score[8] <- tpy %*% pV_theta12 %*% py - tr(temp8)
	# score[9] <- tpy %*% pV_theta13 %*% py - tr(temp9)
	# score[10] <- tpy %*% pV_theta22 %*% py - tr(temp10)
	# score[11] <- tpy %*% pV_theta23 %*% py - tr(temp11)

	
	
	
	
	
	## TEMP 3:11, LOOP VERSION
	temp3 <- matrix(0, m, 9*m)
	for (i in 1:9) {
		temp3[, ((m*(i-1)+1):(m*i))] <- Pst %*% pV[, ((m*(i-1)+1):(m*i))]
	}
		
	# SCORE 3:11, IMPROVED LOOP VERSION, W/ temp3
	for (i in 1:9) {
		score[i+2] <- tpy %*% pV[, ((m*(i-1)+1):(m*i))] %*% py - tr(temp3[, ((m*(i-1)+1):(m*i))])
	}
	

	
	
	

	# # SCORE 3:11, LOOP VERSION
	# for (i in 3:11) {
		# score[i] <- tpy %*% pV[, (m*(i-3)+1): (m*(i-2))] %*% py - tr(Pst %*% pV[, (m*(i-3)+1): (m*(i-2))])
	# }
	
	
	# # SCORE 3:11, NAIVE VERSION
	# score[3] <- tpy %*% pV_phi1 %*% py - tr(Pst %*% pV_phi1)
	# score[4] <- tpy %*% pV_phi2 %*% py - tr(Pst %*% pV_phi2)
	# score[5] <- tpy %*% pV_phi3 %*% py - tr(Pst %*% pV_phi3)
	
	
	# score[6] <- tpy %*% pV_sigma1 %*% py - tr(Pst %*% pV_sigma1)
	# score[7] <- tpy %*% pV_sigma2 %*% py - tr(Pst %*% pV_sigma2)

	
	# score[8] <- tpy %*% pV_theta12 %*% py - tr(Pst %*% pV_theta12)
	# score[9] <- tpy %*% pV_theta13 %*% py - tr(Pst %*% pV_theta13)
	# score[10] <- tpy %*% pV_theta22 %*% py - tr(Pst %*% pV_theta22)
	# score[11] <- tpy %*% pV_theta23 %*% py - tr(Pst %*% pV_theta23)
	
	
	#--------------------------------------------------------------------------
	# FISHER INFO
	fInfo = matrix(0, 11, 11)
	
	# UPPER LEFT 4 ELEMENTS
	fInfo[1, 1] <- tr(temp1 %*% temp1)	

# > system.time(temp1 %*% temp1)
   # user  system elapsed 
  # 4.609   0.029   4.586 
# > system.time(eigenMapMatMult(temp1, temp1))
   # user  system elapsed 
  # 0.854   0.020   0.868 

	fInfo[1, 2] <- tr(temp1 %*% temp2)	
	fInfo[2, 1] <- t(fInfo[1, 2])	
	fInfo[2, 2] <- tr(temp2 %*% temp2)	
	
	
	# ### FISHER INFO, NAIVE VERSION - WON'T WORK FOR I_thetaTheta, AS THERE ARE TOO MANY TERMS.
	
	# # I_tau1Theta 
	# fInfo[1, 3] <- tr(temp1 %*% temp3)
	# fInfo[1, 4] <- tr(temp1 %*% temp4)
	# fInfo[1, 5] <- tr(temp1 %*% temp5)
	# fInfo[1, 6] <- tr(temp1 %*% temp6)
	# fInfo[1, 7] <- tr(temp1 %*% temp7)
	# fInfo[1, 8] <- tr(temp1 %*% temp8)
	# fInfo[1, 9] <- tr(temp1 %*% temp9)
	# fInfo[1, 10] <- tr(temp1 %*% temp10)
	# fInfo[1, 11] <- tr(temp1 %*% temp11)
	# fInfo[3:11, 1] <- t(fInfo[1, 3:11]) # SYMMETRIC PART
	
	
	# # I_tau2Theta
	# fInfo[2, 3] <- tr(temp2 %*% temp3)
	# fInfo[2, 4] <- tr(temp2 %*% temp4)
	# fInfo[2, 5] <- tr(temp2 %*% temp5)
	# fInfo[2, 6] <- tr(temp2 %*% temp6)
	# fInfo[2, 7] <- tr(temp2 %*% temp7)
	# fInfo[2, 8] <- tr(temp2 %*% temp8)
	# fInfo[2, 9] <- tr(temp2 %*% temp9)
	# fInfo[2, 10] <- tr(temp2 %*% temp10)
	# fInfo[2, 11] <- tr(temp2 %*% temp11)
	# fInfo[3:11, 2] <- t(fInfo[2, 3:11]) # SYMMETRIC PART
	
	
	### FISHER INFO, IMPROVED LOOP VERSION, W/ temp3
	# I_tau1Theta
	for (i in 1:9) {
		fInfo[1, i+2] <- tr(temp1 %*% temp3[, ((m*(i-1)+1):(m*i))])
	}
	fInfo[3:11, 1] <- t(fInfo[1, 3:11])


	# I_tau2Theta
	for (i in 1:9) {
		fInfo[2, i+2] <- tr(temp2 %*% temp3[, ((m*(i-1)+1):(m*i))])
	}
	fInfo[3:11, 2] <- t(fInfo[2, 3:11])
	
	
	# I_thetaTheta
	for (i in 1:9) {
		for (j in i:9) {
			fInfo[i+2, j+2] <- tr(temp3[, ((m*(i-1)+1):(m*i))] %*% temp3[, ((m*(j-1)+1):(m*j))])
			fInfo[j+2, i+2] <- t(fInfo[i+2, j+2])
		}
	}
	
	
	
	# # FISHER INFO, LOOP VERSION, VERY SLOW
	# # I_tau1Theta
	# for (i in 3:11) {
		# fInfo[1, i] <- tr(temp1 %*% Pst %*% pV[, (m*(i-3)+1): (m*(i-2))])
	# }
	# fInfo[3:11, 1] <- t(fInfo[1, 3:11])
	
	# # I_tau2Theta
	# for (i in 3:11) {
		# fInfo[2, i] <- tr(temp2 %*% Pst %*% pV[, (m*(i-3)+1): (m*(i-2))])
	# }
	# fInfo[3:11, 2] <- t(fInfo[2, 3:11])

	
	# # I_thetaTheta
	# for (i in 3:11) {
		# for (j in 3:11) {
			# fInfo[i, j] <- tr(Pst %*% pV[, (m*(i-3)+1): (m*(i-2))] %*% Pst %*% pV[, (m*(i-3)+1): (m*(i-2))])
		# }
	# }
	
	
	
	#--------------------------------------------------------------------------
	score <- score/2
	fInfo <- fInfo/2
	Finv <- ginv(fInfo)

	thetaNew <- theta + Finv %*% score	# Parameter estimate after one iteration
	
	
	newlist <- list("thetaNew" = thetaNew, "betaHat" = betaHat, "f1Hat" = f1Hat, "f2Hat" = f2Hat)
	return(newlist)
}

# > system.time(a <- findPar(param, Y))

# Nov. 22
   # user  system elapsed 
# 303.442   2.902 304.679 


# Nov. 23, after improvement, i.e. about 5 mins
   # user  system elapsed 
# 301.739   3.418 302.582 


# Jan. 11, Rcpp
# > system.time(a <- findPar(param, Y))
   # user  system elapsed 
 # 70.517   2.030  73.009 


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
		theta1 <- findPar(theta0, Y)$thetaNew
		
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


# > system.time(b <- fisher.scoring(param, Y, tol = rep(0.1, 11), cap = 10))
# [1] 1
             # [,1]        [,2]      [,3]       [,4]        [,5]      [,6]     [,7]     [,8]
# [1,] 3.533044e-05 8.24524e-06 0.1395518 0.05258852 0.007374105 0.3897276 1.525334 1.610418
         # [,9]     [,10]      [,11]
# [1,] 3.327409 0.0144218 0.01480832
           # [,1]
 # [1,] 1.4999647
 # [2,] 0.5000082
 # [3,] 1.3395518
 # [4,] 0.7474115
 # [5,] 0.6073741
 # [6,] 1.6102724
 # [7,] 4.5253343
 # [8,] 3.6104181
 # [9,] 8.3274094
# [10,] 0.1355782
# [11,] 2.0148083
# [1] 2
             # [,1]         [,2]      [,3]       [,4]        [,5]      [,6]    [,7]     [,8]
# [1,] 3.013006e-05 8.187974e-06 0.1933124 0.05976324 0.009188095 0.9632205 1.89597 2.561902
         # [,9]       [,10]       [,11]
# [1,] 2.925184 0.006995607 0.008411556
            # [,1]
 # [1,]  1.4999345
 # [2,]  0.5000164
 # [3,]  1.5328642
 # [4,]  0.6876482
 # [5,]  0.6165622
 # [6,]  0.6470520
 # [7,]  6.4213039
 # [8,]  6.1723204
 # [9,] 11.2525933
# [10,]  0.1285826
# [11,]  2.0232199
# [1] 3
             # [,1]        [,2]      [,3]       [,4]       [,5]      [,6]     [,7]     [,8]
# [1,] 2.785277e-05 8.58723e-06 0.2758571 0.08451962 0.01374168 0.3989067 2.301631 2.645214
         # [,9]       [,10]      [,11]
# [1,] 2.289587 0.003527085 0.00569751
            # [,1]
 # [1,]  1.4999067
 # [2,]  0.5000250
 # [3,]  1.8087212
 # [4,]  0.6031286
 # [5,]  0.6303039
 # [6,]  0.2481453
 # [7,]  8.7229346
 # [8,]  8.8175340
 # [9,] 13.5421803
# [10,]  0.1250555
# [11,]  2.0289174
# [1] 4
             # [,1]        [,2]       [,3]        [,4]         [,5]      [,6]      [,7]      [,8]
# [1,] 3.548644e-05 1.92644e-05 0.01921302 0.002669503 0.0001846839 0.1858426 0.1390963 0.4382068
          # [,9]        [,10]       [,11]
# [1,] 0.5604225 0.0003856422 0.003413762
             # [,1]
 # [1,]  1.49987120
 # [2,]  0.50004428
 # [3,]  1.82793425
 # [4,]  0.60579812
 # [5,]  0.63011919
 # [6,]  0.06230271
 # [7,]  8.58383830
 # [8,]  9.25574080
 # [9,] 12.98175786
# [10,]  0.12544115
# [11,]  2.03233115
# [1] 5
            # [,1]         [,2]        [,3]         [,4]         [,5]       [,6]       [,7]
# [1,] 3.39048e-05 1.852442e-05 0.008159041 0.0009943117 4.600188e-05 0.04939874 0.05155849
          # [,8]      [,9]        [,10]       [,11]
# [1,] 0.1943071 0.2417612 0.0008198243 0.003190922
             # [,1]
 # [1,]  1.49983730
 # [2,]  0.50006281
 # [3,]  1.83609329
 # [4,]  0.60679243
 # [5,]  0.63007319
 # [6,]  0.01290396
 # [7,]  8.53227981
 # [8,]  9.45004794
 # [9,] 12.73999662
# [10,]  0.12462133
# [11,]  2.03552207
# [1] 6
             # [,1]         [,2]        [,3]         [,4]         [,5]       [,6]       [,7]
# [1,] 3.345302e-05 1.822835e-05 0.003315811 0.0003757442 1.292354e-05 0.01107161 0.01962391
           # [,8]       [,9]       [,10]      [,11]
# [1,] 0.08047409 0.09946831 0.001217783 0.00304923
             # [,1]
 # [1,]  1.49980384
 # [2,]  0.50008104
 # [3,]  1.83940910
 # [4,]  0.60716817
 # [5,]  0.63006027
 # [6,]  0.00183235
 # [7,]  8.51265590
 # [8,]  9.53052203
 # [9,] 12.64052830
# [10,]  0.12340355
# [11,]  2.03857130
   # user  system elapsed 
# 390.424   5.250 392.677 


#==============================================================================
# PART III - SIMULATION 
#==============================================================================


# Store the final simulation estimates into matrices, average of beta and fhat
beta.sim = matrix(0, 2, 1)
f1.sim = matrix(0, tp, 1)
f2.sim = matrix(0, tp, 1)
theta.sim = matrix(0, 11, 1)


count = 0
i = 0 # number of convergences


# n is total number of subjects
# param is the parameter initialization
# sim is the number of convergent simulations
# totalNum is the number of total simulations, including those not convergent
# cap is the maximum number of iteration in fisher-scoring algorithm


simulation = function(n, par, sim, totalNum, tol, cap) {
	
	# Store simulation results, each column corresponds to a simulation result. 
	betaHat = matrix(0, 2, sim)
	f1Hat = matrix(0, tp, sim)
	f2Hat = matrix(0, tp, sim)
	theta = matrix(0, 11, sim)
	
	# For each simulation, obtain the parameter estimates.
	while (i < sim & count < totalNum) {
		
		count = count + 1
		print(count)
			
		
		D <- matrix(c(param[3], param[4], param[4], param[5]), 2, 2) # variance for bi
		# random intercepts;
		b <- mvrnorm(n, c(0, 0), D) 
		
		
		# measurement error - epsilon
		var_eps <- c(param[6], param[7])
		sigma <- diag(var_eps) # variance for epsilon
		
		eps <- mvrnorm(n*tp, c(0, 0), sigma)
		eps1 <- matrix(eps[,1], n, tp)
		eps2 <- matrix(eps[,2], n, tp)
		
		
		# bivariate Gaussian field
		tempU_1 <- rsOU(n*tp, theta=c(0, param[8], param[9])) 
		tempU_2 <- rsOU(n*tp, theta=c(0, param[10], param[11]))
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

	
		
		# Estimates of beta, f and theta from one simulation
		temp = fisher.scoring(param, Y, tol, cap)
		
		# Record the parameter estimates if iternum < 20
		# IF temp = 0, THEN IT MEANS fisher.scoring DIDN'T CONVERGE
		if (all(temp != 0)) {
			
			i = i + 1
			
			theta[, i] = temp
			newPar = findPar(theta[, i], Y)
	 		beta[,i] = newPar$betaHat
	 		f1hat[,i] = newPar$f1Hat
	 		f2hat[,i] = newPar$f2Hat 		
		}		
	}
	
	# Average of simulation results. 
	beta.sim[1] = mean(betaHat[1,])
	beta.sim[2] = mean(betaHat[2,]) 
	
	# f Hat
	for (k in 1:tp) {
		f1.sim[k, ] = mean(f1Hat[k,])
		f2.sim[k, ] = mean(f2Hat[k,])
	}	
	
	# SMOOTHING & VARIANCE PARAMETER ESTIMATES
	for (h in 1:11) {
		theta.sim[h, ] = mean(theta[h,])
	}
	
	newlist <- list("betaHat" = beta.sim, "f1Hat" = f1.sim, "f2Hat" = f2.sim, "theta" = theta.sim, "count" = count)
	return(newlist)

}



system.time(est <- simulation(n=30, param, sim=1, totalNum=2, tol= rep(0.1, 11), cap = 5))






#==============================================================================
# PART IV - PLOTS of SIMULATION RESULTS, beta hat, f hat
#==============================================================================


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
# PART V - BIASES and COVARIANCES
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




