want <- read.table("ageBMI_want.txt", header = T)
data <- read.table("ageBMI.txt", header = T)

library(psych) # to use function tr()

library(Rcpp)
library(microbenchmark)
sourceCpp("matrixmulti_1.cpp") # NEED TO INSTALL Package 'RcppEigen'

#===============================================================================
# PART I - INITIALIZATION, AND SET UP 
# NOTE: DO NOT NEED TO SPECIFY BETAs & f
#===============================================================================

m <- 307 # number of women in the dataset (dataset w/ complete age BMI)

# NUMBER OF OBSERVATIONS FOR EACH SUBJECT
max <- rep(0, m)
for (i in 1:307) {
	max[i] = nrow(data[which(data$womanid == i),])
}

n <- sum(max) # total observations for all subjects

#------------------------------------------------------------------------------
# RESPONSE Y

pdg <- log(data$adjpdg2)
e1c <- log(data$adje1c2)


#------------------------------------------------------------------------------
# COV MATRIX X, FIXED EFFECTS

age <- data$age
BMI <- data$BMI

X <- cbind(age, BMI)


#------------------------------------------------------------------------------
# INCIDENCE MATRIX N 

# FIND NUMBER OF DISTINCT VALUES OF TIME POINTS, r
day <- data$standday
uniqueDay <- sort(unique(day))
r <- length(uniqueDay)


# INCIDENCE MATRIX
N <- matrix(0, n, r)

for (i in 1:m) {
	ni <- max[i]
	
	Ni <- matrix(0, ni, r)
	index <- (sum(max[0:(i-1)]) + 1):sum(max[1:i]) 
	dayi <- day[index]
	for (j in 1:ni){
		for (l in 1:r){
			if(dayi[j] == uniqueDay[l]) {
				Ni[j, l] <- 1
			}
		}
	}
	
	N[index,] <- Ni
}


#------------------------------------------------------------------------------
# MATRIX K

# Define matrix Q, double checked.
Q <- matrix(0, r, r-2)
for(j in 1:(r-2)) {
	Q[j, j] <- 2
	Q[j + 1, j] <- -4
	Q[j + 2, j] <- 2
}

# Define matrix R 
R = matrix(0, r-2, r-2)
for (i in 1:(r-2)) {
	R[i, i] <- (1/3)
}
for (i in 1:(r-3)) {
	R[i, i + 1] <- (1/6)*(1/2)
	R[i + 1, i] <- (1/6)*(1/2)
}

library(MASS) # to use ginv
# K = Q %*% ginv(R) %*% t(Q)
K <- eigenMapMatMult3(Q, ginv(R), t(Q))

# > system.time(Q %*% ginv(R) %*% t(Q))
   # user  system elapsed 
  # 0.008   0.001   0.009 
# > system.time(eigenMapMatMult3(Q, ginv(R), t(Q)))
   # user  system elapsed 
  # 0.001   0.000   0.002 


#------------------------------------------------------------------------------
# B*

# Matrix B
LTranspose = chol(K, pivot = T, LDL = T) # upper triangular matrix
L = t(LTranspose) # lower triangular matrix
# B = L %*% ginv((t(L) %*% L))
B = L %*% ginv(crossprod(L))

# Bstar
Bstar = N %*% B
# Bdou = Bstar %*% t(Bstar)
Bdou <- tcrossprod(Bstar)


B1 <- Bstar
B2 <- t(Bstar)

# > dim(B1)
# [1] 8410   56

# microbenchmark() MAKES BOTH MY LAPTOP AND DESKTOP AT HOME FREEZE.
# microbenchmark(B1%*%B2, eigenMatMult(B1, B2), eigenMapMatMult(B1, B2))

# microbenchmark(Bstar%*%t(Bstar), eigenMatMult(Bstar, t(Bstar)), eigenMapMatMult(Bstar, t(Bstar)), eigenMapMatMult1(Bstar, t(Bstar)), eigenMapMatMult2(Bstar, t(Bstar)))


# CHECK PROCESSING TIME, a = b = c = d, checked using identical(A, B) or is.equal(A, B)
system.time(a <- B1%*%B2)
system.time(b <- eigenMatMult(B1, B2))
system.time(c <- eigenMapMatMult(B1, B2))
system.time(d <- tcrossprod(B1))


# NOTE: 
# User time: how many seconds the computer spent doing your calculations. 
# System time: how much time the operating system spent responding to your program's requests.
# Elapsed time: the sum of those two, plus whatever "waiting around" your program and/or the OS had to do.
# Ref: http://stackoverflow.com/questions/13688840/what-caused-my-elapsed-time-much-longer-than-user-time


# # COMPUTING TIME, 1st TIME, MACBOOK PRO

# > system.time(B1%*%B2)
   # user  system elapsed 
  # 7.852   0.520   8.305 
# > system.time(eigenMatMult(B1, B2))
   # user  system elapsed 
  # 2.075   3.064  12.730 
# > system.time(c <- eigenMapMatMult(B1, B2))
   # user  system elapsed 
  # 1.993   2.665   8.501 



# # COMPUTING TIME, 2nd TIME, MACBOOK PRO

# > system.time(a <- B1%*%B2)
   # user  system elapsed 
  # 8.357   0.622   9.559 
# > system.time(b <- eigenMatMult(B1, B2))
   # user  system elapsed 
  # 3.173   4.361  59.083 
# > system.time(c <- eigenMapMatMult(B1, B2))
   # user  system elapsed 
  # 2.545   3.774  31.323 



# # COMPUTING TIME, 3rd TIME - OFFICE COMPUTER

# > system.time(a <- B1%*%B2)
   # user  system elapsed 
  # 3.453   0.230   3.637 
# > system.time(b <- eigenMatMult(B1, B2))
   # user  system elapsed 
  # 1.180   0.405   1.572 
# > system.time(c <- eigenMapMatMult(B1, B2))
   # user  system elapsed 
  # 1.089   0.393   1.462 
# > system.time(d <- tcrossprod(B1))
   # user  system elapsed 
  # 2.152   0.188   2.315 

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
# Variance initializations
sigma.b = 2 # variance for bi, initially 1.5
sigma =  1 # variance for epsilon, initially 2
theta2 = 2
theta3 = 0.5
lambda = 0.1389352  # based on Zhang et al(1998), lambda = 1/tau

tau = 1/lambda

theta0 = matrix(c(tau, sigma.b, sigma, theta2, theta3), 5, 1)


theta <- theta0 
Y <- as.matrix(pdg)


#===============================================================================
# PART II - PARAMETER/VARIANCE ESTIMATIONS, FISHER SCORING
#===============================================================================

#------------------------------------------------------------------------------
# ONE ITERATION OF FISHER-SCORING - INPUT AN OLD THETA, GIVES A NEW THETA
# GIVE ESTIMATES FOR BETA, F, AND THETA
findPar = function(theta, Y) {
	
	#--------------------------------------------------------------------------
	# VARIANCE MATRIX
	
	W <- matrix(0, n, n)
	
	for (i in 1:m) {
		ni <- max[i]
		
		# Zi
		Zi <- matrix(0, ni, 1) 
		
		# Gamma_i
		Gammai <- matrix(0, ni, ni)
		for (j in 1:ni) {
			for (k in j:ni) {
				Gammai[j, k] <- ouVar(theta[4], theta[5], (j-k))
				Gammai[k, j] <- t(Gammai[j, k])
			}
		}
		
		# index for W
		index <- (sum(max[0:(i-1)]) + 1):sum(max[1:i]) 
 		Vi <- Zi %*% theta[2] %*% t(Zi) + theta[3] * diag(ni) + Gammai
		# Vi <- eigenMapMatMult3(Zi, theta[2], t(Zi)) + theta[3] * diag(ni) + Gammai
		W[index, index] <- solve(Vi)
	}
	
# > microbenchmark(eigenMapMatMult3(Zi, theta[2], t(Zi)), Zi %*% theta[2] %*% t(Zi))
# Unit: microseconds
                                  # expr   min     lq    mean median     uq    max neval
 # eigenMapMatMult3(Zi, theta[2], t(Zi)) 6.895 7.3845 9.39275  7.643 8.0105 59.167   100
             # Zi %*% theta[2] %*% t(Zi) 4.203 4.7395 5.72639  4.988 5.2945 59.590   100
	
	
	#--------------------------------------------------------------------------
	# INVERSE of MATRIX C
	temp1 <- W %*% X # QUICKER THAN C++
	# temp2 <- W %*% N
	temp2 <- eigenMapMatMult(W, N)
	
# > system.time(W %*% X)
   # user  system elapsed 
  # 0.168   0.001   0.167 
# > system.time(eigenMatMult(W, X))
   # user  system elapsed 
  # 0.317   0.213   0.526 
# > system.time(eigenMapMatMult(W, X))
# Error in eigenMapMatMult(W, X) : Wrong R type for mapped matrix
# Timing stopped at: 0.001 0 0.001 



# > system.time(eigenMapMatMult(W, N)) # FASTEST
   # user  system elapsed 
  # 0.830   0.007   0.858 
# > system.time(W %*% N)   
   # user  system elapsed 
  # 3.343   0.020   3.329 
# > system.time(eigenMatMult(W, N))
   # user  system elapsed 
  # 0.955   0.202   1.147 


	
	# The 1st row of block matrix C
	# C11 <- t(X) %*% temp1
	# C12 <- t(X) %*% temp2
	# CRow1 <- cbind(C11, C12)
	# The 2nd row of block matrix C
 	# C21 <- t(N) %*% temp1
	# C22 <- t(N) %*% temp2 + (1/theta[1]) * K
	# CRow2 <- cbind(C21, C22)
	

# crossprod() FASTEST	
# > microbenchmark(t(X) %*% temp1, crossprod(X, temp1), eigenMatMult(t(X), temp1))
# Unit: microseconds
                      # expr     min       lq     mean   median       uq     max neval
            # t(X) %*% temp1 112.509 117.6635 152.4510 136.4355 166.7825 610.402   100
       # crossprod(X, temp1)  50.409  51.1420  74.8155  52.2510  99.9160 135.667   100
 # eigenMatMult(t(X), temp1) 119.316 161.5345 230.8962 222.8080 265.8555 734.697   100
	
	
	# The 1st row of block matrix C
	C11 <- crossprod(X, temp1)
	C12 <- crossprod(X, temp2)
	CRow1 <- cbind(C11, C12)
	# The 2nd row of block matrix C
 	C21 <- crossprod(N, temp1)
	C22 <- crossprod(N, temp2) + (1/theta[1]) * K
	CRow2 <- cbind(C21, C22)
	# Inverse of coefficient matrix
	C <- rbind(CRow1, CRow2)
	invC <- ginv(C)
	
	
	#--------------------------------------------------------------------------
	# ESTIMATES FOR beta AND f
	temp3 <- W %*% Y # FASTER THAN C++
	
# > system.time(W %*% Y)
   # user  system elapsed 
  # 0.106   0.000   0.106 
# > system.time(eigenMapMatMult(W, Y))
   # user  system elapsed 
  # 0.169   0.001   0.169 
# > system.time(eigenMatMult(W, Y))
   # user  system elapsed 
  # 0.289   0.212   0.497 

	
	# temp <- rbind(t(X) %*% temp3, t(N) %*% temp3)
	temp <- rbind(crossprod(X, temp3), crossprod(N, temp3))
	
	estimates <- invC %*% temp
	betaHat <- estimates[1:2] 
	fHat <- matrix(estimates[3:(r+2)])
	
	#--------------------------------------------------------------------------
	# P*
	chi <- cbind(X, N)	
	# temp4 <- W %*% chi
	temp4 <- eigenMapMatMult(W, chi)
	
# FASTEST eigenMapMatMult()
# > system.time(W %*% chi)
   # user  system elapsed 
  # 3.386   0.020   3.378 
# > system.time(eigenMatMult(W, chi))
   # user  system elapsed 
  # 0.986   0.156   1.136 
# > system.time(eigenMapMatMult(W, chi))
   # user  system elapsed 
  # 0.888   0.004   0.887 
  
  	
	# Pst <- W - temp4 %*% invC %*% t(temp4)
	Pst <- W - eigenMapMatMult3(temp4, invC, t(temp4))

# # FASTEST eigenMapMatMult3()
# > system.time(W - temp4 %*% invC %*% t(temp4))
   # user  system elapsed 
  # 3.613   0.224   3.799 
# > system.time(W - eigenMapMatMult(temp4, invC) %*% t(temp4))
   # user  system elapsed 
  # 3.597   0.224   3.784 
# > system.time(W - temp4 %*% tcrossprod(invC, temp4))
   # user  system elapsed 
  # 3.594   0.223   3.780
# > system.time(W - eigenMatMult3(temp4, invC, t(temp4)))
   # user  system elapsed 
  # 1.235   0.413   1.635 
# > system.time(W - eigenMapMatMult3(temp4, invC, t(temp4)))
   # user  system elapsed 
  # 1.224   0.403   1.616 
# > temp4t <- t(temp4)
# > system.time(W - eigenMapMatMult3(temp4, invC, temp4t))
   # user  system elapsed 
  # 1.232   0.534   1.754 
	
	
	#--------------------------------------------------------------------------
	# Values needed in score, py = W * (Y - x*betaHat - N * fHat); tpy 
	py <- W %*% (Y - X %*% betaHat - N %*% fHat)
	tpy <- t(py)

# USING %*% is FASTER

# > system.time(Y - eigenMatMult(X, betaHat) - eigenMatMult(N, fHat))
   # user  system elapsed 
  # 0.002   0.000   0.002 
# > system.time(Y - X %*% betaHat - N %*% fHat)
   # user  system elapsed 
  # 0.001   0.000   0.001 

# > system.time(eigenMapMatMult(W, (Y - X %*% betaHat - N %*% fHat)))
   # user  system elapsed 
  # 0.167   0.000   0.166 
# > system.time(W %*% (Y - X %*% betaHat - N %*% fHat))
   # user  system elapsed 
  # 0.108   0.000   0.108 


	#--------------------------------------------------------------------------
	# PARTIAL DERIVATIVES
	
	# Partial derivative wrt b, epsilon, theta2, theta3
	parV_b <- matrix(0, n, n)
	parV_eps <- matrix(0, n, n)
	parV_theta2 <- matrix(0, n, n)
	parV_theta3 <- matrix(0, n, n)

	for (i in 1:m) {
		ni <- max[i]
		
		Zi <- matrix(0, ni, 1) 
		index <- (sum(max[0:(i-1)]) + 1):sum(max[1:i])
		
		parV_b[index, index] <- tcrossprod(Zi) # Zi %*% t(Zi)
		parV_eps[index, index] <- diag(ni)
		
		
		parV_theta2i <- matrix(0, ni, ni)
		parV_theta3i <- matrix(0, ni, ni)
		for (j in 1:ni) {
			for (k in j:ni) {
				parV_theta2i[j, k] <- ouPv1(theta[4], theta[5], (i-j))
				parV_theta2i[k, j] <- t(parV_theta2i[j, k])
				
				parV_theta3i[j, k] <- ouPv2(theta[4], theta[5], (i-j))
				parV_theta3i[k, j] <- t(parV_theta3i[j, k])			
			}
		}
		parV_theta2[index, index] <- parV_theta2i
		parV_theta3[index, index] <- parV_theta3i
	}
	
	
	pV <- cbind(parV_b, parV_eps, parV_theta2, parV_theta3)

	
	
	#--------------------------------------------------------------------------
	# SCORE
	score <- matrix(0, 5, 1) 
	
	# FOR SIMPLIFICATION OF CODE
	# temp5 <- (Pst %*% Bdou) # Pst, Bdou 8410 * 8410
	temp5 <- eigenMapMatMult(Pst, Bdou)
	
	
	# system.time(a <- Pst%*% Bdou)
	# system.time(b <- eigenMatMult(Pst, Bdou))
	# system.time(c <- eigenMapMatMult(Pst, Bdou))
	
# COMPUTING TIME, OFFICE COMPUTER JAN. 4, 2017
# > 	system.time(a <- Pst%*% Bdou)
   # user  system elapsed 
# 482.721   3.535 481.437 
# > 	system.time(b <- eigenMatMult(Pst, Bdou))
   # user  system elapsed 
# 102.669   1.428 103.107 
# > 	system.time(c <- eigenMapMatMult(Pst, Bdou))
   # user  system elapsed 
# 102.060   0.859 102.181 

	
	
	# temp6, the LOOP TAKES 7 MINS
	temp6 <- matrix(0, n, 5*n)
	for (i in 1:4) {
		ind <- (n*(i-1)+1):(n*i)
		# temp6[,ind] <- Pst %*% pV[, ind]
		temp6[,ind] <- eigenMapMatMult(Pst, pV[, ind])
	}
	
	# TOTAL TIME user 410.196 elapsed 411.2 - QUICKER THAN LOOP
	temp61 <- eigenMapMatMult(Pst, parV_b)
	temp62 <- eigenMapMatMult(Pst, parV_eps)
	temp63 <- eigenMapMatMult(Pst, parV_theta2)
	temp64 <- eigenMapMatMult(Pst, parV_theta3)
	
# > 	temp6 <- matrix(0, n, 5*n)
# > system.time(	for (i in 1:4) {
# + 		ind <- (n*(i-1)+1):(n*i)
# + 		# temp6[,ind] <- Pst %*% pV[, ind]
# + 		temp6[,ind] <- eigenMapMatMult(Pst, pV[, ind])
# + 	})
   # user  system elapsed 
# 425.794   5.976 428.474 	
	
# > system.time(temp61 <- eigenMapMatMult(Pst, parV_b))
   # user  system elapsed 
# 102.484   1.076 102.560 
# > system.time(temp62 <- eigenMapMatMult(Pst, parV_eps))
   # user  system elapsed 
# 102.475   1.160 102.630 
# > system.time(temp63 <- eigenMapMatMult(Pst, parV_theta2))
   # user  system elapsed 
# 102.659   1.076 102.744 
# > system.time(temp64 <- eigenMapMatMult(Pst, parV_theta3))
   # user  system elapsed 
# 102.578   1.356 103.266 

	
	

	# THE FOLLOWING SCORE CODE TAKES ABOUT 1 MIN
	# score[1] <- tpy %*% Bdou %*% py - tr(temp5)
	score[1] <- eigenMapMatMult3(tpy, Bdou, py) - tr(temp5)
	
# > system.time(tpy %*% Bdou %*% py)
   # user  system elapsed 
  # 0.261   0.001   0.261 
# > system.time(eigenMapMatMult3(tpy, Bdou, py))
   # user  system elapsed 
  # 0.201   0.003   0.202 
	
	# SCORE 2:5
	for (i in 1:4) {
		ind <- (n*(i-1)+1):(n*i)
		# score[i+1] <- tpy %*% pV[, ind] %*% py - tr(temp6[, ind])
		score[i+1] <- eigenMapMatMult3(tpy, pV[, ind], py) - tr(temp6[, ind])
	}
	
# > system.time(	for (i in 1:4) {
# + 		ind <- (n*(i-1)+1):(n*i)
# + 		# score[i+1] <- tpy %*% pV[, ind] %*% py - tr(temp6[, ind])
# + 		score[i+1] <- eigenMapMatMult3(tpy, pV[, ind], py) - tr(temp6[, ind])
# + 	}
# + )
   # user  system elapsed 
 # 34.347   2.479  36.592 
 
 	# COMPUTE SCORES 2:5 SEPARATELY
 	score[2] <- eigenMapMatMult3(tpy, parV_b, py) - tr(temp61)
 	score[3] <- eigenMapMatMult3(tpy, parV_eps, py) - tr(temp62)
 	score[4] <- eigenMapMatMult3(tpy, parV_theta2, py) - tr(temp63)
 	score[5] <- eigenMapMatMult3(tpy, parV_theta3, py) - tr(temp64)


# TOTAL TIME: user 0.814, elapsed 0.82, i.e. QUICKER THAN LOOP
# > system.time(score[2] <- eigenMapMatMult3(tpy, parV_b, py) - tr(temp61))
   # user  system elapsed 
  # 0.204   0.004   0.206 
# > system.time(score[3] <- eigenMapMatMult3(tpy, parV_eps, py) - tr(temp62))
   # user  system elapsed 
  # 0.202   0.003   0.204 
# > system.time(score[4] <- eigenMapMatMult3(tpy, parV_theta2, py) - tr(temp63))
   # user  system elapsed 
  # 0.205   0.004   0.206 
# > system.time(score[5] <- eigenMapMatMult3(tpy, parV_theta3, py) - tr(temp64))
   # user  system elapsed 
  # 0.203   0.003   0.204 
	


	
	#--------------------------------------------------------------------------	
	# FISHER INFO
	fInfo <- matrix(0, 5, 5)

	# I_tauTau
	# fInfo[1, 1] <- tr(temp5 %*% temp5)	
	fInfo[1, 1] <- tr(eigenMapMatMult(temp5, temp5))
	
# > 	system.time(tr(eigenMapMatMult(temp5, temp5)))
   # user  system elapsed 
# 105.048   1.644 106.628	

	# I_tauTheta, the LOOP TAKES 7 MINS
	for (i in 1:4) {
		ind <- (n*(i-1)+1):(n*i)
		# fInfo[1, i+1] <- tr(temp5 %*% temp6[, ind])
		fInfo[1, i+1] <- tr(eigenMapMatMult(temp5, temp6[, ind]))		
	}
	fInfo[2:5, 1] <- t(fInfo[1, 2:5])
	
	
	# I_tauTheta, COMPUTE DIRECTLY - FASTER
	fInfo[1, 2] <- tr(eigenMapMatMult(temp5, temp61))
	fInfo[1, 3] <- tr(eigenMapMatMult(temp5, temp62))
	fInfo[1, 4] <- tr(eigenMapMatMult(temp5, temp63))
	fInfo[1, 5] <- tr(eigenMapMatMult(temp5, temp64))
	
	fInfo[2, 1] <- t(fInfo[1, 2])
	fInfo[3, 1] <- t(fInfo[1, 3])
	fInfo[4, 1] <- t(fInfo[1, 4])
	fInfo[5, 1] <- t(fInfo[1, 5])

# > system.time(	for (i in 1:4) {
# + 		ind <- (n*(i-1)+1):(n*i)
# + 		# fInfo[1, i+1] <- tr(temp5 %*% temp6[, ind])
# + 		fInfo[1, i+1] <- tr(eigenMapMatMult(temp5, temp6[, ind]))		
# + 	})
   # user  system elapsed 
# 433.942   6.795 439.334 


# > system.time(fInfo[1, 2] <- tr(eigenMapMatMult(temp5, temp61))
# + )
   # user  system elapsed 
# 102.105   1.090 102.405 

	
	
	
	# I_thetaTheta, FORCE STOPPED THE LOOP AFTER 53 MINS
	for (i in 1:4) {
		for (j in i:4) {
			# fInfo[i+1, j+1] <- tr(temp6[, (n*(i-1)+1):(n*i)] %*% temp6[, (n*(j-1)+1):(n*j)])
			fInfo[i+1, j+1] <- tr(eigenMapMatMult(temp6[, (n*(i-1)+1):(n*i)], temp6[, (n*(j-1)+1):(n*j)]))
			fInfo[j+1, i+1] <- t(fInfo[i+1, j+1])
		}
	}
	
	
	fInfo[2, 2] <- tr(eigenMapMatMult(temp61, temp61))
	fInfo[2, 3] <- tr(eigenMapMatMult(temp61, temp62))
	fInfo[2, 4] <- tr(eigenMapMatMult(temp61, temp63))
	fInfo[2, 5] <- tr(eigenMapMatMult(temp61, temp64))
	
	fInfo[3, 2] <- t(fInfo[2, 3])
	fInfo[3, 3] <- tr(eigenMapMatMult(temp62, temp62))
	fInfo[3, 4] <- tr(eigenMapMatMult(temp62, temp63))
	fInfo[3, 5] <- tr(eigenMapMatMult(temp62, temp64))
	
	fInfo[3, 3] <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_eps))
	fInfo[3, 4] <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_theta2))
	fInfo[3, 5] <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_theta3))


# CONCLUSION: USE eigenMapMatMult FOR fInfo[3, 3:5], fInfo[2, 2:5]
# 			Use eigenMapMatMult4 FOR fInfo[4, 4:5], fInfo[5, 5]

# > system.time(fInfo[3, 3] <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_eps)))
   # user  system elapsed 
# 313.635   2.770 314.322 
# > system.time(fInfo[3, 3] <- tr(eigenMapMatMult(temp62, temp62)))
   # user  system elapsed 
# 106.034   2.411 147.064 


# > system.time(fInfo34 <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_theta2)))
   # user  system elapsed 
# 316.702   3.083 318.086 
# > system.time(fInfo[3, 4] <- tr(eigenMapMatMult(temp62, temp63)))
   # user  system elapsed 
# 106.163   2.177 127.145 


# > system.time(fInfo35 <- tr(eigenMapMatMult4(Pst, parV_eps, Pst, parV_theta3)))
   # user  system elapsed 
# 305.007   2.119 304.331 
# > system.time(fInfo[3, 5] <- tr(eigenMapMatMult(temp62, temp64)))
   # user  system elapsed 
# 103.689   1.490 103.895 

	

	fInfo[4, 2] <- t(fInfo[2, 4])
	fInfo[4, 3] <- t(fInfo[3, 4])
	# fInfo[4, 4] <- tr(eigenMapMatMult(temp63, temp63)) # START FREEZING
	# fInfo[4, 5] <- tr(eigenMapMatMult(temp63, temp64))
	fInfo[4, 4] <- tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta2))
	fInfo[4, 5] <- tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta3))
	
	
# > system.time(tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta2)))
   # user  system elapsed 
# 325.460   3.975 330.485 
# > system.time(fInfo[4, 4] <- tr(eigenMapMatMult(temp63, temp63)))
   # user  system elapsed 
# 796.900   5.984 806.195 


	
# > system.time(fInfo45 <- tr(temp63 %*% temp64))
    # user   system  elapsed 
# 1995.746   18.024 2015.924 
# > system.time(fInfo[4, 5] <- tr(eigenMapMatMult(temp63, temp64)))
   # user  system elapsed 
# 806.893   7.098 815.577
# > system.time(fInfo45C <- tr(eigenMatMult(temp63, temp64)))
   # user  system elapsed 
# 807.639   8.933 829.704 
# > system.time(tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta3)))
   # user  system elapsed 
# 323.599   5.948 345.543 
	
	
	fInfo[5, 2] <- t(fInfo[2, 5])
	fInfo[5, 3] <- t(fInfo[3, 5])
	fInfo[5, 4] <- t(fInfo[3, 5])
	# fInfo[5, 5] <- tr(eigenMapMatMult(temp64, temp64))
	fInfo[5, 5] <- tr(eigenMapMatMult4(Pst, parV_theta3, Pst, parV_theta3))
	
# > system.time(fInfo[5, 5] <- tr(eigenMapMatMult(temp64, temp64)))
   # user  system elapsed 
# 813.827   8.043 852.307 
# > system.time(fInfo55 <- tr(eigenMapMatMult4(Pst, parV_theta3, Pst, parV_theta3)))
   # user  system elapsed 
# 323.828   4.705 380.086 	
	
	score <- score/2
	fInfo <- fInfo/2
	Finv <- ginv(fInfo)

	thetaNew <- theta + Finv %*% score	# Parameter estimate after one iteration
	
	
	newlist <- list("thetaNew" = theta + Finv %*% score, "betaHat" = estimates[1], "fHat" = matrix(estimates[2:16]))
	return(newlist)
}

# > system.time(a <- findPar(theta0, Y))
    # user   system  elapsed 
# 2707.782   28.427 2730.946 
# > 2707.782/60
# [1] 45.1297

# > system.time(a <- findPar(theta0, Y))
    # user   system  elapsed 
# 2751.615   32.714 2776.720 
 

#------------------------------------------------------------------------------
# FISHER-SCORING ALGORITHM

iternum = 0
diff = matrix(1, 1, 5)

# tol is a vector of tolerance levels for different parameters.
# cap is the maximum runs of iteration, after which the iteration will stop.
fisher.scoring = function(par, Y, tol, cap) {
	
	theta1 = par
	
	while (diff[1] > tol[1] | diff[2] > tol[2] | diff[3] > tol[3] | diff[4] > tol[4] | diff[5] > tol[5]) {
		iternum = iternum + 1
		
	    theta0 = theta1
		theta1 = findPar(theta0, Y)$thetaNew
		
		for (i in 1:5){
			while (theta1[i]<0) {
				theta1[i] = (theta1[i] + theta0[i])/2
			}
		}
		
		for (i in 1:5) {
			diff[i] = abs(theta1[i] - theta0[i])
		}
		
		# Abort this iteration if iternum > cap, and set the parameters to be zero
		if (iternum > cap) {
			theta1 = matrix(0, 5)
			return (theta1)
		}
	
		print(iternum)
		print (diff)
		print(theta1)
	}
	
	return(theta1)
}
