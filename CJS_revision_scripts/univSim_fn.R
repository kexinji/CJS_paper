## THIS FILE IS based on univFns_lik.R, USED FOR UNIVARIATE ANALYSIS FOR BOTH SIMULATION AND REAL DATASETS.
## IT OUTPUTS THE PARAMETER ESTIMATES, BIASES AND COVARIANCES.
## No correlations between the Gaussian random fields; no SD provided for theta's.

library(sde) # rsOU()
library(MASS) # mvrnorm()
library(psych) # tr()
library(magic) # adiag()
library(Rcpp) # TO C++
sourceCpp("matrixmulti_1.cpp") # NEED TO INSTALL Package 'RcppEigen'
library(e1071) # rwiener()

#======================================================================================
# SIMULATIONS OF DATASETS
#======================================================================================
#------------------------------------------------------------------------------
# A FUNCTION that RETURNS a UNIV. LONGITUDINAL DATASET by SIMULATION
# m is the number of subjects
# tp is the number of distinct points, assuming each subj. has the same obs'd time pt.
# beta is a vector of fixed effect parameters

newData <- function(theta, beta, m, tp, process) {
	
	n <- m*tp
	#--------------------------------------------------------------------------
	# X
	age <- sample(20:44, size = m, replace = T)
	age <- rep(age, each = tp)
	
	# f 
	f <- sine(tp, 1:tp)
	
	# random intercepts
	bSim <- rnorm(m, 0, theta[2])
	b <- rep(bSim, each = tp)
	
	# measurement error - epsilon
	eps <- rnorm(n, 0, theta[3])
	
	# Gaussian process
	if (process == "OU") {
		u <- rsOU(n, theta=c(0, theta[4], theta[5])) 
	} else if (process == "NOU") {
		t <- 1:tp
		exponent <- log(2*theta[4]) + theta[5] + theta[6]*t + theta[7]*(t^2)
		theta3 <- rep(exp(exponent/2), m)
		u <- sapply(1:n, function(v) return(rsOU(1, theta=c(0, theta[4], theta3[v]))))
	} else if (process == "Wiener"){
		u <- theta[4]*sapply(rep(1, m), function(x) return(rwiener(28, x))) # tp * m matrix
		u <- as.vector(u)
	}
	
	#--------------------------------------------------------------------------
	# univariate Longitudinal dataset
	Y <- beta*age + rep(f, m) + b + u + eps	
	id <- rep(1:m, each = tp)
	day <- rep(1:tp, m)
	simData <- data.frame(id, day, Y, age)
	
	return(simData)	
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS A PERIODIC FUNCTION
sine <- function(tp, x) {
	f1 <- 5*sin((2*pi/tp)*x)
	return(f1)
}

#======================================================================================
# HELPER FUNCTIONS
#======================================================================================
#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS K, N, r, Bdou & niVec from the given dataframe.
helperFn <- function(data, time, id) {
	#--------------------------------------------------------------------------	
	# niVec
	m <- max(id)
	
	niVec <- sapply(1:m, function(v) return(nrow(data[which(id == v),])))
	n <- sum(niVec)
		
	#--------------------------------------------------------------------------
	# r
	uniquet <- sort(unique(time))
	r <- length(uniquet)
	
	#--------------------------------------------------------------------------
	# N
	# colIndex <- sapply(time, function(v) return(which(v == uniquet))) # slower
	colIndex <- sapply(time, function(v) return(match(v, uniquet)))
	N <- matrix(0, n, r)
	N[cbind(1:n, colIndex)] <- 1
	
	#--------------------------------------------------------------------------
	# K
	h <- uniquet[2:r] - uniquet[1:(r-1)]	
	
	rk <- r -2
	# Q 
	Q <- matrix(0, (rk+2), rk)
	diag(Q) <- 1/h[1:rk]
	diag(Q[-1:-2, ]) <- 1/h[2:(rk+1)]
	diag(Q[c(-1, -r), ]) <- -1/h[1:rk] - 1/h[2:(rk+1)]
	
	# R
	dr <- (h[1:rk] + h[2:(rk+1)])/3
	R <- diag(dr, rk)
	diag(R[(2-r), -1]) <- h[2:rk]/6
	diag(R[-1, (2-r)]) <- h[2:rk]/6
	
	K <- eigenMapMatMult3(Q, ginv(R), t(Q))	
	
	#--------------------------------------------------------------------------	
	# B
	B <- Q %*% ginv(crossprod(Q)) %*% t(chol(R))
	
	# LTranspose <- chol(K, pivot=T, LDL = T) # pivot = T to handle positive-semi-definite
	# L <- t(LTranspose)
	# B <- L %*% ginv(crossprod(L))
	
	# Bstar
	Bstar <- N %*% B
	Bdou <- tcrossprod(Bstar)
	
	return(list("niVec" = niVec, "K" = K, "N" = N, "r" = r, "Bdou" = Bdou))
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS Vi 
# ni is the number of obs'ns for subj i
# ti is the obs'n time pts for subj i
matrixVi <- function(theta, ni, ti, process) {
	Zi <- matrix(1, ni)
	if (process == "OU") {
		Gammai <- outer(ti, ti, ouVar, theta[4], theta[5])
	} else if (process == "NOU") {
		Gammai <- outer(ti, ti, nouVar, theta[4], theta[5], theta[6], theta[7])	
	} else if (process == "Wiener") {
		# Gammai <- outer(ti, ti, wVar, theta[4]) # wonder why it doesn't work
		Gammai <- sapply(ti, function(y) return(sapply(ti, function(x) return(wVar(x, y, theta[4])))))
	}
	Vi <- Zi %*% theta[2] %*% t(Zi) + theta[3] * diag(ni) + Gammai
	return("Vi" = Vi)
}

#------------------------------------------------------------------------------
# HELPER FUNCTIONS that RETURN COVARIANCE FUNCTION for PROCESSES
# OU 
ouVar <- function(x, y, theta2, theta3) {
	ouV <- (theta3^2/(2*theta2)) * exp(-theta2*abs(x-y))
	return(ouV)
}

# NOU
nouVar <- function(x, y, rho, a0, a1, a2) {
	logrho <- log(rho)
	avg <- (a1*(x+y) + a2*(x^2 + y^2))/2
	nouV <- exp(logrho*abs(x-y) + a0 + avg)
	return(nouV)
}

# Wiener
wVar <- function(x, y, xi) {
	wV <- xi * min(x, y)
	return(wV)
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS TRUE IF ALL SUBJ HAVE THE SAME NUMBER OF OBSERVATIONS
test1 <- function(id, time) {
	for (i in 1:(max(id)-1)) {
		if (length(time[which(id== i)]) != length(time[which(id == i+1)])) {
			return(FALSE)
		}
	}
	return(TRUE)
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS TRUE IF ALL SUBJ HAVE THE SAME OBSERVATION TIME POINT
test2 <- function(id, time) {
	for (i in 1:(max(id)-1)) {
		if (all(time[which(id == i)] != time[which(id == i+1)])) {
			return(FALSE)
		}
	}
	return(TRUE)
}

#------------------------------------------------------------------------------
# A FUNCTION THAT RETURNS PARAMETER ESTIMATES, GIVEN VARIANCE COMPONENTS
est <- function(theta, X, Y, N, K, niVec, r, time, id, process) {
	m <- max(id)
	
	#--------------------------------------------------------------------------
	# VARIANCE MATRIX
	
	if (test1(id, time) && test2(id, time)){ # IF ALL SUBJ HAVE THE SAME OBSN TIME PT
		Vi <- matrixVi(theta, ni = niVec[1], ti = time[which(id == 1)], process)
        if(any(is.infinite(Vi))) {
            return(NA)
        }
		Wi <- ginv(Vi)
		V <- do.call("adiag", replicate(m, Vi, simplify=F))
		W <- do.call("adiag", replicate(m, Wi, simplify=F))
	} else {
		Vi <- lapply(1:m, function(v) return(matrixVi(theta, niVec[v], time[which(id == v)], process)))
        if(any(is.infinite(Vi))) {
            return(NA)
        }
		Wi <- lapply(1:m, function(v) return(ginv(Vi[[v]])))
		V <- Reduce(adiag, Vi)
		W <- Reduce(adiag, Wi)	
	}
	
	#--------------------------------------------------------------------------
	# INVERSE of MATRIX C

	WX <- W %*% X # QUICKER THAN C++
	WN <- eigenMapMatMult(W, N)
	
	lambda <- 1/theta[1]

	# The 1st row of block matrix C
	C11 <- crossprod(X, WX)
	C12 <- crossprod(X, WN)
	CRow1 <- cbind(C11, C12)
	# The 2nd row of block matrix C
 	C21 <- crossprod(N, WX)
	C22 <- crossprod(N, WN) + lambda * K
	CRow2 <- cbind(C21, C22)
	# Inverse of coefficient matrix
	C <- rbind(CRow1, CRow2)
	invC <- ginv(C)

	WY <- W %*% Y
	temp <- rbind(crossprod(X, WY), crossprod(N, WY))
	
	#--------------------------------------------------------------------------
	# ESTIMATES FOR beta, AND f, TO BE USED IN SCORE IN FISHER-SCORING
	estimates <- invC %*% temp

	nX <- ncol(X)
	betaHat <- estimates[1:nX]
	fHat <- matrix(estimates[(nX+1):(nX+r)])
	
	#--------------------------------------------------------------------------
	# Likelihood function values using new parameter estimates
	residual <- Y - X %*% betaHat - N %*% fHat
	py <- W %*% residual
	lik <- (-log(det(V)) - t(residual) %*% py - lambda * t(fHat) %*% K %*% fHat)/2
	
	newlist <- list("betaHat" = betaHat, "fHat" = fHat, "V" = V, "W" = W, "WX" = WX, "WN" = WN, "C11" = C11, "C22" = C22, "invC" = invC, "py" = py, "lik" = lik)
	return(newlist)
}

#======================================================================================
# FISHER-SCORING
#======================================================================================
#------------------------------------------------------------------------------
# HELPER FUNCTIONS that RETURN PARTIAL DERIVATIVES OF COVARIANCE FUNCTIONS

# OU PROCESS 
ouPv1 <- function(x, y, theta2, theta3) { # THETA2
	ouP1 <- (-1/theta2 - abs(x-y)) * ouVar(x, y, theta2, theta3)
	return(ouP1)
}
ouPv2 <- function(x, y, theta2, theta3) { # THETA3
	ouP2 <- (theta3/theta2) * exp(-theta2*abs(x-y))
	return(ouP2)
}

# NOU PROCESS
nouPv1 <- function(x, y, rho, a0, a1, a2){ # rho
	pv1 <- (abs(x-y)/rho) * nouVar(x, y, rho, a0, a1, a2)
	return(pv1)
}
nouPv2 <- function(x, y, rho, a0, a1, a2){ # a0
	pv2 <- nouVar(x, y, rho, a0, a1, a2)
	return(pv2)
}
nouPv3 <- function(x, y, rho, a0, a1, a2){ # a1
	pv3 <- ((x + y)/2) * nouVar(x, y, rho, a0, a1, a2)
	return(pv3)	
}
nouPv4 <- function(x, y, rho, a0, a1, a2){ # a2
	pv4 <- ((x^2 + y^2)/2) * nouVar(x, y, rho, a0, a1, a2)
	return(pv4)	
}

# Wiener 
wPv <- function(x, y) {
	pv <- min(x, y)
	return(pv)
}


#------------------------------------------------------------------------------
# A HELPER FUNCTION that RETURNS PVi, TO BE USED IN findPar()
PVi <- function(theta, ni, ti, process) {
	
	pV_b <- tcrossprod(matrix(1, ni))
	pV_eps <- diag(ni)
	
	if (process == "OU") {
		pV_theta2i <- outer(ti, ti, ouPv1, theta[4], theta[5])
		pV_theta3i <- outer(ti, ti, ouPv2, theta[4], theta[5])
		pvlist <- list("pV_bi" = pV_b, "pV_epsi" = pV_eps, "pV_theta2i" = pV_theta2i, "pV_theta3i" = pV_theta3i)
		return(pvlist)
	} else if (process == "NOU") {
		pv_rho <- outer(ti, ti, nouPv1, theta[4], theta[5], theta[6], theta[7])	
		pv_a0 <- outer(ti, ti, nouPv2, theta[4], theta[5], theta[6], theta[7])	
		pv_a1 <- outer(ti, ti, nouPv3, theta[4], theta[5], theta[6], theta[7])	
		pv_a2 <- outer(ti, ti, nouPv4, theta[4], theta[5], theta[6], theta[7])	
		pvlist <- list("pV_bi" = pV_b, "pV_epsi" = pV_eps, "pv_rho" = pv_rho, "pv_a0" = pv_a0, "pv_a1" = pv_a1, "pv_a2" = pv_a2)
		return(pvlist)
	} else if (process == "Wiener") {
		pv_xi <- sapply(ti, function(y) return(sapply(ti, function(x) return(wPv(x, y)))))
		pvlist <- list("pV_bi" = pV_b, "pV_epsi" = pV_eps, "pv_xi" = pv_xi)
		return(pvlist)
	}
}

#------------------------------------------------------------------------------
# A HELPER FUNCTION THAT RETURNS A NEW THETA AFTER ONE ITERATION OF FISHER-SCORING
findPar <- function(theta, X, Y, N, K, niVec, r, time, id, Bdou, process, dim) {
	
	res <- est(theta, X, Y, N, K, niVec, r, time, id, process)
    if(any(is.na(res))) {
        return(NA)
    }
	invC <- res$invC
	W <- res$W
	py <- res$py
	tpy <- t(py)
	m <- max(id)
	
	#--------------------------------------------------------------------------	
	# P*
	chi <- cbind(X, N)	
	Wchi <- eigenMapMatMult(W, chi)
	Pst <- W - eigenMapMatMult3(Wchi, invC, t(Wchi))
	
	#--------------------------------------------------------------------------
	# PARTIAL DERIVATIVES
	
	if (test1(id, time) && test2(id, time)){ # IF ALL SUBJ HAVE THE SAME OBSN TIME PT
		pv <- PVi(theta, niVec[1], time[which(id == 1)], process)
		pvlist1 <- lapply(1:(dim-1), function(v) return(do.call("adiag", replicate(m, pv[[v]], simplify=F))))
		pvlist <- c(list(Bdou), pvlist1)
	} else {
		pvAll <- lapply(1:m, function(v) return(PVi(theta, niVec[v], time[which(id == v)], process)))
		pvlist1 <- lapply(1:(dim-1), function(u) return(Reduce(adiag, lapply(1:m, function(v) return(pvAll[[v]][[u]])))))
		pvlist <- c(list(Bdou), pvlist1)
	}
	
	#--------------------------------------------------------------------------	
	# SCORE
	
	# Pst * pv
	Pst_pv <- lapply(1:dim, function(v) return(eigenMatMult(Pst, pvlist[[v]])))
		
	# score
	score <- sapply(1:dim, function(v) return(eigenMatMult3(tpy, pvlist[[v]], py) - tr(Pst_pv[[v]])))
 	
 	#--------------------------------------------------------------------------	
	# FISHER INFO
	
	fInfo <- sapply(1:dim, function(u) return(sapply(1:dim, function(v) return(tr(eigenMapMatMult(Pst_pv[[u]], Pst_pv[[v]]))))))

	score <- score/2
	fInfo <- fInfo/2
	Finv <- ginv(fInfo)

	thetaNew <- theta + Finv %*% score	# Parameter estimate after one iteration
	return(thetaNew)
}

#------------------------------------------------------------------------------
# A FUNCTION that RETURNS THETA, using FISHER-SCORING ALGORITHM
fisher.scoring <- function(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol, cap, process, dim) {
	
	# Initialization
	theta1 <- theta
	res <- est(theta, X, Y, N, K, niVec, r, time, id, process)
    if(any(is.na(res))) {
        return(NA)
    }
	lik1 <- res$lik
	
	iternum <- 0
	diff <- 1
	
	while (diff > tol) {
		iternum <- iternum + 1
		
	    theta0 <- theta1
	    lik0 <- lik1
	    
	    # Update theta0 -> theta1
		theta1 <- findPar(theta0, X, Y, N, K, niVec, r, time, id, Bdou, process, dim)
		if(any(is.na(theta1))) {
			return(NA)
		}
		
		#--------------------------------------------------------------------------	
		# Check whether theta1 is in parameter space
		if (process == "OU") {
			for (i in 1:dim){
				while (theta1[i]<0) {
					theta1[i] <- (theta1[i] + theta0[i])/2
				}
			}
		} else { # NOU case
			for (i in 1:(dim-3)){
				while (theta1[i]<0) {
					theta1[i] <- (theta1[i] + theta0[i])/2
				}
			}	
		}
		
		#--------------------------------------------------------------------------	
		# Find lik1 using updated theta1, which is in parameter space
		res <- est(theta1, X, Y, N, K, niVec, r, time, id, process)
		if(any(is.na(res))) {
			return(NA)
		}
		lik1 <- res$lik
		
		#--------------------------------------------------------------------------	
		# Check whether the lik is increasing 
		while (lik1 <= lik0 || is.infinite(lik1)) {
			theta1 <- (theta1 + theta0)/2
			res <- est(theta1, X, Y, N, K, niVec, r, time, id, process)
			if(any(is.na(res))) {
				return(NA)
			}
			lik1 <- res$lik
		}
		diff <- abs(lik1 - lik0)/abs(lik1)
				
		# Abort this iteration if iternum > cap, and set the parameters to be zero
		if (iternum > cap) {
			theta1 <- matrix(0, dim)
			return(list("theta" = theta1, "capConv" = iternum))
		}
		print(iternum)
		print (diff)
		print(theta1)
	}
	return(list("theta" = theta1, "capConv" = iternum))
}

#======================================================================================
# MAIN FUNCTION - OUTPUTS ESTIMATES, BIASES AND COVARIANCES
#======================================================================================
# data = NAME OF THE DATASET TO BE USED
# response = UNIVARIATE LONGITUDINAL RESPONSE, as.matrix(Y)
# fixed = FIXED EFFECT COVARIATE MATRIX, X<-cbind(X1, ..., Xp)
# random = RANDOM EFFECT COVARIATE MATRIX,
		   # the random effect covariance matrix is assumed to be unstructured.
# process = name of the Gaussian process, a character string. 
			# OU - Ornstein-Uhlenbeck process
			# NOU - Nonhomogeneous OU process with log(v(t)) = a0 + a1*t + a2*t^2
# time = OBSERVATION TIME POINTS OF OBSERVATIONS
# id = IDENTIFICATION OF SUBJECTS
# tol = tolerance level for parameters.
# cap = maximum runs of iterations for FISHER-SCORING
# THE PROGRAM ONLY ALLOWS FOR RANDOM INTERCEPT
results <- function(data, response, fixed, random, process, time, id, tol, cap) {
 	
 	#--------------------------------------------------------------------------	
 	# Determine dimention of parameter theta	
	q <- ncol(as.matrix(random))
	# random effect cov matx unstructured; 2 for smoothing param & measurement error; 
	dim <- q*(q+1)/2 + 2 	
	if (process == "OU") {
		dim <- dim + 2
	} else if (process == "NOU") {
		dim <- dim + 4
	} else if (process == "Wiener") {
		dim <- dim + 1
	}
	
	#--------------------------------------------------------------------------	
 	# Initialize parameter theta	
	theta <- matrix(0, dim)
	s <- 1 # keep track of index position of theta
	theta[s] <- 1 # smoothing parameter
	s <- 2
	if (q > 0) {
		for (i in 1:q){
			for (j in 1:i) {
				if (i == j) {
					theta[s] <- 1
				} else {
					theta[s] <- -0.5
				}
				s <- s + 1
			}
		}
	}
	theta[s] <- 1 # measurement error
	s <- s + 1
	if (process == "OU") { 
		# theta is initialized to ensure that rho = 0.0963, var = 1
		theta[s] <- -log(0.0963) # 2.340287
		theta[s+1] <- sqrt(-log(0.0963)*2) # 2.163463
	} else if (process == "NOU") {
		# theta[s] <- 0.2
		# theta[s+1] <- -0.44
		# theta[s+2] <- 0.3
		# theta[s+3] <- -0.2
		theta[s] <- 0.15
		theta[s+1] <- -1.6
		theta[s+2] <- 0.3
		theta[s+3] <- -0.1
	} else if (process == "Wiener") {
		theta[s] <- 0.7
	}

	#--------------------------------------------------------------------------	
	# Obtain values from helper functions
	helper <- helperFn(data, time, id)
	niVec <- helper$niVec
	K <- helper$K
	N <- helper$N
	r <- helper$r
	Bdou <- helper$Bdou
	
	X <- as.matrix(fixed)
	Y <- as.matrix(response)
	 
	fs <- fisher.scoring(theta, X, Y, N, K, niVec, r, time, id, Bdou, tol, cap, process, dim)
    
    # variance matrix is invertible in est() function
    if(any(is.na(fs))) {
        return(NA)
    }
    
	newTheta <- fs$theta
	
    # fisher-scoring non-convergent
    if (identical(newTheta, matrix(0, dim))) {
        return (fs)
    }

    res <- est(newTheta, X, Y, N, K, niVec, r, time, id, process)
    WX <- res$WX	
    WN <- res$WN
    V <- res$V
    W <- res$W
    
    tX <- t(X)
    tN <- t(N)
            
    # WEIGHT MATRICES, Wx & Wf	
    Wx <- W - eigenMapMatMult3(WN, ginv(res$C22), t(WN))
    Wf <- W - eigenMapMatMult3(WX, ginv(res$C11), t(WX))
    
    #--------------------------------------------------------------------------	
    # COVARIANCE - beta	
    betaInv <- ginv(eigenMatMult3(tX, Wx, X))
    hatB <- eigenMatMult3(betaInv, tX, Wx)
    covBeta <- eigenMapMatMult3(hatB, V, t(hatB))
    
    # COVARIANCE - f
    tNWf <- eigenMapMatMult(tN, Wf)
    fInv <- ginv(eigenMapMatMult(tNWf, N) + (1/theta[1]) * K) 
    hatF <- eigenMapMatMult(fInv, tNWf)
    covF <- eigenMapMatMult3(hatF, V, t(hatF))

    #--------------------------------------------------------------------------	
    # BIAS FOR beta & f, FOR SIMULATION
    if (test1(id, time) && test2(id, time)) {
        
        # BIAS - beta
        f <- sine(niVec[1], time[which(id == 1)])
        Nf <- eigenMapMatMult(N, f)
        biasBeta <- eigenMapMatMult(hatB, Nf)
        
        # BIAS - f
        Kf <- eigenMapMatMult(K, f)
        biasF <- -(1/theta[1]) * eigenMapMatMult(fInv, Kf)
        
        newlist <- list("theta" = newTheta, "betaHat" = res$betaHat, "fHat" = res$fHat, "sdBeta" = sqrt(diag(covBeta)), "sdF" = sqrt(diag(covF)), "biasBeta" = biasBeta, "biasF" = biasF, "conv" = fs$conv, "capConv" = fs$capConv)
        return(newlist)
    }

    newlist <- list("theta" = newTheta, "betaHat" = res$betaHat, "fHat" = res$fHat, "sdBeta" = sqrt(diag(covBeta)), "sdF" = sqrt(diag(covF)), "conv" = fs$conv, "capConv" = fs$capConv)
    return(newlist)
}
