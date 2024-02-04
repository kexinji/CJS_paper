want <- read.table("ageBMI_want.txt", header = T)
data <- read.table("ageBMI.txt", header = T)

library(psych) # to use function tr()
library(MASS) # to use ginv()
library(Rcpp) # TO C++
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

K <- eigenMapMatMult3(Q, ginv(R), t(Q))


#------------------------------------------------------------------------------
# B*

# Matrix B
LTranspose <- chol(K, pivot = T, LDL = T) # upper triangular matrix
L <- t(LTranspose) # lower triangular matrix
B <- L %*% ginv(crossprod(L))

# Bstar
Bstar <- N %*% B
Bdou <- tcrossprod(Bstar)


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

# lambda = 0.1389352  # based on Zhang et al(1998), lambda = 1/tau
# tau = 1/lambda

# theta0 <- matrix(c(tau, sigma.b, sigma, theta2, theta3), 5, 1)
theta0 <- matrix(c(7.2976, 2, 1, 2, 0.5), 5, 1)
Y <- as.matrix(pdg)


#===============================================================================
# PART II - PARAMETER/VARIANCE ESTIMATIONS, FISHER SCORING
#===============================================================================

#------------------------------------------------------------------------------
# GIVEN theta & Y, GIVE ESTIMATES FOR BETA and F.
# p is # of fixed effect parameters; p = 2 here

est = function(theta, Y, p) {
	
	#--------------------------------------------------------------------------
	# VARIANCE MATRIX
	
	W <- matrix(0, n, n)
	V <- matrix(0, n, n)
	
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
 		V[index, index] <- Vi
		W[index, index] <- solve(Vi)
	}

	
	
	#--------------------------------------------------------------------------
	# INVERSE of MATRIX C

	WX <- W %*% X # QUICKER THAN C++
	WN <- eigenMapMatMult(W, N)

	# The 1st row of block matrix C
	C11 <- crossprod(X, WX)
	C12 <- crossprod(X, WN)
	CRow1 <- cbind(C11, C12)
	# The 2nd row of block matrix C
 	C21 <- crossprod(N, WX)
	C22 <- crossprod(N, WN) + (1/theta[1]) * K
	CRow2 <- cbind(C21, C22)
	# Inverse of coefficient matrix
	C <- rbind(CRow1, CRow2)
	invC <- ginv(C)
	
	
	#--------------------------------------------------------------------------
	# ESTIMATES FOR beta AND f
	temp1 <- W %*% Y # FASTER THAN C++	
	temp <- rbind(crossprod(X, temp1), crossprod(N, temp1))
	
	estimates <- invC %*% temp
	
	res <- list("betaHat" = estimates[1:p], "fHat" = matrix(estimates[(p+1):(r+p)]), "V" = V, "W" = W, "WX" = WX, "WN" = WN, "C11" = C11, "C22" = C22, "invC" = invC)
	return(res)
}



#------------------------------------------------------------------------------
# ONE ITERATION OF FISHER-SCORING - INPUT AN OLD THETA, GIVES A NEW THETA

findPar = function(theta, Y) {
	
	
	res <- est(theta, Y, 2)
	betaHat <- res$betaHat
	fHat <- res$fHat
	W <- res$W
	invC <- res$invC


	#--------------------------------------------------------------------------
	# P*
	chi <- cbind(X, N)	
	
	temp4 <- eigenMapMatMult(W, chi)
	Pst <- W - eigenMapMatMult3(temp4, invC, t(temp4))
	
	
	#--------------------------------------------------------------------------
	# Values needed in score, py = W * (Y - x*betaHat - N * fHat); tpy 
	
	py <- W %*% (Y - X %*% betaHat - N %*% fHat)
	tpy <- t(py)

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
	
	
	#--------------------------------------------------------------------------
	# SCORE
	score <- matrix(0, 5, 1) 
	
	# FOR SIMPLIFICATION OF CODE
	temp5 <- eigenMapMatMult(Pst, Bdou)

	# FASTER THAN LOOP
	temp61 <- eigenMapMatMult(Pst, parV_b)
	temp62 <- eigenMapMatMult(Pst, parV_eps)
	temp63 <- eigenMapMatMult(Pst, parV_theta2)
	temp64 <- eigenMapMatMult(Pst, parV_theta3)

	# SCORES, FASTER THAN LOOP
	score[1] <- eigenMapMatMult3(tpy, Bdou, py) - tr(temp5)
	score[2] <- eigenMapMatMult3(tpy, parV_b, py) - tr(temp61)
 	score[3] <- eigenMapMatMult3(tpy, parV_eps, py) - tr(temp62)
 	score[4] <- eigenMapMatMult3(tpy, parV_theta2, py) - tr(temp63)
 	score[5] <- eigenMapMatMult3(tpy, parV_theta3, py) - tr(temp64)

	
	#--------------------------------------------------------------------------	
	# FISHER INFO
	fInfo <- matrix(0, 5, 5)

	# I_tauTau
	fInfo[1, 1] <- tr(eigenMapMatMult(temp5, temp5))

	# I_tauTheta, FASTER THAN LOOP
	fInfo[1, 2] <- tr(eigenMapMatMult(temp5, temp61))
	fInfo[1, 3] <- tr(eigenMapMatMult(temp5, temp62))
	fInfo[1, 4] <- tr(eigenMapMatMult(temp5, temp63))
	fInfo[1, 5] <- tr(eigenMapMatMult(temp5, temp64))
	
	fInfo[2:5, 1] <- t(fInfo[1, 2:5])

	# I_thetaTheta
	fInfo[2, 2] <- tr(eigenMapMatMult(temp61, temp61))
	fInfo[2, 3] <- tr(eigenMapMatMult(temp61, temp62))
	fInfo[2, 4] <- tr(eigenMapMatMult(temp61, temp63))
	fInfo[2, 5] <- tr(eigenMapMatMult(temp61, temp64))
	
	fInfo[3, 2] <- t(fInfo[2, 3])
	fInfo[3, 3] <- tr(eigenMapMatMult(temp62, temp62))
	fInfo[3, 4] <- tr(eigenMapMatMult(temp62, temp63))
	fInfo[3, 5] <- tr(eigenMapMatMult(temp62, temp64))
	
	fInfo[4, 2:3] <- t(fInfo[2:3, 4])
	fInfo[4, 4] <- tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta2))
	fInfo[4, 5] <- tr(eigenMapMatMult4(Pst, parV_theta2, Pst, parV_theta3))
	
	fInfo[5, 2:4] <- t(fInfo[2:4, 5])
	fInfo[5, 5] <- tr(eigenMapMatMult4(Pst, parV_theta3, Pst, parV_theta3))

	score <- score/2
	fInfo <- fInfo/2
	Finv <- ginv(fInfo)

	thetaNew <- theta + Finv %*% score	# Parameter estimate after one iteration
	return(thetaNew)
}

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
		theta1 = findPar(theta0, Y)
		
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

# cpu135 Jan. 17, 2017

# # > system.time(b <- fisher.scoring(theta0, Y, rep(0.1, 5), cap = 10)) 
# [1] 1
            # [,1] [,2]      [,3]         [,4]         [,5]
# [1,] 1.59448e-05    0 0.1748501 3.906539e-05 0.0002036423
          # [,1]
# [1,] 7.1975841
# [2,] 2.0000000
# [3,] 0.8251499
# [4,] 2.0000391
# [5,] 0.4997964
# [1] 2
            # [,1] [,2]        [,3]      [,4]      [,5]
# [1,] 1.59367e-05    0 0.004989596 0.1067546 0.8000246
          # [,1]
# [1,] 7.1975681
# [2,] 2.0000000
# [3,] 0.8201603
# [4,] 1.8932845
# [5,] 1.2998209
# [1] 3
             # [,1] [,2]      [,3]       [,4]      [,5]
# [1,] 1.579727e-05    0 0.4819564 0.09456295 0.2576995
          # [,1]
# [1,] 7.1975523
# [2,] 2.0000000
# [3,] 0.3382039
# [4,] 1.7987215
# [5,] 1.5575204
# [1] 4
             # [,1] [,2]      [,3]       [,4]     [,5]
# [1,] 1.557294e-05    0 0.3196766 0.08510072 0.185646
           # [,1]
# [1,] 7.19753675
# [2,] 2.00000000
# [3,] 0.01852726
# [4,] 1.71362078
# [5,] 1.74316646
# [1] 5
             # [,1] [,2]       [,3]       [,4]      [,5]
# [1,] 1.539049e-05    0 0.01709046 0.07529032 0.1422211
          # [,1]
# [1,] 7.1975214
# [2,] 2.0000000
# [3,] 0.0014368
# [4,] 1.6383305
# [5,] 1.8853876
# [1] 6
             # [,1] [,2]         [,3]     [,4]      [,5]
# [1,] 1.543639e-05    0 0.0008596987 0.069877 0.1158099
             # [,1]
# [1,] 7.1975059214
# [2,] 2.0000000000
# [3,] 0.0005771014
# [4,] 1.5684534645
# [5,] 2.0011974598
# [1] 7
             # [,1]         [,2]         [,3]       [,4]       [,5]
# [1,] 1.550566e-05 2.220446e-16 0.0002957576 0.06447645 0.09554797
             # [,1]
# [1,] 7.1974904157
# [2,] 2.0000000000
# [3,] 0.0002813438
# [4,] 1.5039770105
# [5,] 2.0967454272
    # user   system  elapsed 
# 26474.28    57.64 36648.88 
# > b
# $theta
             # [,1]
# [1,] 7.1974904157
# [2,] 2.0000000000
# [3,] 0.0002813438
# [4,] 1.5039770105
# [5,] 2.0967454272

# $betaHat
# [1] -0.007556841 -0.131601277

# $fHat
              # [,1]
 # [1,]  0.661890432
 # [2,]  0.792447368
 # [3,]  0.524828360
 # [4,]  0.442099244
 # [5,]  0.276281954
 # [6,]  0.230252559
 # [7,]  0.200110719
 # [8,]  0.090416249
 # [9,]  0.022966370
# [10,]  0.008327147
# [11,] -0.006443232
# [12,] -0.078972419
# [13,] -0.115729131
# [14,] -0.140138720
# [15,] -0.067177847
# [16,] -0.164185928
# [17,] -0.092654095
# [18,] -0.086984134
# [19,] -0.073670017
# [20,] -0.117041937
# [21,] -0.047672451
# [22,] -0.080525686
# [23,] -0.021239052
# [24,] -0.014228693
# [25,]  0.086467716
# [26,]  0.129322267
# [27,]  0.196120718
# [28,]  0.279795383
# [29,]  0.299414958
# [30,]  0.501914087
# [31,]  0.548071623
# [32,]  0.662054724
# [33,]  0.897055881
# [34,]  0.935614353
# [35,]  1.110004659
# [36,]  1.204731519
# [37,]  1.355641114
# [38,]  1.473167971
# [39,]  1.648265876
# [40,]  1.765636825
# [41,]  1.819371796
# [42,]  1.858976062
# [43,]  1.947170534
# [44,]  1.896975343
# [45,]  1.937233179
# [46,]  1.963261871
# [47,]  2.118718303
# [48,]  1.994567083
# [49,]  1.871302923
# [50,]  1.877906971
# [51,]  1.879596830
# [52,]  1.738904090
# [53,]  1.604582100
# [54,]  1.520635628
# [55,]  1.559807757
# [56,]  1.251071016


#======================================================================================
# PART III - COVARIANCES for beta, f
#======================================================================================

# INPUT A VECTOR OF PARAMETERS theta, RETURN COVARIANCES for BETA & f.
# theta COMES FROM fisher.scoring()

# system.time(theta <- fisher.scoring(theta0, Y, rep(0.1, 5), cap = 10)) 

Cov <- function(theta){
	
	res <- est(theta, Y, 2) 
	V <- res$V
	W <- res$W
	WN <- res$WN	
	WX <- res$WX
	
	tX <- t(X)
	tN <- t(N)

	# WEIGHT MATRICES
	
	Wx <- W - eigenMapMatMult3(WN, ginv(res$C22), t(WN))
	Wf <- W - eigenMapMatMult3(WX, ginv(res$C11), t(WX))
	
	# COVARIANCES - beta
	
	betaInv <- ginv(eigenMatMult3(tX, Wx, X)) 
	hatB <- eigenMatMult3(betaInv, tX, Wx)
	
	covBeta <- eigenMapMatMult3(hatB, V, t(hatB))
	
	
	# COVARIANCES - f
	
	fInv <- ginv(eigenMapMatMult3(tN, Wf, N) + (1/theta[1]) * K) 
	hatF <- eigenMapMatMult3(fInv, tN, Wf)
	
	covF <- eigenMapMatMult3(hatF, V, t(hatF))
	
	
	newlist <- list("covBeta" = covBeta, "covF" = covF, "betaHat" = res$betaHat, "fHat" = res$fHat)
	return(newlist)
}

write.csv(unlist(Cov(theta)), "Liu_A1_tol0.1_cap10_Cov", row.names = F)