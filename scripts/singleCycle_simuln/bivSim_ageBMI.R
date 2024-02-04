want <- read.table("ageBMI_want.txt", header = T)
data <- read.table("ageBMI.txt", header = T)

m <- 307 # number of women in the dataset (dataset w/ complete age BMI)

# find max in the new dataset after deleting rows with missing data 
max <- rep(0, m)
for (i in 1:307) {
	max[i] = nrow(data[which(data$womanid == i),])
}

n <- sum(max) # total observations for all subjects
n1 <- max[1]

# data visualization for one woman
par(mfrow = c(2, 1))
plot(log(data[1:n1, 3]), type = "l", ylab = "Log progesterone", xlab = "Days in standardized menstrual cycle")
plot(log(data[1:n1, 4]), type = "l", ylab = "Log estrogen", xlab = "Days in standardized menstrual cycle")


# plot(data[1:n1, 3], type = "l", ylab = "Progesterone Level", xlab = "day", ylim = c(1,80))
# par(new = T)
# plot(data[1:n1, 4], type = "l", ylab = "Estrogen Level", xlab = "day", ylim = c(1,80))



# BIVARIATE RESPONSE
Y <- matrix(0, 2*n)
for (i in 1:n) {
	Y[2*i - 1] = log(data$adjpdg2[i])
	Y[2*i] = log(data$adje1c2[i])
}

# COVARIATE X, indexes different from simulation, since ni are different for i
age <- want$age
BMI <- want$BMI
X <- matrix(0, 2*n, 4) # initialization
X[1:(2*n1), 1] = rep(c(age[1], 0), n1) 
X[1:(2*n1), 2] = rep(c(BMI[1], 0), n1) 
X[1:(2*n1), 3] = rep(c(0, age[1]), n1) 
X[1:(2*n1), 4] = rep(c(0, BMI[1]), n1) 
for (i in 2:m) {
	ni = max[i]
	
	# the 1st column of Xi is age[1], 0, age[2], 0...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 1] = rep(c(age[i], 0), ni) 
	# the 2nd column of Xi is BMI[1], 0, BMI[2], 0...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 2] = rep(c(BMI[i], 0), ni) 
	
	# the 3rd column of Xi is 0, age[1], 0, age[2]...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 3] = rep(c(0, age[i]), ni)
	# the 4th column of Xi is 0, BMI[1], 0, BMI[2]...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 4] = rep(c(0, BMI[i]), ni) 
}



# INVERSE VARIANCE MATRIX W
W <- diag(0, 2*n) # initialization of W
V <- diag(0, 2*n) # initialization of V

var_eps <- c(0.8, 0.6) # diagonal variance for epsilon


# Z1
Z1 <- matrix(0, 2*n1, 2) # initialize matrix Zi
Z1[ ,1] <- rep(c(1, 0), n1) # the first column of Zi is 1, 0, 1, 0...
Z1[ ,2] <- rep(c(0, 1), n1) # the second column of Zi is 0, 1, 0, 1...
D = matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2) # variance for bi
# Gamma_1
gamma_1 <- matrix(0, 2*n1, 2*n1)
for (j in 1:n1) {
	gamma_1[2*j - 1, ] <- rep(c(9/4, 0), n1) 
	gamma_1[2*j, ] <- rep(c(25/4, 0), n1)
}

# V1
V1 <- Z1 %*% D %*% t(Z1) + diag(rep(var_eps, n1)) + gamma_1
V[(1:(2*n1)), (1:(2*n1))] <- V1
# W1
W1 <- solve(V1)
W[(1:(2*n1)), (1:(2*n1))] <- W1


# W for subject 2:m, diagonal block matrix
for (i in 2:m) {
	ni <- max[i]
	
	# Zi
	Zi <- matrix(0, 2*ni, 2) # initialize matrix Zi
	Zi[ ,1] <- rep(c(1, 0), ni) # the first column of Zi is 1, 0, 1, 0...
	Zi[ ,2] <- rep(c(0, 1), ni) # the second column of Zi is 0, 1, 0, 1...
	
	# Gamma_i
	gamma_i <- matrix(0, 2*ni, 2*ni)
	for (j in 1:ni) {
		gamma_i[2*j - 1, ] <- rep(c(9/4, 0), ni) 
		gamma_i[2*j, ] <- rep(c(25/4, 0), ni)
	}
	
	# index for W, nrow = ncol
	nrow <- ((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])) 
	Vi <- Zi %*% D %*% t(Zi) + diag(rep(var_eps, ni)) + gamma_i
	V[nrow, nrow] <- Vi
	W[nrow, nrow] <- solve(Vi)
}


# Find ordered distinct values of mod(t_ij, 28)
day <- data$standday
uniqueDay <- sort(unique(day))
r <- length(uniqueDay)

# Incidence matrices N1 & N2
N1 <- matrix(0, 2*n, r)
N2 <- matrix(0, 2*n, r)

N <- matrix(0, n1, r) # equivalent to Ni inside of the loop
day1 <- day[1:n1]
for (j in 1:n1){
	for (l in 1:r) {
		if(day1[j] == uniqueDay[l]){
			N[j, l] <- 1
		}
	}
}
A11 <- matrix(0, 2*n1, n1)
A21 <- matrix(0, 2*n1, n1)
for (j in 1:n1){
	A11[2*j - 1, j] <- 1
	A21[2*j, j] <- 1
}
N1[1:(2*n1),] = A11 %*% N
N2[1:(2*n1),] = A21 %*% N

for (i in 2:m) {
	ni <- max[i]
	
	Ni <- matrix(0, ni, r)
	index <- (sum(max[1:(i-1)]) + 1):sum(max[1:i]) 
	dayi <- day[index]
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
	
	rindex <- ((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])) 
	N1[rindex,] <- A1i %*% Ni
	N2[rindex,] <- A2i %*% Ni	
}



# Matrix K, does not depend on theta

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
K = Q %*% ginv(R) %*% t(Q)


# INVERSE of MATRIX C

lambda1 <- 0.15
lambda2 <- 2

# The 1st row of matrix C
C11 <- t(X) %*% W %*% X
C12 <- t(X) %*% W %*% N1
C13 <- t(X) %*% W %*% N2
CRow1 = cbind(C11, C12, C13)
# The 2nd row of matrix C
C21 <- t(N1) %*% W %*% X
C22 <- t(N1) %*% W %*% N1 + lambda1 * K
C23 <- t(N1) %*% W %*% N2
CRow2 = cbind(C21, C22, C23)
# The 3nd row of matrix C
C31 <- t(N2) %*% W %*% X
C32 <- t(N2) %*% W %*% N1
C33 <- t(N2) %*% W %*% N2 + lambda2 * K
CRow3 <- cbind(C31, C32, C33)
# Inverse of coefficient matrix
C <- rbind(CRow1, CRow2, CRow3)
invC <- ginv(C)

#======================================================================================
#
# ESTIMATES - beta, f1, f2, b, U
#======================================================================================



# ESTIMATES FOR beta, f1 AND f2
temp <- rbind(t(X) %*% W %*% Y, t(N1) %*% W %*% Y, t(N2) %*% W %*% Y)
est <- invC %*% temp
betaHat <- est[1:4] 
f1Hat <- matrix(est[5:60])
f2Hat <- matrix(est[61:116])

par(mfrow = c(2, 2))
plot(log(data[1:n1, 3]), type = "l", ylab = "Log progesterone", xlab = "Days in standardized menstrual cycle")
plot(f1Hat, type = "l", xlab = "56 distinct time points in standardized menstrual cycle", ylab = "Progesterone time effect")
plot(log(data[1:n1, 4]), type = "l", ylab = "Log estrogen", xlab = "Days in standardized menstrual cycle")
plot(f2Hat, type = "l", xlab = "56 distinct time points in standardized menstrual cycle", ylab = "Estrogen time effect")




# ESTIMATES FOR b, U
b <- matrix(0, 2*m, 1)
U <- matrix(0, 2*n, 1)

# b1
ind1 <- 1:(2*n1)
Y1 <- Y[ind1]
X1 <- X[ind1, 1:4]
f11 <- N1[ind1,] %*% f1Hat
f21 <- N2[ind1,] %*% f2Hat
res1 <- Y1 - X1 %*% betaHat - f11 - f21

# b1 and U1
b[1:2] <- D %*% t(Z1) %*% W1 %*% res1
U[ind1] <- gamma_1 %*% W1 %*% res1


# bi for subject 2:m
for (i in 2:m) {
	ni <- max[i]
	
	bi <- matrix(0, 2, 1)
	
	# Zi
	Zi <- matrix(0, 2*ni, 2) # initialize matrix Zi
	Zi[ ,1] <- rep(c(1, 0), ni) # the first column of Zi is 1, 0, 1, 0...
	Zi[ ,2] <- rep(c(0, 1), ni) # the second column of Zi is 0, 1, 0, 1...
	
	# Wi
	index <- ((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])) 
	Wi <- W[index, index]	
	
	# Gamma_i
	gamma_i <- matrix(0, 2*ni, 2*ni)
	for (j in 1:ni) {
		gamma_i[2*j - 1, ] <- rep(c(9/4, 0), ni) 
		gamma_i[2*j, ] <- rep(c(0, 25/4), ni)
	}
	
	# RESIDUAL Yi - Xi\hat\beta - \hat f1i - \hat f2i
	Yi <- Y[index]
	Xi <- X[index, ]
	f1i <- N1[index,] %*% f1Hat
	f2i <- N2[index,] %*% f2Hat
	res <- Yi - Xi %*% betaHat - f1i - f2i
		
	bi <- D %*% t(Zi) %*% Wi %*% res
	b[(2*i - 1):(2*i)] <- bi
	
	Ui <- gamma_i %*% Wi %*% res
	U[index] <- Ui
}

#======================================================================================
#
# COVARIANCES - beta, f1, f2
#======================================================================================

# WEIGHT MATRICES

W1inv <- ginv(C22)
W2inv <- ginv(C33)

a <- W %*% N1
b <- a %*% W1inv
c <- b %*% t(a)
W1 <- W - c
a1 <- W %*% N2
W2 <- W - b %*% W2inv %*% t(b)


invX <- t(N2) %*% W1 %*% N2 + lambda2 * K
WxInv <- ginv(invX)
invF1 <- t(X) %*% W2 %*% X
Wf1Inv <- ginv(invF1)
invF2 <- t(X) %*% W1 %*% X
Wf2Inv <- ginv(invF2)

a2 <- W1 %*% N2
Wx <- W1 - W1 %*% N2 %*% WxInv %*% t(N2) %*% W1
Wf1 <- W2 - W2 %*% X %*% Wf1Inv %*% X %*% W2
Wf2 <- W1 - W1 %*% X %*% Wf2Inv %*% X %*% W1


# COVARIANCES - beta, f1 AND f2

invB <- t(X) %*% Wx %*% X
betaInv <- ginv(invB) 
hatB <-  betaInv %*% t(X) %*% Wx
covBeta <- hatB %*% V %*% t(hatB)

invF1 <- t(N1) %*% Wf1 %*% N1 + lambda1 * K
f1Inv <- solve(invF1)
hatF1 <- f1Inv %*% N1 %*% Wf1
covF1 <- hatF1 %*% V %*% t(hatF1)

invF2 <- t(N2) %*% Wf2 %*% N2 + lambda2 * K
f2Inv <- solve(invF1)
hatF2 <- f2Inv %*% N2 %*% Wf2
covF2 <- hatF2 %*% V %*% t(hatF2)



# biasCov <- function(lambda2){

	# # WEIGHT MATRICES
	
	# W1inv <- ginv(C22)
	# W2inv <- ginv(C33)
	
	# a1 <- W %*% N1
	# W1First <- a1 %*% W1inv
	# % W1Latter <- W1First %*% t(a1)
	
	# W1 <- W - a1 %*% W1inv %*% t(a1)
	
	
	# invX <- t(N2) %*% W1 %*% N2 + lambda2 * K
	# WxInv <- ginv(invX)
	
	
	
	# # BIASES & COVARIANCES - beta
	# invB <- t(X) %*% Wx %*% X
	# betaInv <- ginv(invB) 
	# hatB <-  betaInv %*% t(X) %*% Wx
	
	# biasBeta<- hatB %*% (N1 %*% f1 + N2 %*% f2)
	# covBeta <- hatB %*% V %*% t(hatB)
	
	
	# return(list(biasBeta, covBeta))
# }
