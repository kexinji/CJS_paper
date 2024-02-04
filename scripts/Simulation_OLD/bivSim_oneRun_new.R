# Bivaraite simulations

n = 30 # number of subjects
tp = 28 # 28 time points for one period;

# periodic function with period length equal to 28 days
# use the cyclic spline, alternatively; and see what I can do..
# Q1: use different periodic functions? 
t = seq(1, 28) # 28 for 1 period
f1 = 5*sin((2*pi/28)*t)
f2 = 3*sin((2*pi/28)*t)
#plot(t, f_1, type = "l")

# Simulate bivariate random intercepts. Initiate variance components. 
# sigma.b > sigma; n = 30; 
library(MASS)
D = matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2) # variance for bi
b = mvrnorm(n, c(0, 0), D) 

# noise 
var_eps = c(0.8, 0.6)
sigma = diag(var_eps) # variance for epsilon
eps = mvrnorm(n*tp, c(0, 0), sigma)
eps1 = matrix(eps[,1], n, tp)
eps2 = matrix(eps[,2], n, tp)


# OU process, assuming mean 0 and variance 1^2/2*2 = 1/4.
# now let the variance function depend on t; let xi_0, xi_ = 1; 
library(sde) 
# xi0 = 0.01
# xi1 = 0.02
# xi2 = 0.03
# sigma.ou = exp(xi0 + xi1*t + xi2*t^2)

tempU_1 = rsOU(n*tp, theta=c(0, theta12, theta13)) # theta12 = 2, theta13 = 3
tempU_2 = rsOU(n*tp, theta=c(0, theta22, theta23)) # theta22 = 2, theta23 = 5
u1 = matrix(tempU_1, n, tp)
u2 = matrix(tempU_2, n, tp)



# simulate covariate ages of the women
# the median age of women in the sample is 34 years old.
set.seed(34)
age = sample(20:44, size = n, replace = T)

X = matrix(0, 2*tp*n, 2)
for (i in 1:n) {
	X[((2*tp)*(i-1)+1):(2*tp*i), 1] = rep(c(age[i], 0), tp) # the first column of 1st subject of Xi is age[1], 0, age[1], 0...
	X[((2*tp)*(i-1)+1):(2*tp*i), 2] = rep(c(0, age[i]), tp) # the second column of 1st subject of Xi is 0, age[1], 0, age[1]...
}


id = rep(1:n, each = tp)
day = rep(t, n)

# regression coefficient.
beta1 = -0.5 # initiation of beta1


tempResponse = matrix(rep(0, 2*n*tp), 2*n, tp)
for (i in 1:n) {
	# every two rows are one subject
	tempResponse[(2*i-1),] = beta1*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
	tempResponse[2*i,] = beta1*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
}
Y = as.vector(t(tempResponse))
Y = as.matrix(round(Y))

# # plot of the response 1 and 2 - verify that the responses are indeed cyclic
# for (i in 1:n) {
	# plot(t, tempResponse[(2*i-1),], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	# par(new = TRUE)
# }
# for (i in 1:n) {
	# plot(t, tempResponse[2*i,], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	# par(new = TRUE)
# }


data = data.frame(cbind(id, day, response, age))



#===============================================================================
# PART II 
#===============================================================================



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



# the variance 
Zi = matrix(0, 2*tp, 2) # initialize matrix Zi
Zi[ ,1] = rep(c(1, 0), tp) # the first column of Zi is 1, 0, 1, 0...
Zi[ ,2] = rep(c(0, 1), tp) # the second column of Zi is 0, 1, 0, 1...

sigma_i = diag(rep(var_eps, tp))

gamma_i = matrix(0, 2*tp, 2*tp)
for (i in 1:tp) {
	gamma_i[2*i - 1, ] = rep(c(9/4, 0), tp) 
	gamma_i[2*i, ] = rep(c(25/4, 0), tp)
	}

Vi = Zi %*% D %*% t(Zi) + sigma_i + gamma_i
Wi = solve(Vi)
W = do.call("blockMatrixDiagonal", replicate(n, Wi, simplify = F))


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




# Incidence matrices N1 & N2

NTemp = diag(1, tp, tp) 
z = matrix(0, tp, tp) # matrix of zeros 
N1i = matrix(0, 2*tp, tp) 
for(i in 1:tp) {
	N1i[(2*i-1):(2*i),] = rbind(NTemp[i, ], z[i, ])
}
N2i = matrix(0, 2*tp, tp) 
for(i in 1:tp) {
	N2i[(2*i-1):(2*i),] = rbind(z[i, ], NTemp[i, ])
}
N1  = do.call("rbind", replicate(n, N1i, simplify = F))
N2  = do.call("rbind", replicate(n, N2i, simplify = F))



# INVERSE of MATRIX C

lambda1 = 0.15
lambda2 = 2

# The 1st row of matrix C
C11 = t(X) %*% W %*% X
C12 = t(X) %*% W %*% N1
C13 = t(X) %*% W %*% N2
CRow1 = cbind(C11, C12, C13)
# The 2nd row of matrix C
C21 = t(N1) %*% W %*% X
C22 = t(N1) %*% W %*% N1 + lambda1 * K
C23 = t(N1) %*% W %*% N2
CRow2 = cbind(C21, C22, C23)
# The 2nd row of matrix C
C31 = t(N2) %*% W %*% X
C32 = t(N2) %*% W %*% N1
C33 = t(N2) %*% W %*% N2 + lambda2 * K
CRow3 = cbind(C31, C32, C33)
# Inverse of coefficient matrix
C = rbind(CRow1, CRow2, CRow3)
invC = ginv(C)

# ESTIMATES FOR beta, f1 AND f2
temp = rbind(t(X) %*% W %*% Y, t(N1) %*% W %*% Y, t(N2) %*% W %*% Y)
est = invC %*% temp
betaHat = est[1:2] 
f1Hat = matrix(est[3:30])
f2Hat = matrix(est[31:58])


