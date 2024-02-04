want <- read.table("age_want.txt", header = T)
data <- read.table("ageBMI.txt", header = T)


m = 307 # number of women in the dataset (dataset w/ complete age BMI)

# find max in the new dataset after deleting rows with missing data 
max = rep(0, m)
for (i in 1:307) {
	max[i] = nrow(data[which(data$womanid == i),])
}

n = sum(max) # total observations for all subjects

# BIVARIATE RESPONSE
Y = matrix(0, 2*n)
for (i in 1:n) {
	Y[2*i - 1] = data$adjpdg2[i]
	Y[2*i] = data$adje1c2[i]
}

# COVARIATE X, indexes different from simulation, since ni are different for i
age <- want$age
X <- matrix(0, 2*n, 2) # initialization
n1 <- max[1]
X[1:(2*n1), 1] = rep(c(age[1], 0), n1) 
X[1:(2*n1), 2] = rep(c(0, age[1]), n1) 
for (i in 2:m) {
	ni = max[i]
	
	# the 1st column of 1st subject of Xi is age[1], 0, age[1], 0...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 1] = rep(c(age[i], 0), ni) 
	
	# the 2nd column of 1st subject of Xi is 0, age[1], 0, age[1]...
	X[((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])), 2] = rep(c(0, age[i]), ni) 
}



# INVERSE VARIANCE MATRIX W
W = diag(0, 2*n)
var_eps = c(0.8, 0.6) # diagonal variance for epsilon


# Z1
Z1 = matrix(0, 2*n1, 2) # initialize matrix Zi
Z1[ ,1] = rep(c(1, 0), n1) # the first column of Zi is 1, 0, 1, 0...
Z1[ ,2] = rep(c(0, 1), n1) # the second column of Zi is 0, 1, 0, 1...
D = matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2) # variance for bi
# Gamma_1
gamma_1 = matrix(0, 2*n1, 2*n1)
for (j in 1:n1) {
	gamma_1[2*j - 1, ] = rep(c(9/4, 0), n1) 
	gamma_1[2*j, ] = rep(c(25/4, 0), n1)
}
# W for subject 1
W[(1:(2*n1)), (1:(2*n1))] = solve(Z1 %*% D %*% t(Z1) + diag(rep(var_eps, n1)) + gamma_1)

# W for subject 2:m, diagonal block matrix
for (i in 2:m) {
	ni = max[i]
	
	# Zi
	Zi = matrix(0, 2*ni, 2) # initialize matrix Zi
	Zi[ ,1] = rep(c(1, 0), ni) # the first column of Zi is 1, 0, 1, 0...
	Zi[ ,2] = rep(c(0, 1), ni) # the second column of Zi is 0, 1, 0, 1...
	
	# Gamma_i
	gamma_i = matrix(0, 2*ni, 2*ni)
	for (j in 1:ni) {
		gamma_i[2*j - 1, ] = rep(c(9/4, 0), ni) 
		gamma_i[2*j, ] = rep(c(25/4, 0), ni)
	}
	
	# index for W, nrow = ncol
	nrow = ((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])) 
	W[nrow, nrow] = solve(Zi %*% D %*% t(Zi) + diag(rep(var_eps, ni)) + gamma_i)
}


# Find ordered distinct values of mod(t_ij, 28)
day = data$standday
uniqueDay = sort(unique(day))
r = length(uniqueDay)

# Incidence matrices N1 & N2
N = matrix(0, n, r)
A1 = matrix(0, 2*n, n) # diagonal block matrix
A2 = matrix(0, 2*n, n) # diagonal block matrix

for (j in 1:n) {
	for (l in 1:r) {
		if(day[j] == uniqueDay[l]){
			N[j, l] = 1
		}
	}
}


n1 = max[1]
A11 = matrix(0, 2*n1, n1) # initialize matrix A1i
A21 = matrix(0, 2*n1, n1) # initialize matrix A2i
for (j in 1:n1){
	A1i[2*j-1, j] = 1
	A2i[2*j, j] = 1 
}
A1[(1:(2*n1)), (1:n1)] = A1i
A2[(1:(2*n1)), (1:n1)] = A2i


for (i in 2:m){
	ni = max[i]
	
	A1i = matrix(0, 2*ni, ni) # initialize matrix A1i
	A2i = matrix(0, 2*ni, ni) # initialize matrix A2i
	for (j in 1:ni){
		A1i[2*j-1, j] = 1
		A2i[2*j, j] = 1 
	}
	
	rindex = ((2*sum(max[1:(i-1)])) + 1):(2*sum(max[1:i])) 
	cindex = (sum(max[1:(i-1)]) + 1):sum(max[1:i]) 
	A1[rindex, cindex] = A1i
	A2[rindex, cindex] = A2i
}

N1 = A1 %*% N
N2 = A2 %*% N



# Matrix K, does not depend on theta

# Define matrix Q, double checked.
Q = matrix(0, r, r-2)
for(j in 1:(r-2)) {
	Q[j, j] = 1/2
	Q[j + 1, j] = -1
	Q[j + 2, j] = 1/2
}

# Define matrix R 
R = matrix(0, r-2, r-2)
for (i in 1:(r-2)) {
	R[i, i] = (1/3)*4
}
for (i in 1:(r-3)) {
	R[i, i + 1] = (1/6)*2
	R[i + 1, i] = (1/6)*2
}

library(MASS) # to use ginv
K = Q %*% ginv(R) %*% t(Q)


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











