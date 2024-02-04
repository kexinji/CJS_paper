# INPUT ONLY FIRST 34 SUBJ. TO REDUCE DIMENSION, r is still 56.
# want <- read.table("ageBMI_want.txt", header = T)
# data <- read.table("ageBMI.txt", header = T)

# m = 34
# want <- read.table("ageBMI_want.txt", header = T)[1:34,]
# data <- read.table("ageBMI.txt", header = T)[1:992,]

# m = 50
# want <- read.table("ageBMI_want.txt", header = T)[1:50,]
# data <- read.table("ageBMI.txt", header = T)[1:1463,]

# # m = 100
# want <- read.table("ageBMI_want.txt", header = T)[1:100,]
# data <- read.table("ageBMI.txt", header = T)[1:2810,]

source("biv_fns.R")

# m <- 34 # TO REDUCE THE DIMENSION OF THE DATASET
# m <- 50
m <- 100

# find max in the new dataset after deleting rows with missing data 
max <- rep(0, m)
for (i in 1:m) {
	max[i] = nrow(data[which(data$womanid == i),])
}

# m = 307 is the total number of women in the dataset (dataset w/ complete age BMI)
param <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)
data <- liuData(m, niVec=max)
Y <- data$Y
X <- data$X
N1 <- data$N1
N2 <- data$N2
K <- matrixK(h=1/2, tp=56)


res <- results(param, Y, tol = rep(0.005, 11), cap = 150, m, tp=28, K, X, N1, N2, niVec=max)
write.csv(unlist(res), "liuBiv100_t005c150.txt", row.names = F)	


# # FOR TESTING

#res1 <- est(param, Y, m=307, tp=28, K, X, N1, N2, niVec=max)
# system.time(theta1 <- findPar(param, Y, m=307, tp=28, K, X, N1, N2, niVec=max))

# > system.time(res1 <- est(param, Y, m=307, tp=28, K, X, N1, N2, niVec=max))
   # user  system elapsed 
 # 12.792   1.102  13.811 

# m = 34, tol = 0.1, cap = 50, 5 iterations
# > system.time(source("Liu_biv.R"))
   # user  system elapsed 
# 553.714   0.961 555.780 

# m = 34, tol = 0.01, cap = 50, 7 iterations
# > system.time(source("Liu_biv.R"))
    # user   system  elapsed 
# 1849.639    3.875 3342.786 


# m = 50, tol = 0.01, cap = 50, 14 iterations
# > system.time(source("Liu_biv.R"))
     # user    system   elapsed 
 # 9155.366   125.489 15839.732 
 

# m = 100, tol = 0.01, cap = 50, 27 iterations
# > system.time(source("Liu_biv.R"))
     # user     system    elapsed 
 # 92643.422   1374.549 131636.753 
