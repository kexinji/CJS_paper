## THIS FILE IS USED to TEST univFns.r

source("univFns.R")

#------------------------------------------------------------------------------
# SIMULATION TEST
# theta0 <- c(tau, sigma.b, sigma, theta2, theta3)
theta <- c(rep(1, 3), -log(0.8), sqrt(-log(0.8)*2))
theta <- c(rep(1, 3), -log(0.8), rep(0, 3))
simdata <- newData(theta, -0.5, 30, 28, process = "NOU")

resFinal <- results(simdata, simdata$Y, simdata$age, random = 1, process = "NOU", simdata$day, simdata$id,  tol = 0.1, cap = 5)


	
for(i in 1:100){
	print(paste(i, "th iteration"))
	theta <- c(rep(1, 3), -log(0.8), rep(0, 3))
	set.seed(i)
	simdata <- newData(theta, -0.5, 30, 28, process = "NOU")
	res <- results(simdata, simdata$Y, simdata$age, random = 1, process = "NOU", simdata$day, simdata$id,  tol = 0.001, cap = 100)
	write.csv(unlist(res), paste('t001c100_NOUunivGen_' , i, '.txt', sep=''))
}



#------------------------------------------------------------------------------
# LIU DATASET TEST

liudata34 <- read.table("ageBMI.txt", header = T)[1:992,]
resLiu34 <- results(data = liudata34, response = liudata34$adjpdg2, fixed = cbind(liudata34$age, liudata34$BMI), random = 1, process = "OU", time = liudata34$standday, id = liudata34$womanid, tol = 0.1, cap = 10)


# want <- read.table("ageBMI_want.txt", header = T)
# data <- read.table("ageBMI.txt", header = T)

# m = 34
# want <- read.table("ageBMI_want.txt", header = T)[1:34,]
# liudata34 <- read.table("ageBMI.txt", header = T)[1:992,]

# m = 50
# want <- read.table("ageBMI_want.txt", header = T)[1:50,]
# data <- read.table("ageBMI.txt", header = T)[1:1463,]

# # m = 100
# want <- read.table("ageBMI_want.txt", header = T)[1:100,]
# data <- read.table("ageBMI.txt", header = T)[1:2810,]

