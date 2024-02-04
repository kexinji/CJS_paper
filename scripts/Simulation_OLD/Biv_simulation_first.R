# Bivaraite simulations

n = 30 # number of subjects
tp = 30 # 15 time points for one period; 30 time points for two periods

# periodic function with period length equal to 30 days
# use the cyclic spline, alternatively; and see what I can do..
# Q1: use different periodic functions? 
t = seq(1, 60, by = 2) # 60 for two periods
f1 = 5*sin((2*pi/30)*t)
f2 = 3*sin((2*pi/30)*t)
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

tempU_1 = rsOU(n*tp, theta=c(0,2,3)) # theta2 = 2, theta3 = 3
tempU_2 = rsOU(n*tp, theta=c(0,2,5)) # theta2 = 2, theta3 = 3
u1 = matrix(tempU_1, n, tp)
u2 = matrix(tempU_2, n, tp)



# simulate covariate ages of the women
# the average age of menopause is 51 years old.
set.seed(51)
age = sample(35:60, size = n, replace = T)
#age = rep(age, each = tp)

id = rep(1:n, each = tp)
day = rep(t, n)

# Will not take log here; as the data are already Normally simulated.
beta10 = 27.05 # initiation of beta10
beta20 = 20.05 # initiation of beta20

beta1 = -0.5 # initiation of beta1
tempResponse = matrix(rep(0, 2*n*tp), 2*n, tp)
for(i in 1:n) {
	# every two rows are one subject
	tempResponse[(2*i-1),] = beta10 + beta1*age[i] + f1 + b[i, 1] + u1[i,] + eps1[i,]
	tempResponse[2*i,] = beta20 + beta1*age[i] + f2 + b[i, 2] + u2[i,] + eps2[i,]
}
# response = as.vector(t(tempResponse))
# range(response)
# response = round(response)

# plot of the response 1 and 2 - verify that the responses are indeed cyclic
for (i in 1:n) {
	plot(t, tempResponse[(2*i-1),], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	par(new = TRUE)
}
for (i in 1:n) {
	plot(t, tempResponse[2*i,], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	par(new = TRUE)
}







# the variance 
Zi = matrix(0, 2*tp, 2) # initiate matrix Zi
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






data = data.frame(cbind(id, day, response, age))


export dataset subject into a txt file
write.table(data, "mydata.txt", quote = F, sep = " ", col.names = F)





####################################################################################################
# Fit the data using cyclic splines
####################################################################################################

# This dataset is the same as the one running using SAS. 
data.sas = read.table("mydata.txt", header = F)
data.sas$obs = NULL
head(data.sas)

# cyclic spline	
library(mgcv)
cs = gam(response ~ age + s(day, bs = "cc", k = 20), data = data.sas)




# an R package found online
library(gss)












sigma=matrix(c(
    #y            
    3 ,.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    .5, 3,.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 ,.5, 3,.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 , 0,.5, 3,.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 , 0, 0,.5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 ,0 ,0 ,0 , 0, 3,.5, 0, 0, 0, 0, 0, 0, 0, 0,
    0 ,0 ,0 ,0 ,0 ,.5, 3,.5, 0, 0, 0, 0, 0, 0, 0,
    0 ,0 ,0 ,0 ,0 ,0 ,.5, 3,.5, 0, 0, 0, 0, 0, 0,
    0 ,0 ,0 ,0 ,0 ,0 ,0 ,.5, 3,.5, 0, 0, 0, 0, 0,
    0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,.5, 3, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0 ,0 ,0 ,0 ,0 , 0, 3,.5, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0 ,.5, 3,.5, 0, 0,
    0 , 0, 0,0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,.5, 3,.5, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .5, 3,.5,
    0 ,0 ,0 , 0, 0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,.5, 3

    ),15,15)

mean = rep(0, 15)

library(mvtnorm)
sample = rmvnorm(1000, mean, sigma)
sample = as.matrix(sample)

plot(sample[1, ]+f, type = "l")
plot(sample[2, ]+f, type = "l")


