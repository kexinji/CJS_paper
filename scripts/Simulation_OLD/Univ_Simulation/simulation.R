# Zhang et al 2000

n = 30 # number of subjects
tp = 30 # 15 time points for one period; 30 time points for two periods

# periodic function with period length equal to 30 days
# use the cyclic spline, alternatively; and see what I can do..
t = seq(1, 60, by = 2) # 60 for two periods
f = 5*sin((2*pi/30)*t)
# plot(t, f, type = "l")

# independent random intercepts; sigma.b > sigma; n = 30; random time effect corresponding to different subject; 
sigma.b = 1.5 # variance for bi
bi = rnorm(n, 0, sigma.b) 

# noise 
sigma = 0.8 # variance for epsilon
eps.vec = rnorm(n*tp, 0, sigma)
epsilon = matrix(eps.vec, n, tp)

# OU process, assuming mean 0 and variance 1^2/2*2 = 1/4.
# now let the variance function depend on t; let xi_0, xi_ = 1; 
library(sde) 
xi0 = 0.01
xi1 = 0.02
xi2 = 0.03
sigma.ou = exp(xi0 + xi1*t + xi2*t^2)
tempU = rsOU(n*tp, theta=c(0,2,sigma.ou))
u = matrix(tempU, n, tp)


# simulate covariate ages of the women
# the average age of menopause is 51 years old.
set.seed(51)
ageTemp = sample(35:60, size = n, replace = T)
age = rep(ageTemp, each = tp)

id = rep(1:n, each = tp)
day = rep(t, n)

# Will not take log here; as the data are already Normally simulated.
beta0 = 27.05 # initiation of beta0
beta1 = -0.5 # initiation of beta1
tempResponse = matrix(rep(0, n*tp), n, tp)
for(i in 1:n) {
	# each row is one subject
	tempResponse[i,] = beta0 + beta1*ageTemp[i] + f + bi[i] + u[i,] + epsilon[i,]
}
response = as.vector(t(tempResponse))
range(response)
response = round(response)

# plot of the response - verify that the response is indeed cyclic
for (i in 1:n) {
	plot(t, tempResponse[i,], type = "l", xlab = "day", ylim = c(-10, 20), ylab = "")
	par(new = TRUE)
}


# beta0 = 37.05 # initiation of beta0
# beta1 = -0.5 # initiation of beta1
# tempResponse = matrix(rep(0, n*tp), n, tp)
# for(i in 1:n) {
	# # each row is one subject
	# tempResponse[i,] = beta0 + beta1*ageTemp[i] + f + bi[i] + u[i,] + epsilon[i,]
# }
# range(tempResponse)

# # plot of the log(response) - verify that the log of the response is approximately normal
# for (i in 1:n) {
	# plot(t, log(tempResponse[i,]), type = "l", xlab = "day", ylim = c(-3, 3), ylab = "")
	# par(new = TRUE)
# }
#  
# response = as.vector(t(tempResponse))
# range(response)
# response = round(response)





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


