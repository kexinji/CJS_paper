# Use this script for CJS revision
# Run real data analysis with OU 


#------------------------------------------------------------------------------
# Run 
source("bivFns_lik.R")
source("LiuData_preprocess_fn.R")


data <- preProcessLiuData(50)
system.time(res <- results(data, 
                           response = cbind(log(data$adjpdg2), log(data$adje1c2)), 
                           fixed = cbind(data$age, data$underWeight, data$overWeight), 
                           random = cbind(1, 1), 
                           process = "OU", 
                           time = data$standday, 
                           id = data$womanid, 
                           tol = 0.001, 
                           cap = 50))

# test for equality, with only covariate age
# system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age), random = cbind(1, 1), process = "OU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))


write.csv(unlist(res), "LiuBiv50_OU.txt")
