#------------------------------------------------------------------------------
# create two indicator variables

# dataComp <- read.table("data_3concecutive.txt", header = T)

dataComp <- read.table("data_3cycle.txt", header = T)

dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# m = 10, randomly select 10 subjects with 4 concecutive cycles, but only use 3 cycles
id <- dataComp[,1]
m <- 40

# seed 42 doesn't work for 30 subjects
# Now change it to 60
set.seed(60)
samp <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == samp[v]),])))
data <- do.call(rbind, data)

# Centre the covariate - age
niVec <- sapply(1:m, function(v) return(nrow(data[which(dataComp[,1] == samp[v]),])))
age_index <- sapply(1:(m-1), function(v) return(1+sum(niVec[1:v])))
 
medianAge <- median(data[c(1, age_index), 6]) 
data$age <- (data$age - medianAge)/100


#------------------------------------------------------------------------------
# Run 
source("bivFns_3c.R")

system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age, data$underWeight, data$overWeight), random = cbind(1, 1), process = "NOU", time = data$standday, id = data$womanid, tol = 0.001, cap = 10))

# test for equality, with only covariate age
# system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age), random = cbind(1, 1), process = "NOU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))

write.csv(unlist(res), "Liu3c_40subj_.txt")
