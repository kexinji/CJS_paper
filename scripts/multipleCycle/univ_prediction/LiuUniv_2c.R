#------------------------------------------------------------------------------
# Import data with with 4 concecutive cycles, but only use 2 cycles

dataComp <- read.table("data_2cycle.txt", header = T)

# dataComp <- subset(dataComp, select = c(womanid, standday, adjpdg2, age, BMI))

# create two indicator variables
dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# m = 10, randomly select 10 subjects from the dataset
id <- dataComp[,1]
m <- 10

set.seed(42)
s <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data)

# Centre the covariates, standday and age
data[,2] <- (data[,2] %% 28)/10

medianAge <- median(data[,5][!duplicated(data[,c(1, 5)])])
data[,5] <- (data[,5] - medianAge)/100

niVec <- sapply(1:m, function(v) return(nrow(data[which(id == s[v]),])))

#------------------------------------------------------------------------------
# Run 
source("univFns_2c.R")

system.time(res <- results(data, response = log(data$adjpdg2), fixed = cbind(data$age, data$underWeight, data$overWeight), random = 1, process = "NOU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))

# test for equality, with only covariate age
# system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age), random = cbind(1, 1), process = "NOU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))

write.csv(unlist(res), "LiuUniv_2c.txt")
