setwd("C:/Users/Home/Dropbox/CJS_revision/scripts/singleCycle_LiuData")

#------------------------------------------------------------------------------
# create two indicator variables

dataComp <- read.table("ageBMI.txt", header = T)

dataComp$underWeight <- 0
dataComp$underWeight[which(dataComp$BMI == 1)] <- 1

dataComp$overWeight <- 0
dataComp$overWeight[which(dataComp$BMI == 3)] <- 1

dataComp$BMI <- NULL

#------------------------------------------------------------------------------
# m = 50, randomly select 50 subjects with max 28 in standday
id <- dataComp[,1]
m <- max(id)
niMax <- sapply(1:m, function(v) return(max(dataComp[which(id == v), 2])))
niMax <- cbind(matrix(1:m, m), as.matrix(niMax))

# subject id with 28 standday
niMax28 <- niMax[which(niMax[,2] == 28)]

nSubj <- 100

set.seed(nSubj)
s <- sample(niMax28, nSubj)

data <- lapply(1:nSubj, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data)

# Centre the covariate - age
niVec <- sapply(1:nSubj, function(v) return(nrow(data[which(dataComp[,1] == s[v]),])))
age_index <- sapply(1:(nSubj-1), function(v) return(1+sum(niVec[1:v])))

# Centre the covariates
data[,2] <- (data[,2] - 14)/10

medianAge <- median(data[c(1, age_index), 5]) 
data[,5] <- (data[,5] - medianAge)/100


#------------------------------------------------------------------------------
# Run 
source("bivLiu_sc_fn.R")

system.time(res <- results(data, 
                           response = cbind(log(data$adjpdg2), log(data$adje1c2)), 
                           fixed = cbind(data$age, data$underWeight, data$overWeight), 
                           random = cbind(1, 1), 
                           process = "NOU", 
                           time = data$standday, 
                           id = data$womanid, 
                           tol = 0.001, 
                           cap = 50))

# test for equality, with only covariate age
# system.time(res <- results(data, response = cbind(log(data$adjpdg2), log(data$adje1c2)), fixed = cbind(data$age), random = cbind(1, 1), process = "NOU", time = data$standday, id = data$womanid, tol = 0.001, cap = 50))

write.csv(unlist(res), "bivLiu100_sc_cor.txt")
