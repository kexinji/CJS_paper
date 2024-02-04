## THIS SCIPT PROCESS LIU DATA IN "ageBMI.txt" TO PREPARE FOR BIVARIATE ANALYSIS FOR REAL DATASETS FOR A SINGLE CYCLE.
## IT OUTPUTS DATA READY FOR MODELLING.

#============================================================================================================
preProcessLiuData <- function(nSubj) {
  # Process data from the ageBMI.txt file 
  #
  # Args: 
  #   data.path:  the directory containing the ageBMI.txt file 
  #   nSubj:      number of randomly selected subjects 
  #
  # Returns:
  # 	Output the processed data  
  
  
  dataComp <- read.table("ageBMI.txt", header = T)
  
  #------------------------------------------------------------------------------
  # create two indicator variables
  
  dataComp$underWeight <- 0
  dataComp$underWeight[which(dataComp$BMI == 1)] <- 1
  
  dataComp$overWeight <- 0
  dataComp$overWeight[which(dataComp$BMI == 3)] <- 1
  
  dataComp$BMI <- NULL
  
  #------------------------------------------------------------------------------
  # Randomly select m subjects with max 28 in standday
  id <- dataComp[,1]
  m <- max(id)
  niMax <- sapply(1:m, function(v) return(max(dataComp[which(id == v), 2])))
  niMax <- cbind(matrix(1:m, m), as.matrix(niMax))
  
  # subject id with 28 standday
  niMax28 <- niMax[which(niMax[,2] == 28)]
  
  set.seed(nSubj)
  s <- sample(niMax28, nSubj)
  
  # new dataset
  data <- lapply(1:nSubj, function(v) return(rbind(dataComp[which(id == s[v]),])))
  data <- do.call(rbind, data)
  
  # Centre the covariate - age
  niVec <- sapply(1:nSubj, function(v) return(nrow(data[which(dataComp[,1] == s[v]),])))
  age_index <- sapply(1:(nSubj-1), function(v) return(1+sum(niVec[1:v])))
  
  # Centre the covariates
  data[,2] <- (data[,2] - 14)/10
  
  medianAge <- median(data[c(1, age_index), 5]) 
  data[,5] <- (data[,5] - medianAge)/100
  
  return(data)
}