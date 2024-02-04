#=========================================================================
# Import dataset cycle
#=========================================================================
cycle <- read.csv("cycle.csv", header=TRUE, sep=",")
# drop the last two columns
cycle <- cycle[,1:3]

#-------------------------------------------------------------------
# Keep those subjects that have more than 10 obs'n per cycle

cycle <- cycle[which(cycle[,3] > 10),]

#-------------------------------------------------------------------
# Delete those subjects with less than 4 cycles

id <- unique(cycle[,1])
tol <- length(id) # 328 subjects

# find the number of cycles for each woman
num <- sapply(1:tol, function(v) return(nrow(cycle[which(cycle[,1] == id[v]),])))

# create a variable cycleTol documenting the number of cycles each woman has in record
cycleTol <- unlist(lapply(1:tol, function(v) return(rep(num[v], num[v]))))
cycle$cycleTol <- cycleTol

# Keep those subjs with > 3 cycles
cycle <- cycle[which(cycle[,4] > 3),]

#-------------------------------------------------------------------
# Delete those subjects that don't have consecutive cycles

idnew <- unique(cycle[,1])
tolnew <- length(idnew) # 144 subjects
# return true if the subject has concecutive cycles
concecutive <- sapply(1:tolnew, function(v) return(all(diff(cycle[which(cycle[,1] == idnew[v]), 2]) == 1)))

# create a dataset indicating whether subjects have concecutive cycles
concecutiveTemp <- cbind(idnew, concecutive)

# subject id that have at least 4 concecutive cycles
cId <-  concecutiveTemp[which(concecutiveTemp[,2] == 1), 1] # concecutiveId
tolCon <- length(cId) # 86 subjects

# Keep those that have at least 3 concecutive cycles
concecutivelist <- lapply(1:tolCon, function(v) return(subset(cycle, cycle[,1] == cId[v])))
concecutiveCycle <- do.call("rbind", concecutivelist)

#-------------------------------------------------------------------
# Keep first 3 cycles for subjects
# Subj. with 4 concecutive cycles

concecutiveCycle4 <- do.call("rbind", by(concecutiveCycle, concecutiveCycle[,1], head, n = 4))
cc4 <- concecutiveCycle4[,1:3] # concecutiveCycle4

#=========================================================================
# Add responses - import dataset after
#=========================================================================
after <- read.csv("after.csv", header=TRUE, sep=",")
after <- after[,1:5] # 31945 * 5

# Go through each row and determine if a value is zero
row_sub = apply(after, 1, function(row) all(row !=0 ))
# Subset as usual
after <- after[row_sub,] # 30084 * 5

# Merge datasets after and concecutiveCycle3 cc3
afterlist <- lapply(1:tolCon, function(v) return(after[which(after[,1] == cId[v]),][which(after[which(after[,1] == cId[v]),][,2] %in% cc4[which(cc4[,1] == cId[v]), 2]),]))
afterC <- do.call("rbind", afterlist) # 9695 * 5

#=========================================================================
# Add age and BMI into dataset afterC
#=========================================================================
ready_cycle <- read.csv("ready_cycle.csv", header=TRUE, sep=",")
ready_cycle <- ready_cycle[,c(1, 29, 31)] # 309 subjects

# merge with concecutiveId
rc <- subset(ready_cycle, ready_cycle[,1] %in% cId) # 85 subj.
rc <- rc[-which(rc[,1] == 156), 1:3] # 156 doesn't have BMI, 84 subj. 

# delete subj. that don't have age BMI data
temp <- cbind(cId, as.numeric(cId %in% ready_cycle[,1])) # 1019, 156
afterC <- subset(afterC, !(HSN4 %in% c(156, 1019))) # 9481 * 5

# find the number of cycles for each woman
fid <- rc[,1] # final id
fidlength <- length(fid)
numObs <- sapply(1:fidlength, function(v) return(nrow(afterC[which(afterC[,1] == fid[v]),])))

age <- sapply(1:fidlength, function(v) return(rep(rc[v,2], numObs[v])))
BMI <- sapply(1:fidlength, function(v) return(rep(rc[v,3], numObs[v])))

# Add age and BMI to afterC, 9592 * 7
afterC$age <- unlist(age)
afterC$BMI <- unlist(BMI)

#=========================================================================
# Final dataset
#=========================================================================

# obtain index for subj id and cycle, obtain new cc4
cc4N <- subset(cc4, !(HSN4 %in% c(156, 1019))) #  336 * 3

# cycleLength for each subj each cycle
# clength is the actual length of the cycle
# cMAx is the max of the nominal length
# e.g. a cycle is 1, 3, 5, 6, clength = 4, cMax = 6
clength <- sapply(1:dim(cc4N)[1], function(v) return(nrow(afterC[which(afterC[,1] == cc4N[v,1]),][which(afterC[which(afterC[,1] == cc4N[v,1]),][,2] == cc4N[v,2]),])))
cMax <- sapply(1:dim(cc4N)[1], function(v) return(max(afterC[which(afterC[,1] == cc4N[v,1]),][which(afterC[which(afterC[,1] == cc4N[v,1]),][,2] == cc4N[v,2]),3])))

#/---------------------------------------------
# testing for clength
# v <- 3
# a <- afterC[which(afterC[,1] == cc3N[v,1]),]
# nrow(a[which(a[,2]== cc3N[v,2]),])
# 
# max(a[which(a[,2]== cc3N[v,2]),3])
#----------------------------------------------/

# add standday
afterC$cycleMax <- unlist(sapply(1:(fidlength*4), function(v) return(rep(cMax[v], clength[v])))) 
afterC$standard <- afterC$DAY*28/afterC$cycleMax
afterC$standday <- round(afterC$standard, 1)

# add womanid
womanid <- sapply(1:fidlength, function(v) return(rep(v, numObs[v])))
afterC$womanid <- unlist(womanid)

# Export afterC for reference
write.table(afterC, "raw_4cData.txt", sep="\t", row.names=F)

# final data set that has 3 concecutive cycles with age and BMI.
# 112 subj. in total, 9592*6
# delete cycleLength and HSN4
# data <- afterC[,-c(1,3,8,9)]
# data <- data[,c(7,1,6,2:5)]
# 
# data[,2] <- unlist(sapply(1:fidlength, function(v) return(c(rep(1, clength[1+3*(v-1)]), rep(2, clength[2+3*(v-1)]), rep(3, clength[3+3*(v-1)]), rep(4, clength[4+4*(v-1)])))))

data <- afterC[,-c(1:3,8,9)]
data <- data[,c(6,5,1:4)]

# Export data to txt file
write.table(data, "data_4concecutive.txt", sep="\t", row.names=F)