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
# Delete those subjects with less than 3 cycles

id <- unique(cycle[,1])
tol <- length(id) # 328 subjects

# find the number of cycles for each woman
num <- sapply(1:tol, function(v) return(nrow(cycle[which(cycle[,1] == id[v]),])))

# create a variable cycleTol documenting the number of cycles each woman has in record
cycleTol <- unlist(lapply(1:tol, function(v) return(rep(num[v], num[v]))))
cycle$cycleTol <- cycleTol

# Keep those subjs with > 2 cycles
cycle <- cycle[which(cycle[,4] > 2),]

#-------------------------------------------------------------------
# Delete those subjects that don't have consecutive cycles

idnew <- unique(cycle[,1])
tolnew <- length(idnew) # 204 subjects
# return true if the subject has concecutive cycles
concecutive <- sapply(1:tolnew, function(v) return(all(diff(cycle[which(cycle[,1] == idnew[v]), 2]) == 1)))

# create a dataset indicating whether subjects have concecutive cycles
concecutiveTemp <- cbind(idnew, concecutive)

# subject id that have at least 3 concecutive cycles, 116 subjects
cId <-  concecutiveTemp[which(concecutiveTemp[,2] == 1), 1] # concecutiveId
tolCon <- length(cId)

# Keep those that have at least 3 concecutive cycles
concecutivelist <- lapply(1:tolCon, function(v) return(subset(cycle, cycle[,1] == cId[v])))
concecutiveCycle <- do.call("rbind", concecutivelist)

#-------------------------------------------------------------------
# Keep first 3 cycles for subjects
# Subj. with 3 concecutive cycles

concecutiveCycle3 <- do.call("rbind", by(concecutiveCycle, concecutiveCycle[,1], head, n = 3))
cc3 <- concecutiveCycle3[,1:3] # concecutiveCycle3

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
afterlist <- lapply(1:tolCon, function(v) return(after[which(after[,1] == cId[v]),][which(after[which(after[,1] == cId[v]),][,2] %in% cc3[which(cc3[,1] == cId[v]), 2]),]))
afterC <- do.call("rbind", afterlist) # 9921 * 5

#=========================================================================
# Add age and BMI into dataset afterC
#=========================================================================
ready_cycle <- read.csv("ready_cycle.csv", header=TRUE, sep=",")
ready_cycle <- ready_cycle[,c(1, 29, 31)] # 309 subjects

# merge with concecutiveId
rc <- subset(ready_cycle, ready_cycle[,1] %in% cId) # 113 subj.
rc <- rc[-which(rc[,1] == 156), 1:3] # 112 subj. 

# delete subj. that don't have age BMI data
temp <- cbind(cId, as.numeric(cId %in% ready_cycle[,1])) # 118, 335, 1019, 156
afterC <- subset(afterC, !(HSN4 %in% c(118, 156, 335, 1019))) # 9592 * 5

# find the number of cycles for each woman
fid <- rc[,1] # final id
fidlength <- length(fid)
numObs <- sapply(1:fidlength, function(v) return(nrow(afterC[which(afterC[,1] == fid[v]),])))

age <- sapply(1:fidlength, function(v) return(rep(rc[v,2], numObs[v])))
BMI <- sapply(1:fidlength, function(v) return(rep(rc[v,3], numObs[v])))

# Add age and BMI to afterC, 9592 * 7, 112 subjects
afterC$age <- unlist(age)
afterC$BMI <- unlist(BMI)

#=========================================================================
# Delete missing data
#=========================================================================

# subset afterC, those with 4 concecutive cycles
# fid_4c <- c(2, 4, 6, 9, 10, 11, 12, 14, 15, 16, 17, 20, 22, 25, 28, 30, 32, 34, 38, 39, 47, 60, 69, 81, 82, 89, 91, 92, 93, 99, 
#             103, 107, 108, 123, 126, 128, 134, 138, 143, 144, 150, 166, 168, 170, 171, 174, 192, 193, 194, 198, 
#             200, 201, 202, 210, 213, 217, 225, 234, 239, 251, 256, 274, 275, 284, 286, 287, 294, 
#             302, 316, 1004, 1010, 1022, 1039, 1040, 1053, 1057, 1066, 1071, 1074, 1083, 1094, 1098, 1104, 1127)
# 
# afterC <- subset(afterC, afterC[,1] %in% fid_4c) # 7112 * 10

# obtain index for subj id and cycle, obtain new cc3
cc3N <- subset(cc3, !(HSN4 %in% c(118, 156, 335, 1019))) # 9592 * 5

# cc3N <- subset(cc3, (HSN4 %in% fid_4c)) # 252*3

# cycleLength for each subj each cycle
# clength is the actual length of the cycle
# cMAx is the max of the nominal length
# e.g. a cycle is 1, 3, 5, 6, clength = 4, cMax = 6
clength <- sapply(1:dim(cc3N)[1], function(v) return(nrow(afterC[which(afterC[,1] == cc3N[v,1]),][which(afterC[which(afterC[,1] == cc3N[v,1]),][,2] == cc3N[v,2]),])))
cMax <- sapply(1:dim(cc3N)[1], function(v) return(max(afterC[which(afterC[,1] == cc3N[v,1]),][which(afterC[which(afterC[,1] == cc3N[v,1]),][,2] == cc3N[v,2]),3])))

# range of cycle lengths for each woman
cRange <- sapply(0:111, function(i) return(max(cMax[(1+3*i):(3+3*i)]) - min(cMax[(1+3*i):(3+3*i)])))

length(cRange[which(cRange == 1)]) # 7/112
length(cRange[which(cRange <= 5 & cRange > 1)]) # 66/112
length(cRange[which(cRange <= 10 & cRange > 5)]) # 21/112
length(cRange[which(cRange <= 20 & cRange > 10)]) # 106/112


hist(cRange, main = "Within-Woman Between-Cycle Range of Days", xlab = "within-woman between-cycle range (days)")

# missing data
temp <- 1-clength/cMax  # 112*3 = 336
hist(temp, main = "1 - num obs'n/max \n for each woman for each cycle", xlab = "1-clength/cMax")
length(temp[which(temp == 0)]) # 126
length(temp[which(temp <= 0.1)]) # 285


# Only keep those with less than or equal to 10 cRange, down to 94
cc3N$range <- rep(cRange, each = 3)
cc3N$missing <- temp
cc3N <- cc3N[cc3N[,4] <= 10,] 


#/---------------------------------------------
# testing for clength
# v <- 3
# a <- afterC[which(afterC[,1] == cc3N[v,1]),]
# nrow(a[which(a[,2]== cc3N[v,2]),])
# 
# max(a[which(a[,2]== cc3N[v,2]),3])
#----------------------------------------------/
#=========================================================================
# Final dataset
#=========================================================================

#-----------------------------------------------------------
# Subset afterC according to the final woman id
fid_4c <- unique(cc3N[,1]) # final id
afterC <- subset(afterC, afterC[,1] %in% fid_4c) 


#/-----------------------------------------------------------
# add standday

# new cMax clength, it's after the trimming for 94 women
cMax <- sapply(1:dim(cc3N)[1], function(v) return(max(afterC[which(afterC[,1] == cc3N[v,1]),][which(afterC[which(afterC[,1] == cc3N[v,1]),][,2] == cc3N[v,2]),3])))
clength <- sapply(1:dim(cc3N)[1], function(v) return(nrow(afterC[which(afterC[,1] == cc3N[v,1]),][which(afterC[which(afterC[,1] == cc3N[v,1]),][,2] == cc3N[v,2]),])))

afterC$cycleMax <- unlist(sapply(1:(length(fid_4c)*3), function(v) return(rep(cMax[v], clength[v])))) 
afterC$standard <- afterC$DAY*28/afterC$cycleMax
afterC$standday <- round(afterC$standard)

# update standday
# afterC$day <- afterC$standday + 28*(afterC$CYCLE-1)

#-----------------------------------------------------------
# change cycle index 
afterC[,2] <- unlist(sapply(1:length(fid_4c), function(v) return(c(rep(1, clength[1+3*(v-1)]), rep(2, clength[2+3*(v-1)]), rep(3, clength[3+3*(v-1)])))))


#-----------------------------------------------------------
# add womanid
numObs <- sapply(1:length(fid_4c), function(v) return(nrow(afterC[which(afterC[,1] == fid_4c[v]),])))
womanid <- sapply(1:length(fid_4c), function(v) return(rep(v, numObs[v])))
afterC$womanid <- unlist(womanid)

# final data set that has 3 concecutive cycles with age and BMI.
# 112 subj. in total, 9592*6
# delete cycleLength and HSN4

# data <- afterC[,-c(1:3,8,9)]
# data <- data[,c(6,5,1:4)]

# data_3c <- afterC[,c(12,2,11,4:7)]

#-----------------------------------------------------------
# final dataset
data_3c <- subset(afterC, select = c(womanid, CYCLE, standday, adjpdg2, adje1c2, age, BMI, day))


write.table(data_3c, "data_3cycle.txt", sep="\t", row.names=F)

#=========================================================================
# data with 3 concecutive cycles

# # add womanid
# womanid <- sapply(1:fidlength, function(v) return(rep(v, numObs[v])))
# afterC$womanid <- unlist(womanid)
# 
# # Export afterC for reference
# write.table(afterC, "raw_3cData.txt", sep="\t", row.names=F)
# 
# # final data set that has 3 concecutive cycles with age and BMI.
# # 112 subj. in total, 9592*6
# # delete cycleLength and HSN4
# data <- afterC[,-c(1,3,8,9)]
# data <- data[,c(7,1,6,2:5)]
# 
# data[,2] <- unlist(sapply(1:fidlength, function(v) return(c(rep(1, clength[1+3*(v-1)]), rep(2, clength[2+3*(v-1)]), rep(3, clength[3+3*(v-1)])))))
# 
# # Export data to txt file
# write.table(data, "data_3concecutive.txt", sep="\t", row.names=F)


