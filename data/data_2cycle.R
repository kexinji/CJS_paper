#=============================================================================
# Import dataset cycle
#=============================================================================
cycle <- read.csv("cycle.csv", header=TRUE, sep=",") # 341 subjects

# drop the last two columns
cycle <- cycle[,1:3]
# Keep those subjects that have more than 10 obs'n per cycle
cycle <- cycle[which(cycle[,3] > 10),]

#-----------------------------------------------------------------------
# Subset dataset cycle and keep those subjects with at least 4 cycles

id <- unique(cycle[,1])
tol <- length(id) # 328 subjects

# find the number of cycles for each woman
num <- sapply(1:tol, function(v) return(nrow(cycle[which(cycle[,1] == id[v]),])))

# create a variable cycleTol documenting the number of cycles each woman has in record
cycleTol <- unlist(lapply(1:tol, function(v) return(rep(num[v], num[v]))))
cycle$cycleTol <- cycleTol

# Keep those subjs with > 3 cycles
cycle <- cycle[which(cycle[,4] > 3),] # want to use the same id with dataset with more than 4 concecutive cycles
# cycle <- cycle[which(cycle[,4] > 1),] 

#-----------------------------------------------------------------------
# Delete those subjects that don't have consecutive cycles

idnew <- unique(cycle[,1])
tolnew <- length(idnew) # 272 subjects 2 cycles, 144 subjects 4 cycles

# return true if the subject has concecutive cycles
concecutive <- sapply(1:tolnew, function(v) return(all(diff(cycle[which(cycle[,1] == idnew[v]), 2]) == 1)))
# create a dataset indicating whether subjects have concecutive cycles
concecutiveTemp <- cbind(idnew, concecutive)

# Extract subject id that have at least 4 concecutive cycles
cId <-  concecutiveTemp[which(concecutiveTemp[,2] == 1), 1] # concecutiveId
tolCon <- length(cId) # 160 subjects 2 cycles, 86 subjects 4 cycles

#-----------------------------------------------------------------------
# Create dataset concecutiveCycle that have at least 4 concecutive cycles
concecutivelist <- lapply(1:tolCon, function(v) return(subset(cycle, cycle[,1] == cId[v])))
concecutiveCycle <- do.call("rbind", concecutivelist)

# Keep first 2 cycles for subjects
concecutiveCycle2 <- do.call("rbind", by(concecutiveCycle, concecutiveCycle[,1], head, n = 2))
cc2 <- concecutiveCycle2[,1:3] # concecutiveCycle2

#=============================================================================
# Add responses - import dataset after
#=============================================================================
after <- read.csv("after.csv", header=TRUE, sep=",")
after <- after[,1:5] # 31945 * 5

# Go through each row and determine if a value is zero
row_sub = apply(after, 1, function(row) all(row !=0 ))
# Subset as usual
after <- after[row_sub,] # 30084 * 5

# Merge datasets after and concecutiveCycle2 cc2
afterlist <- lapply(1:tolCon, function(v) return(after[which(after[,1] == cId[v]),][which(after[which(after[,1] == cId[v]),][,2] %in% cc2[which(cc2[,1] == cId[v]), 2]),]))
afterC <- do.call("rbind", afterlist) # 4831 * 5

#=============================================================================
# Add age and BMI into dataset afterC
#=============================================================================
ready_cycle <- read.csv("ready_cycle.csv", header=TRUE, sep=",")
ready_cycle <- ready_cycle[,c(1, 29, 31)] # 309 subjects

# Go through each row and determine if a value is zero
row_sub_ageBMI = apply(ready_cycle, 1, function(row) all(row !=-99 ))
# Subset as usual
ready_cycle <- ready_cycle[row_sub_ageBMI,] # 307 subjects

# Subset ready_cycle according to concecutiveId
rc <- subset(ready_cycle, ready_cycle[,1] %in% cId) # 84 subj.
fid <- rc[,1] # final id
fidlength <- length(fid)

# Subset dataset afterC and delete those subjects that don't have age BMI 
afterC <- subset(afterC, HSN4 %in% fid) # 4722 * 5

# find the number of observations for each woman
numObs <- sapply(1:fidlength, function(v) return(nrow(afterC[which(afterC[,1] == fid[v]),])))

age <- sapply(1:fidlength, function(v) return(rep(rc[v,2], numObs[v])))
BMI <- sapply(1:fidlength, function(v) return(rep(rc[v,3], numObs[v])))
womanid <- sapply(1:fidlength, function(v) return(rep(v, numObs[v])))

# Add age, BMI and womanid to afterC, 4722 * 7
afterC$age <- unlist(age)
afterC$BMI <- unlist(BMI)
afterC$womanid <- unlist(womanid)

#=============================================================================
# Add standardized days to dataset afterC
#=============================================================================
# Update dataset cc2 by the final id
cc2 <- subset(cc2, HSN4 %in% fid) # 4722 * 5

# cycleLength for each subj each cycle
# clength is the actual length of the cycle; cMAx is the max of the nominal length
# e.g. a cycle is 1, 3, 5, 6, clength = 4, cMax = 6
clength <- sapply(1:dim(cc2)[1], function(v) return(nrow(afterC[which(afterC[,1] == cc2[v,1]),][which(afterC[which(afterC[,1] == cc2[v,1]),][,2] == cc2[v,2]),])))
cMax <- sapply(1:dim(cc2)[1], function(v) return(max(afterC[which(afterC[,1] == cc2[v,1]),][which(afterC[which(afterC[,1] == cc2[v,1]),][,2] == cc2[v,2]),3])))

#/---------------------------------------------
# testing for clength
# v <- 3
# a <- afterC[which(afterC[,1] == cc2[v,1]),]
# nrow(a[which(a[,2]== cc2[v,2]),])
# 
# max(a[which(a[,2]== cc2[v,2]),3])
#----------------------------------------------/

# change cycle index 
afterC[,2] <- unlist(sapply(1:fidlength, function(v) return(c(rep(1, clength[1+2*(v-1)]), rep(2, clength[2+2*(v-1)])))))

# add standday
afterC$cycleMax <- unlist(sapply(1:(fidlength*2), function(v) return(rep(cMax[v], clength[v])))) 
afterC$standard <- afterC$DAY*28/afterC$cycleMax
afterC$standday <- round(afterC$standard)

# update standday
# afterC$standday <- afterC$standardR + 28*(afterC$CYCLE-1)


# Final dataset that has 4 concecutive cycles with age and BMI,
# outputing only 2 of the 4 cycles, with 84 subjects in total, 4722*6.

data <- subset(afterC, select = c(womanid, CYCLE, standday, adjpdg2, adje1c2, age, BMI))

write.table(data, "data_2cycle.txt", sep="\t", row.names=F)

