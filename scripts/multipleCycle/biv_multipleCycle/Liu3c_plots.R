# Plots for data with 3 concecutive cycles

dataComp <- read.table("data_3cycle.txt", header = T)

#------------------------------------------------------------------------------
# randomly select m subjects with 3 concecutive cycles
id <- dataComp[,1]
m <- 20

set.seed(42)
s <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data) # 906*6

data$day <- data$standday + 28*(data$CYCLE-1)

#------------------------------------------------------------------------------
# f hat
# liuResult <- read.table("Liu3c_test_v2.txt", sep = ",", header = T)
liuResult <- read.table("Liu3c_20subj.txt", sep = ",", header = T)

p <- 15
r <- 28 # number of distinct time points
# r <- 353
b <- 6 # number of fixed effect parameters, only age

f1 <- liuResult[(15+b+1):(15+b+r), 2]
f2 <- liuResult[(15+b+r+1):(15+b+2*r), 2]

plot(f1, type = "l")
plot(f2, type = "l")


sdF1 <- liuResult[(p+2*b+2*r+1):(p+2*b+3*r), 2]
sdF2 <- liuResult[(p+2*b+3*r+1):(p+2*b+4*r), 2]

lowerF1 <- f1 - 1.96*sdF1
upperF1 <- f1 + 1.96*sdF1

lowerF2 <- f2 - 1.96*sdF2
upperF2 <- f2 + 1.96*sdF2


#------------------------------------------------------------------------------
# plot 1

# par(mfrow=c(2,1))

dev.new(width=8, height=4)

plot(data[which(data[,1] == s[1]), 8], log(data[which(data[,1] == s[1]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone")
par(new = T)
plot(data[which(data[,1] == s[2]), 8], log(data[which(data[,1] == s[2]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[3]), 8], log(data[which(data[,1] == s[3]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[4]), 8], log(data[which(data[,1] == s[4]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[5]), 8], log(data[which(data[,1] == s[5]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[6]), 8], log(data[which(data[,1] == s[6]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[7]), 8], log(data[which(data[,1] == s[7]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[8]), 8], log(data[which(data[,1] == s[8]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[9]), 8], log(data[which(data[,1] == s[9]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[10]), 8], log(data[which(data[,1] == s[10]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[11]), 8], log(data[which(data[,1] == s[11]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone")
par(new = T)
plot(data[which(data[,1] == s[12]), 8], log(data[which(data[,1] == s[12]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[13]), 8], log(data[which(data[,1] == s[13]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[14]), 8], log(data[which(data[,1] == s[14]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[15]), 8], log(data[which(data[,1] == s[15]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[16]), 8], log(data[which(data[,1] == s[16]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[17]), 8], log(data[which(data[,1] == s[17]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[18]), 8], log(data[which(data[,1] == s[18]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[19]), 8], log(data[which(data[,1] == s[19]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[20]), 8], log(data[which(data[,1] == s[20]), 4]), xlim = c(1,84), ylim = c(-3, 4), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(f1, 3), type = "l", lwd=2, col = "red", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(lowerF1, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(upperF1, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")

#------------------------------------------------------------------------------
# estrogen

dev.new(width=8, height=4)

plot(data[which(data[,1] == s[1]), 8], log(data[which(data[,1] == s[1]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen")
par(new = T)
plot(data[which(data[,1] == s[2]), 8], log(data[which(data[,1] == s[2]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[3]), 8], log(data[which(data[,1] == s[3]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[4]), 8], log(data[which(data[,1] == s[4]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[5]), 8], log(data[which(data[,1] == s[5]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[6]), 8], log(data[which(data[,1] == s[6]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[7]), 8], log(data[which(data[,1] == s[7]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[8]), 8], log(data[which(data[,1] == s[8]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[9]), 8], log(data[which(data[,1] == s[9]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[10]), 8], log(data[which(data[,1] == s[10]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[11]), 8], log(data[which(data[,1] == s[11]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen")
par(new = T)
plot(data[which(data[,1] == s[12]), 8], log(data[which(data[,1] == s[12]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[13]), 8], log(data[which(data[,1] == s[13]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[14]), 8], log(data[which(data[,1] == s[14]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[15]), 8], log(data[which(data[,1] == s[15]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[16]), 8], log(data[which(data[,1] == s[16]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[17]), 8], log(data[which(data[,1] == s[17]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[18]), 8], log(data[which(data[,1] == s[18]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[19]), 8], log(data[which(data[,1] == s[19]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[20]), 8], log(data[which(data[,1] == s[20]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(f2, 3), type = "l", lwd=2, col = "red", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(lowerF2, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(upperF2, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")

