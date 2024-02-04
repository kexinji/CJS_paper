# Plots for data with 3 concecutive cycles

dataComp <- read.table("data_3cycle.txt", header = T)

#======================================================================================
# Panel plots, using raw longitudinal data
#======================================================================================
# plot the resonses for the first subj. 

data1 <- subset(data, data[,1] == 1)

logpro <- log(data1[,3])
logetg <- log(data1[,4])

plot(data1[,2], logpro, type = "l")
plot(data1[,2], logetg, type = "l")

data2 <- subset(data, data[,1] == 2)
data3 <- subset(data, data[,1] == 3)
data4 <- subset(data, data[,1] == 4)
data5 <- subset(data, data[,1] == 5)

pdf('plots_3c_5subj_1pg.pdf')
par(mfrow=c(5,3))
#------------------------------------------------------------------------------
# subj. 1
# cycle 1, subj 1
plot(log(subset(data1, data1[,2] == 1)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.1")
par(new = T)
plot(log(subset(data1, data1[,2] == 1)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.1")
# cycle 2, subj. 1
plot(log(subset(data1, data1[,2] == 2)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.1")
par(new = T)
plot(log(subset(data1, data1[,2] == 2)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.1")
# cycle 3, subj. 1
plot(log(subset(data1, data1[,2] == 3)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.1")
par(new = T)
plot(log(subset(data1, data1[,2] == 3)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.1")

#------------------------------------------------------------------------------
# subj. 2
# cycle 1
plot(log(subset(data2, data2[,2] == 1)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.2")
par(new = T)
plot(log(subset(data2, data2[,2] == 1)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.2")
# cycle 2
plot(log(subset(data2, data2[,2] == 2)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.2")
par(new = T)
plot(log(subset(data2, data2[,2] == 2)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.2")
# cycle 3
plot(log(subset(data2, data2[,2] == 3)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.2")
par(new = T)
plot(log(subset(data2, data2[,2] == 3)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.2")

#------------------------------------------------------------------------------
# subj. 3
# cycle 1
plot(log(subset(data3, data3[,2] == 1)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.3")
par(new = T)
plot(log(subset(data3, data3[,2] == 1)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.3")
# cycle 2
plot(log(subset(data3, data3[,2] == 2)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.3")
par(new = T)
plot(log(subset(data3, data3[,2] == 2)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.3")
# cycle 3
plot(log(subset(data3, data3[,2] == 3)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.3")
par(new = T)
plot(log(subset(data3, data3[,2] == 3)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.3")
# dev.off()

#------------------------------------------------------------------------------
# subj. 4
# cycle 1
plot(log(subset(data4, data4[,2] == 1)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.4")
par(new = T)
plot(log(subset(data4, data4[,2] == 1)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.4")
# cycle 2
plot(log(subset(data4, data4[,2] == 2)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.4")
par(new = T)
plot(log(subset(data4, data4[,2] == 2)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.4")
# cycle 3
plot(log(subset(data4, data4[,2] == 3)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.4")
par(new = T)
plot(log(subset(data4, data4[,2] == 3)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.4")


#------------------------------------------------------------------------------
# subj. 5
# cycle 1
plot(log(subset(data5, data5[,2] == 1)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.5")
par(new = T)
plot(log(subset(data5, data5[,2] == 1)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c1 subj.5")
# cycle 2
plot(log(subset(data5, data5[,2] == 2)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.5")
par(new = T)
plot(log(subset(data5, data5[,2] == 2)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c2 subj.5")
# cycle 3
plot(log(subset(data5, data5[,2] == 3)[,4]), type = "l", xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.5")
par(new = T)
plot(log(subset(data5, data5[,2] == 3)[,5]), type = "l", lty = 2, xlim = c(0, 28), ylim = c(-3, 5), xlab = "Days", ylab = "Log responses c3 subj.5")

dev.off()


#======================================================================================
# Plots, using data of m subjects of random selections
#======================================================================================

#------------------------------------------------------------------------------
# m = 10, randomly select 10 subjects with 3 concecutive cycles
id <- dataComp[,1]
m <- 20

set.seed(42)
s <- sample(1:max(id), m)

data <- lapply(1:m, function(v) return(rbind(dataComp[which(id == s[v]),])))
data <- do.call(rbind, data) # 906*6

data$day <- data$standday + 28*(data$CYCLE-1)

#------------------------------------------------------------------------------
# Organize the dataset according to observed timepoints

tp <- sort(unique(data[,3]))
toltp <- length(tp)
logPro <- lapply(1:toltp, function(v) return(log(data[which(data$standday == tp[v]), 4])))
logEtg <- lapply(1:toltp, function(v) return(log(data[which(data$standday == tp[v]), 5])))
day <- lapply(1:toltp, function(v) return(rep(tp[v], length(logPro[[v]]))))

pro <- cbind(unlist(day), unlist(logPro))
etg <- cbind(unlist(day), unlist(logEtg))

f1 <- res$f1Hat
f2 <- res$f2Hat

sdF1 <- res$sdF1
sdF2 <- res$sdF2

lowerF1 <- f1 - 1.96*sdF1
upperF1 <- f1 + 1.96*sdF1

lowerF2 <- f2 - 1.96*sdF2
upperF2 <- f2 + 1.96*sdF2

#------------------------------------------------------------------------------
# plots
pdf("LiuBiv_3c_fnew.pdf")
par(mfrow=c(2,1))

# progesterone
plot(pro, pch = 39, ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n")
par(new = T)
plot(tp, rep(f1, 3), type = "l",ylim = c(-3, 4), lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(lowerF1, 3), type = "l", col = "red", ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(upperF1, 3), type = "l", col = "green", ylim = c(-3, 4), xlab = "", ylab = "", xaxt="n", yaxt="n")
title(xlab="Days in standardized menstrual cycle", ylab="Log progesterone")
axis(1, at = seq(0, 84, by = 7), labels = c(seq(0, 28, by = 7), rep(seq(7, 28, by = 7),2)))

# estrogen
plot(etg, pch = 39, ylim = c(1, 5.5), xlab = "", ylab = "", xaxt="n")
par(new = T)
plot(tp, rep(f2, 3), type = "l",ylim = c(1, 5.5), lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(lowerF2, 3), type = "l", col = "red", ylim = c(1, 5.5), xlab = "", ylab = "", xaxt="n", yaxt="n")
par(new = T)
plot(tp, rep(upperF2, 3), type = "l", col = "green", ylim = c(1, 5.5), xlab = "", ylab = "", xaxt="n", yaxt="n")
title(xlab="Days in standardized menstrual cycle", ylab="Log estrogen")
axis(1, at = seq(0, 84, by = 7), labels = c(seq(0, 28, by = 7), rep(seq(7, 28, by = 7),2)))
dev.off()

plot(etg, pch = 39, ylim = c(1.5, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n")

scatter.smooth(pro, span = 0.1)
scatter.smooth(etg, span = 0.1, xlab = "Days in standardized menstrual cycle", ylab = "Variance of log estrogen")


# progesterone
pdf("plots_3cTest.pdf")
par(mfrow=c(2,1))

plot(data[which(data[,1] == s[1]), 9], log(data[which(data[,1] == s[1]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone")
par(new = T)
plot(data[which(data[,1] == s[2]), 9], log(data[which(data[,1] == s[2]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[3]), 9], log(data[which(data[,1] == s[3]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[4]), 9], log(data[which(data[,1] == s[4]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[5]), 9], log(data[which(data[,1] == s[5]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[6]), 9], log(data[which(data[,1] == s[6]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[7]), 9], log(data[which(data[,1] == s[7]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[8]), 9], log(data[which(data[,1] == s[8]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[9]), 9], log(data[which(data[,1] == s[9]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[10]), 9], log(data[which(data[,1] == s[10]), 4]), xlim = c(1,84), ylim = c(-3.5,3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(f1, 3), type = "l", lwd=2, col = "red", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(lowerF1, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")
par(new = T)
plot(rep(upperF1, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(-3.5, 3.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", yaxt="n", xaxt="n")

#------------------------------------------------------------------------------
# estrogen

plot(data[which(data[,1] == s[1]), 9], log(data[which(data[,1] == s[1]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen")
par(new = T)
plot(data[which(data[,1] == s[2]), 9], log(data[which(data[,1] == s[2]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[3]), 9], log(data[which(data[,1] == s[3]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[4]), 9], log(data[which(data[,1] == s[4]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[5]), 9], log(data[which(data[,1] == s[5]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[6]), 9], log(data[which(data[,1] == s[6]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[7]), 9], log(data[which(data[,1] == s[7]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[8]), 9], log(data[which(data[,1] == s[8]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[9]), 9], log(data[which(data[,1] == s[9]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(data[which(data[,1] == s[10]), 9], log(data[which(data[,1] == s[10]), 5]), xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(f2, 3), type = "l", lwd=2, col = "red", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(lowerF2, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
par(new = T)
plot(rep(upperF2, 3), type = "l", col = "green", xlim = c(1,84), ylim = c(1, 6), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", yaxt="n", xaxt="n")
dev.off()

#------------------------------------------------------------------------------
# f hat
# liuResult <- read.table("Liu3c_test_v2.txt", sep = ",", header = T)
liuResult <- read.table("Liu3c_test.txt", sep = ",", header = T)

p <- 15
r <- 28 # number of distinct time points
# r <- 353
b <- 6 # number of fixed effect parameters, only age

f1 <- liuResult[(15+b+1):(15+b+r), 2]
f2 <- liuResult[(15+b+r+1):(15+b+2*r), 2]

plot(tp, f1, type = "l")

pdf("plots_3c_1cycle.pdf")
par(mfrow=c(2,1))
plot(f1, type = "l")
plot(f2, type = "l")
dev.off()

pdf("plots_3c.pdf")
par(mfrow=c(2,1))
plot(rep(f1,3), type = "l")
plot(rep(f2,3), type = "l")
dev.off()



sdF1 <- liuResult[(p+2*b+2*r+1):(p+2*b+3*r), 2]
sdF2 <- liuResult[(p+2*b+3*r+1):(p+2*b+4*r), 2]

lowerF1 <- f1 - 1.96*sdF1
upperF1 <- f1 + 1.96*sdF1

lowerF2 <- f2 - 1.96*sdF2
upperF2 <- f2 + 1.96*sdF2
