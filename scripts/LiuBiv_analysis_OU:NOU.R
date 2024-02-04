# Liu data biv analysis


#======================================================================================
#
# PLOTS, using data of 50 subjects of random selection
#======================================================================================


pOU <- 11
pNOU <- 15
r <- 56 # number of distinct time points
b <- 4 # number of fixed effect parameters

ou <- read.table("LiuBiv50_OUeq_cat.txt", sep = ",", header = T)
nou <- read.table("LiuBiv50_NOUeq_cat.txt", sep = ",", header = T)


#------------------------------------------------------------------------------
# beta hat

# OU

betaOU <- ou[(pOU+1):(pOU+b), 2]
sdBetaOU <- ou[(pOU+b+2*r+1):(pOU+2*b+2*r), 2]

ouLower <- tableOU[,1] - 1.96* tableOU[,2]
ouUpper <- tableOU[,1] + 1.96* tableOU[,2]

tableOU <- round(cbind(betaOU, sdBetaOU, ouLower, ouUpper), 4)


# NOU
betaNOU <- nou[(pNOU+1):(pNOU+b), 2]
sdBetaNOU <- nou[(pNOU+b+2*r+1):(pNOU+2*b+2*r), 2]

tableNOU <- round(cbind(betaNOU, sdBetaNOU, CIou), 4)
#------------------------------------------------------------------------------
# f hat

# OU
f1ou <- ou[(pOU+b+1):(pOU+b+r), 2]
f2ou <- ou[(pOU+b+r+1):(pOU+b+2*r), 2]

sdF1ou <- ou[(pOU+2*b+2*r+1):(pOU+2*b+3*r), 2]
sdF2ou <- ou[(pOU+2*b+3*r+1):(pOU+2*b+4*r), 2]

ouLower1 <- f1ou - 1.96*sdF1ou
ouUpper1 <- f1ou + 1.96*sdF1ou

ouLower2 <- f2ou - 1.96*sdF2ou
ouUpper2 <- f2ou + 1.96*sdF2ou

# NOU 

f1nou <- nou[(15+b+1):(15+b+r), 2]
f2nou <- nou[(15+b+r+1):(15+b+2*r), 2]

sdF1nou <- nou[(pNOU+2*b+2*r+1):(pNOU+2*b+3*r), 2]
sdF2nou <- nou[(pNOU+2*b+3*r+1):(pNOU+2*b+4*r), 2]

nouLower1 <- f1nou - 1.96*sdF1nou
nouUpper1 <- f1nou + 1.96*sdF1nou

nouLower2 <- f2nou - 1.96*sdF2nou
nouUpper2 <- f2nou + 1.96*sdF2nou

# plots OU vs NOU

par(mfrow = c(2, 1))

plot(f1ou, type = "l", col = "red", ylim = c(-0.5, 2.1), ylab = "f1 hat")
par(new = T)
plot(f1nou, type = "l", col = "blue", ylim = c(-0.5, 2.1), ylab = "f1 hat")

plot(f2ou, type = "l", col = "red", ylim = c(3.2, 4.5), ylab = "f2 hat")
par(new = T)
plot(f2nou, type = "l", col = "blue", ylim = c(3.2, 4.5), ylab = "f2 hat")


#------------------------------------------------------------------------------
# fitted 


fitted_ou <- ou[3062:5775, 2]
fitted1_ou <- fitted_ou[seq(1, 2714, by = 2)]
fitted2_ou <- fitted_ou[seq(2, 2714, by = 2)]


fitted_nou <- ou[3066:5779, 2]
fitted1_nou <- fitted_ou[seq(1, 2714, by = 2)]
fitted2_nou <- fitted_ou[seq(2, 2714, by = 2)]

# subject 1
par(mfrow = c(2, 1))
plot(fitted1_ou[1:29], type = "l", col = "blue", ylim = c(-1.7, 2.7), ylab = "f1 hat for subject 98")
par(new = T)
plot(log(data[1:29,3]), type = "l", col = "red", ylim = c(-1.7, 2.7), ylab = "f1 hat for subject 98")
par(new = T)
plot(fitted1_nou[1:29], type = "l", col = "green", ylim = c(-1.7, 2.7), ylab = "f1 hat for subject 98")

plot(fitted2_ou[1:29], type = "l", col = "blue", ylim = c(3.15, 5), ylab = "f2 hat for subject 98")
par(new = T)
plot(log(data[1:29,4]), type = "l", col = "red", ylim = c(3.15, 5), ylab = "f2 hat for subject 98")
par(new = T)
plot(fitted2_nou[1:29], type = "l", col = "green", ylim = c(3.15, 5), ylab = "f2 hat for subject 98")



#------------------------------------------------------------------------------
# random intercept
bOU <- ou[(pOU+2*b+4*r+5):(pOU+2*b+4*r+105), 2]

#======================================================================================
#
# PLOTS, using data of 50 subjects of random selection
#======================================================================================
#------------------------------------------------------------------------------
# Organize the dataset according to observed timepoints

tp <- seq(0.5, 28, by = 0.5)
logPro <- lapply(1:56, function(v) return(log(data[which(data$standday == tp[v]), 3])))
logEtg <- lapply(1:56, function(v) return(log(data[which(data$standday == tp[v]), 4])))
day <- lapply(1:56, function(v) return(rep(tp[v], length(logPro[[v]]))))

pro <- cbind(unlist(day), unlist(logPro))
etg <- cbind(unlist(day), unlist(logEtg))

#------------------------------------------------------------------------------
# Figure 1 in Zhang et al 1998

# OU
pdf('bivLiuFig1_OU.pdf')
par(mfrow = c(2, 1))
# PROGESTERONE
plot(pro, ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n")
par(new = T)
plot(f1ou, type = "l", lwd=2, ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n")
par(new = T)
plot(ouLower1, type = "l", col = "red", ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n", yaxt="n")
par(new = T)
plot(ouUpper1, type = "l", col = "green", ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n", yaxt="n")
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2))

# ESTROGEN
plot(etg, ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n")
par(new = T)
plot(f2ou, type = "l", lwd=2, ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen",  xaxt="n")
par(new = T)
plot(ouLower2, type = "l", col = "red", ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n", yaxt="n")
par(new = T)
plot(ouUpper2, type = "l", col = "green", ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n", yaxt="n")
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2))
dev.off()

# NOU
pdf('bivLiuFig1_NOU.pdf')
par(mfrow = c(2, 1))
# PROGESTERONE
plot(pro, ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n")
par(new = T)
plot(f1nou, type = "l", lwd=2, ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n")
par(new = T)
plot(nouLower1, type = "l", col = "red", ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n", yaxt="n")
par(new = T)
plot(nouUpper1, type = "l", col = "green", ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n", yaxt="n")
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2))

# ESTROGEN
plot(etg, ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n")
par(new = T)
plot(f2nou, type = "l", lwd=2, ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen",  xaxt="n")
par(new = T)
plot(nouLower2, type = "l", col = "red", ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n", yaxt="n")
par(new = T)
plot(nouUpper2, type = "l", col = "green", ylim = c(1.5, 5.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n", yaxt="n")
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2))
dev.off()




# FOR COMPARISON
# plot(pro, xaxt="n", xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", ylim = c(-3, 4.5))
# lines(lowess(pro, f = 0.4),  lwd=3)
# axis(1, at = seq(0, 28, by = 2), labels = seq(0, 28, by = 2), lwd = 0.5, xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone")
# par(new = T)
# plot(f1ou, type = "l", lwd=2, ylim = c(-3, 4.5), xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone")

# plot(etg, ylim = c(1.5, 5.9), xaxt="n", xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen")
# lines(lowess(etg, f = 0.4),  lwd=3)
# axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2), lwd = 0.5)
# par(new = T)
# plot(f2ou, type = "l", lwd=2, ylim = c(1.5, 6.9), xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen")


# scatter.smooth(proPlot, span = 0.4, xlab = "Days in standardized menstrual cycle", ylab = "Log progesterone", xaxt="n", lwd = 0.7)
# axis(1, at = seq(0, 28, by = 2), labels = seq(0, 28, by = 2), lwd = 0.5)
# scatter.smooth(etgPlot, span = 0.4, xlab = "Days in standardized menstrual cycle", ylab = "Log estrogen", xaxt="n")
# axis(1, at = seq(0, 28, by = 2), labels = seq(0, 28, by = 2), lwd = 0.5)


#------------------------------------------------------------------------------
# Figure 2 in Zhang et al 1998

# Obtain sample variances of the log progesterone/estrogen values at each time points
tp <- seq(0.5, 28, by = 0.5)
sv1 <- sapply(1:56, function(v) return(var(log(data[which(data$standday == tp[v]), 3]))))
sv2 <- sapply(1:56, function(v) return(var(log(data[which(data$standday == tp[v]), 4]))))

tp_plot <- seq(0, 27.5, by = 0.5)
pdf('bivLiuFig2.pdf')
par(mfrow = c(2, 1))
# Adding the xaxt="n" says to not use the default number scheme for the x-axis
scatter.smooth(sv1, span = 0.4, xlab = "Days in standardized menstrual cycle", ylab = "Variance of log progesterone", xaxt="n")
# The axis function places the custom labels.
# at/labels = a vector containing the position/label of the ticks; at and labels must be of the same length
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2), lwd = 0.5)
scatter.smooth(sv2, span = 0.4, xlab = "Days in standardized menstrual cycle", ylab = "Variance of log estrogen", xaxt="n")
axis(1, at = seq(0, 56, by = 4), labels = seq(0, 28, by = 2), lwd = 0.5)
dev.off()

#------------------------------------------------------------------------------
# 
