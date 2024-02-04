pred_results <- read.table("pred_univSim50.txt")

pse <- round(pred_results[1, 1:7],4)

dev.new(width=8, height=2)
#pdf("RP_simUnivPred50.pdf")
par(mfrow=c(2,2))
matplot(1:28, cbind(pred_results[2:29,c(1:3,5:6)],CI_lower1, CI_upper1), type="l", lty = c(1:5,6,6), col=c("black", "blue", "red", "green", "orange", "pink", "pink"), xlab = "time points", ylab = "response", xaxt="n", main = "Smallest PSE ")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(1:28, cbind(pred_results[30:57,c(1:3,5:6)],CI_lower2, CI_upper2), type="l", lty = c(1:5,6,6), col=c("black", "blue", "red", "green", "orange", "pink", "pink"), xlab = "time points", ylab = "response", xaxt="n", main = "Second Smallest PSE")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(1:28, cbind(pred_results[58:85,c(1:3,5:6)],CI_lower3, CI_upper3), type="l", lty = c(1:5,6,6), col=c("black", "blue", "red", "green", "orange", "pink", "pink"), xlab = "time points", ylab = "response", xaxt="n", main = "Median PSE")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))

matplot(1:28, cbind(pred_results[86:113,c(1:3,5:6)],CI_lower4, CI_upper4), type="l", lty = c(1:5,6,6), col=c("black", "blue", "red", "green", "orange", "pink", "pink"), xlab = "time points", ylab = "response", xaxt="n", main = "Largest PSE")
axis(1, at = seq(0, 28, by = 4), labels = seq(0, 28, by = 4))
#dev.off()

#--------------------------------------------------------------------------
t <- (1:56 %% 28)/10
ij <- 29:56
ij_1 <- 28:55

theta <- c(rep(1, 3), 0.1, -3, -5, -2)
alpha <- exp(abs(t[ij] - t[ij_1])*log(theta[4]))

# A HELPER FUNCTION that RETURNS coefficients of periodic cubic splines
s <- function(knot, t) {
  T <- 28
  a <- -T*(T - knot)/2 + 3*(T - knot)^2/2 - (T - knot)^3/T
  b <- 3*(T - knot)/2 - 3*(T - knot)^2/(2*T)
  c <- -(T - knot)/T
  sj <- a*t + b*t^2 + c*t^3 + (t - knot)^3*((t - knot)^3>0)
  return(sj)
}

# Variance
expo <- theta[5] + (theta[6]*(s(28/3, t[ij_1]) + s(28/3, t[ij])) + theta[7]*(s(56/3, t[ij_1]) + s(56/3, t[ij])))/2
variance <- (1 + alpha^2)*1 + (1 + alpha^2)*1 + exp(expo)*(1 - alpha^2)

# Confidence Intervals
sd <- sqrt(variance)
CI_lower1 <- as.matrix(pred_results[2:29,2] - 1.96*sd)
CI_upper1 <- as.matrix(pred_results[2:29,2] + 1.96*sd)

CI_lower2 <- as.matrix(pred_results[30:57,2] - 1.96*sd)
CI_upper2 <- as.matrix(pred_results[30:57,2] + 1.96*sd)

CI_lower3 <- as.matrix(pred_results[58:85,2] - 1.96*sd)
CI_upper3 <- as.matrix(pred_results[58:85,2] + 1.96*sd)

CI_lower4 <- as.matrix(pred_results[86:113,2] - 1.96*sd)
CI_upper4 <- as.matrix(pred_results[86:113,2] + 1.96*sd)

# matplot(1:28, cbind(rowMeans(resp), rowMeans(mu),  rowMeans(mu_pc), rowMeans(mu_pcv), rowMeans(mu_pre), CI_lower, CI_upper), type="l", lty = c(1:5, 6, 6), col=c("black", "orange", "blue", "red", "green", "pink", "pink"))
