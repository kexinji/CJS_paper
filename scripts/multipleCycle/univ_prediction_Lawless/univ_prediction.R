library(sde) # rsOU()

#--------------------------------------------------------------
# Parameter Initialization for OU process, a = theta[2], sigma = theta[3]
thetaOU <- c(0, 2, 2)

# obtain true alpha, and var
alpha <- exp(-thetaOU[2])
vari <- thetaOU[3]^2/(2*thetaOU[2])

# Initialize the first observation for OU
u0 <- rnorm(1, 0, vari)

#--------------------------------------------------------------
# Generate obs'ns from OU

# Generate 10 observations from OU, m iterations
n <- 10 
m <- 20
ugen <- replicate(m, rsOU(n, theta=thetaOU)) # each column is a new iteration

# combine u0 and u
u <- rbind(u0, ugen)

#-------------------------------------------------------------
# obtain alpha_star, a_star, and var_star from generated data

# obtain alpha star & varStar using u, for m iterations
alphaNu <- colSums(u[1:n,]*ugen)
alphaDe <- colSums(u^2)
alphaStar <- alphaNu/alphaDe

varStar <- colSums((ugen - t(alphaStar * t(u[1:n,])))^2/n)

# obtain alpha star, variance using u, for one iteration
# 
# alphaNu <- sum(u[1:n]*ugen)
# alphaDe <- sum(u^2)
# alphaStar <- alphaNu/alphaDe
  
# obtain a star
# aStar <- -log(alphaStar)

# varStar <- sum((ugen - alphaStar*u[1:n])^2/n)

#-------------------------------------------------------------
# predictive distribution function for U

# for one iteration
pred <- function(w){
  num <- (alphaStar - alpha)*ugen[n] + w*sqrt(varStar * (1-alphaStar^2))
  den <- sqrt(vari*(1 - alpha^2))
  pnorm(num/den)
}
plot(pred(-5:5), type = "l", ylim = c(0, 1), xlab = "w", ylab = "G(w)", xaxt="n")
axis(1, at = seq(2, 11, by = 2), labels = seq(-4, 4, by = 2))

# for the i^th iterations
pred <- function(w, i) {
  # num <- (alphaStar - alpha) * ugen[n,] + w*sqrt(varStar * (1-alphaStar^2))
  den <- sqrt(vari*(1 - alpha^2))
  # return(num/den)
  
  num <- ((alphaStar - alpha) * ugen[n,])[i] + w*(sqrt(varStar * (1-alphaStar^2)))[i]
  predf <- pnorm(num[i]/den)
  
  return(mean(predf))
}

# for m iterations
pred <- function(w) {
  # num <- (alphaStar - alpha) * ugen[n,] + w*sqrt(varStar * (1-alphaStar^2))
  den <- sqrt(vari*(1 - alpha^2))
  # return(num/den)
  predf <- matrix(0, m, 1)
  num <- matrix(0, m, 1)
  for (i in 1:n) {
    num[i] <- ((alphaStar - alpha) * ugen[n,])[i] + w*(sqrt(varStar * (1-alphaStar^2)))[i]
    predf[i] <- pnorm(num[i]/den)
  }
  return(mean(predf))
}
