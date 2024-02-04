library(sde) # rsOU()

#--------------------------------------------------------------
# Parameter Initialization for OU process, a = theta[2], sigma = theta[3]
thetaOU <- c(0, 2, 2)

# obtain true alpha, and var
alpha <- exp(-thetaOU[2])
vari <- thetaOU[3]^2/(2*thetaOU[2])

# Initialize the first observation for OU
u0 <- rnorm(1, 0, vari)

#-------------------------------------------------------------
# predictive distribution function for U

# for one iteration
pred <- function(w, simulationData, n){
  
  # obtain alpha star, variance from generated data u, for one iteration
  alphaNu <- sum(u[1:n]*ugen)
  alphaDe <- sum(u^2)
  alphaStar <- alphaNu/alphaDe
  
  varStar <- sum((ugen - alphaStar*u[1:n])^2/n)
  
  # obtain the distribution function from u
  num <- (alphaStar - alpha)*ugen[n] + w*sqrt(varStar * (1-alphaStar^2))
  den <- sqrt(vari*(1 - alpha^2))
  return(pnorm(num/den))
}


#--------------------------------------------------------------
# Generate obs'ns from OU

n <- 10

# predf <- function(m){

  m <- 10
  store <- matrix(0, 11, m) # store predictive distribution function values
  
  for (i in 1:m) {
    
    # Generate 10 observations from OU, m iterations, each column is a new iteration
    ugen <- rsOU(n, theta=thetaOU)
    
    # combine u0 and u
    u <- c(u0, ugen)
    
    store[,i] <- pred(-5:5, u, n)
    storeMean <- rowMeans(store)
  }
#   return(storeMean)
# }


plot(storeMean, type = "l", ylim = c(0, 1), xlab = "w", ylab = "G(w)", xaxt="n") # n = 100 iterations
par(new = T)
plot(storeMean, col="red", type = "l", ylim = c(0, 1), xlab = "w", ylab = "G(w)", xaxt="n") # n = 10 iterations
par(new = T)
plot(storeMean, col="blue", type = "l", ylim = c(0, 1), xlab = "w", ylab = "G(w)", xaxt="n") # 500 iterations
axis(1, at = seq(2, 11, by = 2), labels = seq(-4, 4, by = 2))
