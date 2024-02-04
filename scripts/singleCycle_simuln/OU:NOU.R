ornstein_uhlenbeck <- function(T,n,nu,lambda,sigma,x0){
  dw  <- rnorm(n, 0, sqrt(T/n))
  dt  <- T/n
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + lambda*(nu-x[i-1])*dt + sigma*dw[i-1]
  }
  return(x);
}


library(sde)
> rsOU
function (n = 1, theta) 
{
    checkOU(theta)
    rnorm(n, mean = theta[1]/theta[2], sd = theta[3]/sqrt(2 * 
        theta[2]))
}

theta <- c(0, 1, 1)
# rnorm and rsOU returns the same values. 
set.seed(1)
rnorm(10, mean = theta[1]/theta[2], sd = theta[3]/sqrt(2 * theta[2]))

set.seed(1)
rsOU(10, c(0, 1, 1))

set.seed(1)
ornstein_uhlenbeck(10, 10, 0, 1, 1, -0.4429697)