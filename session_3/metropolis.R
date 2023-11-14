metropolis <- function(n){
  theta0 = runif(1)
  theta = c()
  theta[1] = theta0
  for(k in 2:n){
    thetastar <- rnorm(1, theta[k-1], sd=0.05)
    lpstar = logposterior(survived, thetastar)
    lpprev= logposterior(survived, theta[k-1])
    lr = lpstar-lpprev
    r = exp(lr)
    accept = rbinom(1,1,min(1,r))
    if (accept == 1){
      theta[k] = thetastar}
    else {
      theta[k] = theta[k-1]}
  }
  return(theta)
}
theta = metropolis(10000)
plot(theta)
density(theta)
plot(density(theta))