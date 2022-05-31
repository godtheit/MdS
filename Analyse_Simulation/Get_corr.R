get_cor <- function(theta, observes, p, u = 0.1, v = 0.3){
  diag(theta) = 0
  omega = theta * v
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  
  
  omega = solve(sigma)
  x = mvrnorm(observes, rep(0, p), sigma)
  if(observes > 1){
    sigmahat = cor(x) 
    return(list(sigma,x, sigmahat))
  }
  else return(list(sigma,x))
}

