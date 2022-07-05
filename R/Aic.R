mds_aic <- function(Theta, S ,n, lambda1, lambda2) {

  M <- length(Theta)
   #Number of non zero Elements
  crit = 0
  for (m in 1:M) {
    E <- length(which(Theta[[m]] != 0))
    crit = crit +( n[m]* sum(diag((S[[m]] %*% Theta[[m]])) - n[m]*log(det(Theta[[m]])) + 2*E ))
  }


  return(crit)
}
