mds_aic <- function(Theta, S ,n) {

  M <- length(Theta)
   #Number of non zero Elements
  crit = 0
  for (m in 1:M) {
    E <- length(which(Theta[[m]] != 0))

    ldet <- log(det(Theta[[m]]))
    tr <- sum(diag(S[[m]] %*% Theta[[m]]))

    crit = crit +( (n[m]* (-ldet + tr))   + 2*E )
  }


  return(crit)
}


mds_bic <- function(Theta, S, n , p){
  M <- length(Theta)
  E <- ((p^2)-p) /2
  crit = 0

  for (m in 1:M) {

    ldet <- log(det(Theta[[m]]))
    tr <- sum(diag(S[[m]] %*% Theta[[m]]))
    crit = crit +( (2*n[m]* (-ldet + tr))   + (E*log(n[[m]]))   )
  }
  return(crit)
}
