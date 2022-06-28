lambda_selector <- function(Theta, S ,n, lambda1, lambda2) {

  M <- length(Theta)

  crit = 0
  for (m in 1:M) {
    crit = crit + n[m]* trace()
  }


  return(crit, best_lam1,best_lam2)
}
