# MdS
library(MdS2)

MdS_wrapper <- function(data, job, instance, ...) {

  current_network <- instance[[1]]
  best_mds <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "MdS", w = w)

return(best_mds)
}
#jgl
JGL_wrapper <- function(data, job, instance, ...) {

  current_network <- instance[[1]]
  best_jgl <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "JGL")

  return(best_mds)
}
#glasso
glasso_wrapper <- function(data, job, instance, ...) {

  current_network <- isntance[[1]]
  best_mds <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "MdS")

  return(best_mds)
}
