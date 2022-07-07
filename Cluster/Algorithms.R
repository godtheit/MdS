# MdS
library(MdS2)
#Lambdas for MdS




MdS_wrapper <- function(data, job, instance, a, b, c, d, L1, L2, weight_matrix , ...) {
  data <- instance$samples
  M <- length(data)

  if (weight_matrix == "full"){
    w <- matrix(1 , M, M)
    w[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0
  }

  if (weight_matrix == "horizontal-vertikal") {
    th11 <- c(0,0,0,0,0,0,0,0,0) #1
    th12 <- c(1,0,0,0,0,0,0,0,0) #2
    th13 <- c(1,1,0,0,0,0,0,0,0) #3
    th21 <- c(1,0,0,0,0,0,0,0,0) #4
    th22 <- c(0,1,0,1,0,0,0,0,0) #5
    th23 <- c(0,0,1,1,1,0,0,0,0) #6
    th31 <- c(1,0,0,1,0,0,0,0,0) #7
    th32 <- c(0,1,0,0,1,0,1,0,0) #8
    th33 <- c(0,0,1,0,0,1,1,1,0) #9

    w <- rbind(th11,th12, th13, th21, th22, th23, th31, th32, th33)
    colnames(w) <- rownames(w)
  }

  lam1 <- seq(a,b, length.out = L1)
  lam2 <- seq(c,d, length.out = L2)



  best_mds <- procedure(Y = data, lam1 = lam1, lam2 = lam2, method = "MdS", w = w)

return(best_mds)
}
#jgl
JGL_wrapper <- function(data, a, b, c, d, L1, L2, job, instance, ...) {
  #Lambdas for JGL
  lam1 <- seq(a,b, length.out = L1)
  lam2 <- seq(c,d, length.out = L2)
  data <- instance$samples
  best_jgl <- procedure(Y = data, lam1 = lam1, lam2 = lam2, method = "fused")

  return(best_jgl)
}
#glasso
glasso_wrapper <- function(data, job, instance, ...) {
  lam1 <-  1  # these
  lam2 <-  2  # are not used
  data <- instance$samples
  best_glasso <- procedure(Y = data, lam1 = lam1, lam2 = lam2, method = "glasso")

  return(best_glasso)
}
