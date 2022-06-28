MdS <- function(Y, rho = 1, lambda1 = 0.1, lambda2= 0.1, w , dims, epsilon = 1e-3, cores = 1,  maxIter = 1000){

  library(DescTools)
  library(parallel)

  p = dim(Y[[1]])[2] # number of variables
  M = length(Y)      # number of networks
  n = rep(0, M)      # number observations for each network
  for (m in 1:M) {
    n[m] = dim(Y[[m]])[1]
  }
  if (length(dimnames(Y[[1]])[[2]]) == 0) {
    for (m in 1:M) {
      dimnames(Y[[m]])[[2]] = paste("V", 1:p, sep = "")
    }
  }

  #Daten zentralisieren
  for (m in 1:M) {
    for (j in 1:p) {
      Y[[m]][, j] = Y[[m]][, j] - mean(Y[[m]][, j])
    }
  }


  #Zum Test der Funktion

  #library(CVTuningCov)
  # M <- 3
  # load("allsamples.RData")
  # data <- BigListof_samples[[1]]
  # #load("allcovars.RData")
  # #data <- BigListof_covs[[1]]
  #
  # for (m in 1:M) {
  #   dat <- round(data[[m]][1:20,1:20], digits = 3)
  #   data[[m]] <- cov(dat)
  # }
  # psi <- "L1"
  # dims <- c(3,1)
  # n <- 5
  # p <- 20
  # rho <- 2
  # lambda <- 0.5
  # epsilon = 1e-3
  # w <- matrix(0.5,prod(dims),prod(dims))
  # w[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0
  # maxIter <- 570
  #
  #
  #





  Empty <- matrix(0,p,p)

  Theta <- Z  <- U <- betas <- list()



  for (m in 1:M) {
    Theta[[m]] <- diag(p)
    Z[[m]] <- U[[m]] <- betas[[m]] <- Empty    #Start the algorithm with a empty network

  }

  k <- 0






  diff_value = 10
  while((k==0) || ((k < maxIter) && diff_value > epsilon))
  {

    Theta_old <- Theta
    Z_old <- Z
    U_old <- U




  # repeat{
  #   for (m in 1:M) {                            # Zu diesem Zeitpunkt ist U == U_old, da U erst am Ende geupdated wird
  #     if((all(abs(U0_old[[m]]-U0_old2[[m]]) < Eps) == TRUE) &
  #        (all(abs(U1_old[[m]]-U1_old2[[m]]) < Eps) == TRUE) &
  #        (all(abs(U2_old[[m]]-U2_old2[[m]]) < Eps))
  #        == TRUE){
  #       exitStatement[m] <- TRUE
  #     } else {
  #       exitStatement[m] <- FALSE
  #       }
  #   }
  #   if (k == maxIter) {
  #     break("MaxIter has been reached")
  #   }
  #   if(all(exitStatement == TRUE)){
  #     break("Converged")
  #     cat("Hurray")
  #   }





  # M1 <- 1:M
  # M2 <- 1:M
  #
  # if(verschiebung == (M-1))
  # {
  #   verschiebung <- 0
  # } else {
  #   verschiebung <- verschiebung + 1
  # }
  #
  #
  #   if (verschiebung == 0){
  #     M2 <- M1
  #   } else {
  #
  #     M2[1:verschiebung] <-  M1[(M-verschiebung+1):M]
  #     M2[(verschiebung+1):M] <- M1[1:(M-verschiebung)]
  #  }

  #Update Step Theta
    S <- list()
    eta <- list()
    A <- list()

  for (m in 1:M) {
    S[[m]] <- (t(Y[[m]]) %*% Y[[m]]) / n[m]   #empirical covariance
    A[[m]] <- ((Z_old[[m]] - U_old[[m]]) )
    eta[[m]] <- n[m] / rho

  }


    for (m in 1:M) {

      #formula for theta from danaher/hallac
    G <- (eta[[m]]^(-1)) * ((A[[m]] + t(A[[m]]))/2 ) - S[[m]]
    G_Eigen <- eigen(G) #-1 for danaher formula

    Q <- G_Eigen$vectors
    D <- diag(G_Eigen$values)

    V <- D + sqrt(D^2 + ((4*(eta[[m]]^(-1)))*diag(p)))

    Theta_k <- (1/(2*(eta[[m]]^(-1)))) * (Q%*%  V  %*% t(Q))
    Theta[[m]] <- Theta_k
  }

    #Update Step Z


    y <- list()

      for (m in 1:M) {
       y[[m]] <- Theta[[m]] + U_old[[m]]
      }


    #y <- mapply(function(x,y){x+y}, Theta, U_old)

   eta_1 <- lambda1/rho
   eta_2 <- lambda2/rho

## Create Matrix D
     signum <- c(c(1,1), vector(mode = "integer",length = M-2))
     B <-  Permn(signum)
     Flipp <- FALSE
     for (i in 1:length(B[,1])) {
       for (j in 1:length(B[1,])) {
          if(Flipp == TRUE & B[i,j] == 1){
             B[i,j] <- -1
             Flipp <- FALSE
          } else if (B[i,j] == 1){
         Flipp <- TRUE
           }

        }
     }
     B_final <- matrix(0,length(B[,1]),length(B[1,]))


     for (i in 1:length(B[,1])) {
       weights <- which( abs(B[i,])== 1)
       B_final[i,] <- eta_2 * w[weights[2],weights[1]] * B[i,]
     }
     C <- eta_1 * diag(M)
     D <- rbind(C, B_final)



     row_with_only_zeros <- apply(D, 1, function(row) all(row == 0))
     D <- D[!row_with_only_zeros, ]


     for (m in 1:M) {
       diag(Z[[m]]) <- diag(y[[m]])
     }


     ### Use generalized Lasso on all possible Edges
combinations <- combn(1:p, 2, simplify = FALSE)

mclapply(combinations, function(combs){
  yps <- sapply(y, function(M) M[combs[1], combs[2]])



  #genL <- dualpath(y = yps, D = D, verbose = FALSE)
  #out <- genL$beta

  genL <- genlasso(y = yps, diag(1, m), D, minlam = 1, verbose = F)
  out <- coef(genL, lambda = 1)$beta

  out[which(abs(out) <= 10^-12)] <- 0

  for (m in 1:M) {
    Z[[m]][combs[1],combs[2]] <<- Z[[m]][combs[2],combs[1]]  <<- out[m]
  }

}, mc.cores = cores )

#also fill in the diagonal



     # for (i in 1:p) {
     #   for (j in 1:p) {
     #     for (m in 1:M) {
     #       yps[m] <- -y[[m]][i,j]
     #     }
     #     genL <- dualpath(y = yps,  D = D, verbose = F)
     #       #genL<- genlasso::genlasso (y = yps, X = I, D =D, verbose = T, minlam = 1) # Hier fehlt irgendwie das lambda
     #      for (m in 1:M) {
     #        Z[[m]][i,j] <- genL$beta[,1][m]
     #      }
     #   }
     # }

     # mclapply(combinations, function(combination) {
     #   # obtain the vector y for the generalized LASSO
     #   y <- sapply(B, function(M) M[combination[1], combination[2]])
     #
     #   # apply the generalized LASSO
     #   out <- genlasso(y, diag(1, m), D, minlam = 1)
     #   beta <- coef(out, lambda = 1)$beta
     #
     #   beta[which(abs(beta) <= 10^-12)] <- 0
     #
     #   # update the matrix Z (use that it is symmetric)
     #   for (m in 1:M) {
     #     Z[[m]][combination[1], combination[2]] <<- beta[i]
     #     Z[[m]][combination[2], combination[1]] <<- beta[i]
     #   }
     # })






  # #Update Step Z0
  #  for (m in 1:M) {
  #   H <- Theta[[m]] + U0_old[[m]]
  #   Z0[[m]] <- soft.thresholding(H, c = (lambda/rho))
  # }
  #
  #   #Update Step Z1, Z2
  #
  #
  #
  #
  #   E <- matrix(0,p,p)
  #
  #   for (m in 1:M) {
  #     m1 <- M1[m]
  #     m2 <- M2[m]
  #
  #       delta <-(2*w[m1,m2]/rho)
  #
  #
  #       B <- Theta[[m1]] - Theta[[m2]] + U1[[m1]] - U2[[m2]]
  #       C <- Theta[[m1]] + Theta[[m2]] + U1[[m1]] + U2[[m2]]
  #
  #
  #
  #
  #         if(psi == "L1") {
  #           E <- soft.thresholding(B,c = delta)
  #
  #         }else if(psi == "Lagrangian"){
  #           for (i in 1:p) {
  #             for (j in 1:p) {
  #               E[i,j] <- 1/(1+(2*delta)) * B[i,j]
  #             }
  #           }
  #         }
  #
  #
  #     Z1[[m1]]  <- 0.5 * C + 0.5 * -E
  #     Z2[[m2]]  <- 0.5 * C + 0.5 * E
  #
  #   }


    for (m in 1:M) {
      U[[m]] <- U_old[[m]]  + (Theta[[m]] - Z[[m]])
      #U1[[m]] <- U1_old[[m]]      + (abs(Theta[[m]]) - abs(Z1[[m]]))
      #U2[[m]] <- U2_old[[m]]      + (abs(Theta[[m]]) - abs(Z2[[m]]))
    }






    #save them for later and compare them


     # Stopping Criterea from Danaher
    k <- k + 1
    diff_value = 0
    for(m in 1:M) {diff_value = diff_value + sum(abs(Theta[[m]] - Theta_old[[m]])) / sum(abs(Theta_old[[m]]))}

    cat(sprintf("iteration %d  |  %f\n", k, diff_value))


   }#Ende While, bzw Repeat

#aus danaher, absoluter unterschied von Theta und Z
  diff = 0; for(m in 1:M){diff = diff + sum(abs(Theta[[m]]-Z[[m]]))}

  adj_matrices <- list()


  for (m in 1:M) {
      adj <- Empty
      adj[which(Z[[m]]!=0)] <- 1
      diag(adj) <- 0
      adj[upper.tri(adj, diag = T) == TRUE] <- 0
      #adj <- adj[which(lower.tri(adj))==TRUE]
      adj_matrices[[m]]<- adj
    }




  res <- list(
    Theta = Z,
    adj_matrices = adj_matrices,
    S = S,
    n = n,
    lambda1 = lambda1,
    lambda2 = lambda2,
    rho = rho
  )

  class(res) <- "Multismoothed_Func"

  return(res)

}

