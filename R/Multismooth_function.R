MdS <- function(Y, rho = 1, lambda1 = 0.1,lambda2= 0.1, w ,psi = "L1", dims, epsilon = 1e-3, maxIter = 1000){

  library(DescTools)

  p = dim(Y[[1]])[2]
  M = length(Y)
  n = rep(0, M)
  for (m in 1:M) {
    n[m] = dim(Y[[m]])[1]
  }
  if (length(dimnames(Y[[1]])[[2]]) == 0) {
    for (m in 1:M) {
      dimnames(Y[[m]])[[2]] = paste("V", 1:p, sep = "")
    }
  }

  #Daten Normalisieren
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
  Ones <- matrix (1,p,p)
  Eps <- matrix(epsilon,p,p)
  Theta <- Z  <- U <- betas <- list()



  for (m in 1:M) {
    Theta[[m]] <- Z[[m]] <- U[[m]] <- betas[[m]] <- Empty

  }

  k <- 0
  verschiebung <- 0


  #_old ist immer die k-te Iteration, _old2 die k-1 te
  # ohne index k+1
  #exitStatement <- matrix(FALSE, M, 1)


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



  cat(" Iterationstep ", k)

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
    S[[m]] <- (t(Y[[m]]) %*% Y[[m]]) / n[m]
    A[[m]] <- ((Z_old[[m]] - U_old[[m]]) )
    eta[[m]] <- n[m] / rho    #n probably not the same for all m, maybe change this later

  }


    for (m in 1:M) {
    G <- (eta[[m]]^(-1)) * ((A[[m]] + t(A[[m]]))/2 ) - S[[m]]
    G_Eigen <- eigen(G)

    Q <- G_Eigen$vectors
    D <- diag(G_Eigen$values)

    V <- D + sqrt(D^2 + ((4*(eta[[m]]^(-1)))*diag(p)))

    Theta_k <- (1/(2*(eta[[m]]^(-1)))) * (Q%*%  V  %*% t(Q))
    Theta[[m]] <- Theta_k
  }

    y <- list()
     #Update Step Z
      for (m in 1:M) {
       y[[m]] <- Theta[[m]] + U_old[[m]]
     }
   eta_1 <- lambda1/rho
   eta_2 <- lambda2/rho


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



     yps <- betas <- vector(mode = "double", length = M)
     I <- diag(M)
     counter <- 1
     for (i in 1:p) {
       for (j in 1:p) {
         for (m in 1:M) {
           yps[m] <- -y[[m]][i,j]
           betas[m] <- Z_old[[m]][i,j]

         }
         betas <- as.matrix(betas)



         genL <- dualpath(y = yps,  D = D, verbose = TRUE)
          #genL <- genlasso(y = yps, X = I, D =D, verbose = T, minlam = 1) # Hier fehlt irgendwie das lambda
          for (m in 1:M) {
            Z[[m]][i,j] <- genL$beta[,1][m]
          }
       }
     }






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

    k <- k + 1
    diff_value = 0
    for(m in 1:M) {diff_value = diff_value + sum(abs(Theta[[m]] - Theta_old[[m]])) / sum(abs(Theta_old[[m]]))}


   }#Ende While, bzw Repeat





  return(Theta)
}

