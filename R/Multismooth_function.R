MdS <- function(Y, rho = 1, lambda1 = 2, lambda2= 0.5, w , epsilon = 1e-3, cores = 1,  maxIter = 500){

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

 ################### #Daten zentralisieren
  for (m in 1:M) {
    for (j in 1:p) {
      Y[[m]][, j] = Y[[m]][, j] - mean(Y[[m]][, j])
    }
  }





############## Prepare Constants #####################

  Empty <- matrix(0,p,p)
  Theta <- Z  <- U <- betas <- list()

    S <- list()
    eta <- list()

  for (m in 1:M) {
    S[[m]] <- (t(Y[[m]]) %*% Y[[m]]) / n[m]
    Theta[[m]] <- diag(1/diag(S[[m]]))
    Z[[m]] <- U[[m]] <- Empty    #Start the algorithm with an empty network
    eta[[m]] <- n[m] / rho

  }


    ####### Prepare weight matrix ##########################
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



    row_with_only_zeros <- apply(D, 1, function(row) all(row == 0))
    matrix_D <- D[!row_with_only_zeros, ]



    ### Emergency exit for non converging cases
    diff_value1 <- 0
    diff_value2 <- 0
    diff_value3 <- 0
    diff_value4 <- 0



####################initiate iterationssteps ###########
  k <- 0
  diff_value = 10

  while((k==0) || ((k < maxIter) && diff_value > epsilon))
  {

    Theta_old <- Theta
    Z_old <- Z
    U_old <- U











 ############# #Update Step Theta ###################
    for (m in 1:M) {
      #copied from Danaher
      G = eigen(S[[m]] - rho*Z[[m]]/n[m] + rho*U[[m]]/n[m])
      D = G$values
      Q = G$vectors
      D2 = n[m]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[m]) )
      Theta[[m]] = Q %*% diag(D2) %*% t(Q)


      #formula for theta from hallac
      #maybe this performs worse, idk

    # G <- (eta[[m]]^(-1)) * (A[[m]]) - S[[m]]
    # G_Eigen <- eigen(G)
    # Q <- G_Eigen$vectors
    # D <- diag(G_Eigen$values)
    # V <- D + sqrt(D^2 + ((4*(eta[[m]]^(-1)))*diag(p)))
    # Theta_k <- (1/(2*(eta[[m]]^(-1)))) * (Q%*%  V  %*% t(Q))
    # Theta[[m]] <- Theta_k
  }

   #########################Update Step Z #########################


    y <- list()

      for (m in 1:M) {
       y[[m]] <- Theta[[m]] + U_old[[m]]
      }


     for (m in 1:M) {
       diag(Z[[m]]) <- diag(y[[m]])
     }


     ### Use generalized Lasso on all possible Edges
combinations <- combn(1:p, 2, simplify = FALSE)

for (all_edges in 1:length(combinations)) {
  i <- combinations[[all_edges]][1]
  j <- combinations[[all_edges]][2]

yps <- c()
  for (m in 1:M) {
  yps[m] <- y[[m]][[i,j]]
  }

  genL <- genlasso(y = yps, diag(1, m), matrix_D, minlam = 1, verbose = F)
  #genL <- genlasso(y = yps, diag(1, m), D = NewWeights, minlam = 1, verbose = F)
  out <- coef(genL, lambda = 1)$beta

  #out[which(abs(out) <= 10^-12)] <- 0

  for (m in 1:M) {
    Z[[m]][i,j] <- Z[[m]][j,i]  <- out[m]
  }

  }






    for (m in 1:M) {
      U[[m]] <- U_old[[m]]  + (Theta[[m]] - Z[[m]])

    }









     # Stopping Criterea from Danaher
    k <- k + 1
    diff_value = 0

    for(m in 1:M) {
      diff_value = diff_value + sum(abs(Theta[[m]] - Theta_old[[m]])) / sum(abs(Theta_old[[m]]))
      }

    cat(sprintf("iteration %d  |  %f\n", k, diff_value))

    diff_value2 <- diff_value1
    diff_value1 <- diff_value


exiter <- diff_value2 - diff_value1

    if (k == maxIter) {
      warning("We reached the max Iteration, without meeting coversion criteria, last value was ",
              diff_value)
    }


    if (abs(exiter) < 1e-7 & diff_value > epsilon) {
      k <- maxIter
      warning("Did Not Converge for lambda1 = ", lambda1, " lambda2 = " ,
              lambda2 , "||   last Iterationvalue at ", diff_value)
    }

   }#Ende While, bzw Repeat


#Store the adjacency matrices
  adj_matrices <- list()


  for (m in 1:M) {
     Z[[m]][which(abs(Z[[m]]) < 1e-12)] <- 0

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

