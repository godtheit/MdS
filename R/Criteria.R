procedure <- function(Y, lam1, lam2, method, w = FALSE, criteria = "Aic"){

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
  obs <- Y
  #Daten zentralisieren
  for (m in 1:M) {
    for (j in 1:p) {
      Y[[m]][, j] = Y[[m]][, j] - mean(Y[[m]][, j])
    }
  }

  S <- list()
  for (m in 1:M) {
    S[[m]] <- (t(Y[[m]]) %*% Y[[m]]) / n[m]
  }




  if (method == "MdS") {
    if(is.logical(w) == TRUE){
      stop("You need to add a weight matrix if you use MdS")
    }

    #Lists and matrices to store stuff later
    MdS_tmp1 <- list()
    MdS_crit1_list <- matrix(0,length(lam1),1)
    MdS_tmp2 <- list()
    MdS_crit2_list <- matrix(0,length(lam2),1)

    #Do MdS for all lamba1
    counter <- 1
    for (a in lam1) {
      MdS_tmp1[[counter]] <-  MdS(Y = obs, lambda1 = a, lambda2 = 0.2 , w = w)
      counter <- counter + 1
    }

if(criteria == "Aic"){
    for (c in 1:length(lam1)) {

      MdS_crit1_list[c] <- mds_aic(Theta = MdS_tmp1[[c]]$Theta,
                                   S = S,
                                   n = n)
    }
} else if(criteria == "Bic"){
  for (c in 1:length(lam1)) {


  MdS_crit1_list[c] <- mds_bic(Theta = MdS_tmp1[[c]]$Theta,
                               S = S,
                               n = n,
                               p = p)
  }
  }
    lam1_MdS <- lam1[which(MdS_crit1_list == min(MdS_crit1_list))]

    L1 <- length(lam1_MdS)
    if(length(lam1_MdS) > 1){
      lam1_MdS <- lam1_MdS[1]
      multiple_lambda1 <- TRUE
    } else {multiple_lambda1 <- FALSE}




    counter <- 1
    for (b in lam2) {
      MdS_tmp2[[counter]] <-  MdS(Y = obs, lambda1 = lam1_MdS, lambda2 = b , w = w)
      counter <- counter + 1
    }

if(criteria == "Aic"){
    for (c in 1:length(lam2)) {
      MdS_crit2_list[c] <- mds_aic(Theta = MdS_tmp2[[c]]$Theta,
                                   S = S,
                                   n = n)
    }
} else if(criteria == "Bic"){
  for (c  in 1:length(lam2)) {
    MdS_crit2_list[c] <- mds_bic(Theta = MdS_tmp2[[c]]$Theta,
                                 S = S,
                                 n = n,
                                 p = p)
  }
}
    lam2_MdS <- lam2[which((MdS_crit2_list) == min(MdS_crit2_list))]


    L2 <- length(lam2_MdS)


    if(length(lam2_MdS) > 1){
    lam2_MdS <- lam2_MdS[1]
    multiple_lambda2 <- TRUE
  } else {multiple_lambda2 <- FALSE}


    estimation <- MdS_tmp2[[which(lam2 == lam2_MdS)]]


    estimation$multiples <- list(multiple_lambda1, multiple_lambda2, L1, L2)

    return(list(estimation = estimation,
                list_lam1 = MdS_crit1_list,
                list_lam2 = MdS_crit2_list)
          )
  }













  if(method == "fused"){
    estimation <- list()
    JGL_tmp1 <- list()
    JGL_crit1_list <- matrix(0,length(lam1),1)
    JGL_tmp2 <- list()
    JGL_crit2_list <- matrix(0,length(lam2),1)


    counter <- 1
    for (a in lam1) {
      JGL_tmp1[[counter]] <- JGL(Y = obs, lambda1 = a, lambda2 = 0.2, penalty = "fused", return.whole.theta = T)
      counter <- counter + 1
    }






    #Suche für alle Ergebnisse das Kriterium heraus
if (criteria == "Aic"){
    for (c in 1:length(lam1)) {
      JGL_crit1_list[c] <-       mds_aic(Theta = JGL_tmp1[[c]]$theta,
                                         S = S,
                                         n = n)
    }
} else if (criteria == "Bic"){
  for (c in 1:length(lam1)) {
    JGL_crit1_list[c] <-       mds_bic(Theta = JGL_tmp1[[c]]$theta,
                                       S = S,
                                       n = n,
                                       p = p)
  }
}

    #Suche das größte Kriterium
    lam1_JGL <- lam1[which(JGL_crit1_list == min(JGL_crit1_list))]
    L1 <- length(lam1_JGL)
    if(length(lam1_JGL) > 1){
      lam1_JGL <- lam1_JGL[1]
      multiple_lambda1 <- TRUE
    } else {multiple_lambda1 <- FALSE}


    #Analog für lam2 aber dieses mal ist lam1 bekannt
    counter <- 1
    for (b in lam2) {
      JGL_tmp2[[counter]] <- JGL(Y = obs, lambda1 = lam1_JGL, lambda2 = b, penalty = "fused", return.whole.theta = T)
      counter <- counter + 1
    }


if (criteria == "Aic") {

    for (c in 1:length(lam2)) {
      JGL_crit2_list[c] <- mds_aic(Theta = JGL_tmp2[[c]]$theta,
                                   S = S,
                                   n = n)

    }
} else if (criteria == "Bic"){

  for (c in 1:length(lam2)) {
    JGL_crit2_list[c] <- mds_bic(Theta = JGL_tmp2[[c]]$theta,
                                 S = S,
                                 n = n,
                                 p = p)

  }
}

    lam2_JGL <- lam2[which((JGL_crit2_list) == min(JGL_crit2_list))]
    L2 <- length(lam2_JGL)
    if(length(lam2_JGL) > 1){
      lam2_JGL <- lam2_JGL[1]
      multiple_lambda2 <- TRUE
    } else {multiple_lambda2 <- FALSE}

    estimation[[1]] <- JGL_tmp2[[which(lam2 == lam2_JGL)]]
    estimation[[2]] <- c(lam1_JGL,lam2_JGL)

    estimation$multiples <- list(multiple_lambda1, multiple_lambda2, L1, L2)
    return(list(estimation = estimation,
                list_lam1 = JGL_crit1_list,
                list_lam2 = JGL_crit2_list)
            )
  }





  if (method == "glasso"){
    estimation <- list()
    estimation_tmp <- list()
    for (z in 1:M) {
      sample <- obs[[z]]
      x <- huge(sample,  method = "glasso")
      y <- huge.select(x, criterion = "ebic")
      estimation_tmp[[1]] <- y$opt.icov
      estimation_tmp[[2]] <- y$opt.lambda
      estimation[[z]] <- estimation_tmp
    }
    return(estimation)
  }


}
