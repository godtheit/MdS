rm(list = ls())

library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

#setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Analyse_simulation")
setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#setwd("H:/Documents/MdS/Analyse_simulation")

load("allsamples.RData")
load("alladjacs.RData")





filefortest <- readRDS("test-dataset.RDS")

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

lam1 <- seq(0.1,2, length.out = 15)
lam2 <- seq(0.1,2, length.out = 8)

List_Estimations_glasso <- list()
List_Estimations_FusedLasso <- list()
List__Estimations_Mds <- list()
dims <- c(1,2)

for (t in 1:1) {#:length(BigListof_samples)) {
  estimation_tmp <- list()


  for (z in 1:prod(dims)) {
    sample <- filefortest[[2]][[z]]
    #sample <- BigListof_samples[[t]][[z]]
    #covariance <- cov(sample)
    #n <- length(sample[,1])

    #glasso(s = covariance, rho = lam1)
    x <- huge(sample,  method = "glasso")
    y <- huge.select(x, criterion = "ebic")

    estimation_tmp[[z]] <- y$opt.icov

  }
List_Estimations_glasso[[t]] <- estimation_tmp


#Setup the network for each t

current_network <- filefortest[[2]]

p = dim(current_network[[1]])[2]
K = length(current_network)
n = rep(0, K)
for (k in 1:prod(dims)) {
  n[k] = dim(current_network[[k]])[1]
}

for (k in 1:K) {
  for (j in 1:p) {
    current_network[[k]][, j] = current_network[[k]][, j] - mean(current_network[[k]][, j])
  }
}
S <- list()
for (k in 1:K) {
  S[[k]] <- (t(current_network[[k]]) %*% current_network[[k]]) / n[k]
}


counter <- 1
JGL_tmp <- list()
MdS_tmp <- list()
#Mache für alle lambda1 die 2 Methoden und speichere diese in Liste
for (a in lam1) {
  #JGL_tmp[[counter]] <- JGL(BigListof_samples[[t]], lambda1 = a, lambda2 = 1, penalty = "fused", return.whole.theta = T)
  #MdS_tmp[[counter]] <-  MdS(BigListof_samples[[t]], lambda1 = a, lambda2 = 1 , w = w, dims = c(3,3))
  JGL_tmp[[counter]] <- JGL(Y = current_network, lambda1 = 1, lambda2 = 1, penalty = "fused", return.whole.theta = T)
  MdS_tmp[[counter]] <-  MdS(Y = current_network, lambda1 = a, lambda2 = 1 , w = w, dims = c(1,2), maxIter = 20)
  counter <- counter + 1
}

#Suche für alle Ergebnisse das Kriterium heraus
JGL_crit1_list <- matrix(0,length(lam1),1)
MdS_crit1_list <- matrix(0,length(lam1),1)
  for (c in 1:length(lam1)) {
    JGL_crit1_list[c] <-       mds_aic(Theta = JGL_tmp[[c]]$theta,
                                       S = S,
                                       n = n,
                                       lambda1 = lam1[c],
                                       lambda2 = 1)

    MdS_crit1_list[c] <- mds_aic(Theta = MdS_tmp[[c]]$Theta,
                                 S = S,
                                 n = n,
                                 lambda1 = lam1[c],
                                 lambda2 = 1)
  }

#Suche das größte Kriterium
lam1_JGL <- lam1[which(JGL_crit1_list == max(JGL_crit1_list))]
lam1_JGL <- lam1_JGL[1]
lam1_MdS <- lam1[which(abs(MdS_crit1_list) == max(abs(MdS_crit1_list)))] #too small?




#Analog für lam2 aber dieses mal ist lam1 bekannt

counter <- 1
JGL_tmp <- list()
MdS_tmp <- list()
for (b in lam2) {
  #JGL_tmp[[counter]] <- JGL(current_network, lambda1 = lam1_JGL, lambda2 = b, penalty = "fused", return.whole.theta = T)
  #MdS_tmp[[counter]] <-  MdS(current_network, lambda1 = lam1_MdS, lambda2 = b , w = w, dims = c(3,3))
  JGL_tmp[[counter]] <- JGL(Y = current_network, lambda1 = lam1_JGL, lambda2 = b, penalty = "fused", return.whole.theta = T)
  MdS_tmp[[counter]] <-  MdS(Y = current_network, lambda1 = lam1_MdS, lambda2 = b , w = w, dims = c(1,2), maxIter = 20)
  counter <- counter + 1
}




JGL_crit2_list <- matrix(0,length(lam2),1)
MdS_crit2_list <- matrix(0,length(lam2),1)
for (c in 1:length(lam2)) {
  JGL_crit2_list[c] <- mds_aic(Theta = JGL_tmp[[c]]$theta,
                               S=S,
                               n = n,
                               lambda1 = lam1_JGL,
                               lambda2 = lam2[c])

  MdS_crit2_list[c] <- mds_aic(Theta = MdS_tmp[[c]]$Theta,
                              S = S,
                              n = n,
                              lambda1 = lam1_MdS,
                              lambda2 = lam2[c])
}


lam2_JGL <- lam2[which(abs(JGL_crit2_list) == max(abs(JGL_crit2_list)))]
lam2_JGL <- lam2_JGL[1]
lam2_MdS <- lam1[which(abs(MdS_crit2_list) == max(abs(MdS_crit2_list)))] #too small?


List_Estimations_FusedLasso[[t]]<- JGL_tmp[[which(lam2 == lam2_JGL)]]
List__Estimations_Mds[[t]] <- MdS_tmp[[which(lam2 == lam2_MdS)]]
}


