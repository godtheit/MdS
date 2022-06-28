rm(list = ls())

#library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

#setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_simulation")
setwd("H:/Documents/MdS/Analyse_simulation")

load("allsamples.RData")
load("alladjacs.RData")





#filefortest <- readRDS("test-dataset.RDS")

#dims <-c(3,3)

#w <- matrix(1 ,prod(dims),prod(dims))
#w[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0

#wat2 <- MdS(Y = filefortest[[2]],w = w, dims = c(1,2), lambda1 = 10, lambda2 = 1, rho = 1, maxIter = 20, epsilon = 1e-4)



#bigtestsample <- BigListof_samples[[1]]

#wat <- MdS(Y = bigtestsample, w = w, dims = c(3,3), lambda1 = 1,lambda2 = 1, rho = 1, maxIter = 2000)

lam1 <- seq(0.1,2, length.out = 15)
lam2 <- seq(0.1,2, length.out = 8)

List_Estimations_glasso <- list()
List_Estimations_FusedLasso <- list()
List__Estimations_Mds <- list()

for (t in 1:length(BigListof_samples)) {
  estimation_tmp <- list()


  for (z in 1:sum(dims)) {
    #sample <- BigListof_samples[[t]][[z]]
    #covariance <- cov(sample)
    #n <- length(sample[,1])

    #glasso(s = covariance, rho = lam1)
    x <- huge(sample,  method = "glasso")
    y <- huge.select(x, criterion = "ebic")

    estimation_tmp[[z]] <- y$opt.icov

  }
List_Estimations_glasso[[t]] <- estimation_tmp



counter <- 1
JGL_tmp <- list()
MdS_tmp <- list()
#Mache für alle lambda1 die 2 Methoden und speichere diese in Liste
for (a in lam1) {
  JGL_tmp[[counter]] <- JGL(BigListof_samples[[t]], lambda1 = a, lambda2 = 1, penalty = "fused", return.whole.theta = T)
  MdS_tmp[[counter]] <-  MdS(BigListof_samples[[t]], lambda1 = a, lambda2 = 1 , w = w, dims = c(3,3))
  counter <- counter + 1
}


#Suche für alle Ergebnisse das Kriterium heraus
JGL_crit1_list <- matrix(0,length(lam1),1)
MdS_crit1_list <- matrix(0,length(lam1),1)
  for (c in length(lam1)) {
    JGl_crit1_list[c] <- lambda_selector(JGL_tmp[[c]])
    MdS_crit1_list[c] <- lambda_selector(Theta = MdS_tmp[[c]]$Theta,
                                           n = MdS_tmp[[c]]$n,
                                           lam1 = MdS_tmp[[c]]$lamda1,
                                           lam2 =MdS_tmp[[c]]$lambda2)
  }

#Suche das größte Kriterium
lam1_JGL <- lam1[which(max(JGL_crit1_list))]
lam1_MdS <- lam1[which(max(MdS_crit1_list))]




#Analog für lam2 aber dieses mal ist lam1 bekannt

counter <- 1
JGL_tmp <- list()
MdS_tmp <- list()
for (b in lam2) {
  JGL_tmp[[counter]] <- JGL(BigListof_samples[[t]], lambda1 = lam1_JGL, lambda2 = b, penalty = "fused", return.whole.theta = T)
  MdS_tmp[[counter]] <-  MdS(BigListof_samples[[t]], lambda1 = lam1_MdS, lambda2 = b , w = w, dims = c(3,3))
  counter <- counter + 1
}




JGL_crit2_list <- matrix(0,length(lam2),1)
MdS_crit2_list <- matrix(0,length(lam2),1)
for (c in length(lam2)) {
  JGl_crit2_list[c] <- lambda_selector(JGL_tmp[[c]])
  MdS_crit2_list[c] <- lambda_selector(Theta = MdS_tmp[[c]]$Theta,
                                       n = MdS_tmp[[c]]$n,
                                       lam1 = MdS_tmp[[c]]$lamda1,
                                       lam2 =MdS_tmp[[c]]$lambda2)
}


lam2_JGL <- lam2[which(max(JGL_crit1_list))]
lam2_MdS <- lam2[which(max(MdS_crit1_list))]

List_Estimations_FusedLasso[[t]]<-
List__Estimations_Mds[[t]][[counter]] <-

}


