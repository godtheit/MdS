rm(list = ls())

library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

#setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Analyse_simulation")
setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#load("alladjacs.RData")
#load("allcovars.RData")
#load("allgraphs.RData")
#source("Berechnung_Z.R")

load("allsamples.RData")
load("alladjacs.RData")
#testsample1 <- list()
#testsample2 <- list()
set.seed(20)
#for (m in 1:2) {
# testsample <- huge.generator(n = 1000, d = 5, graph ="scale-free")
# testsample1[[m]] <- testsample$data
# testsample2[[m]] <- testsample$omega
#}

#fgl.results = JGL(Y=testsample1, penalty="fused",lambda1=.25,lambda2=.1)
filefortest <- readRDS("test-dataset.RDS")
#source("Multismooth_function.R")
dims <-c(3,3)

w <- matrix(1 ,prod(dims),prod(dims))
w[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0

wat2 <- MdS(Y = filefortest[[2]],w = w, dims = c(1,2), lambda1 = 10, lambda2 = 1, rho = 1, maxIter = 20, epsilon = 1e-4)


wat2$adj_matrices[[1]]





bigtestsample <- BigListof_samples[[1]]

wat <- MdS(Y = bigtestsample, w = w, dims = c(3,3), lambda1 = 1,lambda2 = 1, rho = 1, maxIter = 2000)



List_Estimations_glasso <- list()
List_Estimations_FusedLasso <- list()
List__Estimations_Mds <- list()

for (t in 1){            ##length(BigListof_samples)) {
  estimation_tmp <- list()

  for (z in 1:sum(dims)) {
    tmp <- BigListof_samples[[t]][[z]]
    x <- huge(tmp,  method = "glasso")
    y <- huge.select(x, criterion = "ric")

    estimation_tmp[[z]] <- y$opt.icov

  }
List_Estimations_glasso[[t]] <- estimation_tmp
List_Estimations_FusedLasso[[t]] <- JGL(BigListof_samples[[t]], lambda1 = 5, lambda2 = 1, penalty = "fused", return.whole.theta = T)
List__Estimations_Mds[[t]] <- MdS(BigListof_samples[[t]], lambda1 = 5, lambda2 = 2, w = w, dims = c(3,3))
}

mds1 <- MdS(filefortest[[2]],dims = c(1,2), lambda1 = 3, lambda2 = 1, rho = 1, w = w)
jgl1 <- JGL(filefortest[[2]], lambda1 = 0.1, lambda2 = 0.1, penalty = "fused", return.whole.theta = T)
