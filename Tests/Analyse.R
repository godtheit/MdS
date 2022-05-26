rm(list = ls())

library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

#setwd("C:/Users/Arne/Dropbox/Masterarbeit/Multidimensional-Smoothing/R")
#setwd("C:/Users/arne2/Dropbox/Masterarbeit/Multidimensional-Smoothing/R")
#load("alladjacs.RData")
#load("allcovars.RData")
#load("allgraphs.RData")
#source("Berechnung_Z.R")

#load("allsamples.RData")

testsample1 <- list()
testsample2 <- list()

for (m in 1:9) {
 testsample <- huge.generator(n = 10, d = 50, graph ="scale-free")
 testsample1[[m]] <- testsample$data
 testsample2[[m]] <- testsample$omega
}

#fgl.results = JGL(Y=testsample1, penalty="fused",lambda1=.25,lambda2=.1)

#source("Multismooth_function.R")
dims <-c(3,3)

w <- matrix(0.5,prod(dims),prod(dims))
w[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0





wat <- MdS(Y = testsample1, w = w, dims = c(3,3), lambda1 = 0.2,lambda2 = 0.3, rho = 1, psi = "L1", maxIter = 200)

wat2 <- JGL(testsample1, lambda1 = 0.1, lambda2 = 0.2)
y <- vector(mode = "double", length = 9)
for (m in 1:9) {
  y[m] <- testsample1[[m]][1,2]
}
wat3 <- genlasso(y = y, D = diag(9), minlam = 1)





t <- 1
z <- 1
List_Estimations <- list()

for (t in 1:length(BigListof_samples)) {
  estimation_tmp <- list()

  for (z in 1:dims) {
    tmp <- BigListof_samples[[t]][[z]]
    x <- huge(tmp, method = "glasso")
    y <- huge.select(x, criterion = "stars")
    test <- y$refit
  estimation_tmp[[z]] <- test
  }
List_Estimations[[t]] <- estimation_tmp
}









Samples <- BigListof_samples[[1]]
Samp_11 <- Samples[[1]]
glasso()
help(glasso)
H <- huge.glasso(Samp_11)
