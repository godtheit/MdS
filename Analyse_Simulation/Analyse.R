rm(list = ls())

library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

#setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_simulation")
setwd("H:/Documents/MdS/Analyse_simulation")
path <- "H:/Documents/MdS/Analyse_simulation"


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

lam1 <- seq(0.02,2, length.out = 60)
lam2 <- seq(0.1,2, length.out = 8)



List_Estimations_glasso <- list()
List_Estimations_FusedLasso <- list()
List__Estimations_Mds <- list()


for (t in 1:1) {#:length(BigListof_samples)) {
current_network <- filefortest[[2]]
List_Estimations_glasso[[t]] <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "glasso")
List_Estimations_FusedLasso[[t]] <-  procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "fused")
List__Estimations_Mds[[t]]  <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "MdS", w = w)
}
#hurray


test <- sim_graphs(
p = 150, perc1 = 0.05,
perc2 = 0.05,
length_var1 = 3,
length_var2 = 3,
observes = 52)


