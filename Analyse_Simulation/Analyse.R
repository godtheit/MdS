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

a11 <- c(0,0,0,0,0,0,0,0,0) #1
a12 <- c(8,0,0,0,0,0,0,0,0) #2
a13 <- c(6,8,0,0,0,0,0,0,0) #3
a21 <- c(8,6,4,0,0,0,0,0,0) #4
a22 <- c(6,8,6,8,0,0,0,0,0) #5
a23 <- c(4,6,8,6,8,0,0,0,0) #6
a31 <- c(6,4,2,8,6,4,0,0,0) #7
a32 <- c(4,6,4,6,8,6,8,0,0) #8
a33 <- c(2,4,6,4,6,8,6,8,0) #9

w_distance2 <- rbind(a11,a12, a13, a21, a22, a23, a31, a32, a33)
w_distance <- w_distance2/10
colnames(w_distance) <- rownames(w_distance)


g <- matrix(1 , 9, 9)
g[which(upper.tri(w, diag = TRUE)==TRUE)] <- 0

lam1 <- seq(0.1,30, length.out = 20)
lam2 <- seq(0.01,2, length.out = 8)



List_Estimations_glasso <- list()
List_Estimations_FusedLasso <- list()
List__Estimations_Mds <- list()


for (t in 1:1) {#:length(BigListof_samples)) {
current_network <- filefortest[[2]]
List_Estimations_glasso[[t]] <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2, method = "glasso")

List_Estimations_FusedLasso[[1]] <-  procedure(Y = current_network,
                                               lam1 = lam1, lam2 = lam2, method = "fused")
List_Estimations_FusedLasso[[2]] <- procedure(Y = current_network, lam1 = lam1, lam2 = lam2,
                                              method = "fused", criteria = "Bic")
List__Estimations_Mds[[1]]  <- procedure(Y = current_network, lam1 = lam1,
                                         lam2 = lam2, method = "MdS", w = w)
List__Estimations_Mds[[2]]  <- procedure(Y = current_network, lam1 = lam1,
                                         lam2 = lam2, method = "MdS", w = w, criteria = "Bic")
}
#hurray

vergleich1 <- MdS(Y = current_network, lambda1 = 1, lambda2 = 1, w = w)
vergleich2 <- JGL(Y = current_network, lambda1 = 1, lambda2 = 1, penalty = "fused",
                  return.whole.theta = T)

data(example.data)
str(example.data)
## run FGL:
fgl.results = JGL(Y=example.data,penalty="fused",lambda1=3,lambda2=10)
str(fgl.results)
print.jgl(fgl.results)



u <- matrix(1 , 2, 2)
u[which(upper.tri(u, diag = TRUE)==TRUE)] <- 0

MdS.results = MdS(Y = example.data, lambda1 = 0.26, lambda2 = 0.1, w = u)
