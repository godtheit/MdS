rm(list = ls())

library(huge)
library(glasso)
library(JGL)
library(genlasso)

library(MdS2)

setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_simulation")
#setwd("H:/Documents/MdS/Analyse_simulation")
path <- "H:/Documents/MdS/Analyse_simulation"


load("allsamples.RData")
load("alladjacs.RData")




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

data <- readRDS("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Daten/Final_Gennetwork2.RDS")

Adjs <- data$adj_matrices

for (m in 1:9) {
  Adjs[[m]][which(upper.tri( Adjs[[m]] ) == TRUE)] <- t(Adjs[[m]][which(lower.tri( Adjs[[m]] ) == TRUE)])
}


