rm(list = ls())
library(MdS2)
#setwd("H:/Documents/MdS/Analyse_Simulation/")
setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/")
source("get_rates.R")





JGL_list  <- list()
Results_as_frame <- c()
Durchschnitte <- as.data.frame(matrix(0,320,7))
colnames(Durchschnitte) <- c("Hamming-Distance", "Precision", "Recall",  "lambda1", "lambda2",  "id")


############## plotting JGL
for(k in 1:320){
  #Daten <- readRDS(paste0("H:/Documents/MdS/Analyse_Simulation/MdS_Simulation/",k, ".RDS")  )
  Daten <- readRDS(paste0("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/JGL_Simulation2/", k, ".RDS"))
  TrueAdjs <- Daten$trueAdj
  MdS_est    <- Daten$adj_matrices

  theta <- Daten$Theta

  JGL  <- as.data.frame(matrix(0,9,3))
  colnames(JGL<- c("Hamming-Distance", "Precision", "Recall"))


  for (m in 1:9) {

    JGL[m,1:3] <- get_rates(TrueAdjs[[m]], JGL_est[[m]])

    #JGL[m,] <- get_rates(TrAdj, JGL_est[[m]])
    #Glasso[m,] <- get_rates(TrAdj, Glasso_est[[m]])


  }
  JGL$lambda1 <- Daten$lambda1
  JGL$lambda2 <- Daten$lambda2
  JGL$id <- k

  for (r in 1:7) {
    Durchschnitte[k,r] <- mean(MdS[,r])
  }


  Results_as_frame <- rbind(Results_as_frame, JGL)
  JGL_list[[k]] <- JGL



}
