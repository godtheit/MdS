rm(list = ls())
setwd("H:/Documents/MdS/Analyse_Simulation/")
source("get_rates.R")






MdS_list  <- JGL_list <- Glasso_list <- list()
Results_as_frame <- c()
Durchschnitte <- as.data.frame(matrix(0,320,7))
colnames(Durchschnitte) <- c("Hamming-Distance", "Precision", "Recall", "Aic", "lambda1", "lambda2",  "id")


############## plotting MdS
for(k in 1:320){
  Daten <- readRDS(paste0("H:/Documents/MdS/Analyse_Simulation/MdS_Simulation/",k, ".RDS")  )
  TrueAdjs <- Daten$trueAdj
  MdS_est            <- Daten$adj_matrices

  theta <- Daten$Theta
  S <- Daten$S
  MdS  <- JGL <- Glasso <- as.data.frame(matrix(0,9,3))
  colnames(MdS) <- names(JGL) <- names(Glasso) <- c("Hamming-Distance", "Precision", "Recall")
  MdS$Aic <- 0

  for (m in 1:9) {

    MdS[m,1:3] <- get_rates(TrueAdjs[[m]], MdS_est[[m]])

    #JGL[m,] <- get_rates(TrAdj, JGL_est[[m]])
    #Glasso[m,] <- get_rates(TrAdj, Glasso_est[[m]])

    MdS$Aic <- mds_aic(Theta = Daten$Theta, S = Daten$S, n = Daten$n)
  }
 MdS$lambda1 <- Daten$lambda1
 MdS$lambda2 <- Daten$lambda2
 MdS$id <- k

for (r in 1:7) {
  Durchschnitte[k,r] <- mean(MdS[,r])
}


 Results_as_frame <- rbind(Results_as_frame, MdS)
  MdS_list[[k]] <- MdS



  #JGL_list[[k]] <- JGL
  #Glasso_list[[k]] <- Glasso
}

############## plots
path <- "H:/Documents/MdS/Ergebnisse_Sim/"
for (a in 11:14) {
  for (b in c(1.5,2.25,3,3.75)) {
    plotter <- subset(Results_as_frame, Results_as_frame$lambda1 == a & Results_as_frame$lambda2 == b)
    plotter2 <- subset(Durchschnitte, Durchschnitte$lambda1 == a & Durchschnitte$lambda2 == b)

    png(filename = paste0(path, "_lam1_", a ,"_lam2_", b, ".png"))
    plot(y = plotter$Precision, x = plotter$Recall)
    dev.off()

    plot(y = plotter2$Precision, x = plotter2$Recall)

  }
}






