rm(list = ls())
library(MdS2)
#setwd("H:/Documents/MdS/Analyse_Simulation/")
setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/")
source("get_rates.R")





MdS_list  <- list()
Results_as_frame <- c()
Durchschnitte <- as.data.frame(matrix(0,320,7))
colnames(Durchschnitte) <- c("Hamming-Distance", "Precision", "Recall", "Aic", "lambda1", "lambda2",  "id")


############## plotting MdS
for(k in 1:320){
  #Daten <- readRDS(paste0("H:/Documents/MdS/Analyse_Simulation/MdS_Simulation/",k, ".RDS")  )
  Daten <- readRDS(paste0("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/MdS_Simulation/", k, ".RDS"))
  TrueAdjs <- Daten$trueAdj
  MdS_est    <- Daten$adj_matrices

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



}

## create glasso's fits

source("get_rates.R")
Glasso_list  <- list()
Results_as_frame_glasso <- c()
Results_as_frame_glasso_oracle <- c()
Durchschnitte_glasso <- as.data.frame(matrix(0,20,5))
Durchschnitte_oracle <- as.data.frame(matrix(0,180,5))
colnames(Durchschnitte_glasso) <- colnames(Durchschnitte_oracle) <- c("Hamming-Distance", "Precision", "Recall",  "lambda", "id")


lambdas <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

counter <- 1
for (k in 1:20) {
  daten_glasso <- readRDS(paste0("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Cluster/registries/MdS_registry/results/",k, ".rds"))
  #automatic_glasso <- list()

  #Glasso <- as.data.frame(matrix(0,9,4))
  #colnames(Glasso) <- c("Hamming-Distance", "Precision", "Recall", "lambda")
  GesamtesGlasso <- c()
  #colnames(GesamtesGlasso) <- c("Hamming-Distance", "Precision", "Recall", "lambda", "GraphId")
  for (m in 1:9) {
    #automatic_glasso <- huge.select(daten_glasso[[m]], criterion =  "ebic")
    #opt_glasso <- automatic_glasso$opt.icov


    #opt_glasso[which(upper.tri(opt_glasso) == TRUE)]
    #opt_glasso[which(abs(opt_glasso) > 0)] <- 1

    #Glasso[m,1:3] <- get_rates(daten_glasso$trueAdj[[m]], opt_glasso)
    #Glasso$lambda[m] <- automatic_glasso$opt.lambda


    Glasso_Oracle <- as.data.frame(matrix(0,9,5))
    colnames(Glasso_Oracle) <- c("Hamming-Distance", "Precision", "Recall", "lambda", "GraphId")

    for (x in 1:9) {
      glasso_tmp <- daten_glasso[[m]]$icov[[x+1]]
      glasso_tmp[which(upper.tri(glasso_tmp) == TRUE)]
      glasso_tmp[which(abs(glasso_tmp) > 0)] <- 1

      Glasso_Oracle[x,1:3] <- get_rates(daten_glasso$trueAdj[[m]], glasso_tmp) # hier sind die 9 lambdavalue gespeichert
      Glasso_Oracle$lambda[x] <- daten_glasso[[m]]$lambda[[x+1]]
    }
    Glasso_Oracle$GraphId <- m
    GesamtesGlasso <- rbind(GesamtesGlasso, Glasso_Oracle)
  }

  GesamtesGlasso$id <- k
  #Glasso$id <- k

  for (L in lambdas) {
    subdata_glasso <- subset(GesamtesGlasso, GesamtesGlasso$lambda == L )
    for (r in 1:4) {
      Durchschnitte_oracle[counter,r] <- mean(subdata_glasso[,r])
    }
    Durchschnitte_oracle$id[counter] <- subdata_glasso$id[1]
    counter <- counter + 1
  }





 # for (r in 1:5) {
 #   Durchschnitte_glasso[k,r] <- mean(Glasso[,r])
 # }


  Results_as_frame_glasso <- rbind(Results_as_frame_glasso, GesamtesGlasso)

}

Results_as_frame_JGL <- c()
Durchschnitte_JGL <- as.data.frame(matrix(0,320,7))
colnames(Durchschnitte_JGL) <- c("Hamming-Distance", "Precision", "Recall", "lambda1", "lambda2",  "id")

for (k in 41:60) {
  data_jgl <- readRDS(file = paste0("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation/JGL_Simulation/", k, ".RDS"))



  TrueAdjs <- data_jgl$trueAdj
  JGL_est    <- data_jgl[[1]]$theta


  JGL <- as.data.frame(matrix(0,9,3))
  colnames(JGL) <- c("Hamming-Distance", "Precision", "Recall")


  for (m in 1:9) {
    dat <- JGL_est[[m]]

    dat[which(upper.tri(dat,diag = TRUE) == TRUE)] <- 0
    dat[which(abs(dat) > 0)] <- 1
    JGL[m,1:3] <- get_rates(TrueAdjs[[m]], dat)

    #JGL[m,] <- get_rates(TrAdj, JGL_est[[m]])
    #Glasso[m,] <- get_rates(TrAdj, Glasso_est[[m]])


  }
  JGL$lambda1 <- data_jgl[[2]][[1]]
  JGL$lambda2 <- data_jgl[[2]][[2]]
  JGL$id <- k

  for (r in 1:6) {
    Durchschnitte_JGL[k,r] <- mean(JGL[,r])
  }


  Results_as_frame_JGL <- rbind(Results_as_frame_JGL, JGL)



}




























ForPlots_MdS_Oracle <- subset(Durchschnitte,Durchschnitte$lambda1 == 11 & Durchschnitte$lambda2 == 3.75)
ForPlots_MdS_Aic <- subset(Durchschnitte, Durchschnitte$lambda1 == 14 & Durchschnitte$lambda2 == 3.75)
ForPlots_glasso <- subset(Durchschnitte_oracle, Durchschnitte_oracle$lambda == 0.5)
ForPlots_JGL <- subset(Durchschnitte_JGL)


save(ForPlots_MdS_Oracle, ForPlots_MdS_Aic, ForPlots_glasso, ForPlots_JGL, file = "PlottingData.RData")









