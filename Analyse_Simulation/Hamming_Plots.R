setwd("H:/Documents/MdS/Analyse_Simulation/")
source("get_rates.R")






funfzig <- readRDS("test500.RDS")

zwoelf <- readRDS("Test12.RDS")

eins <- readRDS("TestMdS1.RDS")
zwei <- readRDS("TestMdS2.RDS")
drei <- readRDS("TestMdS3.RDS")
vier <- readRDS("TestMdS4.RDS")

load("alladjacs.RData")


get_rates(BigListof_adjencies[[1]][[2]],zwoelf$adj_matrices[[4]])
#TrueAdjs <- BigListof_adjencies[[1]]
#sharedAdjs         <- eins$shared_adj
#   zeros_shared <- matrix(0,67,67)
#shared_adjacency <- rbind(cbind(zeros_shared, zeros_shared), cbind(sharedAdjs,zeros_shared)    )




#JGL_est            <- Daten$placehodler3
#Glasso_est         <- Daten$placehodler4


MdS_list  <- JGL_list <- Glasso_list <- list()

for(k in 1:20){
  Daten <- readRDS(paste0("H:/Documents/MdS/Cluster/Ergebnisse/",k, ".RDS")  )
  TrueAdjs <- Daten$trueAdj
  MdS_est            <- Daten$adj_matrices
  MdS  <- JGL <- Glasso <- as.data.frame(matrix(0,9,3))
  colnames(MdS) <- names(JGL) <- names(Glasso) <- c("Hamming-Distance", "Precision", "Recall")
  for (m in 1:9) {





    MdS[m,] <- get_rates(TrueAdjs[[m]], MdS_est[[m]])
    #JGL[m,] <- get_rates(TrAdj, JGL_est[[m]])
    #Glasso[m,] <- get_rates(TrAdj, Glasso_est[[m]])
  }

  MdS_list[[k]] <- MdS
  JGL_list[[k]] <- JGL
  Glasso_list[[k]] <- Glasso
}

############ Start plotting ######################



for (k in 1:20) {
  for (m in 1:10) {




    if(m == 10)
    {
      #plot summe
    }
  }
}


testmatrix1 <-matrix(sample(0:1,64,replace=TRUE),ncol=8)
testmatrix2 <-matrix(sample(0:1,64,replace=TRUE),ncol=8)

