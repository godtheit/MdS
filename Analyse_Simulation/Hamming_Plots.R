setwd("H:/Documents/MdS/Analyse_Simulation")
source("get_rates.R")


Daten <- readRDS("ImportantDataOmg_3fortest.RDS")




TrueAdjs           <- Daten$trueAdj
sharedAdjs         <- Daten$sharedAdj
MdS_decreasing_est <- Daten$placehodler1
MdS_hv_est         <- Daten$placehodler2         #horizontal/vertikal
JGL_est            <- Daten$placehodler3
Glasso_est         <- Daten$placehodler4


MdS1_list <- MdS2_list <- JGL_list <- Glasso_list <- list()

for(k in 1:length(Daten)){ #20
  for (m in 1:9) {
    TrAdj <- TrueAdj[[m]]

    MdS1 <- MdS2 <- JGL <- Glasso <- as.data.frame(matrix(0,9,3))
    colnames(MdS1) <- names(MdS2) <- names(JGL) <- names(Glasso) <- c("Hamming-Distance", "Precision", "Recall")

    MdS1[m,] <- get_rates(TrAdj, MdS_decreasing_est[[m]])
    MdS2[m,] <- get_rates(TrAdj, MdS_hv_est[[m]])
    JGL[m,] <- get_rates(TrAdj, JGL_est[[m]])
    Glasso[m,] <- get_rates(TrAdj, Glasso_est[[m]])
  }

  MdS1_list[[k]] <- MdS1
  MdS2_list[[k]] <- MdS2
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

