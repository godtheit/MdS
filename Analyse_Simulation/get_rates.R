   get_rates <- function(TrueAdj, SimAdj){


   tet <- (TrueAdj == SimAdj)
   tet[which(upper.tri(tet, diag = TRUE))] <- "irrelevant"

   TP <- sum(testmatrix2[which(tet == TRUE & testmatrix2 == 1)])
   TN <- sum(testmatrix2[which(tet == TRUE & testmatrix2 == 0)])
   FP <- sum(testmatrix2[which(tet == FALSE & testmatrix2 == 1)])
   FN <- sum(testmatrix2[which(tet == FALSE & testmatrix2 == 0)])

   hamming <- FP + FN
   Prec    <- TP/(TP + FP)
   Rec     <- TP/(TP + FN)

   rate <- c(hamming, Prec, Rec)
   names(rate) <- c("Hamming-distance", "Precision", "Recall")
   return(rate)
   }
