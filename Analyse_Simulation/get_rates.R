   get_rates <- function(TrueAdj, SimAdj){

   TrueAdj[which(upper.tri(TrueAdj,diag = TRUE))] <- 5
   SimAdj[which(upper.tri(SimAdj,diag = TRUE))] <- 6

   tet <- (TrueAdj == SimAdj)
   tet[which(upper.tri(tet, diag = TRUE))] <- "irrelevant"

   TP <- length(which(tet == TRUE & SimAdj == 1))
   TN <- length(which(tet == TRUE & SimAdj == 0))
   FP <- length(which(tet == FALSE & SimAdj == 1))
   FN <- length(which(tet == FALSE & SimAdj == 0))

   hamming <- FP + FN
   Prec    <- TP/(TP + FP)
   Rec     <- TP/(TP + FN)

   rate <- c(hamming, Prec, Rec)
   names(rate) <- c("Hamming-distance", "Precision", "Recall")
   return(rate)
   }
