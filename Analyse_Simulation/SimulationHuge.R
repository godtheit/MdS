rm(list = ls())

library(mvtnorm)
library(igraph)
library(matlib)
library(huge)
library(MASS)
library(sjmisc)
library(matrixcalc)
library(ggplot2)


#path <-"C:/Users/Arne/Dropbox/Masterarbeit/Multidimensional-Smoothing/R"
#path <- "C:/Users/arne2/Dropbox/Masterarbeit/Multidimensional-Smoothing/R"
source("Emulate_graph.R")
source("Get_corr.R")
set.seed(20)
N <- 100 #Anzahl Durchl?ufe
samples <- 1 #Anzahl samples gezogen aus simulierte Kovarianz
observes <- 52
length_var1 = 3
length_var2 = 3

p <- 150 #Anzahl variablen
perc1 <- 0.05
perc2 <- 0.05

BigListof_adjencies <- list()
BigListof_covs <- list()
BigListof_samples <- list()
BigListof_graphs <-  list()
BigCount <- 1
for (t in 1:N) {
  tmp <- Emulate_graph(p)
  listof_adjencies <- tmp[[1]]
  listof_covs <- list()
  listof_samples <- list()
  listof_graphs <- tmp[[2]]

  for (x in 1:length(listof_adjencies)) {
   test <- get_cor(theta = listof_adjencies[[x]], p = p, observes = observes)

   listof_covs[[x]] <- test[[1]]
   listof_samples[[x]] <- test [[2]]
  }

  BigListof_adjencies[[BigCount]] <- listof_adjencies
  BigListof_covs[[BigCount]] <- listof_covs
  BigListof_samples[[BigCount]] <- listof_samples
  BigListof_graphs[[BigCount]] <- listof_graphs


  BigCount <- BigCount + 1


}

save(BigListof_adjencies, file = "alladjacs.RData")
save(BigListof_covs, file = "allcovars.RData")
save(BigListof_graphs, file = "allgraphs.RData")
save(BigListof_samples, file = "allsamples.RData")






##Plot degree distributions
for (t in 1:N) {
  path <- paste("C:/Users/arne2/Dropbox/Masterarbeit/Simulation/Distributions2/Graph_", t, sep = "")

#par(mfrow=c(length_var1,length_var2))
maxMat <- matrix(0,1,length_var1*length_var2)
  for (r in 1:(length_var1*length_var2)) {
    maxMat[r] <- max(degree(BigListof_graphs[[t]][[r]]))
  }

count <- 1
plotter <- matrix(0,(length_var1*length_var2), max(maxMat)+1)
colnames(plotter) <- 0:max(maxMat)
for (i in 1:length_var1) {
  for(j in 1:length_var2){
    temp <- degree.distribution(BigListof_graphs[[t]][[count]])
    plotter[count,1:length(temp)] <- temp
    names(plotter) <- 0:(length(plotter)-1)
    test <- as.data.frame(plotter)



    png(file=paste(path, "/Theta_",i,j , ".png", sep = ""))
    par(mar = c(5.1, 3, 4.1, 10))
    barplot(plotter[count,],
            xlab = paste("G_", i,j, sep = ""), xlim = c(0, max(maxMat)),
            ylim = c(0,0.5),  col = "skyblue")
    dev.off()
    count <- count + 1
  }
}
}

# for (t in 1:N) {
#   dir.create(paste("C:/Users/arne2/Dropbox/Masterarbeit/Simulation/Distributions2/Graph_", t, sep = ""))
# }




# zahl <- 1
#
#   for (i in 1:length_var1) {
#     for (j in 1:length_var2) {
#
#    somegraph <- tmp[[2]][[zahl]]
#
#       png(file=paste(path, "/Pictures/Graphs_",t,"/Theta_",
#                      i, j, ".png", sep = ""))
#
#       plot.igraph(somegraph)
#       dev.off()
#
#     }
#     zahl <- zahl + 1
#   }


