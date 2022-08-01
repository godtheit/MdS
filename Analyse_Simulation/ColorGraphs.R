#setwd("H:/Documents/MdS/Daten_der_Analyse")
setwd("C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse")
library(igraph)
library(visNetwork)
library(gridExtra)



Gen_network_MdS <- readRDS("Daten/Final_Gennetwork2.RDS")
Adjs_MdS <- Gen_network_MdS$adj_matrices

graphs_MdS <- list()
subgraphs <- list()
add_all_adjs <- matrix(0,134,134)
shared_network <- matrix(0,134,134)
for (m in 1:9) {
  graphs_MdS[[m]] <- graph_from_adjacency_matrix(Adjs_MdS[[m]], mode= "undirected")
  V(graphs_MdS[[m]])$name <- as.character(1:134)

  add_all_adjs <- add_all_adjs + Adjs_MdS[[m]]

}
###Diese sind überall ohne Kanten
nicht_null  <- which((degree(graphs_MdS[[1]]) != 0) | (degree(graphs_MdS[[2]]) != 0) |
                   (degree(graphs_MdS[[3]]) != 0) | (degree(graphs_MdS[[4]]) != 0) |
                   (degree(graphs_MdS[[5]]) != 0) | (degree(graphs_MdS[[6]]) != 0) |
                   (degree(graphs_MdS[[7]]) != 0) | (degree(graphs_MdS[[8]]) != 0) |
                   (degree(graphs_MdS[[1]]) != 0)
                 )



for (m in 1:9) {
  subgraphs[[m]] <- induced.subgraph(graphs_MdS[[m]], vids = nicht_null)
  E(subgraphs[[m]])$color <- "red"
}
shared_network[which(add_all_adjs == 9)] <- 1

shared_network2 <- intersection(graphs_MdS[[1]], graphs_MdS[[2]], graphs_MdS[[3]],
                                graphs_MdS[[4]], graphs_MdS[[5]], graphs_MdS[[6]],
                                graphs_MdS[[7]], graphs_MdS[[8]], graphs_MdS[[9]], keep.all.vertices = FALSE)

graphs_shared <- graph_from_adjacency_matrix(shared_network, mode= "undirected")

shared <- induced.subgraph(shared_network2, vids = nicht_null)
E(shared)$color <- "blue"

graphs <- list()
for (m in 1:9) {
  graphs[[m]] <- union(shared,subgraphs[[m]])

  for (i in 1:length(E(graphs[[m]]))) {
    if (is.na(E(graphs[[m]])$color_1[i]) == TRUE) {
      E(graphs[[m]])$color[i] <- "red"
    } else { E(graphs[[m]])$color[i] <- "blue"}

  }
}

Gruppe <- c("GruppeA", "GruppeB", "GruppeC")
Strahlung <- c("0gy", "0.05gy", "2gy")


# path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe"
# counter <- 1
# for (i in 1:3) {
#   for (j in 1:3) {
#     png(file=paste(path, "/Theta", i,j, ".png", sep = ""), width= 960, height = 960, res = 120)
#     plot(graphs[[counter]], layout = layout_on_grid ,vertex.size = 0.1)
#     legend("topright", legend=c("Kernnetzwerk", "Änderung"),
#            col=c("blue", "red"), lty=1, cex=0.8,
#            box.lty=0)
#     title(paste(Gruppe[j],Strahlung[i]))
#     dev.off()
#     counter <- counter + 1
#   }
# }



library(bnstruct)
#bnstruct::shd(as.matrix(Adjs_MdS[[1]]), as.matrix())
#bnstruct::shd(Adjs_MdS[[1]], Adjs_MdS[[3]])

diff_matrix <- matrix(0, nrow=9, ncol = 9)

for (i in 1:9) {
  for (j in 1:9) {
    diff_matrix[i,j] <- bnstruct::shd(Adjs_MdS[[i]], Adjs_MdS[[j]])
    diff_matrix[j,i] <- diff_matrix[i,j]
  }
}
colnames(diff_matrix) <- rownames(diff_matrix) <- c("GrpA_0gy", "GrpB_0gy", "GrpC_0gy",
                                                    "GrpA_0.05gy", "GrpB_0.05gy", "GrpC_0.05gy",
                                                    "GrpA_2gy", "GrpB_2gy", "GrpC_2gy")
heatmap(diff_matrix)
library(reshape2)
melted_cormat <- melt(diff_matrix)
head(melted_cormat)

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile()

diffmatrix2 <- diff_matrix/100

ggcorrplot::ggcorrplot(diffmatrix2)
heatmap(diff_matrix)













#################### Unterschiede pro Gruppe




Subgraphs_0gy <- list(subgraphs[[1]], subgraphs[[2]], subgraphs[[3]])
Subgraphs_005gy <- list(subgraphs[[4]], subgraphs[[5]], subgraphs[[6]])
Subgraphs_2gy <- list(subgraphs[[7]], subgraphs[[8]], subgraphs[[9]])

Subgraphs_GrpA <- list(subgraphs[[1]], subgraphs[[4]], subgraphs[[7]])
Subgraphs_GrpB <- list(subgraphs[[2]], subgraphs[[5]], subgraphs[[8]])
Subgraphs_GrpC <- list(subgraphs[[3]], subgraphs[[6]], subgraphs[[9]])



###################### 0 gy
final_0gy <- list()
diff_0gy <- list()
shared_network_0gy <- intersection(Subgraphs_0gy[[1]],Subgraphs_0gy[[2]],Subgraphs_0gy[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_0gy[[m]] <- difference(Subgraphs_0gy[[m]], shared_network_0gy)
}


nicht_null_0gy <- which((degree(diff_0gy[[1]]) != 0) | (degree(diff_0gy[[2]]) != 0) |
                          (degree(diff_0gy[[3]]) != 0))

for (m in 1:3) {
final_0gy[[m]] <- induced.subgraph(diff_0gy[[m]], vids = nicht_null_0gy)
E(final_0gy[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/Ohne0gy"

for (j in 1:3) {
    png(file=paste(path, "/Theta", 1,j, ".png", sep = ""), width= 960, height = 960, res = 120)
    plot(final_0gy[[j]], layout = layout_on_grid ,vertex.size = 0.1)
    dev.off()


}

###################### 0.05 gy
final_005gy <- list()
diff_005gy <- list()
shared_network_005gy <- intersection(Subgraphs_005gy[[1]],Subgraphs_005gy[[2]],Subgraphs_005gy[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_005gy[[m]] <- difference(Subgraphs_005gy[[m]], shared_network_005gy)
}


nicht_null_005gy <- which((degree(diff_005gy[[1]]) != 0) | (degree(diff_005gy[[2]]) != 0) |
                          (degree(diff_005gy[[3]]) != 0))

for (m in 1:3) {
  final_005gy[[m]] <- induced.subgraph(diff_005gy[[m]], vids = nicht_null_005gy)
  E(final_005gy[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/Ohne005gy"

for (j in 1:3) {
  png(file=paste(path, "/Theta", 2,j, ".png", sep = ""), width= 960, height = 960, res = 120)
  plot(final_005gy[[j]], layout = layout_on_grid ,vertex.size = 0.1)
  dev.off()


}

###################### 2 gy
final_2gy <- list()
diff_2gy <- list()
shared_network_2gy <- intersection(Subgraphs_2gy[[1]],Subgraphs_2gy[[2]],Subgraphs_2gy[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_2gy[[m]] <- difference(Subgraphs_2gy[[m]], shared_network_2gy)
}


nicht_null_2gy <- which((degree(diff_2gy[[1]]) != 0) | (degree(diff_2gy[[2]]) != 0) |
                          (degree(diff_2gy[[3]]) != 0))

for (m in 1:3) {
  final_2gy[[m]] <- induced.subgraph(diff_2gy[[m]], vids = nicht_null_2gy)
  E(final_2gy[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/Ohne2gy"

for (j in 1:3) {
  png(file=paste(path, "/Theta", 3,j, ".png", sep = ""), width= 960, height = 960, res = 120)
  plot(final_2gy[[j]], layout = layout_on_grid ,vertex.size = 0.1)
  dev.off()


}
###################### GruppeA
final_GrpA <- list()
diff_GrpA <- list()
shared_network_GrpA <- intersection(Subgraphs_GrpA[[1]],Subgraphs_GrpA[[2]],Subgraphs_GrpA[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_GrpA[[m]] <- difference(Subgraphs_GrpA[[m]], shared_network_GrpA)
}


nicht_null_GrpA <- which((degree(diff_GrpA[[1]]) != 0) | (degree(diff_GrpA[[2]]) != 0) |
                          (degree(diff_GrpA[[3]]) != 0))

for (m in 1:3) {
  final_GrpA[[m]] <- induced.subgraph(diff_GrpA[[m]], vids = nicht_null_GrpA)
  E(final_GrpA[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/OhneGruppeA"

for (j in 1:3) {
  png(file=paste(path, "/Theta", j,1, ".png", sep = ""), width= 960, height = 960, res = 120)
  plot(final_GrpA[[j]], layout = layout_on_grid ,vertex.size = 0.1)
  dev.off()


}
###################### GruppeB
final_GrpB <- list()
diff_GrpB <- list()
shared_network_GrpB <- intersection(Subgraphs_GrpB[[1]],Subgraphs_GrpB[[2]],Subgraphs_GrpB[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_GrpB[[m]] <- difference(Subgraphs_GrpB[[m]], shared_network_GrpB)
}


nicht_null_GrpB <- which((degree(diff_GrpB[[1]]) != 0) | (degree(diff_GrpB[[2]]) != 0) |
                           (degree(diff_GrpB[[3]]) != 0))

for (m in 1:3) {
  final_GrpB[[m]] <- induced.subgraph(diff_GrpB[[m]], vids = nicht_null_GrpB)
  E(final_GrpB[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/OhneGruppeB"

for (j in 1:3) {
  png(file=paste(path, "/Theta", j,2, ".png", sep = ""), width= 960, height = 960, res = 120)
  plot(final_GrpB[[j]], layout = layout_on_grid ,vertex.size = 0.1)
  dev.off()


}

###################### GruppeC
final_GrpC <- list()
diff_GrpC <- list()
shared_network_GrpC <- intersection(Subgraphs_GrpC[[1]],Subgraphs_GrpC[[2]],Subgraphs_GrpC[[3]], keep.all.vertices = FALSE)
for (m in 1:3) {
  diff_GrpC[[m]] <- difference(Subgraphs_GrpC[[m]], shared_network_GrpC)
}


nicht_null_GrpC <- which((degree(diff_GrpC[[1]]) != 0) | (degree(diff_GrpC[[2]]) != 0) |
                           (degree(diff_GrpC[[3]]) != 0))

for (m in 1:3) {
  final_GrpC[[m]] <- induced.subgraph(diff_GrpC[[m]], vids = nicht_null_GrpC)
  E(final_GrpC[[m]])$color <- "blue"
}

path <- "C:/Users/Arne/Dropbox/Masterarbeit/MdS/Daten_der_Analyse/Graphs/Graphs_mit_Farbe/OhneGruppeC"

for (j in 1:3) {
  png(file=paste(path, "/Theta", j,3, ".png", sep = ""), width= 960, height = 960, res = 120)
  plot(final_GrpC[[j]], layout = layout_on_grid ,vertex.size = 0.1)
  dev.off()


}
