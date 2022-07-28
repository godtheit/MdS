setwd("H:/Documents/MdS/Daten_der_Analyse")

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
  V(graphs_MdS[[m]])$name <- 1:134

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
}
shared_network[which(add_all_adjs == 9)] <- 1

shared_network2 <- intersection(graphs_MdS[[1]], graphs_MdS[[2]], graphs_MdS[[3]],
                                graphs_MdS[[4]], graphs_MdS[[5]], graphs_MdS[[6]],
                                graphs_MdS[[7]], graphs_MdS[[8]], graphs_MdS[[9]], keep.all.vertices = FALSE)

graphs_shared <- graph_from_adjacency_matrix(shared_network, mode= "undirected")



Edges_gemeinsam <- E(shared_network2)




for (m in 1:9) {
  for (k in 1:length(E(subgraphs[[m]])) ) {
    if(E(subgraphs[[m]])[k] %in% Edges_gemeinsam == TRUE){
      E(subgraphs[[m]])$color[k] <- 'blue'
    } else if (E(subgraphs[[m]])[k] %in% Edges_gemeinsam == FALSE) {
      cat(k, "   ")
      E(subgraphs[[m]])$color[k] <- 'red'
      }
  }
}


plot(subgraphs[[1]], layout = layout_on_grid ,vertex.size = 0.1)
legend("topright", legend=c("Kernnetzwerk", "Andere"),
       col=c("red", "blue"), lty=1, cex=0.8,
       box.lty=0)
title("GruppeB, 0gy")

Edges_gemeinsam[33]
E(subgraphs[[1]])[143]

