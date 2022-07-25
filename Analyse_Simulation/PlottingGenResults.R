setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation")

library(igraph)
library(visNetwork)
library(gridExtra)
Gen_network_MdS <- readRDS("Gen_network_MdS.RDS")
Gen_network_Glasso <- readRDS("Gen_network_glasso.RDS")


Adjs_MdS <- Gen_network_MdS$adj_matrices

Adjs_Glasso <- list()
Glasso_lambdas <- c()
graphs_Glasso <- list()
graphs_MdS <- list()
graphs_MdS_sub <- list()

for (m in 1:9) {
  Adjs_Glasso[[m]] <- Gen_network_Glasso[[m]][[1]]
  #diag(Adjs_Glasso[[m]]) <- 0
  #Adjs_Glasso[[m]][which(abs(Adjs_Glasso[[m]])> 0)] <-1
  Glasso_lambdas[m] <- Gen_network_Glasso[[m]][[2]]
}

for (m in 1:9) {
  graphs_Glasso[[m]] <- graph_from_adjacency_matrix(Adjs_Glasso[[m]],mode = "undirected")
  graphs_MdS[[m]] <- graph_from_adjacency_matrix(Adjs_MdS[[m]], mode= "undirected")

  V(graphs_MdS[[m]])$name <- 1:134

  nix <- which(degree(graphs_MdS[[m]]) != 0)
  graphs_MdS_sub[[m]] <- induced.subgraph(graphs_MdS[[m]],
                                          vids = nix)
}



######################### Graph Plots

path3 <- "C:/Users/arne2/Dropbox/Masterarbeit/MdS/Gen_graphs"
counter <- 1
for (i in 1:3) {
  for (j in 1:3) {


data <- toVisNetworkData(graphs_MdS[[counter]])
data$nodes$label <- as.character(data$nodes$label)
network1 <- visNetwork(nodes = data$nodes, edges = data$edges)%>%
  #visIgraphLayout(layout = "layout_in_circle") %>%
  #visIgraphLayout(layout = "layout_with_sugiyama") %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
             nodesIdSelection = TRUE, width = "100%", height = "100%")

data2 <- toVisNetworkData(graphs_MdS_sub[[counter]])
data2$nodes$label <- as.character(data2$nodes$label)
network2 <- visNetwork(nodes = data2$nodes, edges = data2$edges)%>%
  #visIgraphLayout(layout = "layout_in_circle") %>%
  #visIgraphLayout(layout = "layout_with_sugiyama") %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
             nodesIdSelection = TRUE, width = "100%", height = "100%")

visSave(network1, file = paste(path3, "/Theta_FULL_",i,j , ".html", sep = ""))
visSave(network2, file = paste(path3, "/Theta_Subnetwork_",i,j , ".html", sep = ""))

    counter <- counter + 1
  }
}


















####################### plot tables
path2 <- "C:/Users/arne2/Dropbox/Masterarbeit/MdS/Gen_tables"
counter <- 1
for (i in 1:3) {
  for (j in 1:3) {
    degrees <- degree(graphs_MdS[[counter]])
    sorted <- as.data.frame(sort(degrees, decreasing = TRUE))
    names(sorted) <- "Number of Edges"
    x <- head(sorted)

    png(file=paste(path2, "/Theta_",i,j , ".png", sep = ""))
    #png(file = paste(path,"table_for_graph",i,j, ".png" sep = ""), height=6, width=4)
    grid.table(x)
    dev.off()



    counter <- counter + 1
  }
}












######################### Degree Plots

#for (t in 1:N) {
  path <-  "C:/Users/arne2/Dropbox/Masterarbeit/MdS/DistributionsAnalyse"
  #par(mfrow=c(length_var1,length_var2))
  maxMat <- matrix(0,1,9)
  for (r in 1:9) {
    maxMat[r] <- max(degree(graphs_MdS[[r]]))
  }

  count <- 1
  plotter <- matrix(0,9, max(maxMat)+1)
  colnames(plotter) <- 0:max(maxMat)
  for (i in 1:3) {
    for(j in 1:3){
      temp <- degree.distribution(graphs_MdS[[count]])
      plotter[count,1:length(temp)] <- temp
      names(plotter) <- 0:(length(plotter)-1)
      test <- as.data.frame(plotter)



      png(file=paste(path, "/Theta_",i,j , ".png", sep = ""))
      par(mar = c(5.1, 3, 4.1, 10))
      barplot(plotter[count,],
              xlab = paste("G_", i,j, sep = ""), xlim = c(0, max(maxMat)),
              ylim = c(0,0.5),  col = "skyblue", legend = test)
      dev.off()
      count <- count + 1
    }
  }
#}

# for (t in 1:N) {
 #  dir.create(paste("C:/Users/arne2/Dropbox/Masterarbeit/Simulation/Distributions2/Graph_", t, sep = ""))
# }
#dir.create("C:/Users/arne2/Dropbox/Masterarbeit/MdS/DistributionsAnalyse")
