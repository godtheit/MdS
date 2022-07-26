#setwd("C:/Users/arne2/Dropbox/Masterarbeit/MdS/Analyse_Simulation")

setwd("H:/Documents/MdS/Daten_der_Analyse")

library(igraph)
library(visNetwork)
library(gridExtra)

#Gen_network_MdS <- readRDS("Daten/Gen_network_MdS.RDS")
Gen_network_Glasso <- readRDS("Daten/Gen_network_glasso.RDS")
Gen_network_JGL <- readRDS("Daten/Gen_network_JGL.RDS")

Gen_network_MdS <- readRDS("Daten/Final_Gennetwork2.RDS")



Adjs_MdS <- Gen_network_MdS$adj_matrices
Adjs_JGL <- list()
Adjs_Glasso <- list()
Glasso_lambdas <- c()
graphs_Glasso <- list()
graphs_MdS <- list()
graphs_MdS_sub <- list()
graphs_JGL <- list()
graphs_JGL_sub <- list()

for (m in 1:9) {
  Adjs_Glasso[[m]] <- Gen_network_Glasso[[m]][[1]]
  Adjs_JGL[[m]] <- Gen_network_JGL[[1]]$theta[[m]]

  diag(Adjs_Glasso[[m]]) <- 0
  diag(Adjs_JGL[[m]]) <- 0

  Adjs_Glasso[[m]][which(abs(Adjs_Glasso[[m]])> 0)] <-1
  Adjs_JGL[[m]][which(abs(Adjs_JGL[[m]])> 0)] <-1

  Glasso_lambdas[m] <- Gen_network_Glasso[[m]][[2]]
}

for (m in 1:9) {
  graphs_Glasso[[m]] <- graph_from_adjacency_matrix(Adjs_Glasso[[m]],mode = "undirected")
  graphs_JGL[[m]] <- graph_from_adjacency_matrix(Adjs_JGL[[m]],mode = "undirected")


  graphs_MdS[[m]] <- graph_from_adjacency_matrix(Adjs_MdS[[m]], mode= "undirected")

  V(graphs_MdS[[m]])$name <- 1:134

  V(graphs_JGL[[m]])$name <- as.character(1:134)
  nix2 <- which(degree(graphs_JGL[[m]]) != 0)
  graphs_JGL_sub[[m]] <- induced.subgraph(graphs_JGL[[m]],
                                      vids = nix2)

V(graphs_JGL[[m]])$name <- as.character(1:134)
  nix <- which(degree(graphs_MdS[[m]]) != 0)
  graphs_MdS_sub[[m]] <- induced.subgraph(graphs_MdS[[m]],
                                          vids = nix)
}



######################### Graph Plots

path3 <- "H:/Documents/MdS/Daten_der_Analyse/Graphs/Graphs_JGL"
counter <- 1
for (i in 1:3) {
  for (j in 1:3) {

#
# data <- toVisNetworkData(graphs_MdS[[counter]])
# data$nodes$label <- as.character(data$nodes$label)
# network1 <- visNetwork(nodes = data$nodes, edges = data$edges)%>%
#   #visIgraphLayout(layout = "layout_in_circle") %>%
#   #visIgraphLayout(layout = "layout_with_sugiyama") %>%
#   visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
#              width = "100%", height = "100%")
#
# data2 <- toVisNetworkData(graphs_MdS_sub[[counter]])
# data2$nodes$label <- as.character(data2$nodes$label)
# network2 <- visNetwork(nodes = data2$nodes, edges = data2$edges)%>%
#   #visIgraphLayout(layout = "layout_in_circle") %>%
#   #visIgraphLayout(layout = "layout_with_sugiyama") %>%
#   visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
#               width = "100%", height = "100%")
#
#
# visSave(network1, file = paste(path3, "/Theta_FULL_",i,j , ".html", sep = ""))
# visSave(network2, file = paste(path3, "/Theta_Subnetwork_",i,j , ".html", sep = ""))
g <- graphs_JGL[[counter]]
g2 <- graphs_JGL_sub[[counter]]
max_deg <- max(degree(g)) # with a length of X vertices :
max_deg2 <- max(degree(g2))
colors <- (max_deg+20):20

gray_colors <- as.data.frame(paste0("gray", colors))
colnames(gray_colors)  <- "farbe"

for (n in 1:length(V(g))) {
  V(g)[n]$color <- gray_colors[degree(g)[n] + 1,]
}
for (n in 1:length(V(g2))) {
  V(g2)[n]$color <- gray_colors[degree(g2)[n] + 1,]
}

png(file=paste(path3, "/Theta_",i,j , ".png", sep = ""), width= 960, height = 960, res = 120)
plot(g, layout=layout_with_fr,vertex.label.dist=0.85, vertex.label.degree=pi/2, vertex.size=5)
dev.off()


png(file=paste(path3, "/Theta_SUB_",i,j , ".png", sep = ""), width= 960, height = 960, res = 120)
plot(g2, layout=  layout_on_grid,vertex.label.dist=0.85, vertex.label.degree=pi/2,vertex.size = 5)
dev.off()


    counter <- counter + 1
  }
}
layout_as_tree



















####################### plot tables
path2 <- "H:/Documents/MdS/Daten_der_Analyse/Tables/Tables_Final2"
counter <- 1
for (i in 1:3) {
  for (j in 1:3) {
    degrees <- degree(graphs_MdS[[counter]])
    sorted <- as.data.frame(sort(degrees, decreasing = TRUE))



    names(sorted) <- "Anzahl Knoten"
x <- head(sorted, n = 12)

    png(file=paste(path2, "/Theta_",i,j , ".png", sep = ""))
    #png(file = paste(path,"table_for_graph",i,j, ".png" sep = ""), height=6, width=4)
    grid.table(x)
    dev.off()



    counter <- counter + 1
  }
}












######################### Degree Plots

#for (t in 1:N) {
  path <-  "H:/Documents/MdS/Daten_der_Analyse/Distributions/Distributions_Final2"
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
      temp <- degree.distribution(graphs_JGL[[count]])
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
#}

################################# For all Datasets (n = 32)

  lambda1 <- c("0.5", "2", "5", "10")
  lambda2 <- c("0.5", "2" ,"5", "10")


   for (w in 1:2) {
     for (l1 in lambda1) {
       for (l2 in lambda2) {

         ## Prepare the data
         datapath <- paste("H:/Documents/MdS/Daten_der_Analyse/Daten/lambda1_",
                           l1, "_lambda2_", l2,  "_weightmatrix_", w ,"_Gennetwork.RDS", sep = "")


         Gen_network <- readRDS(datapath)
         Adjs <- Gen_network$adj_matrices

         graphs <- list()
         subgraphs <- list()


         for (m in 1:9) {
           graphs[[m]] <- graph_from_adjacency_matrix(Adjs[[m]], mode= "undirected")
           V(graphs[[m]])$name <- as.character(1:134)

           nix <- which(degree(graphs[[m]]) != 0)
           subgraphs[[m]] <- induced.subgraph(graphs[[m]],
                                                   vids = nix)
         }

         ######################### Degree Plots


         path <-  paste("H:/Documents/MdS/Daten_der_Analyse/Distributions/Distributions_MdS_", l1, "_", l2,  "_", "w", w , sep = "")

         maxMat <- matrix(0,1,9)
         for (r in 1:9) {
           maxMat[r] <- max(degree(graphs[[r]]))
         }

         count <- 1
         plotter <- matrix(0,9, max(maxMat)+1)
         colnames(plotter) <- 0:max(maxMat)
         for (i in 1:3) {
           for(j in 1:3){
             temp <- degree.distribution(graphs[[count]])
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


         ####################### plot tables
         path2 <- paste("H:/Documents/MdS/Daten_der_Analyse/Tables/Tables_MdS_", l1, "_", l2,  "_", "w", w , sep = "")
         counter <- 1
         for (i in 1:3) {
           for (j in 1:3) {
             degrees <- degree(graphs[[counter]])
             sorted <- as.data.frame(sort(degrees, decreasing = TRUE))



             names(sorted) <- "Anzahl Knoten"
             x <- head(sorted, n = 12)

             png(file=paste(path2, "/Theta_",i,j , ".png", sep = ""))
             #png(file = paste(path,"table_for_graph",i,j, ".png" sep = ""), height=6, width=4)
             grid.table(x)
             dev.off()



             counter <- counter + 1
           }
         }


         ######################### Graph Plots

         path3 <- paste("H:/Documents/MdS/Daten_der_Analyse/Graphs/Graphs_MdS_", l1, "_", l2,  "_", "w", w , sep = "")
         counter <- 1
         for (i in 1:3) {
           for (j in 1:3) {

             #
             # data <- toVisNetworkData(graphs_MdS[[counter]])
             # data$nodes$label <- as.character(data$nodes$label)
             # network1 <- visNetwork(nodes = data$nodes, edges = data$edges)%>%
             #   #visIgraphLayout(layout = "layout_in_circle") %>%
             #   #visIgraphLayout(layout = "layout_with_sugiyama") %>%
             #   visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
             #              width = "100%", height = "100%")
             #
             # data2 <- toVisNetworkData(graphs_MdS_sub[[counter]])
             # data2$nodes$label <- as.character(data2$nodes$label)
             # network2 <- visNetwork(nodes = data2$nodes, edges = data2$edges)%>%
             #   #visIgraphLayout(layout = "layout_in_circle") %>%
             #   #visIgraphLayout(layout = "layout_with_sugiyama") %>%
             #   visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
             #               width = "100%", height = "100%")
             #
             #
             # visSave(network1, file = paste(path3, "/Theta_FULL_",i,j , ".html", sep = ""))
             # visSave(network2, file = paste(path3, "/Theta_Subnetwork_",i,j , ".html", sep = ""))
             g <- graphs[[counter]]
             g2 <- subgraphs[[counter]]
             max_deg <- max(degree(g)) # with a length of X vertices :
             max_deg2 <- max(degree(g2))

             if (max_deg < 80 ) {
               colors <- (max_deg+20):20
             } else if (max_deg >= 80){
               x <- 100 - max_deg
               colors <- 100:x
             }


             gray_colors <- as.data.frame(paste0("gray", colors))
             colnames(gray_colors)  <- "farbe"

             for (n in 1:length(V(g))) {
               V(g)[n]$color <- gray_colors[degree(g)[n] + 1,]
             }
             for (n in 1:length(V(g2))) {
               V(g2)[n]$color <- gray_colors[degree(g2)[n] + 1,]
             }

             png(file=paste(path3, "/Theta_",i,j , ".png", sep = ""), width= 960, height = 960, res = 120)
             plot(g, layout=layout_with_fr, vertex.size=5)
             dev.off()


             png(file=paste(path3, "/Theta_SUB_",i,j , ".png", sep = ""), width= 960, height = 960, res = 120)
             plot(g2, layout=layout_with_fr, vertex.size=5)
             dev.off()


             counter <- counter + 1
           }
         }



       }
     }
   }

