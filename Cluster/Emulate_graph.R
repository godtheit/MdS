Emulate_graph <- function(p, perc1 = 0.05, perc2 = 0.05, length_var1 = 3, length_var2 = 3){
  library(sjmisc)
  p=150
  # length_var1 = 3
  # length_var2 = 3
  # perc1 <- perc2 <- 0.05
  if(!is_even(p)){stop("Number of Variables must be an even number!")}
  Adjecency_matrices <- list()
  graphs <- list()



  # graph_2 <- huge.generator(d = p, graph = "scale-free")
  # test <- as.matrix(graph_2$theta)
  # number_edges <- table(test)["1"]/2
  # a <- round(number_edges*perc1, digits = 0)  #Anzahl Ver?nderungen Schritt 1
  # b <- round(number_edges*perc2, digits = 0) #Anzahl Ver?nderungen Schritt 2
  # adjmat <- lower.triangle(test)
  #graph_1 <- sample_pa(n=p, power = 2.5, directed = FALSE)

  graph_2 <- barabasi.game(p, directed = FALSE)
  graph_1 <-    permute.vertices(graph_2, sample(1:p,p))

  a <- round(gsize(graph_1)*perc1, digits = 0)  #Anzahl Ver?nderungen Schritt 1
  b <- round(gsize(graph_1)*perc2, digits = 0) #Anzahl Ver?nderungen Schritt 2
  adjmat <- as.matrix(as_adj(graph_1, type = "lower")) #get the adjacency matrice of the scale free network





  Shared_lowerM <- adjmat[((p/2)+1):p , 1:(p/2)] #Save the shared adjacency
  I <- matrix(0,p/2,p/2)




for (i in 1:length_var1) {
  assign(paste("A_",i,sep = ""), adjmat[1:(p/2) , 1:(p/2)])
}
for (j in 1:length_var2) {
  assign(paste("B_",j,sep = ""), adjmat[((p/2)+1):p,((p/2)+1):p])
}




  t <- which(lower.tri(A_1)==TRUE) #Lower Triangular of p/2 x p/2 matrice, generic


  t0_A <- t[which(A_1[t]==0)]
  t1_A <- t[which(A_1[t]==1)]


  for (i in 2:length_var1) {
    tmp <- eval(as.name(paste("A_",i, sep = "")))

    tmp[sample(t0_A, a, replace = FALSE)]  <- 1
    tmp[sample(t1_A, a, replace = FALSE)]  <- 0

    assign(paste("A_", i, sep = ""), tmp)

    t0_A <- t[which(tmp[t] == 0)]
    t1_A <- t[which(tmp[t] == 1)]
  }



  t0_B <- t[which(B_1[t]==0)]
  t1_B <- t[which(B_1[t]==1)]

  for (j in 2:length_var2) {
    tmp <- eval(as.name(paste("B_",j, sep = "")))

    tmp[sample(t0_B, b, replace = FALSE)]  <- 1
    tmp[sample(t1_B, b, replace = FALSE)]  <- 0

    assign(paste("B_", j, sep = ""), tmp)

    t0_B <- t[which(tmp[t] == 0)]
    t1_B <- t[which(tmp[t] == 1)]
  }


  U_11 <- adjmat
  rownames(U_11) <- colnames(U_11) <- 1:p

counter <- 1

 for (i in 1:length_var1) {
   for (j in 1:length_var2) {
     A <- eval(as.name(paste("A_", i, sep = "")))
     B <- eval(as.name(paste("B_", j, sep = "")))



     U <- rbind((cbind(A,I)), cbind(Shared_lowerM,B))
        assign(paste("U_", i,j, sep = ""),  U)




     g <- as.undirected(graph_from_adjacency_matrix(eval(as.name(paste("U_", i,j, sep = "")))))

     V(g)[1:(p/2)]$color <- "lightblue"
     V(g)[((p/2)+1):p]$color <- "orange"
        assign(paste("g_",i,j, sep = ""), g)

        #Bisher betrachten wir nur untere adjecency matrix
      U <- U+t(U)

        Adjecency_matrices[[counter]] <- U
        graphs[[counter]] <- g
        counter <- counter + 1

   }
 }



return(list(Adjecency_matrices,graphs))


}
