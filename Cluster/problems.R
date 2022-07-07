
source("Emulate_graph.R")
source("Get_Corr.R")
library(MASS)

sim_graphs <- function(data, job, observes, length_var1, length_var2 , p, perc1, perc2) {
  graph <- Emulate_graph(p, perc1, perc2, length_var1, length_var2)

  samples <- list()
  for (m in 1:length(graph[[1]])) {
    cor <- get_cor(theta = graph[[1]][[m]], p = p, observes = observes)
    samples[[m]] <- cor[[2]]
  }



  return(list(samples = samples,
       adjs = graph[[1]]))
}

