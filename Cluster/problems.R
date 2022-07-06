
source("Emulate_graph.R")
source("Get_Corr.R")


sim_graphs <- function(data, job, observes, length_var1, length_var2 , p, perc1, perc2) {
  graph <- Emulate_graph(p, perc1, perc2, length_var1, length_var2)

  samples <- list()
  for (m in 1:length(graph)) {
    cor <- get_cor(theta = graph[[m]], p = p, observes = observes)
    samples[[m]] <- cor[[2]]
  }



  return(list(samples = samples,
       adjs = graph[[1]]))
}

