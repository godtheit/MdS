
library(batchtools)
setwd("H:/Documents/MdS/Cluster")

set.seed <- 42

# Registry ----------------------------------------------------------------

### Setting up the repository
start_from_scratch <- T # if true, removes all repository and creates a new one

reg_name <- "MdS_registry"
reg_dir <- sprintf("%s/registries/%s", getwd(), reg_name)

if (start_from_scratch) {
  dir.create("registries", showWarnings = FALSE)
  unlink(reg_dir, recursive = TRUE)
  reg <- makeExperimentRegistry(file.dir = reg_dir)
} else {
  reg <- loadRegistry(file.dir = reg_dir, writeable = TRUE)
}


# Problems -----------------------------------------------------------
source("problems.R")
addProblem(name = "sim_graphs", fun = sim_graphs, seed = 43)


# Algorithms -----------------------------------------------------------
source("Algorithms.R")
addAlgorithm(name = "MdS", fun = MdS_wrapper)
addAlgorithm(name = "JGL", fun = JGL_wrapper)
addAlgorithm(name = "glasso", fun = glasso_wrapper)

# Experiments -----------------------------------------------------------
prob_design <- list(sim_graphs = expand.grid(p = 150,
                                           perc1 = 0.05,
                                           perc2 = 0.05,
                                           length_var1 = 3,
                                           length_var2 = 3,
                                           observes = 52))


#Setup the weight matrix


###############################


# a min lambda1
# b max lambda1
# c min lambda2
# d max lambda2

# L1 step count for lambda1
# L2 step count for lambda2

algo_design <- list(MdS = expand.grid(a = 0.1, b = 30, c = 0.1, d = 2, L1 = 20, L2 = 8 ,
                                      weight_matrix = c("horizontal-vertikal",  "full")),
                    JGL = expand.grid(a = 0.1, b = 2, c = 0.1, d = 2, L1 = 20, L2 = 8 ),
                    glasso = expand.grid())

addExperiments(prob_design, algo_design, repls = 2)
summarizeExperiments()

# Test jobs -----------------------------------------------------------

testJob(id = 6)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = reg_name, chunks.as.arrayjobs = TRUE,
  				              ncpus = 1, memory = 6000, walltime = 10*24*3600,
  							        max.concurrent.jobs = 200))
} else {
  submitJobs()
}
waitForJobs()

getStatus()

# Get results -------------------------------------------------------------
res <- ijoin(reduceResultsDataTable(), getJobPars())

saveRDS(res, file = "ImportantDataOmg.RDS")

