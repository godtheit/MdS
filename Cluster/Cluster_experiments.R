
library(batchtools)


set.seed <- 42

# Registry ----------------------------------------------------------------
reg_name <- "example_experiments"
reg_dir <- file.path("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir)

# Problems -----------------------------------------------------------
source("problems.R")
addProblem(name = "sim_graphs", fun = sim_graphs, seed = 43)


# Algorithms -----------------------------------------------------------
source("Algorithms.R")
addAlgorithm(name = "MdS", fun = MdS_wrapper)
addAlgorithm(name = "JGL", fun = jgl_wrapper)
addAlgorithm(name = "glasso", fun = glasso_wrapper)

# Experiments -----------------------------------------------------------
prob_design <- list(sim_graphs = expand.grid(p = 150,
                                           perc1 = 0.05,
                                           perc2 = 0.05,
                                           length_var1 = 3,
                                           length_var2 = 3,
                                           observes = 52))


#Setup the weight matrix
th11 <- c(0,0,0,0,0,0,0,0,0) #1
th12 <- c(1,0,0,0,0,0,0,0,0) #2
th13 <- c(1,1,0,0,0,0,0,0,0) #3
th21 <- c(1,0,0,0,0,0,0,0,0) #4
th22 <- c(0,1,0,1,0,0,0,0,0) #5
th23 <- c(0,0,1,1,1,0,0,0,0) #6
th31 <- c(1,0,0,1,0,0,0,0,0) #7
th32 <- c(0,1,0,0,1,0,1,0,0) #8
th33 <- c(0,0,1,0,0,1,1,1,0) #9

w <- rbind(th11,th12, th13, th21, th22, th23, th31, th32, th33)
colnames(w) <- rownames(w)

###############################

#Lambdas for MdS
lam1_MdS <- seq(0.01,15, length.out = 20)
lam2_MdS <- seq(0.1,2, length.out = 8)

#Lambdas for JGL
lam1_JGL <- seq(0.02,2, length.out = 20)
lam2_JGL <- seq(0.1,2, length.out = 8)


algo_design <- list(MdS = expand.grid(lam1 = lam1_MdS, lam2 = lam2_MdS, w = w),
                    JGL = expand.grid(lam1 = lam1_JGL, lam2 = lam2_JGL),
                    glasso = expand.grid(lam1 = 1, lam2 = 1))

addExperiments(prob_design, algo_design, repls = 5) #100
summarizeExperiments()

# Test jobs -----------------------------------------------------------
testJob(id = 1)

# Submit -----------------------------------------------------------
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = reg_name, chunks.as.arrayjobs = TRUE,
  				              ncpus = 1, memory = 6000, walltime = 10*24*3600,
  							        max.concurrent.jobs = 40))
} else {
  submitJobs()
}
waitForJobs()

# Get results -------------------------------------------------------------
res <-  flatten(ijoin(reduceResultsDataTable(), getJobPars()))
res

# Plot results -------------------------------------------------------------
ggplot(res, aes(x = algorithm, y = error)) +
  facet_grid(problem ~ p) +
  geom_boxplot()
