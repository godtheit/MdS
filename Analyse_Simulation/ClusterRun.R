#Run this Script on Cluster
rm(list = ls())
library(batchtools)
library(MdS2)
library(huge)
library(glasso)
library(JGL)
library(genlasso)



load("allsamples.RData")


###############################

#Lambdas for MdS
lam1_MdS <- seq(0.01,15, length.out = 20)
lam2_MdS <- seq(0.1,2, length.out = 8)

#Lambdas for JGL
lam1_JGL <- seq(0.02,2, length.out = 20)
lam2_JGL <- seq(0.1,2, length.out = 8)

### do search over the 28 Lambda Values for each of the 100 networks


#Registry
reg_name <- "Reg_MdS"
reg_dir <- file.path("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeRegistry(file.dir = reg_dir)


batchMap(fun = procedure,  args = list(Y = BigListof_samples, lam1 = lam1_MdS, lam2 =lam2_MdS, method = "MdS", w = w))
submitJobs()
batchMap(fun = procedure, args = list(Y = BigListof_samples, lam1 = lam1_JGL, lam2 = lam2_JGL, method = "JGL"))
