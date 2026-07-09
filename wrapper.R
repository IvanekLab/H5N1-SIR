# Wrapper for pipeline

## Package dependencies
library(deSolve)
library(ggplot2)
library(reshape2)
library(rBayesianOptimization)
library(dplyr)
library(tidyr)

## Run baseline
source("params_baseline.R")
source("model.R")
source("calculate_results.R")
write.csv(output, "baseline_results.csv")
baseline_vector <- vector

## Run Bayesian optimization
outbreak_data <- read.csv("outbreak_data.csv")
set.seed(11)
source("bayesian_optimization.R")

## Run model with optimized parameters
parameters["kappa_1"] <- best_params["kappa_1"] 
parameters["kappa_2"] <- best_params["kappa_2"]
parameters["kappa_3"] <- best_params["kappa_3"]
parameters["I_Ds0"] <- best_params["I_Ds0"]
source("model.R")
source("calculate_results.R")
write.csv(output, "optimized_results.csv")
optimized_vector <- vector

bayesian_table <- data.frame(baseline_vector,target_vector,optimized_vector)
