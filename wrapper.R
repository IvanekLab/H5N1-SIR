# Wrapper for pipeline

## Package dependencies
library(deSolve)
library(ggplot2)
library(reshape2)
library(rBayesianOptimization)
library(dplyr)
library(tidyr)

## Perform baseline
source("params_baseline.R")
source("model.R")
# source()

## Run Bayesian optimization
outbreak_data <- read.csv("outbreak_data.csv")
set.seed(11)
source("bayesian_optimization.R")
