library(cowplot)
library(xtable)
library(ggplot2)
library(rethinking)
library(cmdstanr)
library(rstan)

seed <- 1312

# make output folders 

dir.create("stanfits", recursive = TRUE)

dir.create("output", recursive = TRUE)
dir.create("output/trace", recursive = TRUE)
dir.create("output/tables", recursive = TRUE)
dir.create("output/figures", recursive = TRUE)

# fit models

source("R/fit_father.R")
source("R/fit_mother.R")

# plot results

source("R/plot_survival.R")
source("R/plot_height.R")
source("R/plot_weight.R")
source("R/plot_education.R")
