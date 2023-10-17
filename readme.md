Pimbwe parent marriage status and child health outcomes
============

This repository contains code to reproduce results from:

Age-specific impacts of time-varying family structures on on children's well-being in Mpimbwe, Tanzania --- Riana Minocher, Monique Borgerhoff Mulder, Cody T. Ross

The code was written in R (v4.0.4) and Stan. Statistical models are fit using the Stan engine in R, implemented with the `cmdstanr` package (0.5.2), which requires a C++ compiler. Installation instructions are available at https://mc-stan.org/cmdstanr/. 

A number of packages are required to process the model outputs. These include: rstan, cmdstanr, rethinking, xtable, ggplot2, cowplot. 

## How to run

Make sure you have all packages and dependencies installed. Navigate to the directory, and source the script "run_all.R".

## About the data

We provide anonymized, processed data to reproduce all analyses reported in the paper. The Stan code is general and can be extended to other data of similar structure. Raw data (i.e., longitudinal anthropometric measures, educational achievements, personal family information) and scripts including identifying information about individuals will not be shared online. 

```
