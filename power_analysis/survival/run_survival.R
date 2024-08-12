# conduct power analysis 

library(cmdstanr)
library(rethinking)
library(ggplot2)
library(parallel)

# create data as in analysis
source("R/load_base_data.R")

# set parameter values
source("R/load_base_parameters.R")

# function to update parameter values 
source("R/modify_parameters_function.R")

# function to simulate data
source("R/simulate_data_function.R")

# define simulation parameters

S = 9 # number of variables
Q = 11 # set of values to iterate over
A = 19 # age categories

results_list = vector("list", Q*S)
model_data_list = vector("list", Q*S)
model_parameters_list = vector("list", Q*S)
slot_list = vector("list", Q*S)

# which age parameter to vary
for (slot in 2:10) {

# set the range of effect sizes [-2 to 2]
mes = seq(-2, 2, length.out = Q)

# create the simulated sets of data
for (q in 1:Q) {
 
 model_parameters_list[[q + (slot-2)*Q]] = modify_age_parameter_set(model_parameters, seq(mes[q], 0, length.out = A), slot)
 
 model_data_list[[q + (slot-2)*Q]] = simulate_data(model_parameters_list[[q + (slot-2)*Q]])

 slot_list[[q + (slot-2)*Q]] = slot

 }

}

# compile model
m = cmdstan_model("stan/survival_father.stan")

# fit in parallel
fit = mclapply(1:(Q*S), 
               function(d) m$sample(data = model_data_list[[d]],
                                    chains = 1,
                                    iter_sampling = 1000,
                                    parallel_chains = 1,
                                    adapt_delta = 0.99,
                                    max_treedepth = 15), 
               mc.cores = (Q*S))

# save output
for (q in 1:(Q*S)) {
  
  stanfit = posterior::as_draws_rvars(fit[[q]]$draws())
  post = posterior::draws_of(stanfit$"a_age")
  slot = slot_list[[q]]
  
  true_values = model_parameters_list[[q]]$a_age[slot, ]
  median_values = apply(post[, slot, ], 2, median)
  L_values = apply(post[, slot, ], 2, function(x) HPDI(x)[1])
  H_values = apply(post[, slot, ], 2, function(x) HPDI(x)[2])
  
  results_list[[q]] = data.frame(true_values = true_values,
                                 median_values = median_values,
                                 L_values = L_values,
                                 H_values = H_values,
                                 Variable = slot,
                                 Age = c(1:19),
                                 MES = true_values[1])
  
}

full_res = do.call(rbind, results_list)

save(full_res, file = "output/full_res.Rda")

