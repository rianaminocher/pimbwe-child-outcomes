# simulate education data 

 zerotruncated_poisson_rng = function(lambda){
    u = runif(1, exp(-lambda), 1)
    t = -log(u)
    k = 1 + rpois(1, lambda - t)
    return(k) 
  }

  the_sneaky_coin = function(theta){
    x = rbinom(1, size=1, prob=theta)
    y = 1 - x 
    return(x)
  }

  zerotruncated_poisson_rng(2)

simulate_data = function(model_parameters) {
  
  source("R/load_base_data.R") # re-load base to replace
  
  # simulate outcome
  theta = matrix(-99, nrow = data$N, ncol = data$A)
  lambda = matrix(-99, nrow = data$N, ncol = data$A)
  p_edu = matrix(-99, nrow = data$N, ncol = data$A)
    
  # add "data" objects to the environment
  for (i in 1:length(data)) {
    assign(names(data)[i], unlist(data[[i]]))
  }

  # add "parameters" objects to the environment
  for (i in 1:length(model_parameters)) {
    assign(names(model_parameters)[i], unlist(model_parameters[[i]]))
  }
  
  # impute male
  male[male == -99] = rbinom(length(male[male == -99]), 1, 0.5)
  
  for (i in 1:N) {
    for (a in 1:alive[i]) { # for ages until aoc[i]
      
      if (dob[i] != -99) { # for children born after 1976
        if (edu[i,a] != -99) { # for kids with actual data

        # if father is known
        if (skip[i, a] == 0) {
          
          # p of survival of ind i at age a
          # is a function of age and time-varying parent marriage states
          # and birth-order, birth-year
          theta[i, a] = inv_logit(
              alpha_t + 
              a_bo_t[birthorder[i]] +
              a_year_t[dob[i]] + 
              a_mother_t[mother_id[i]] * mother_sigma_t +
              a_father_t[father_id[i]] * father_sigma_t +
              a_age_t[1, a] + 
              a_age_t[2, a] * male[i] + 
              a_age_t[3, a] * twin[i] + 
              a_age_t[4, a] * mother_dead[i, a] +
              a_age_t[5, a] * father_dead[i, a] +
              a_age_t[6, a] * father_unmarried[i, a] +
              a_age_t[7, a] * father_married_to_notmother_monogamy[i, a] +
              a_age_t[8, a] * father_married_to_notmother_polygyny[i, a] +
              a_age_t[9, a] * father_married_to_mother_polygyny[i, a]
          )
          
          lambda[i, a] = exp(
            alpha_l + 
              a_bo_l[birthorder[i]] +
              a_year_l[dob[i]] + 
              a_mother_l[mother_id[i]] * mother_sigma_l +
              a_father_l[father_id[i]] * father_sigma_l +
              a_age_l[1, a] + 
              a_age_l[2, a] * male[i] + 
              a_age_l[3, a] * twin[i] + 
              a_age_l[4, a] * mother_dead[i, a] +
              a_age_l[5, a] * father_dead[i, a] +
              a_age_l[6, a] * father_unmarried[i, a] +
              a_age_l[7, a] * father_married_to_notmother_monogamy[i, a] +
              a_age_l[8, a] * father_married_to_notmother_polygyny[i, a] +
              a_age_l[9, a] * father_married_to_mother_polygyny[i, a]
          )  
          
        }
        
        # if father is unknown
        if (skip[i, a] == 1) {
          
          theta[i, a] = inv_logit(
              alpha_t + 
              a_bo_t[birthorder[i]] +
              a_year_t[dob[i]] + 
              a_mother_t[mother_id[i]] * mother_sigma_t +
              a_father_t[father_id[i]] * father_sigma_t +
              a_age_t[1, a] + 
              a_age_t[2, a] * male[i] + 
              a_age_t[3, a] * twin[i] + 
              a_age_t[10, a]
          )
          
          lambda[i, a] = exp(
              alpha_l + 
              a_bo_l[birthorder[i]] +
              a_year_l[dob[i]] + 
              a_mother_l[mother_id[i]] * mother_sigma_l +
              a_father_l[father_id[i]] * father_sigma_l +
              a_age_l[1, a] + 
              a_age_l[2, a] * male[i] + 
              a_age_l[3, a] * twin[i] + 
              a_age_l[10, a]
          )
          
        }
        
      }}
    }
  } # simulate lambda/theta
      
  for (i in 1:N) {
    for (a in 1:alive[i]) {
      
      if (theta[i, a] != -99) {
        if (edu[i,a] != -99) {
        
        p_edu[i, a] = the_sneaky_coin(theta[i, a])*zerotruncated_poisson_rng(lambda[i, a])
        
      }}
      
    }
  }
  
  data$edu = p_edu
  data$male = male
  return(data)
  
}
