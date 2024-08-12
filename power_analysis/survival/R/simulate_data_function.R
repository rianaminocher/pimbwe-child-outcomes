# simulate survival data 

simulate_data = function(model_parameters) {
  
  source("R/load_base_data.R") # re-load base to replace
  
  # simulate outcome
  p_survival = matrix(-99, nrow = data$N, ncol = data$A)
  survival = matrix(-99, nrow = data$N, ncol = data$A)
  
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
      
      # if father is known
      if (skip[i, a] == 0) {
        
        # p of survival of ind i at age a
        # is a function of age and time-varying parent marriage states
        # and birth-order, birth-year
        p_survival[i, a] = logistic(
            alpha + 
            a_bo[birthorder[i]] +
            a_year[dob[i]] + 
            a_mother[mother_id[i]] * mother_sigma +
            a_father[father_id[i]] * father_sigma +
            a_age[1, a] + 
            a_age[2, a] * male[i] + 
            a_age[3, a] * twin[i] + 
            a_age[4, a] * mother_dead[i, a] +
            a_age[5, a] * father_dead[i, a] +
            a_age[6, a] * father_unmarried[i, a] +
            a_age[7, a] * father_married_to_notmother_monogamy[i, a] +
            a_age[8, a] * father_married_to_notmother_polygyny[i, a] +
            a_age[9, a] * father_married_to_mother_polygyny[i, a]
        )
        
      }
      
      # if father is unknown
      if (skip[i, a] == 1) {
        
        p_survival[i, a] = logistic(
          alpha + 
            a_bo[birthorder[i]] +
            a_year[dob[i]] + 
            a_mother[mother_id[i]] * mother_sigma +
            a_father[father_id[i]] * father_sigma +
            a_age[1, a] + 
            a_age[2, a] * male[i] + 
            a_age[3, a] * twin[i] + 
            a_age[10, a]
        )
        
      }
      
    }
  }
  
  for (i in 1:N) {
    for (a in 1:alive[i]) {
      
      # if its in the first year of life
      if (a == 1) {
        survival[i, a] = rbinom(1, size = 1, prob = p_survival[i, a])
      } 
      
      else {
        
        # if later, conditional on survival to the prev year
        if (survival[i, a-1] == 1) {
          survival[i, a] = rbinom(1, size = 1, prob = p_survival[i, a])
        } else {
          survival[i, a] = -99
        }
        
      }        
    }
  }
  
  data$alive = survival
  data$male = male
  return(data)
  
}
