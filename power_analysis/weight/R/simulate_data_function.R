# simulate survival data 

simulate_data = function(model_parameters) {
  
  source("R/load_base_data.R") # re-load base to replace
  
  # simulate outcome
  m_weight = matrix(-99, nrow = data$N, ncol = data$A)
  p_weight = matrix(-99, nrow = data$N, ncol = data$A)
  
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
        if (weight[i,a] != -99) { # for children with observed weight
          if (height[i,a] != -99) { # for children with observed height
        
        # if father is known
        if (skip[i, a] == 0) {
          
          # p of survival of ind i at age a
          # is a function of age and time-varying parent marriage states
          # and birth-order, birth-year
          m_weight[i, a] = (
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
              a_age[9, a] * father_married_to_mother_polygyny[i, a] +
              a_age[11, a] * log(height[i,a])
          )
          
        }
        
        # if father is unknown
        if (skip[i, a] == 1) {
          
          m_weight[i, a] = (
            alpha + 
              a_bo[birthorder[i]] +
              a_year[dob[i]] + 
              a_mother[mother_id[i]] * mother_sigma +
              a_father[father_id[i]] * father_sigma +
              a_age[1, a] + 
              a_age[2, a] * male[i] + 
              a_age[3, a] * twin[i] + 
              a_age[10, a] +
              a_age[11, a] * log(height[i,a])
          )
          
        }
        
      }}}
    }
  }
      
  for (i in 1:N) {
    for (a in 1:alive[i]) {
      
      if (m_weight[i, a] != -99) {
        
        p_weight[i, a] = exp(rnorm(1, m_weight[i, a], scale)) 
        
      }
      
    }
  }
  
  data$weight = p_weight
  data$male = male
  return(data)
  
}
