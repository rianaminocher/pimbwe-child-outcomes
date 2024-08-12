# update parameter values

# slot = which a_age vector it is
# new vector = the slope

modify_age_parameter_set = function (model_parameters, new_vector, slot) {
  
  model_parameters$a_age[slot, ] = new_vector
  return(model_parameters)
  
}
