# set initial param values from a fitted model 
real_fit = readRDS("stanfits/father_height.rds")
post = as_draws_rvars(real_fit)

# extract param values
alpha = median(draws_of(post$alpha))

a_bo = apply(draws_of(post$a_bo), 2, median)
a_year = apply(draws_of(post$a_year), 2, median)

a_age = matrix(NA, nrow = 10, ncol = 19)

# age-offset
a_age[1, ] = apply(draws_of(post$a_age)[ , 1, ], 2, median)
# male effect
a_age[2, ] = apply(draws_of(post$a_age)[ , 2, ], 2, median) 
# twin effect
a_age[3, ] = apply(draws_of(post$a_age)[ , 3, ], 2, median)
# mother_dead
a_age[4, ] = apply(draws_of(post$a_age)[ , 4, ], 2, median)
# father_dead
a_age[5, ] = apply(draws_of(post$a_age)[ , 5, ], 2, median)
# father_unmarried
a_age[6, ] = apply(draws_of(post$a_age)[ , 6, ], 2, median)
# father_married_to_notmother_monogamy
a_age[7, ] = apply(draws_of(post$a_age)[ , 7, ], 2, median)
# father_married_to_notmother_polygyny 
a_age[8, ] = apply(draws_of(post$a_age)[ , 8, ], 2, median)
# father_married_to_mother_polygyny
a_age[9, ] = apply(draws_of(post$a_age)[ , 9, ], 2, median)
# father_unknown
a_age[10, ] = apply(draws_of(post$a_age)[ , 10, ], 2, median)


a_mother = apply(draws_of(post$a_mother), 2, median)
a_father = apply(draws_of(post$a_father), 2, median)

mother_sigma = median(draws_of(post$mother_sigma))
father_sigma = median(draws_of(post$father_sigma))

scale = median(draws_of(post$scale))

# make model_param list

model_parameters = list(alpha = alpha,
                         a_bo = a_bo,
                         a_year = a_year,
                         a_age = a_age,
                         a_mother = a_mother,
                         a_father = a_father,
                         mother_sigma = mother_sigma,
                         father_sigma = father_sigma, 
                         scale = scale)

