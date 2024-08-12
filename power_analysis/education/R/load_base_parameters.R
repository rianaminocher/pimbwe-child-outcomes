# set initial param values from a fitted model 

real_fit = readRDS("stanfits/father_education.rds")
post = as_draws_rvars(real_fit)

# extract param values
alpha_l = median(draws_of(post$alpha_l))
alpha_t = median(draws_of(post$alpha_t))

a_bo_l = apply(draws_of(post$a_bo_l), 2, median)
a_year_l = apply(draws_of(post$a_year_l), 2, median)

a_bo_t = apply(draws_of(post$a_bo_t), 2, median)
a_year_t = apply(draws_of(post$a_year_t), 2, median)

a_age_l = matrix(NA, nrow = 10, ncol = 15)
a_age_t = matrix(NA, nrow = 10, ncol = 15)

# age-offset
a_age_l[1, ] = apply(draws_of(post$a_age_l)[ , 1, ], 2, median)
# male effect
a_age_l[2, ] = apply(draws_of(post$a_age_l)[ , 2, ], 2, median) 
# twin effect
a_age_l[3, ] = apply(draws_of(post$a_age_l)[ , 3, ], 2, median)
# mother_dead
a_age_l[4, ] = apply(draws_of(post$a_age_l)[ , 4, ], 2, median)
# father_dead
a_age_l[5, ] = apply(draws_of(post$a_age_l)[ , 5, ], 2, median)
# father_unmarried
a_age_l[6, ] = apply(draws_of(post$a_age_l)[ , 6, ], 2, median)
# father_married_to_notmother_monogamy
a_age_l[7, ] = apply(draws_of(post$a_age_l)[ , 7, ], 2, median)
# father_married_to_notmother_polygyny 
a_age_l[8, ] = apply(draws_of(post$a_age_l)[ , 8, ], 2, median)
# father_married_to_mother_polygyny
a_age_l[9, ] = apply(draws_of(post$a_age_l)[ , 9, ], 2, median)
# father_unknown
a_age_l[10, ] = apply(draws_of(post$a_age_l)[ , 10, ], 2, median)

# age-offset
a_age_t[1, ] = apply(draws_of(post$a_age_t)[ , 1, ], 2, median)
# male effect
a_age_t[2, ] = apply(draws_of(post$a_age_t)[ , 2, ], 2, median) 
# twin effect
a_age_t[3, ] = apply(draws_of(post$a_age_t)[ , 3, ], 2, median)
# mother_dead
a_age_t[4, ] = apply(draws_of(post$a_age_t)[ , 4, ], 2, median)
# father_dead
a_age_t[5, ] = apply(draws_of(post$a_age_t)[ , 5, ], 2, median)
# father_unmarried
a_age_t[6, ] = apply(draws_of(post$a_age_t)[ , 6, ], 2, median)
# father_married_to_notmother_monogamy
a_age_t[7, ] = apply(draws_of(post$a_age_t)[ , 7, ], 2, median)
# father_married_to_notmother_polygyny 
a_age_t[8, ] = apply(draws_of(post$a_age_t)[ , 8, ], 2, median)
# father_married_to_mother_polygyny
a_age_t[9, ] = apply(draws_of(post$a_age_t)[ , 9, ], 2, median)
# father_unknown
a_age_t[10, ] = apply(draws_of(post$a_age_t)[ , 10, ], 2, median)

a_mother_l = apply(draws_of(post$a_mother_l), 2, median)
a_father_l = apply(draws_of(post$a_father_l), 2, median)

a_mother_t = apply(draws_of(post$a_mother_t), 2, median)
a_father_t = apply(draws_of(post$a_father_t), 2, median)

mother_sigma_l = median(draws_of(post$mother_sigma_l))
father_sigma_l = median(draws_of(post$father_sigma_l))

mother_sigma_t = median(draws_of(post$mother_sigma_t))
father_sigma_t = median(draws_of(post$father_sigma_t))

# make model_param list

model_parameters = list(alpha_l = alpha_l,
                         alpha_t = alpha_t,
                         a_bo_l = a_bo_l,
                         a_bo_t = a_bo_t,
                         a_year_l = a_year_l,
                         a_year_t = a_year_t,
                         a_age_l = a_age_l,
                         a_age_t = a_age_t,
                         a_mother_l = a_mother_l,
                         a_mother_t = a_mother_t,
                         a_father_l = a_father_l,
                         a_father_t = a_father_t,
                         mother_sigma_l = mother_sigma_l,
                         mother_sigma_t = mother_sigma_t,
                         father_sigma_l = father_sigma_l,
                         father_sigma_t = father_sigma_t)
