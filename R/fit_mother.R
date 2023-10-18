# fit all mother models

# read processed data objects

load("processed_data/ppl.robj")
load("processed_data/obs.robj")

load("processed_data/father_dead.robj")
load("processed_data/mother_dead.robj")
load("processed_data/mother_unmarried.robj")
load("processed_data/mother_married_to_notfather.robj")
load("processed_data/mother_married_to_father_with_cowife.robj")

# make data list for stan

n <- nrow(ppl)

# 1: fit survival

any(is.na(ppl$date_of_censor) & is.na(ppl$date_of_death))

alive <- matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
    dob <- ppl$dob[i]
    dod <- ppl$date_of_death[i]
    doc <- ppl$date_of_censor[i]
    
    if (!is.na(dod)) {
      aod <- (dod - dob) + 1
      
      if (aod > 19) {
        
        alive[i, ] <- 1
        
      }
      
      else {
        
          alive[i, 1:aod] <- 1
          alive[i, aod] <- 0
          
          if (aod < 19) {
            alive[i, (aod+1):19] <- -99 
          }
          
        }
        
    }
    
    if (!is.na(doc)) {
      
      aoc <- (doc - dob) + 1
      
      if (aoc > 19) {
        
        alive[i, ] <- 1
        
      }
      
      else {
        
        alive[i, 1:aoc] <- 1
        alive[i, aoc] <- 1
        
        if (aoc < 19) {
          alive[i, (aoc+1):19] <- -99 
        }
        
      }
      
    }

}

years <- 1930:2015
n_years <- length(years)

dob <- ppl$dob
dob <- match(dob, years)

male <- ppl$male
male[is.na(male)] <- -99

birthorder <- ppl$birth_order
birthorder[is.na(birthorder)] <- 16

twin <- ppl$twin

father_id <- ppl$father_id
mother_id <- ppl$mother_id

# assign parent ids for all the NA values
# these are the "external" parents
# each gets their own unique ID to be conservative

father_id <- as.integer(as.factor(father_id))
start_unk_id <- max(father_id, na.rm = TRUE) + 1 # 461 is the max ID, so start from there + 1
unk_id_n <- length(which(is.na(father_id)))
father_id[is.na(father_id)] <- start_unk_id:(start_unk_id + unk_id_n - 1)

mother_id <- as.integer(as.factor(mother_id))
start_unk_id <- max(mother_id, na.rm = TRUE) + 1 # 461 is the max ID, so start from there + 1
unk_id_n <- length(which(is.na(mother_id)))
mother_id[is.na(mother_id)] <- start_unk_id:(start_unk_id + unk_id_n - 1)

I <- length(unique(mother_id))
L <- length(unique(father_id))

data <- list(N = n, 
             A = 19, 
             B = 16, 
             Y = 86, 
             I = I,
             L = L,
             alive = alive, 
             birthorder = birthorder, 
             dob = dob, 
             male = male, 
             twin = twin, 
             mother_id = mother_id,
             father_id = father_id,
             father_dead = father_dead, 
             mother_dead = mother_dead, 
             mother_unmarried = mother_unmarried, 
             mother_married_to_notfather = mother_married_to_notfather,
             mother_married_to_father_with_cowife = mother_married_to_father_with_cowife)

skip <- matrix(0, nrow = n, ncol = 19)

for (i in 1:n) {
  for (a in 1:19) {
    
    if (data$mother_dead[i, a] == -99) {
      skip[i, a] <- 1
    }
    
    
    if (data$father_dead[i, a] == -99) {
      skip[i, a] <- 1
    }
  }
}

data$skip <- skip

# compile model

m <- cmdstan_model("stan/survival_mother_twin.stan")

# fit model

fit <- m$sample(data = data, 
                chains = 1, 
                refresh=1,
                iter_warmup=1000,
                iter_sampling=1000,
                parallel_chains = 1, 
                adapt_delta = 0.95,
                max_treedepth = 12,
                init = 0)

# save fit

fit <- rstan::read_stan_csv(fit$output_files())
saveRDS(fit, "stanfits/mother_survival.rds")

# 2: fit height model 

height <- matrix(nrow = n, ncol = 19)

id <- ppl$id

for (i in 1:n) {
  
  dob <- ppl$dob[i]
  tmp <- obs[obs$id == id[i], ]
  tmp <- tmp[which(!is.na(tmp$height)), ]
  tmp$age <- (tmp$year - dob) + 1
  
  if (nrow(tmp) > 0) {
    
    for (j in 1:nrow(tmp)) {
      
      if (tmp$age[j] <= 19) {
        
        height[i, tmp$age[j]] <- tmp$height[j]
        
      }
    }
  }
}

height[is.na(height)] <- -99

data$alive <- height
names(data)[names(data) == "alive"] <- "height"

# are any height measures non NA for missing sex kids?
all(data$height[which(data$male == -99), ] == -99) # all are -99

# revise Y for year-effects estimated
# earliest year of obs in the data is earliest yob for somebody who could have received a height measure
# between 1995-2014 and is under 18 during this period
# i.e., 1976
years_short <- 1976:2014
data$Y <- length(years_short)

# also redo dobs so that 1976 is index 1
data$dob <- ppl$dob
data$dob <- (data$dob - 1976) + 1
table(data$dob)
data$dob[data$dob < 1] <- -99

# compile model

m <- cmdstan_model("stan/height_mother.stan")

# fit model

fit <- m$sample(data = data, 
                chains = 4, 
                parallel_chains = 4, 
                adapt_delta = 0.95,
                max_treedepth = 13,
                init = 0)

# save fit

fit <- rstan::read_stan_csv(fit$output_files())
saveRDS(fit, "stanfits/mother_height.rds")

# 3: fit weight model

weight <- matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
  dob <- ppl$dob[i]
  tmp <- obs[obs$id == id[i], ]
  tmp <- tmp[which(!is.na(tmp$weight)), ]
  tmp$age <- (tmp$year - dob) + 1
  
  if (nrow(tmp) > 0) {
    
    for (j in 1:nrow(tmp)) {
      
      if (tmp$age[j] <= 19) {
        
        weight[i, tmp$age[j]] <- tmp$weight[j]
        
      }
    }
  }
}

# convert weight to BMI

# height -99 to NA for the calculation
height[height == -99] <- NA

# height from cm to m
height <- height / 100

# bmi -> kg/m^2
bmi <- weight / (height*height)

# there are definitely some height outliers which may explain bmis < 12
# one bmi > 28

data$height <- bmi
names(data)[names(data) == "height"] <- "weight"

range(data$weight, na.rm = TRUE) 
# assume 12-30 is a reasonable range based on WHO growth standards

# flag outliers:
which(data$weight < 11, arr.ind = TRUE) # child 13 at age 6
which(data$weight > 30, arr.ind = TRUE) # child 3294 at age 2

# remove for now
data$weight[which(data$weight < 11)] <- NA
data$weight[which(data$weight > 30)] <- NA

data$weight[is.na(data$weight)] <- -99

# compile model

m <- cmdstan_model("stan/weight_mother.stan")

# fit model

fit <- m$sample(data = data, 
                chains = 4, 
                parallel_chains = 4, 
                adapt_delta = 0.95,
                max_treedepth = 13,
                init = 0)

# save fit

fit <- rstan::read_stan_csv(fit$output_files())
saveRDS(fit, "stanfits/mother_weight.rds")

# 4: fit education model

edu <- matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
  dob <- ppl$dob[i]
  tmp <- obs[obs$id == id[i], ]
  tmp <- tmp[which(!is.na(tmp$edu)), ]
  tmp$age <- (tmp$year - dob) + 1
  
  if (nrow(tmp) > 0) {
    
    for (j in 1:nrow(tmp)) {
      
      if (tmp$age[j] <= 19) {
        
        edu[i, tmp$age[j]] <- tmp$edu[j]
        
      }
    }
  }
}

edu[is.na(edu)] <- -99

data$weight <- edu
names(data)[names(data) == "weight"] <- "edu"
data$height <- NULL

# are any edu measures non NA for missing sex kids?
all(data$edu[which(data$male == -99), ] == -99) # all are -99

# trim edu data to ages which have data
data$A <- 15
data$edu <- data$edu[, 5:19]
data$mother_dead <- data$mother_dead[, 5:19]
data$father_dead <- data$father_dead[, 5:19]
data$mother_unmarried <- data$mother_unmarried[, 5:19]
data$mother_married_to_notfather <- data$mother_married_to_notfather[, 5:19]
data$mother_married_to_father_with_cowife <- data$mother_married_to_father_with_cowife[, 5:19]
data$skip <- data$skip[, 5:19]
str(data)

# compile model

m <- cmdstan_model("stan/education_hurdle_mother.stan")

# fit model

fit <- m$sample(data = data, 
                chains = 4, 
                parallel_chains = 4, 
                adapt_delta = 0.95,
                max_treedepth = 13,
                init = 0)

# save fit

fit <- rstan::read_stan_csv(fit$output_files())
saveRDS(fit, "stanfits/mother_education.rds")
