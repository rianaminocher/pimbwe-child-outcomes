# read processed data objects

load("processed_data/ppl.robj")
load("processed_data/obs.robj")

load("processed_data/mother_dead.robj")
load("processed_data/father_dead.robj")
load("processed_data/father_unmarried.robj")
load("processed_data/father_married_to_notmother_monogamy.robj")
load("processed_data/father_married_to_notmother_polygyny.robj")
load("processed_data/father_married_to_mother_polygyny.robj")

# get covariate data

n = nrow(ppl)

years = 1930:2015
n_years = length(years)

alive = matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
  # get child dob, dod, doc
  dob = ppl$dob[i]
  dod = ppl$date_of_death[i]
  doc = ppl$date_of_censor[i]
  
  # if dod is known
  if (!is.na(dod)) {
    # age at death is yod - yob + 1 # so if it dies in the same year as born = 1
    aod = (dod - dob) + 1
    
    if (aod > 19) {
      
      alive[i, ] = 1
      
    }
    
    else {
      
      alive[i, 1:aod] = 1
      alive[i, aod] = 0
      
      if (aod < 19) {
        alive[i, (aod+1):19] = -99 
      }
      
    }
    
  }
  
  if (!is.na(doc)) {
    
    aoc = (doc - dob) + 1
    
    if (aoc > 19) {
      
      alive[i, ] = 1
      
    }
    
    else {
      
      alive[i, 1:aoc] = 1
      alive[i, aoc] = 1
      
      if (aoc < 19) {
        alive[i, (aoc+1):19] = -99 
      }
      
    }
    
  }
  
}

dob = ppl$dob
dob = match(dob, years)

male = ppl$male
male[is.na(male)] = -99

birthorder = ppl$birth_order
birthorder[is.na(birthorder)] = 16

twin = ppl$twin

father_id = ppl$father_id
mother_id = ppl$mother_id

# assign parent ids for all the NA values
# these are the "external" parents
# each gets their own unique ID to be conservative

father_id = as.integer(as.factor(father_id))
start_unk_id = max(father_id, na.rm = TRUE) + 1 # 461 is the max ID, so start from there + 1
unk_id_n = length(which(is.na(father_id)))
father_id[is.na(father_id)] = start_unk_id:(start_unk_id + unk_id_n - 1)

mother_id = as.integer(as.factor(mother_id))
start_unk_id = max(mother_id, na.rm = TRUE) + 1 # 461 is the max ID, so start from there + 1
unk_id_n = length(which(is.na(mother_id)))
mother_id[is.na(mother_id)] = start_unk_id:(start_unk_id + unk_id_n - 1)

I = length(unique(mother_id))
L = length(unique(father_id))

data = list(N = n, 
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
             father_id = father_id,
             mother_id = mother_id,
             mother_dead = mother_dead, 
             father_dead = father_dead, 
             father_unmarried = father_unmarried, 
             father_married_to_notmother_monogamy = father_married_to_notmother_monogamy,
             father_married_to_notmother_polygyny = father_married_to_notmother_polygyny,
             father_married_to_mother_polygyny = father_married_to_mother_polygyny)

skip = matrix(0, nrow = n, ncol = 19)

for (i in 1:n) {
  for (a in 1:19) {
    
    if (data$mother_dead[i, a] == -99) {
      skip[i, a] = 1
    }
    
    
    if (data$father_dead[i, a] == -99) {
      skip[i, a] = 1
    }
  }
}

data$skip = skip

# store the aoc based on alive
data$alive = apply(data$alive, 1, function(x) {which(x == -99)[1]-1})
data$alive[is.na(data$alive)] = 19

