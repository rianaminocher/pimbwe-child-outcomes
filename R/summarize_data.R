# produce some descriptive figures/tables

library(ggplot2)
library(rethinking)
library(openxlsx)
library(xtable)

load("processed_data/ppl.robj")
load("processed_data/obs.robj")

load("processed_data/mother_dead.robj")
load("processed_data/father_dead.robj")

load("processed_data/father_unmarried.robj")
load("processed_data/father_married_to_notmother_monogamy.robj")
load("processed_data/father_married_to_notmother_polygyny.robj")
load("processed_data/father_married_to_mother_polygyny.robj")

load("processed_data/mother_unmarried.robj")
load("processed_data/mother_married_to_notfather.robj")
load("processed_data/mother_married_to_father_with_cowife.robj")

n <- nrow(ppl)
n == 3693

# check if kids with missing parents have any missing predictor data

any(is.na(ppl[is.na(ppl$father_id), ]$dob)) # all kids have birth-years
which(is.na(ppl[is.na(ppl$father_id), ]$male)) # missing sex
any(is.na(ppl[is.na(ppl$father_id), ]$birth_order)) # none are missing birth-order

ppl$missing_kid <- is.na(ppl$mother_id) | is.na(ppl$father_id)
obs$missing_kid <- ppl$missing_kid[match(obs$id, ppl$id)]

min(ppl$dob) == 1931
max(ppl$dob) == 2014

png("output/figures/hist_birthyear.png", 
    res = 250, 
    height = 1400, 
    width = 1800)

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE))

plot(table(ppl$dob),
     ylab = "number of children", 
     main = "full sample")

plot(table(ppl$dob[ppl$missing_kid]),
     ylab = "number of children", 
     main = "unknown parent")

plot(table(ppl$dob[!ppl$missing_kid]),
     ylab = "number of children", 
     main = "parental status known")

dev.off()

# how many children are observed retrospectively?
# prospective = at least one height/weight/edu obs

height_ids <- unique(obs[complete.cases(obs$height), ]$id)
weight_ids <- unique(obs[complete.cases(obs$weight), ]$id)
edu_ids <-  unique(obs[complete.cases(obs$edu), ]$id)

prosp_ids <- unique(c(height_ids, weight_ids, edu_ids))
retro_ids <- setdiff(ppl$id, prosp_ids)
length(retro_ids) + length(prosp_ids) == nrow(ppl)

any(!is.na(obs[obs$id %in% retro_ids, ]$height)) # none are observed

# when are these retro kids born?

table(ppl[ppl$id %in% retro_ids, ]$dob) # they are also born after 1995

# how many died before 1995, these are certainly retro
sum(ppl[ppl$id %in% retro_ids & !is.na(ppl$date_of_death), ]$date_of_death < 1995) # 515
# how many censored before 1995?
sum(ppl[ppl$id %in% retro_ids & !is.na(ppl$date_of_censor), ]$date_of_censor < 1995) # none

# make observations by height, weight, edu histogram

# height recorded by child 
length(unique(obs[!is.na(obs$height), ]$id)) # 1139 kids measured at least once
table(table(obs[!is.na(obs$height), ]$id)) # 501 children measured just once

# weight recorded by child
length(unique(obs[!is.na(obs$weight), ]$id)) # 1298 kids measured at least once
table(table(obs[!is.na(obs$weight), ]$id)) # 526 children measured just once

# education recorded by child 
length(unique(obs[!is.na(obs$edu), ]$id)) # 1377 kids measured at least once
table(table(obs[!is.na(obs$edu), ]$id)) # 439 children measured just once

png("output/figures/long_obs_by_child.png",
    res = 250,
    height = 1000, 
    width = 3000)

par(mfrow = c(1, 3))

var_list <- c("height", "weight", "edu")

for (i in 1:3) {
  
  hist_by_child <- tapply(obs[, var_list[i]], obs$id, function(x) sum(!is.na(x)))
  
  plot(table(hist_by_child[hist_by_child > 0]), 
       col = "slateblue", 
       xlab = "number of observations", 
       ylab = "frequency", 
       lwd = 6, 
       main = var_list[i], 
       cex.lab = 2, 
       cex.axis = 2,
       cex.main = 3)
  
}

dev.off()

# sex of child unknown

missing_sex <- ppl[is.na(ppl$male), ]
nrow(missing_sex) # 278 unknown

missing_sex$age <- NA

for (i in 1:nrow(missing_sex)) {
  
  if (!is.na(missing_sex$date_of_death[i])) {
    missing_sex$age[i] <- missing_sex$date_of_death[i] - missing_sex$dob[i] 
  }
  
  else missing_sex$age[i] <- missing_sex$date_of_censor[i] - missing_sex$dob[i]
  
}

missing_sex$baby <- ifelse(missing_sex$age <= 1, TRUE, FALSE)
sum(!missing_sex$baby) == 120
sum(missing_sex$baby) == 158

missing_sex$born_before_95 <- ifelse(missing_sex$dob < 1995, TRUE, FALSE)

table(missing_sex$baby & missing_sex$born_before_95)

any(missing_sex$id %in% obs[!is.na(obs$edu), ]$id)
any(missing_sex$id %in% obs[!is.na(obs$height), ]$id)
any(missing_sex$id %in% obs[!is.na(obs$weight), ]$id)
any(missing_sex$id %in% obs$id)

sum(missing_sex[!missing_sex$baby, ]$missing_kid) == 83
sum(!missing_sex[!missing_sex$baby, ]$missing_kid) == 37

table(missing_sex[!missing_sex$missing_kid & !missing_sex$baby, ]$born_before_95)

write.csv(missing_sex, 
          "for_monique/kids_missing_sex.csv",
          row.names = FALSE)

# to make the new fig
# make survival data as per model fit
# two routes to censorship: 
# child doesn't turn 18 before 2015 (died or unobserved)
# parent censored before 2015 

alive <- matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
  # get child dob, dod, doc
  dob <- ppl$dob[i]
  dod <- ppl$date_of_death[i]
  doc <- ppl$date_of_censor[i]
  
  # if dod is known
  if (!is.na(dod)) {
    # age at death is yod - yob + 1 # so if it dies in the same year as born = 1
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

skip <- matrix(0, nrow = n, ncol = 19)

for (i in 1:n) {
  for (a in 1:19) {
    
    if (mother_dead[i, a] == -99) {
      skip[i, a] <- 1
    }
    
    
    if (father_dead[i, a] == -99) {
      skip[i, a] <- 1
    }
  }
}

# we should distinguish:
# child censored
# parent censored

mother_state <- matrix(NA, nrow = n, ncol = 19)

for (i in 1:n) {
      
  for (a in 1:19) {
    
    # correct child censored or died
    if (alive[i, a] == -99) mother_state[i, a] <- "child cens"
    
    if (alive[i, a] != -99) {
      
      # skips are censored + unknown
      if (is.na(ppl$mother_id[i])) mother_state[i, a] <- "unknown parent"
      # if (skip[i, a] == 1) mother_state[i, a] <- "unknown parent"
      
      # if mother is dead
      if (mother_dead[i, a] == 1) mother_state[i, a] <- "dead"
      
      if (mother_dead[i, a] == 0) {
        
        if (mother_unmarried[i, a] == 1) mother_state[i, a] <- "unmarried"
        
        if (mother_unmarried[i, a] == 0) {
          
          if (mother_married_to_notfather[i, a] == 1) mother_state[i, a] <- "married to step-father"
          
          if (mother_married_to_notfather[i, a] == 0) {
            
            if (mother_married_to_father_with_cowife[i, a] == 1) mother_state[i, a] <- "married to father with co-wife"
            
            if (mother_married_to_father_with_cowife[i, a] == 0) mother_state[i, a] <- "married to father"
            
          }
          
        }
        
      }
      
    }

  }
  
}

# where when mothers are censored in 2011, they are getting 0's 
# for married/unmarried and so on until 2015
# but this shouldn't be an issue, as long as we're skipping over -99 in the mother_dead
# so check that all mother state NA's are also -99 in mother_dead
# same for fathers. essentially not a prob but check where I'm filling in 0 to 2015

all(mother_dead[is.na(mother_state)] == -99)

# then fill in as parent cens
mother_state[is.na(mother_state)] <- "parent cens"

mother_state[which(mother_state == "dead")] <- 1
mother_state[which(mother_state == "unmarried")] <- 2
mother_state[which(mother_state == "married to father")] <- 3
mother_state[which(mother_state == "married to step-father")] <- 4
mother_state[which(mother_state == "married to father with co-wife")] <- 5
mother_state[which(mother_state == "unknown parent")] <- 6
mother_state[which(mother_state == "child cens")] <- 7
mother_state[which(mother_state == "parent cens")] <- 8

mother_state <- matrix(as.numeric(mother_state), ncol = 19)

father_state <- matrix(NA, nrow = n, ncol = 19)

for (i in 1:n) {
  
  for (a in 1:19) {
    
    # correct child censored
    if (alive[i, a] == -99) father_state[i, a] <- "child cens"
    
    if (alive[i, a] != -99) {
      
      # skips are censored + unknown
      if (is.na(ppl$father_id[i])) father_state[i, a] <- "unknown parent"
      
      if (father_dead[i, a] == 1) father_state[i, a] <- "dead"
      
      if (father_dead[i, a] == 0) {
        
        if (father_unmarried[i, a] == 1) father_state[i, a] <- "unmarried"
        
        if (father_unmarried[i, a] == 0) {
          
          if (father_married_to_mother_polygyny[i, a] == 0) father_state[i, a] <- "to mother (m)"
          if (father_married_to_notmother_monogamy[i, a] == 1) father_state[i, a] <- "to step-mother (m)"
          if (father_married_to_notmother_polygyny[i, a] == 1) father_state[i, a] <- "to step-mother (p)"
          if (father_married_to_mother_polygyny[i, a] == 1) father_state[i, a] <- "to mother (p)"
          
        }
        
      }
      
    }
    
  }
}

table(father_state)

# check which are NAs, are all parent cens?
table(father_dead[is.na(father_state)])

father_state[which(is.na(father_state))] <- "parent cens"

father_state[which(father_state == "dead")] <- 1
father_state[which(father_state == "unmarried")] <- 2
father_state[which(father_state == "to step-mother (m)")] <- 3
father_state[which(father_state == "to step-mother (p)")] <- 4
father_state[which(father_state == "to mother (m)")] <- 5
father_state[which(father_state == "to mother (p)")] <- 6
father_state[which(father_state == "unknown parent")] <- 7
father_state[which(father_state == "child cens")] <- 8
father_state[which(father_state == "parent cens")] <- 9

father_state <- matrix(as.numeric(father_state), ncol = 19)

# plot mother status

kids <- sample(1:nrow(ppl), 20)

png("output/figures/mother_status.png", 
    res = 250, 
    height = 2000, 
    width = 1600)

par(mfrow = c(1, 1))

state <- mother_state[kids, ]

cols <- c("navy",
          "goldenrod",
          "firebrick",
          "cyan4", 
          "mediumorchid", 
          "chocolate3",
          "gray", 
          "darkgray")

plot(NULL,
     xlab = "", 
     ylab = "", 
     ylim = c(0, 20),
     xlim = c(0, 18), 
     xaxt = "n", 
     yaxt = "n", 
     cex.lab = 1.5,
     xaxs = "i",
     yaxs = "i")

for (i in 1:20) {
  for (a in 1:19) {
    
    rect(xleft = 0 + (a-1), 
         xright = 1 + (a-1), 
         ybottom = 0 + (i-1), 
         ytop = 1 + (i-1), 
         col = col.alpha(cols[state[i, a]], 0.8), 
         border = col.alpha(cols[state[i, a]], 1))
    
  }
}

axis(side = 1, 
     at = 0:18, 
     labels = FALSE,
     tick = TRUE, 
     tcl = -0.2)

axis(side = 1, 
     at = c(1, 10, 18), 
     labels = TRUE, 
     cex.axis = 2)

axis(side = 2, 
     at = 0:20, 
     labels = FALSE,
     tick = TRUE, 
     tcl = -0.2)

dev.off()

# plot father status

png("output/figures/father_status.png", 
    res = 250, 
    height = 2000, 
    width = 1600)

par(mfrow = c(1, 1))

state <- father_state[kids, ]

cols <- c("navy", 
          "goldenrod", 
          "cyan4",
          "sienna4",
          "firebrick", 
          "mediumorchid",
          "chocolate3",
          "gray",
          "darkgray")

plot(NULL,
     xlab = "", 
     ylab = "", 
     ylim = c(0, 20),
     xlim = c(0, 18), 
     xaxt = "n", 
     yaxt = "n", 
     cex.lab = 2.5,
     xaxs = "i",
     yaxs = "i")

for (i in 1:20) {
  for (a in 1:19) {
    
    rect(xleft = 0 + (a-1), 
         xright = 1 + (a-1), 
         ybottom = 0 + (i-1), 
         ytop = 1 + (i-1), 
         col = col.alpha(cols[state[i, a]], 0.8), 
         border = col.alpha(cols[state[i, a]], 1))
    
  }
}

axis(side = 1, 
     at = 0:18, 
     labels = FALSE,
     tick = TRUE, 
     tcl = -0.2)

axis(side = 1, 
     at = c(1, 10, 18), 
     labels = TRUE, 
     cex.axis = 2)

axis(side = 2, 
     at = 0:20, 
     labels = FALSE,
     tick = TRUE, 
     tcl = -0.2)

dev.off()

# table summarizing child years in different parent states

# mother states

# mother: dead = 1, unmarried = 2, married to f = 3, married to sf = 4, married to f with cw = 5, unknown = 6

mother <- data.frame("dead" = apply(mother_state, 2, function(x) length(which(x == 1))),
                     "unmarried" = apply(mother_state, 2, function(x) length(which(x == 2))),
                     "married_to_stepfather" = apply(mother_state, 2, function(x) length(which(x == 4))),
                     "married_to_father" = apply(mother_state, 2, function(x) length(which(x == 3))),
                     "married_to_father_with_cowife" = apply(mother_state, 2, function(x) length(which(x == 5))), 
                     "missing" = apply(mother_state, 2, function(x) length(which(x == 6| x == 7 | x == 8))))

all(apply(mother, 1, sum) == 3693) 

# group ages 1-5, 6-10, 11-15, 16-19

mother <- data.frame("0-5" = apply(mother[1:5, ], 2, sum),
                     "6-10" = apply(mother[6:10, ], 2, sum),
                     "11-15" = apply(mother[11:15, ], 2, sum),
                     "16-18" = apply(mother[16:19, ], 2, sum))

mother <- t(mother)
mother <- mother[, 1:5]

print(xtable(mother), 
      file = "output/tables/mother_age_bin_sums.txt", 
      only.contents = TRUE,
      sanitize.text.function = function(x) {x})

# father states

# dead = 1, unm = 2, step m = 3, step p = 4, mother m = 5, mother p = 6, unk = 7

father <- data.frame("dead" = apply(father_state, 2, function(x) length(which(x == 1))),
                     "unmarried" = apply(father_state, 2, function(x) length(which(x == 2))),
                     "married_to_stepmother_m" = apply(father_state, 2, function(x) length(which(x == 3))),
                     "married_to_stepmother_p" = apply(father_state, 2, function(x) length(which(x == 4))),
                     "married_to_mother_m" = apply(father_state, 2, function(x) length(which(x == 5))), 
                     "married_to_mother_p" = apply(father_state, 2, function(x) length(which(x == 6))), 
                     "dead/censored/missing" = apply(father_state, 2, function(x) length(which(x == 7 | x == 8 |  x == 9))))

# group ages 1-5, 6-10, 11-15, 16-19

all(apply(father, 1, sum) == 3693) # check they sum to 3693

father <- data.frame("0-5" = apply(father[1:5, ], 2, sum),
                     "6-10" = apply(father[6:10, ], 2, sum),
                     "11-15" = apply(father[11:15, ], 2, sum),
                     "16-18" = apply(father[16:19, ], 2, sum))

father <- t(father)
father <- father[, 1:6]

print(xtable(father), 
      file = "output/tables/father_age_bin_sums.txt", 
      only.contents = TRUE,
      sanitize.text.function = function(x) {x})

# get data descriptives

length(unique(ppl$mother_id[!is.na(ppl$mother_id)]))
length(unique(ppl$father_id[!is.na(ppl$father_id)]))

# do children have both parent missing

length(which(is.na(ppl$mother_id)))
length(which(is.na(ppl$father_id)))

table(ppl$male)
length(which(is.na(ppl$male)))

table(ppl$twin)

table(ppl$birth_order)

# reproductive stat data
# check how pregnant/breastfeeding relates to height weight data

# first plot height by age data

data <- obs[, c("id", "height", "year")]
data$male <- ppl[match(obs$id, ppl$id), ]$male
data$dob <- ppl[match(obs$id, ppl$id), ]$dob
data$age <- (data$year - data$dob) + 1

# check whether all negative ages are missing data values
all(which(data$age < 0) %in% which(is.na(data$height)))

# drop missing values
data <- data[!is.na(data$height), ]

# no kids of unknown sex have height measures
any(which(is.na(data$male)))

data$male[data$male == "0"] <- "female"
data$male[data$male == "1"] <- "male"

# plot age x height for all ids

out <- ggplot(data,
              aes(x = age,
                  y = height,
                  color = male),
              fill = male) +
  
  geom_point() + 
  
  theme_bw() +
  
  labs(y = "height (cm)", 
       x = "age") +
  
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.position = "right", 
        legend.title = element_blank(),
        plot.title = element_text(size = 13, face = "italic")) +
  
  scale_color_manual(values = c("male" = col.alpha("darkblue", 0.8),
                                "female" = col.alpha("lightblue", 0.5)))

png("output/figures/heightforage.png",
    res = 250,
    height = 1200, 
    width = 2000)

print(out)

dev.off()

# tables of outcome data

load("processed_data/obs.robj")

obs$dob <- ppl$dob[match(obs$id, ppl$id)]
obs$age <- (obs$year - obs$dob) + 1

obs <- obs[obs$age < 20 & obs$age > 0, ]

tab <- data.frame(age = c("0---5", "6---10", "11---15", "16---18"))

tab$height[1] <- nrow(obs[obs$age > 0 & obs$age <= 5 & 
                            !is.na(obs$height), ])

tab$height[2] <- nrow(obs[obs$age > 5 & obs$age <= 10
                          & !is.na(obs$height), ])

tab$height[3] <- nrow(obs[obs$age > 10 & obs$age <= 15
                          & !is.na(obs$height), ])

tab$height[4] <- nrow(obs[obs$age > 15 & obs$age <= 19
                          & !is.na(obs$height), ])

# check
tapply(obs$height, obs$age, function(x) length(which(!is.na(x))))

tab$weight[1] <- nrow(obs[obs$age > 0 & obs$age <= 5 & 
                            !is.na(obs$weight), ])

tab$weight[2] <- nrow(obs[obs$age > 5 & obs$age <= 10
                          & !is.na(obs$weight), ])

tab$weight[3] <- nrow(obs[obs$age > 10 & obs$age <= 15
                          & !is.na(obs$weight), ])

tab$weight[4] <- nrow(obs[obs$age > 15 & obs$age <= 19
                          & !is.na(obs$weight), ])

# check weight
tapply(obs$weight, obs$age, function(x) length(which(!is.na(x))))

tab$bmi <- tab$height

tab$education[1] <- nrow(obs[obs$age > 0 & obs$age <= 5 & 
                            !is.na(obs$edu), ])

tab$education[2] <- nrow(obs[obs$age > 5 & obs$age <= 10
                          & !is.na(obs$edu), ])

tab$education[3] <- nrow(obs[obs$age > 10 & obs$age <= 15
                          & !is.na(obs$edu), ])

tab$education[4] <- nrow(obs[obs$age > 15 & obs$age <= 19
                          & !is.na(obs$edu), ])
# check edu
tapply(obs$edu, obs$age, function(x) length(which(!is.na(x))))

print(xtable(tab), 
      file = "output/tables/outcome_age_bins.txt", 
      only.contents = TRUE,
      sanitize.text.function = function(x) {x})

# make figures/tables on the structure of the dataset

# birth-death by birth year
# only for 18 years of life

d <- ppl[, c("id", "dob", "date_of_death", "date_of_censor")]
d <- d[order(d$dob), ]
d$year_at_18 <- d$dob + 18
d$died_or_censor_before_18 <- ifelse(is.na(d$date_of_death), d$date_of_censor, d$date_of_death)
d$died_or_censor_before_18 <- d$died_or_censor_before_18 < d$year_at_18

d$date_of_death_or_censor <- ifelse(is.na(d$date_of_death), d$date_of_censor, d$date_of_death)

pdf("output/figures/lifelines.pdf",
    height = 25,
    width = 20)

par(mar = c(6, 6, 2, 2))

plot(NULL,
     xlim = c(1930, 2014), 
     ylim = c(1, 3693),
     xlab = "year",
     ylab = "child",
     cex.lab = 2.5,
     cex.axis = 2.5)

for (i in 1:nrow(d)) {
  
  arrows(y0 = i,
         y1 = i, 
         x0 = d$dob[i], 
         x1 = ifelse(d$died_or_censor_before_18[i] == TRUE, d$date_of_death_or_censor[i], d$year_at_18[i]),
         length = 0,
         col = col.alpha("grey30", 0.4))
  
}

abline(v = 1995, col = "red", lwd = 2)

dev.off()

# outcome 1 survival

alive <- matrix(nrow = n, ncol = 19)

for (i in 1:n) {
  
  # get child dob, dod, doc
  dob <- ppl$dob[i]
  dod <- ppl$date_of_death[i]
  doc <- ppl$date_of_censor[i]
  
  # if dod is known
  if (!is.na(dod)) {
    # age at death is yod - yob + 1 # so if it dies in the same year as born = 1
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

alive[alive == -99] <- NA

apply(alive, 2, function(x) sum(x, na.rm = TRUE))

height <- matrix(nrow = n, ncol = 19)

id <- ppl$id

for (i in 1:n) {
  
  # get child dob
  dob <- ppl$dob[i]
  # subset to obs for that child
  tmp <- obs[obs$id == id[i], ]
  # subset to non-NA obs
  tmp <- tmp[which(!is.na(tmp$height)), ]
  # fill in the "age from 1:19"
  tmp$age <- (tmp$year - dob) + 1
  
  if (nrow(tmp) > 0) {
    
    for (j in 1:nrow(tmp)) {
      
      if (tmp$age[j] <= 19) {
        
        height[i, tmp$age[j]] <- tmp$height[j]
        
      }
    }
  }
}

# how many total height measures?
sum(apply(height, 2, function(x) length(which(!is.na(x)))))

# how many inds have an edu measurement?
height_measures <- apply(height, 1, function(x) length(which(!is.na(x))))
sum(height_measures > 0)
table(height_measures)

# BMI is the same as height, because we need height to calc BMI

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

# how many total measurements?

sum(apply(edu, 2, function(x) length(which(!is.na(x)))))

# how many inds have an edu measurement?
edu_measures <- apply(edu, 1, function(x) length(which(!is.na(x))))
sum(edu_measures > 0)

# print summary table

tab <- data.frame(survival = apply(alive, 2, function(x) length(which(!is.na(x)))),
                  height = apply(height, 2, function(x) length(which(!is.na(x)))),
                  education = apply(edu, 2, function(x) length(which(!is.na(x)))))

print(xtable(tab), 
      file = "output/tables/sample_by_age.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# calculate mother/father/male sums for height & edu

ids_with_height <- which(apply(height, 1, function(x) length(which(!is.na(x))) > 0)) 

length(unique(ppl[ppl$id %in% ids_with_height, ]$father_id)) - 1 # 245 fathers
length(unique(ppl[ppl$id %in% ids_with_height, ]$mother_id)) - 1 # 324 mothers

length(which(is.na(ppl[ppl$id %in% ids_with_height, ]$father_id))) # 140
length(which(is.na(ppl[ppl$id %in% ids_with_height, ]$mother_id))) # 13

sum(ppl[ppl$id %in% ids_with_height, ]$male) # 412
sum(!ppl[ppl$id %in% ids_with_height, ]$male) # 469

ids_with_edu <- which(apply(edu, 1, function(x) length(which(!is.na(x))) > 0)) 

length(unique(ppl[ppl$id %in% ids_with_edu, ]$father_id)) - 1 # 314 fathers
length(unique(ppl[ppl$id %in% ids_with_edu, ]$mother_id)) - 1 # 423 mothers

length(which(is.na(ppl[ppl$id %in% ids_with_edu, ]$father_id))) # 279
length(which(is.na(ppl[ppl$id %in% ids_with_edu, ]$mother_id))) # 53

sum(ppl[ppl$id %in% ids_with_edu, ]$male) # 667
sum(!ppl[ppl$id %in% ids_with_edu, ]$male) # 703

# why do we have more measures in the obs table than when we split by age
# is this only because of duplicated measurements for a single age?
# my loop is always using the second measure made

length(which(!is.na(obs[obs$id %in% ids_with_edu & obs$age < 20, ]$edu)))
length(which(!is.na(obs[obs$id %in% ids_with_height & obs$age < 20, ]$height)))

sum(apply(edu, 2, function(x) length(which(!is.na(x)))))
sum(apply(height, 2, function(x) length(which(!is.na(x)))))

# check how many obs have a duplicated id/year combo
obs$idyear <- paste(obs$id, obs$year)
sum(duplicated(obs[obs$id %in% ids_with_edu, ]$idyear))
3845 - 3693

sum(duplicated(obs[obs$id %in% ids_with_height, ]$idyear))
1922 - 1744
