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

# check if kids with missing parents have any missing predictor data

any(is.na(ppl[is.na(ppl$father_id), ]$dob))
which(is.na(ppl[is.na(ppl$father_id), ]$male))
any(is.na(ppl[is.na(ppl$father_id), ]$birth_order))

ppl$missing_kid <- is.na(ppl$mother_id) | is.na(ppl$father_id)
obs$missing_kid <- ppl$missing_kid[match(obs$id, ppl$id)]

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

# to make the new fig
# make survival data as per model fit
# two routes to censorship: 
# child doesn't turn 18 before 2015
# parent censored before 2015 
# select out censor child-years

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

mother_state <- matrix(NA, nrow = n, ncol = 19)

for (i in 1:n) {
      
  for (a in 1:19) {
    
    # skips are censored + unknown
    if (is.na(ppl$mother_id[i])) mother_state[i, a] <- "unknown parent"
    # if (skip[i, a] == 1) mother_state[i, a] <- "unknown parent"
    
    # correct child censored
    if (alive[i, a] == -99) mother_state[i, a] <- "child cens"
    
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
mother_state[which(mother_state == "parent cens")] <- 7

mother_state <- matrix(as.numeric(mother_state), ncol = 19)

father_state <- matrix(NA, nrow = n, ncol = 19)

for (i in 1:n) {
  
  for (a in 1:19) {
    
    # skips are censored + unknown
    if (is.na(ppl$father_id[i])) father_state[i, a] <- "unknown parent"
    # if (skip[i, a] == 1) mother_state[i, a] <- "unknown parent"
    
    # correct child censored
    if (alive[i, a] == -99) father_state[i, a] <- "child cens"
    
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
father_state[which(father_state == "parent cens")] <- 8

father_state <- matrix(as.numeric(father_state), ncol = 19)

# plot mother status

png("output/figures/mother_status.png", 
    res = 250, 
    height = 2000, 
    width = 1600)

par(mfrow = c(1, 1))

state <- mother_state[sample(1:nrow(ppl), 20), ]

cols <- c("navy",
          "goldenrod",
          "firebrick",
          "cyan4", 
          "mediumorchid", 
          "chocolate3",
          "gray", 
          "gray")

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

state <- father_state[sample(1:nrow(ppl), 20), ]

cols <- c("navy", 
          "goldenrod", 
          "cyan4",
          "sienna4",
          "firebrick", 
          "mediumorchid",
          "chocolate3",
          "gray",
          "gray")

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
                     "missing" = apply(mother_state, 2, function(x) length(which(x == 6| x == 7))))

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
                     "missing" = apply(father_state, 2, function(x) length(which(x == 7 | x == 8))))

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
