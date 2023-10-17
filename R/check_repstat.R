# load in repstat, ppl, obs tables 

load("processed_data/ppl.robj")
load("processed_data/obs.robj")

repstat <- read.xlsx("raw_data/repstat2.xlsx")

# load raw ppl table to get id's
d <- read.csv("raw_data/Copyofgoodkidsn3693forcodyJan202015.csv")

# scramble them back for repstat/anthrop
ppl$id_true <- d$code
obs$id_true <- ppl$id_true[match(obs$id, ppl$id)]

ppl$id <- ppl$id_true
obs$id <- obs$id_true

ppl$id_true <- NULL
obs$id_true <- NULL

data <- ppl[, c("id", "male", "dob")]

# merge year to round
anthrop <- read.xlsx("raw_data/ASP with R16(May 2018 5th corrected) best anthrops_UI_sentRiana.xlsx")

year_by_round <- unique(anthrop[, c("yrwindow", "round")])
year_by_round <- year_by_round[order(year_by_round$yrwindow), ]

# add year to repstat
repstat$year <- year_by_round[match(repstat$ROUND, year_by_round$round), ]$yrwindow

# just one NA
repstat[which(is.na(repstat$year)), ]

data <- merge(data, 
              repstat[, c("CODE", "ROUND", "repstat", "year")], 
              by.x = c("id"), 
              by.y = "CODE", 
              all.x = TRUE)

# check if any male ids have P/L/NL 
data[which(is.na(data$male)), ]$repstat
table(data$repstat, data$male)

data[which(!is.na(data$repstat) & data$male == 1), ]
# drop these
data <- data[-which(!is.na(data$repstat) & data$male == 1), ]

# calculate age
data$age <- (data$year - data$dob) + 1
# trim to only the non-NA repstat values
data <- data[!is.na(data$repstat), ]

# check out age range
table(data$age)
data[data$age < 12, ]
# 2 infants, less than 6 year old
data[data$age < 12, ]$id # ids 16405 20604 23103 28301 have L values
# and I checked their birthdays, they seem correct

data <- data[!data$age < 12, ]
data[data$age > 45, ]

# merge repstat with obs table: 

obs <- merge(obs, 
             data[, c("id", "repstat", "year")], 
             by.x = c("id","year"),
             by.y = c("id", "year"), 
             all.x = TRUE)

# calculate BMI for all obs as done for models/plots

obs$dob <- ppl$dob[match(obs$id, ppl$id)]
obs$male <- ppl$male[match(obs$id, ppl$id)]
obs$age <- (obs$year - obs$dob) + 1

obs$height <- obs$height / 100
obs$bmi <- obs$weight / (obs$height * obs$height)

# drop the ones we drop for analysis (i.e., outside 12-30 range)
obs <- obs[-which(obs$bmi > 30), ]
obs <- obs[-which(obs$bmi < 11), ]

# make the BMI plot, highlight the repstat values

obs$repstat_binary <- obs$repstat
obs$repstat_binary[which(obs$repstat_binary == "NL")] <- NA
obs$repstat_binary[which(obs$repstat_binary == "L")] <- "lactating"
obs$repstat_binary[which(obs$repstat_binary == "P")] <- "pregnant"
obs$repstat_binary[is.na(obs$repstat_binary)] <- "not pregnant or lactating"

obs <- obs[obs$male == 0, ]
obs <- obs[obs$age < 20, ]
obs <- obs[!is.na(obs$bmi), ]

# add mother dead data

load("processed_data/mother_dead.robj")

# add mother_dead to obs

obs$mother_dead <- NA
id <- ppl$id

for (i in 1:nrow(obs)) {
  
  tmp_ind <- which(ppl$id == obs$id[i])
  tmp_age <- obs$age[i]
  obs$mother_dead[i] <- mother_dead[tmp_ind, tmp_age]
  
}

mother_dead_check <- obs[obs$mother_dead == 1, ]

out <- ggplot() +
  
  geom_point(data = obs, 
             aes(x = age, 
                 y = bmi, 
                 col = repstat_binary), 
             alpha = 0.4, 
             size = 2) + 
  
  theme_bw() +
  
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  
  xlim(0, 19)  + 
  
  scale_color_manual(values = c("not pregnant or lactating" = "orange", 
                                "pregnant" = "purple", 
                                "lactating" = "darkgreen"), 
                     name = NULL) +
  
  geom_point(data = mother_dead_check, 
             aes(x = age, 
                 y = bmi),
             col = "red",
             shape = 1,
             size = 3.5)

png("output/figures/pregnancy_bmi.png", 
    res = 250, 
    height = 1100,
    width = 1900)

print(out)

dev.off()
