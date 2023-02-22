# process WHO data so we have a single df throughout childhood

library(openxlsx)

# height-for-age, boys

d1 <- read.xlsx("raw_data_WHO/lhfa-boys-percentiles-expanded-tables.xlsx")
d2 <- read.xlsx("raw_data_WHO/hfa-boys-perc-who2007-exp.xlsx")

# convert day to year

b1 <- data.frame(age = 1:6) # 1 represents first year of life

d1$age <- floor(d1$Day/365.25) + 1

b1$mean_height <- tapply(d1$M, d1$age, mean)
b1$sex <- "boys"

# convert month to year

d2$age <- floor(d2$Month/12) + 1

b2 <- data.frame(age = 6:19)
b2$mean_height <- tapply(d2$M, d2$age, mean)[1:14]
b2$sex <- "boys"

who_height <- rbind(b1[1:5, ], b2)

# height-for-age, girls

d1 <- read.xlsx("raw_data_WHO/lhfa-girls-percentiles-expanded-tables.xlsx")
d2 <- read.xlsx("raw_data_WHO/hfa-girls-perc-who2007-exp.xlsx")

# convert day to year

g1 <- data.frame(age = 1:6) # 1 represents first year of life

d1$age <- floor(d1$Day/365.25) + 1

g1$mean_height <- tapply(d1$M, d1$age, mean)
g1$sex <- "girls"

# convert month to year

d2$age <- floor(d2$Month/12) + 1

g2 <- data.frame(age = 6:19)
g2$mean_height <- tapply(d2$M, d2$age, mean)[1:14]
g2$sex <- "girls"

who_height <- rbind(who_height, g1[1:5, ], g2)

write.csv(who_height, 
          "processed_data/who_height.csv",
          row.names = FALSE)

# BMI boys

d1 <- read.xlsx("raw_data_WHO/bfa-boys-percentiles-expanded-tables.xlsx")
d2 <- read.xlsx("raw_data_WHO/bmi-boys-perc-who2007-exp.xlsx")

# convert day to year

b1 <- data.frame(age = 1:6) # 1 represents first year of life

d1$age <- floor(d1$Age/365.25) + 1

b1$mean_weight <- tapply(d1$M, d1$age, mean)
b1$sex <- "boys"

# convert month to year

d2$age <- floor(d2$Month/12) + 1

b2 <- data.frame(age = 6:20)
b2$mean_weight <- tapply(d2$M, d2$age, mean)
b2$sex <- "boys"

who_weight <- rbind(b1[1:5, ], b2)

# BMI girls

d1 <- read.xlsx("raw_data_WHO/bfa-girls-percentiles-expanded-tables.xlsx")
d2 <- read.xlsx("raw_data_WHO/bmi-girls-perc-who2007-exp.xlsx")

# convert day to year

g1 <- data.frame(age = 1:6) # 1 represents first year of life

d1$age <- floor(d1$Age/365.25) + 1

g1$mean_weight <- tapply(d1$M, d1$age, mean)
g1$sex <- "girls"

# convert month to year

d2$age <- floor(d2$Month/12) + 1

g2 <- data.frame(age = 6:20)
g2$mean_weight <- tapply(d2$M, d2$age, mean)
g2$sex <- "girls"

who_weight <- rbind(who_weight, g1[1:5, ], g2)

write.csv(who_weight, 
          "processed_data/who_weight.csv",
          row.names = FALSE)
