library(cowplot)
library(xtable)
library(ggplot2)
library(rethinking)

# plot full results for child height models

fit <- readRDS("stanfits/father_height.rds")

# check fit

png("output/trace/father_height.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(traceplot(fit, pars = c("alpha",  
                              "a_bo_tau", 
                              "a_bo_kappa",
                              "a_bo_delta",
                              "a_year_tau", 
                              "a_year_kappa",
                              "a_year_delta",
                              "a_age_tau",
                              "a_age_delta",
                              "a_age_kappa", 
                              "mother_sigma",
                              "father_sigma",
                              "sum_parent_sigma")))
dev.off()

# print summary table

tab <- precis(fit, 3, pars = c("alpha",  
                               "a_bo_tau", 
                               "a_bo_kappa",
                               "a_bo_delta",
                               "a_year_tau", 
                               "a_year_kappa",
                               "a_year_delta",
                               "a_age_tau",
                               "a_age_delta",
                               "a_age_kappa",
                               "mother_sigma",
                               "father_sigma",
                               "sum_parent_sigma"))

rownames(tab) <- c("$\\alpha$",
                   "$\\gamma_{\\tau}$", 
                   "$\\gamma_{\\kappa}$",
                   "$\\gamma_{\\delta}$",
                   "$\\epsilon_{\\tau}$",
                   "$\\epsilon_{\\kappa}$",
                   "$\\epsilon_{\\delta}$",
                   "$\\beta_{\\tau_1}$",
                   "$\\beta_{\\tau_2}$",
                   "$\\beta_{\\tau_3}$",
                   "$\\beta_{\\tau_4}$",
                   "$\\beta_{\\tau_5}$",
                   "$\\beta_{\\tau_6}$",
                   "$\\beta_{\\tau_7}$",
                   "$\\beta_{\\tau_8}$",
                   "$\\beta_{\\tau_9}$",
                   "$\\beta_{\\tau_{10}}$",
                   "$\\beta_{\\kappa_1}$",
                   "$\\beta_{\\kappa_2}$",
                   "$\\beta_{\\kappa_3}$",
                   "$\\beta_{\\kappa_4}$",
                   "$\\beta_{\\kappa_5}$",
                   "$\\beta_{\\kappa_6}$",
                   "$\\beta_{\\kappa_7}$",
                   "$\\beta_{\\kappa_8}$",
                   "$\\beta_{\\kappa_9}$",
                   "$\\beta_{\\kappa_{10}}$",
                   "$\\beta_{\\delta_1}$",
                   "$\\beta_{\\delta_2}$",
                   "$\\beta_{\\delta_3}$",
                   "$\\beta_{\\delta_4}$",
                   "$\\beta_{\\delta_5}$",
                   "$\\beta_{\\delta_6}$",
                   "$\\beta_{\\delta_7}$",
                   "$\\beta_{\\delta_8}$",
                   "$\\beta_{\\delta_9}$",
                   "$\\beta_{\\delta_{10}}$",
                   "$\\kappa_{\\sigma}$", 
                   "$\\eta_{\\sigma}$", 
                   "$\\pi_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_father_height.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# read processed data objects

load("processed_data/ppl.robj")
load("processed_data/obs.robj")
load("processed_data/father_dead.robj")
load("processed_data/mother_dead.robj")
load("processed_data/mother_unmarried.robj")
load("processed_data/mother_married_to_notfather.robj")
load("processed_data/mother_married_to_father_with_cowife.robj")

# prep data for plotting

n <- nrow(ppl)

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

data <- list(N = n, 
             A = 19, 
             B = 16, 
             Y = 86, 
             alive = alive, 
             birthorder = birthorder, 
             dob = dob, 
             male = male, 
             twin = twin, 
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

# turn -99 to NA for mean
data$height[data$height == -99] <- NA

data$male[data$male == -99] <- NA

boys <- which(data$male == 1)
girls <- which(data$male == 0)

height_boys <- data$height[boys, ]
height_girls <- data$height[girls, ]
  
height_boys <- reshape2::melt(height_boys, na.rm = TRUE)
height_boys <- height_boys[, 2:3]
colnames(height_boys) <- c("age", "height")
height_boys$sex <- "boys"

height_girls <- reshape2::melt(height_girls, na.rm = TRUE)
height_girls <- height_girls[, 2:3]
colnames(height_girls) <- c("age", "height")
height_girls$sex <- "girls"

raw_height <- rbind(height_girls, height_boys)

# extract samples

post <- extract.samples(fit)

# look at father effects

png("output/figures/height_father_random_effects.png", 
    res = 250,
    height = 1000,
    width = 2800)

par(mfrow = c(1, 3))

plot(apply(post$a_father, 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(post$a_mother, 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(post$father_sigma + post$mother_sigma, 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot predictions for base case

plot_data <- list()

sex <- c("girls", "boys")

for (s in 1:2) {
  
  p <- post$m_base[ , s, ]
  
  plot_data[[s]] <- data.frame(age = 1:19, 
                               sex = sex[s],
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, prob = 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, prob = 0.9))[2, ],
                               type = "biological parents alive and monogamously married") 
  
}

plot_data <- do.call(rbind, plot_data)

# load who growth curves

who_height <- read.csv("processed_data/who_height.csv")

a <-
  
ggplot() +
  
  labs(y = "height (cm)", 
       x = "child age") +
  
  theme_linedraw() +
  
  geom_line(data = plot_data, 
            aes(x = age, 
                y = offset,
                color = sex),
            size = 1.5) +
  
  geom_ribbon(data = plot_data,
              aes(ymin = low,
                  ymax = upp,
                  x = age,
                  fill = sex), 
              alpha = 0.2, 
              linetype = 2) +
  
  facet_grid(. ~ type, 
             scales = "free") +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(1, "cm"), 
        legend.position = c(0.8, 0.2),
        legend.background = element_rect(linetype = "solid", 
                                         color = "black"),
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic")) +
  
  scale_color_manual(values = c("girls" = "goldenrod", 
                                "boys" = "navy")) +
  
  scale_fill_manual(values = c("girls" = "goldenrod", 
                               "boys" = "navy")) + 
  
  ggtitle("A. Predicted child height (cm), across childhood") +
  
  geom_line(data = who_height, 
            aes(x = age, 
                y = mean_height, 
                col = sex), 
            lty = 2, 
            lwd = 1) +
  
  geom_point(data = raw_height, 
             aes(x = age, 
                 y = height, 
                 col = sex), 
             alpha = 0.15, 
             size = 1.5, 
             position = "jitter")

# plot age-specific parameters

par_names <- c("intercept", 
               "male", 
               "twin")

plot_data <- list()

for (i in 1:3) {
 
  p <- post$a_age[, i, ]
  
  plot_data[[i]] <- data.frame(age = 1:19, 
                               cat = par_names[i], 
                               mean = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ],
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ])
  
}

plot_data <- do.call(rbind, plot_data)

b <- 
  
ggplot(plot_data, 
       aes(x = age, 
           y = mean)) +
  
  theme_linedraw() +
  
  facet_wrap(. ~ cat, scales = "free") +

  geom_pointrange(aes(ymin = low, ymax = upp)) +
  
  ylab("estimate") +
  
  xlab("child age") +
  
  geom_hline(yintercept = 0, col = col.alpha("indianred", 0.6), size = 1) +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic")) + 
  
  ggtitle("B. Age-specific effects")

# plot birth-order parameters

p <- post$a_bo

plot_data <- data.frame(bo = 1:15, 
                        cat = "birth-order", 
                        mean = apply(p, 2, mean)[1:15], 
                        upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, 1:15],
                        low = apply(p, 2, function(x) HPDI(x, 0.9))[2, 1:15])

c <- 
  
  ggplot(plot_data, 
         aes(x = bo, 
             y = mean)) +
  
  theme_linedraw() +
  
  facet_wrap(. ~ cat, scales = "free") +
  
  geom_pointrange(aes(ymin = low, ymax = upp)) +
  
  ylab("estimate") +
  
  xlab("birth-order") +
  
  geom_hline(yintercept = 0, col = col.alpha("indianred", 0.6), size = 1) +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic")) + 
  
  ggtitle("C. Birth-order effects")

# plot birth-year parameters

p <- post$a_year

plot_data <- data.frame(year = 1930:2015, 
                        cat = "birth-year", 
                        mean = apply(p, 2, mean), 
                        upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ],
                        low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ])

d <- 
  
  ggplot(plot_data, 
         aes(x = year, 
             y = mean)) +
  
  theme_linedraw() +
  
  facet_wrap(. ~ cat, scales = "free") +
  
  geom_pointrange(aes(ymin = low, ymax = upp)) +
  
  ylab("estimate") +
  
  xlab("birth-year") +
  
  geom_hline(yintercept = 0, col = col.alpha("indianred", 0.6), size = 1) +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic")) + 
  
  geom_rect(aes(xmin = 1930,
                xmax = 1975,
                ymin = -0.5,
                ymax = 0.8),
            fill = col.alpha("grey", 0.006)) +
  
  ggtitle("D. Birth-year effects")

# plot deviations from base-case on prediction scale

post_list <- list(post$m_father_dead,
                  post$m_father_unmarried,
                  post$m_father_married_to_notmother_monogamy,
                  post$m_father_married_to_notmother_polygyny,
                  post$m_father_married_to_mother_polygyny)

type <- c("father dead", 
          "father unmarried",
          "father married \nto step-mother (monogamy)",
          "father married \nto step-mother (polygyny)",
          "father married \nto mother (polygyny)")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - post$m_base
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:19, 
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ],
                               type = type[z]) 
  
}

plot_data <- do.call(rbind, plot_data)

plot_data$type <- factor(plot_data$type, levels = c("father dead",
                                                    "father unmarried",
                                                    "father married \nto step-mother (monogamy)", 
                                                    "father married \nto step-mother (polygyny)",
                                                    "father married \nto mother (polygyny)"))
e <- 
  
ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           ymin = low, 
           ymax = upp, 
           color = type)) +
  
  theme_linedraw() +
  
  geom_pointrange(size = 0.6) +
  
  facet_grid(. ~ type) +
  
  labs(y = "contrast (cm)", 
       x = "child age") +
  
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 10, angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        plot.margin = unit(c(5.5, 0, 5.5, 0), "points")) +
  
  geom_hline(yintercept = 0, size = 1.1, color = col.alpha("firebrick", 0.5)) +
  
  scale_color_manual(values = c("father dead" = "navy",
                                "father unmarried" = "goldenrod",
                                "father married \nto step-mother (monogamy)" = "cyan4",
                                "father married \nto step-mother (polygyny)" = "sienna4",
                                "father married \nto mother (polygyny)" = "purple4")) +
  
  ggtitle("E. Age-specific contrasts to children of different father-states")

# produce estimates reported in text

mean(post$m_father_unmarried[, 2, 7])
HPDI(post$m_father_unmarried[, 2, 7], 0.9)

mean(post$m_base[, 2, 7] - post$m_father_unmarried[, 2, 7])
HPDI(post$m_base[, 2, 7] - post$m_father_unmarried[, 2, 7], 0.9)

mean(post$m_base[, 2, 7])
HPDI(post$m_base[, 2, 7], 0.9)

# mother model

fit <- readRDS("stanfits/mother_height.rds")

# check fit

png("output/trace/mother_height.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(traceplot(fit, pars = c("alpha",  
                              "a_bo_tau", 
                              "a_bo_kappa",
                              "a_bo_delta",
                              "a_year_tau", 
                              "a_year_kappa",
                              "a_year_delta",
                              "a_age_tau",
                              "a_age_delta",
                              "a_age_kappa",
                              "father_sigma",
                              "mother_sigma",
                              "sum_parent_sigma")))
dev.off()

# print summary table

tab <- precis(fit, 3, pars = c("alpha",  
                               "a_bo_tau", 
                               "a_bo_kappa",
                               "a_bo_delta",
                               "a_year_tau", 
                               "a_year_kappa",
                               "a_year_delta",
                               "a_age_tau",
                               "a_age_delta",
                               "a_age_kappa",
                               "father_sigma",
                               "mother_sigma",
                               "sum_parent_sigma"))

rownames(tab) <- c("$\\alpha$",
                   "$\\gamma_{\\tau}$", 
                   "$\\gamma_{\\kappa}$",
                   "$\\gamma_{\\delta}$",
                   "$\\epsilon_{\\tau}$",
                   "$\\epsilon_{\\kappa}$",
                   "$\\epsilon_{\\delta}$",
                   "$\\beta_{\\tau_1}$",
                   "$\\beta_{\\tau_2}$",
                   "$\\beta_{\\tau_3}$",
                   "$\\beta_{\\tau_4}$",
                   "$\\beta_{\\tau_5}$",
                   "$\\beta_{\\tau_6}$",
                   "$\\beta_{\\tau_7}$",
                   "$\\beta_{\\tau_8}$",
                   "$\\beta_{\\tau_9}$",
                   "$\\beta_{\\kappa_1}$",
                   "$\\beta_{\\kappa_2}$",
                   "$\\beta_{\\kappa_3}$",
                   "$\\beta_{\\kappa_4}$",
                   "$\\beta_{\\kappa_5}$",
                   "$\\beta_{\\kappa_6}$",
                   "$\\beta_{\\kappa_7}$",
                   "$\\beta_{\\kappa_8}$",
                   "$\\beta_{\\kappa_9}$",
                   "$\\beta_{\\delta_1}$",
                   "$\\beta_{\\delta_2}$",
                   "$\\beta_{\\delta_3}$",
                   "$\\beta_{\\delta_4}$",
                   "$\\beta_{\\delta_5}$",
                   "$\\beta_{\\delta_6}$",
                   "$\\beta_{\\delta_7}$",
                   "$\\beta_{\\delta_8}$",
                   "$\\beta_{\\delta_9}$",
                   "$\\kappa_{\\sigma}$", 
                   "$\\eta_{\\sigma}$", 
                   "$\\pi_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_mother_height.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# extract samples

post <- extract.samples(fit)

# check mother effects

png("output/figures/height_mother_random_effects.png", 
    res = 250,
    height = 1000,
    width = 2800)

par(mfrow = c(1, 3))

plot(apply(post$a_father, 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(post$a_mother, 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(post$father_sigma + post$mother_sigma, 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot deviations from base-case on prediction scale

post_list <- list(post$m_mother_dead,
                  post$m_mother_unmarried,
                  post$m_mother_married_to_notfather,
                  post$m_mother_married_to_father_with_cowife,
                  post$m_unknown_parent)

type <- c("mother dead", 
          "mother unmarried",
          "mother married \nto step-father",
          "mother married \nto father (with co-wife)",
          "mother or father unknown")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - post$m_base
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:19, 
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ],
                               type = type[z]) 
  
}

plot_data <- do.call(rbind, plot_data)

plot_data$type <- factor(plot_data$type, levels = c("mother dead", 
                                                    "mother unmarried",
                                                    "mother married \nto step-father",
                                                    "mother married \nto father (with co-wife)",
                                                    "mother or father unknown"))
f <-

ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           ymin = low, 
           ymax = upp, 
           color = type)) +
  
  theme_linedraw() +
  
  geom_pointrange(size = 0.6) +
  
  facet_wrap(. ~ type, scales = "free", ncol = 5) +
  
  labs(y = "contrast (cm)", 
       x = "child age") +
  
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 10, angle = 0), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        plot.margin = unit(c(5.5, 0, 5.5, 0), "points")) +
  
  geom_hline(yintercept = 0, size = 1.1, color = col.alpha("firebrick", 0.5)) +
  
  scale_color_manual(values = c("mother dead" = "navy",
                                "mother unmarried" = "goldenrod",
                                "mother married \nto step-father" = "cyan4",
                                "mother married \nto father (with co-wife)" = "purple4",
                                "mother or father unknown" = "chocolate3")) +

  ggtitle("F. Age-specific contrasts to children of different mother-states")

# produce estimates reported in text

mean(post$m_base[, 2, 7] - post$m_mother_dead[, 2, 7])
HPDI(post$m_base[, 2, 7] - post$m_mother_dead[, 2, 7], 0.9)

mean(post$m_base[, 2, 7] - post$m_mother_unmarried[, 2, 7])
HPDI(post$m_base[, 2, 7] - post$m_mother_unmarried[, 2, 7], 0.9)

# plot them all

tmp <- plot_grid(c,
                 d, 
                 rel_widths = c(1.5, 2.5))

tmp2 <- plot_grid(b,
                  tmp, 
                  nrow = 2)

tmp3 <- plot_grid(a, 
                  tmp2, 
                  rel_widths = c(1.5, 2))

tmp4 <- plot_grid(tmp3, 
                  e, 
                  f, 
                  nrow = 3, 
                  rel_heights = c(1.1, 0.6, 0.6))

png("output/figures/height.png",
    res = 250, 
    height = 2600, 
    width = 3000)

print(tmp4)

dev.off()
