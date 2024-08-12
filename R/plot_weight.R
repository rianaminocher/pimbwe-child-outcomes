# plot full results for weight-for-height models

# read fit

fit <- readRDS("stanfits/mother_weight.rds")

# extract samples

posterior_samples <- fit$draws()
mcmc_samples <- as_draws_df(posterior_samples)

# check fit
# plot trace of all main parameters

parameter_names <- names(mcmc_samples)
pattern <- "alpha|a_bo_tau|a_bo_kappa|a_bo_delta|a_year_tau|a_year_kappa|a_year_delta|a_age_tau\\[|a_age_delta\\[|a_age_kappa\\[|mother_sigma|father_sigma"
matching_parameters <- grep(pattern, parameter_names, value = TRUE)

# check fit

png("output/trace/mother_weight.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(bayesplot::mcmc_trace(mcmc_samples, pars = matching_parameters))

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
                   "$\\beta_{\\tau_10}$",                  
                   "$\\beta_{\\kappa_1}$",
                   "$\\beta_{\\kappa_2}$",
                   "$\\beta_{\\kappa_3}$",
                   "$\\beta_{\\kappa_4}$",
                   "$\\beta_{\\kappa_5}$",
                   "$\\beta_{\\kappa_6}$",
                   "$\\beta_{\\kappa_7}$",
                   "$\\beta_{\\kappa_8}$",
                   "$\\beta_{\\kappa_9}$",
                   "$\\beta_{\\kappa_10}$",                   
                   "$\\beta_{\\delta_1}$",
                   "$\\beta_{\\delta_2}$",
                   "$\\beta_{\\delta_3}$",
                   "$\\beta_{\\delta_4}$",
                   "$\\beta_{\\delta_5}$",
                   "$\\beta_{\\delta_6}$",
                   "$\\beta_{\\delta_7}$",
                   "$\\beta_{\\delta_8}$",
                   "$\\beta_{\\delta_9}$",
                   "$\\beta_{\\delta_10}$",
                   "$\\kappa_{\\sigma}$", 
                   "$\\eta_{\\sigma}$", 
                   "$\\pi_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_mother_weight.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# read processed data objects

load("processed_data/ppl.robj")
load("processed_data/obs.robj")
load("processed_data/father_dead.robj")
load("processed_data/father_unmarried.robj")
load("processed_data/mother_dead.robj")
load("processed_data/mother_unmarried.robj")
load("processed_data/mother_married_to_notfather.robj")
load("processed_data/mother_married_to_father_with_cowife.robj")

# prep data for plotting

n <- nrow(ppl)

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
             birthorder = birthorder, 
             dob = dob, 
             male = male, 
             twin = twin, 
             father_dead = father_dead, 
             father_unmarried = father_unmarried,
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

id <- ppl$id

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

data$weight <- weight

boys <- which(data$male == 1)
girls <- which(data$male == 0)

weight_boys <- data$weight[boys, ]
weight_girls <- data$weight[girls, ]

weight_boys <- reshape2::melt(weight_boys, na.rm = TRUE)
weight_boys <- weight_boys[, 2:3]
colnames(weight_boys) <- c("age", "bmi")
weight_boys$sex <- "boys"

weight_girls <- reshape2::melt(weight_girls, na.rm = TRUE)
weight_girls <- weight_girls[, 2:3]
colnames(weight_girls) <- c("age", "bmi")
weight_girls$sex <- "girls"

raw_weight <- rbind(weight_girls, weight_boys)

age <- matrix(NA, nrow = nrow(data$weight), ncol = ncol(data$weight))

for (i in 1:nrow(age)) {
  age[i, ] <- 1:19
}

mother_dead_check <- data.frame(weight = c(data$weight[which(data$mother_dead == 1)]),
                                age = c(age[which(data$mother_dead == 1)]))
mother_dead_check <- mother_dead_check[which(complete.cases(mother_dead_check)), ]

# extract samples

post <- as_draws_rvars(fit)

# plot parent random effects

png("output/figures/weight_mother_random_effects.png",
    res = 250,
    height = 1000,
    width = 2800)

par(mfrow = c(1, 3))

plot(apply(draws_of(post$a_father), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma) + draws_of(post$mother_sigma), 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot predictions for base case

plot_data <- list()

sex <- c("girls", "boys")

for (s in 1:2) {
  
  p <- draws_of(post$m_base)[ , s, ]
  
  plot_data[[s]] <- data.frame(age = 1:19, 
                               sex = sex[s],
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, prob = 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, prob = 0.9))[2, ],
                               type = "biological parents alive and monogamously married") 
  
}

plot_data <- do.call(rbind, plot_data)

# load who growth curves

who_weight <- read.csv("processed_data/who_weight.csv")

a <-
  
ggplot() +
  
  labs(y = expression(paste("weight ", (kg))), 
       x = "child age") +
  
  theme_linedraw() +
  
  geom_line(data = plot_data, 
            aes(x = age, 
                y = offset,
                color = sex),
            lwd = 1.5) +
  
  geom_ribbon(data = plot_data,
              aes(ymin = low, 
                  ymax = upp, 
                  x = age,
                  fill = sex), 
              alpha = 0.2, 
              linetype = 2) +
  
  facet_grid(. ~ type, 
             scales = "free") +
  
  ggtitle(expression(paste("A. Child weight ", (kg)))) +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(1, "cm"), 
        legend.position.inside = c(0.1, 0.85),
        legend.background = element_rect(linetype = "solid", 
                                         color = "black"),
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1),
        panel.grid.minor = element_line(colour = "grey70", linewidth = 0.05)) +
  
  scale_color_manual(values = c("girls" = "goldenrod", 
                                "boys" = "navy")) +
  
  scale_fill_manual(values = c("girls" = "goldenrod", 
                               "boys" = "navy")) + 
  
  geom_point(data = raw_weight, 
             aes(x = age, 
                 y = bmi, 
                 col = sex), 
             alpha = 0.15, 
             size = 1.5, 
             position = "jitter") +
  
  geom_point(data = mother_dead_check, 
             aes(x = age, 
                 y = weight),
             size = 2,
             col = col.alpha("red", 0.5))

# geom_line(data = who_weight, 
         # aes(x = age, 
         #     y = mean_weight, 
         #     col = sex), 
         #  lty = 2, 
         # lwd = 1) +
  
# plot age-specific parameters

par_names <- c("intercept", 
               "male", 
               "twin")

plot_data <- list()

for (i in 1:3) {
 
  p <- draws_of(post$a_age)[, i, ]
  
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
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) + 
  
  ggtitle("Age-specific effects")

# plot birth-order parameters

p <- draws_of(post$a_bo)

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
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) + 
  
  ggtitle("Birth-order effects")

# plot birth-year parameters

p <- draws_of(post$a_year)

plot_data <- data.frame(year = 1976:2014, 
                        cat = "year", 
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
  
  xlab("year of measurement") +
  
  geom_hline(yintercept = 0, col = col.alpha("indianred", 0.6), size = 1) +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) +
  
  ggtitle("Year-specific effects")

# plot deviations from base-case on prediction scale

post_list <- list(draws_of(post$m_unknown_parent),
                  draws_of(post$m_mother_dead),
                  draws_of(post$m_mother_unmarried),
                  draws_of(post$m_mother_married_to_notfather),
                  draws_of(post$m_mother_married_to_father_with_cowife))

type <- c("either parent external", 
          "mother deceased",
          "mother unmarried",
          "mother married to \nstep-father",
          "mother married to \nbio-father (with co-wife)")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - draws_of(post$m_base)
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:19, 
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ],
                               type = type[z]) 
  
}

plot_data <- do.call(rbind, plot_data)

plot_data$type <- factor(plot_data$type, levels = c("either parent external", 
                                                    "mother deceased",
                                                    "mother unmarried",
                                                    "mother married to \nstep-father",
                                                    "mother married to \nbio-father (with co-wife)"))

e <- 
  
ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           ymin = low, 
           ymax = upp, 
           color = type)) +
  
  theme_linedraw() +
  
  geom_pointrange(size = 0.6) +
  
  ylim(-8, 18) +
  
  facet_grid(. ~ type) +
  
  labs(y = expression(paste("contrast ", (kg))), 
       x = "child age") +
  
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 10, angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        plot.margin = unit(c(5.5, 0, 5.5, 0), "points"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) +
  
  geom_hline(yintercept = 0, size = 1.1, color = col.alpha("firebrick", 0.5)) +
  
  scale_color_manual(values = c("either parent external" = "chocolate3",
                                "mother deceased" = "navy",
                                "mother unmarried" = "goldenrod",
                                "mother married to \nstep-father" = "cyan4",
                                "mother married to \nbio-father (with co-wife)" = "purple4")) +
  
  ggtitle("B. Age-specific contrasts to children of different mother-states")

# in text predictions

mean(draws_of(post$m_mother_dead)[, 1, 16])
HPDI(draws_of(post$m_mother_dead)[, 1, 16], 0.9)

mean(draws_of(post$m_mother_dead)[, 1, 16] - draws_of(post$m_base)[, 1, 16])
HPDI(draws_of(post$m_mother_dead)[, 1, 16] - draws_of(post$m_base)[, 1, 16], 0.9)

# father model

fit <- readRDS("stanfits/father_weight.rds")

# extract samples

posterior_samples <- fit$draws()
mcmc_samples <- as_draws_df(posterior_samples)

# check fit
# plot trace of all main parameters

parameter_names <- names(mcmc_samples)
pattern <- "alpha|a_bo_tau|a_bo_kappa|a_bo_delta|a_year_tau|a_year_kappa|a_year_delta|a_age_tau\\[|a_age_delta\\[|a_age_kappa\\[|mother_sigma|father_sigma"
matching_parameters <- grep(pattern, parameter_names, value = TRUE)

# check fit

png("output/trace/father_weight.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(bayesplot::mcmc_trace(mcmc_samples, pars = matching_parameters))

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
                   "$\\beta_{\\tau_{10}}$",
                   "$\\beta_{\\tau_{11}}$",                
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
                   "$\\beta_{\\kappa_{11}}$",
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
                   "$\\beta_{\\delta_{11}}$",
                   "$\\kappa_{\\sigma}$", 
                   "$\\eta_{\\sigma}$", 
                   "$\\pi_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_father_weight.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# extract samples

post <- as_draws_rvars(fit)

# look at mother effects

png("output/figures/weight_father_random_effects.png",
    res = 250,
    height = 1000,
    width = 2800)

par(mfrow = c(1, 3))

plot(apply(draws_of(post$a_father), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma) + draws_of(post$mother_sigma), 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot deviations from base-case on prediction scale

post_list <- list(draws_of(post$m_father_dead),
                  draws_of(post$m_father_unmarried),
                  draws_of(post$m_father_married_to_notmother_monogamy),
                  draws_of(post$m_father_married_to_notmother_polygyny),
                  draws_of(post$m_father_married_to_mother_polygyny))

type <- c("father deceased", 
          "father unmarried",
          "father married to one \nstep-mother",
          "father polygynously married \n(w/o bio-mother)",
          "father polygynously married \n(with bio-mother)")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - draws_of(post$m_base)
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:19, 
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ],
                               type = type[z]) 
  
}

plot_data <- do.call(rbind, plot_data)

plot_data$type <- factor(plot_data$type, levels = c("father deceased", 
                                                    "father unmarried",
                                                    "father married to one \nstep-mother",
                                                    "father polygynously married \n(w/o bio-mother)",
                                                    "father polygynously married \n(with bio-mother)"))

f <-

ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           ymin = low, 
           ymax = upp, 
           color = type)) +
  
  theme_linedraw() +
  
  ylim(-8, 18) +
  
  geom_pointrange(size = 0.6) +
  
  facet_wrap(. ~ type, ncol = 5) +
  
  labs(y = expression(paste("contrast ", (kg))), 
       x = "child age") +
  
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 10, angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 10), 
        legend.position = "none", 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"), 
        plot.margin = unit(c(5.5, 0, 5.5, 0), "points"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) +
  
  geom_hline(yintercept = 0, size = 1.1, color = col.alpha("firebrick", 0.5)) +
  
  scale_color_manual(values = c("father deceased" = "navy",
                                "father unmarried" = "goldenrod",
                                "father married to one \nstep-mother" = "cyan4",
                                "father polygynously married \n(w/o bio-mother)" = "sienna4",
                                "father polygynously married \n(with bio-mother)" = "purple4")) +

  ggtitle("C. Age-specific contrasts to children of different father-states")

# produce estimates reported in text

mean(draws_of(post$m_father_dead)[, 1, 19] - draws_of(post$m_base)[, 1, 19])
HPDI(draws_of(post$m_father_dead)[, 1, 19] - draws_of(post$m_base)[, 1, 19], 0.9)

mean(draws_of(post$m_father_married_to_notmother_monogamy)[, 1, 19] - draws_of(post$m_base)[, 1, 19])
HPDI(draws_of(post$m_father_married_to_notmother_monogamy)[, 1, 19] - draws_of(post$m_base)[, 1, 19], 0.9)

mean(draws_of(post$m_base)[, 1, 19] - draws_of(post$m_father_married_to_mother_polygyny)[, 1, 19])
HPDI(draws_of(post$m_base)[, 1, 19] - draws_of(post$m_father_married_to_mother_polygyny)[, 1, 19], 0.9)

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

png("output/figures/weight.png",
    res = 250, 
    height = 2600, 
    width = 3000)

print(tmp4)

dev.off()

tmp <- plot_grid(e,
                 f, 
                 nrow = 2)

tmp2 <- plot_grid(a, 
                  tmp,
                  rel_widths = c(1.7, 3, 3))

pdf("output/figures/weight.pdf",
    height = 5.5, 
    width = 16)

print(tmp2)

dev.off()

pdf("output/figures/weight_age_effects.pdf",
    height = 3.5, 
    width = 10)

print(b)

dev.off()

tmp <- plot_grid(c,
                 d, 
                 rel_widths = c(1.5, 2.5))

pdf("output/figures/weight_birth_effects.pdf",
    height = 3.5, 
    width = 10.5)

print(tmp)

dev.off()

