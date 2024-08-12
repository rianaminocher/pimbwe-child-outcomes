# plot full results for education models

# read fit

fit <- readRDS("stanfits/mother_education.rds")

# extract samples

posterior_samples <- fit$draws()
mcmc_samples <- as_draws_df(posterior_samples)

# check fit
# plot trace of all main parameters

parameter_names <- names(mcmc_samples)
pattern <- "alpha|a_bo_tau|a_bo_kappa|a_bo_delta|a_year_tau|a_year_kappa|a_year_delta|a_age_tau\\[|a_age_delta\\[|a_age_kappa\\[|mother_sigma|father_sigma"
matching_parameters <- grep(pattern, parameter_names, value = TRUE)

# check fit

png("output/trace/mother_education.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(bayesplot::mcmc_trace(mcmc_samples, pars = matching_parameters))

dev.off()

# write summary tab

tab <- precis(fit, 3, pars = c("alpha_l",
                               "a_bo_tau_l", 
                               "a_bo_kappa_l",
                               "a_bo_delta_l",
                               "a_year_tau_l", 
                               "a_year_kappa_l",
                               "a_year_delta_l",
                               "a_age_tau_l",
                               "a_age_delta_l",
                               "a_age_kappa_l",
                               "father_sigma_l",
                               "mother_sigma_l"))

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
                   "$\\eta_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_mother_education_1.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

tab <- precis(fit, 3, pars = c("alpha_t",
                               "a_bo_tau_t", 
                               "a_bo_kappa_t",
                               "a_bo_delta_t",
                               "a_year_tau_t", 
                               "a_year_kappa_t",
                               "a_year_delta_t",
                               "a_age_tau_t",
                               "a_age_delta_t",
                               "a_age_kappa_t",
                               "father_sigma_t",
                               "mother_sigma_t"))

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
                   "$\\eta_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_mother_education_2.txt", 
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

n <- nrow(ppl)
id <- ppl$id

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
             mother_dead = mother_dead, 
             mother_unmarried = mother_unmarried, 
             mother_married_to_notfather = mother_married_to_notfather,
             mother_married_to_father_with_cowife = mother_married_to_father_with_cowife)

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

data$edu <- edu

data$male[data$male == -99] <- NA

boys <- which(data$male == 1)
girls <- which(data$male == 0)

edu_boys <- data$edu[boys, ]
edu_girls <- data$edu[girls, ]

edu_boys <- reshape2::melt(edu_boys, na.rm = TRUE)
edu_boys <- edu_boys[, 2:3]
colnames(edu_boys) <- c("age", "edu")
edu_boys$sex <- "boys"

edu_girls <- reshape2::melt(edu_girls, na.rm = TRUE)
edu_girls <- edu_girls[, 2:3]
colnames(edu_girls) <- c("age", "edu")
edu_girls$sex <- "girls"

raw_edu <- rbind(edu_girls, edu_boys)

# extract samples

post <- as_draws_rvars(fit)

# plot parent random effects

png("output/figures/education_mother_random_effects.png", 
    res = 250,
    height = 2000,
    width = 2800)

par(mfrow = c(2, 3))

plot(apply(draws_of(post$a_father_t), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother_t), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma_t) + draws_of(post$mother_sigma_t), 
     xlab = "sum of variance on mother/father random effects")

plot(apply(draws_of(post$a_father_l), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother_l), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma_l) + draws_of(post$mother_sigma_l), 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot predictions for base case

plot_data <- list()

sex <- c("girls", "boys")

for (s in 1:2) {
  
  p <- draws_of(post$m_base_l)[ , s, ]
  
  plot_data[[s]] <- data.frame(age = 4+(1:15), 
                               sex = sex[s],
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, prob = 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, prob = 0.9))[2, ],
                               type = "biological parents alive and monogamously married") 
  
}

plot_data <- do.call(rbind, plot_data)

a <-
  
  ggplot() +
  
  labs(y = "education (years)", 
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
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(1, "cm"), 
        legend.position.inside = c(0.2, 0.8),
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
  
  ggtitle("A. Childhood education (years)") +
  
  geom_point(data = raw_edu, 
             aes(x = age,
                 y = edu,
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
 
  p <- draws_of(post$a_age_l)[, i, ]
  
  plot_data[[i]] <- data.frame(age = 1:15, 
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
  
  facet_wrap(. ~ cat, 
             scales = "free") +

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
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1),
        panel.grid.minor = element_line(colour = "grey70", linewidth = 0.05)) + 
  
  ggtitle("Age-specific effects")

# plot birth-order parameters

p <- draws_of(post$a_bo_l)

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

p <- draws_of(post$a_year_l)

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

post_list <- list(draws_of(post$m_unknown_parent_l),
                  draws_of(post$m_mother_dead_l),
                  draws_of(post$m_mother_unmarried_l),
                  draws_of(post$m_mother_married_to_notfather_l),
                  draws_of(post$m_mother_married_to_father_with_cowife_l))

type <- c("either parent external", 
          "mother deceased",
          "mother unmarried",
          "mother married to \nstep-father",
          "mother married to  \nbio-father (with co-wife)")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - draws_of(post$m_base_l)
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:15, 
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
                                                    "mother married to  \nbio-father (with co-wife)"))

e <- 
  
ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           ymin = low, 
           ymax = upp, 
           color = type)) +
  
  theme_linedraw() +
  
  ylim(-2, 4) +
  
  geom_pointrange(size = 0.6) +
  
  facet_grid(. ~ type) +
  
  labs(y = "contrast (years)", 
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
                                "mother married to  \nbio-father (with co-wife)" = "purple4")) +
  
  ggtitle("B. Age-specific contrasts to children of different mother-states")

# produce estimates reported in text

mean(draws_of(post$m_mother_married_to_notfather_l)[, 1, 15])
HPDI(draws_of(post$m_mother_married_to_notfather_l)[, 1, 15], 0.9)

mean(draws_of(post$m_base_l)[, 1, 15])
mean(draws_of(post$m_mother_married_to_notfather_l)[, 1, 15])
mean(draws_of(post$m_mother_unmarried_l)[, 1, 15])
mean(draws_of(post$m_mother_dead_l)[, 1, 15])
mean(draws_of(post$m_mother_married_to_father_with_cowife_l)[, 1, 15])

# father model

fit <- readRDS("stanfits/father_education.rds")

# extract samples

posterior_samples <- fit$draws()
mcmc_samples <- as_draws_df(posterior_samples)

# check fit
# plot trace of all main parameters

parameter_names <- names(mcmc_samples)
pattern <- "alpha|a_bo_tau|a_bo_kappa|a_bo_delta|a_year_tau|a_year_kappa|a_year_delta|a_age_tau\\[|a_age_delta\\[|a_age_kappa\\[|mother_sigma|father_sigma"
matching_parameters <- grep(pattern, parameter_names, value = TRUE)

# check fit

png("output/trace/father_education.png", 
    res = 250, 
    height = 3000, 
    width = 3000)

print(bayesplot::mcmc_trace(mcmc_samples, pars = matching_parameters))

dev.off()

# print summary table

tab <- precis(fit, 3, pars = c("alpha_l",
                               "a_bo_tau_l", 
                               "a_bo_kappa_l",
                               "a_bo_delta_l",
                               "a_year_tau_l", 
                               "a_year_kappa_l",
                               "a_year_delta_l",
                               "a_age_tau_l",
                               "a_age_delta_l",
                               "a_age_kappa_l",
                               "father_sigma_l",
                               "mother_sigma_l"))

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
                   "$\\eta_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_father_education_1.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

tab <- precis(fit, 3, pars = c("alpha_t",
                               "a_bo_tau_t", 
                               "a_bo_kappa_t",
                               "a_bo_delta_t",
                               "a_year_tau_t", 
                               "a_year_kappa_t",
                               "a_year_delta_t",
                               "a_age_tau_t",
                               "a_age_delta_t",
                               "a_age_kappa_t",
                               "father_sigma_t",
                               "mother_sigma_t"))

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
                   "$\\eta_{\\sigma}$")

print(xtable(tab), 
      file = "output/tables/summary_father_education_2.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# extract samples

post <- as_draws_rvars(fit)

# check father sigma

png("output/figures/education_father_random_effects.png", 
    res = 250,
    height = 2000,
    width = 2800)

par(mfrow = c(2, 3))

plot(apply(draws_of(post$a_father_l), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother_l), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma_l) + draws_of(post$mother_sigma_l), 
     xlab = "sum of variance on mother/father random effects")

plot(apply(draws_of(post$a_father_t), 2, mean), 
     ylab = "father random effects",
     xlab = "")

plot(apply(draws_of(post$a_mother_t), 2, mean), 
     ylab = "mother random effects",
     xlab = "")

dens(draws_of(post$father_sigma_t) + draws_of(post$mother_sigma_t), 
     xlab = "sum of variance on mother/father random effects")

dev.off()

# plot deviations from base-case on prediction scale

post_list <- list(draws_of(post$m_father_dead_l),
                  draws_of(post$m_father_unmarried_l),
                  draws_of(post$m_father_married_to_notmother_monogamy_l),
                  draws_of(post$m_father_married_to_notmother_polygyny_l),
                  draws_of(post$m_father_married_to_mother_polygyny_l))

type <- c("father deceased", 
          "father unmarried",
          "father married to one \nstep-mother",
          "father polygynously married \n(w/o bio-mother)",
          "father polygynously married \n(with bio-mother)")

plot_data <- list()

for (z in 1:5) {
  
  p <- post_list[[z]] - draws_of(post$m_base_l)
  # plot for boys
  p <- p[ , 2, ]
  
  plot_data[[z]] <- data.frame(age = 1:15, 
                               offset = apply(p, 2, mean), 
                               upp = apply(p, 2, function(x) HPDI(x, 0.9))[1, ], 
                               low = apply(p, 2, function(x) HPDI(x, 0.9))[2, ],
                               type = type[z]) 
  
}

plot_data <- do.call(rbind, plot_data)

plot_data$age <- plot_data$age + 4

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
  
  ylim(-2, 4) +
  
  geom_pointrange(size = 0.6) +
  
  facet_grid(. ~ type) +
  
  labs(y = "contrast (years)", 
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

png("output/figures/education.png",
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

pdf("output/figures/education.pdf",
    height = 5.5, 
    width = 16)

print(tmp2)

dev.off()

pdf("output/figures/edu_age_effects.pdf",
    height = 3.5, 
    width = 10)

print(b)

dev.off()

tmp <- plot_grid(c,
                 d, 
                 rel_widths = c(1.5, 2.5))

pdf("output/figures/edu_birth_effects.pdf",
    height = 3.5, 
    width = 10.5)

print(tmp)

dev.off()
