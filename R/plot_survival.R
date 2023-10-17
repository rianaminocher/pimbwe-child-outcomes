# plot full results for child survival models

# mother model

# read model fit

fit <- readRDS("stanfits/mother_survival.rds")

# check fit
# plot trace of all main parameters

png("output/trace/mother_survival.png", 
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

# print summary table for all main parameters

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
      file = "output/tables/summary_mother_survival.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# extract samples

post <- extract.samples(fit)

# plot parent random effects

png("output/figures/survival_mother_random_effects.png", 
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

# plot predictions for base case (biological 2-parent)

plot_data <- list()

# loop over sex
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

a <-
  
ggplot(plot_data, 
       aes(x = age, 
           y = offset,
           color = sex)) +
  
  ylim(c(0.5, 1)) +
  
  labs(y = "cumulative probability of survival", 
       x = "child age") +
  
  theme_linedraw() +
  
  geom_line(size = 1.5) +
  
  geom_ribbon(aes(ymin = low, ymax = upp, fill = sex), 
              alpha = 0.2, 
              linetype = 0) +
  
  facet_grid(. ~ type, scales = "free") +
  
  theme(strip.text.x = element_text(size = 10, color = "white"), 
        strip.text.y = element_text(size = 10, color = "white", angle = 0), 
        axis.text = element_text(size = 10.5), 
        axis.title = element_text(size = 10.5),
        legend.key.size = unit(1, "cm"), 
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 13, face = "italic"),
        panel.grid.major = element_line(colour = "grey70", size = 0.1),
        panel.grid.minor = element_line(colour = "grey70", size = 0.05)) +
  
  scale_color_manual(values = c("girls" = "goldenrod", 
                                "boys" = "navy")) +
  
  scale_fill_manual(values = c("girls" = "goldenrod", 
                               "boys" = "navy")) + 
  
  ggtitle("A. Child survival")

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

# plot year parameters

p <- post$a_year

plot_data <- data.frame(year = 1930:2015, 
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
# i.e., diff between a family state and the bio 2-parent

post_list <- list(post$m_unknown_parent,
                  post$m_mother_dead,
                  post$m_mother_unmarried,
                  post$m_mother_married_to_notfather,
                  post$m_mother_married_to_father_with_cowife)

type <- c("either parent external", 
          "mother deceased",
          "mother unmarried",
          "mother married to \nstep-father",
          "mother married to \nbio-father (with co-wife)")

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
  
  facet_grid(. ~ type) +
  
  ylim(c(-0.5, 0.2)) +
  
  labs(y = "contrast", 
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
  
  ggtitle("B. Age-specific contrasts in survival for children experiencing different mother-states")

# produce estimates reported in-text

# male vs. female

mean(post$m_base[, 2, 2] - post$m_base[, 1, 2])
HPDI(post$m_base[, 2, 2] - post$m_base[, 1, 2], prob = 0.9)

# mother/father unknown

mean(post$m_base[, 1, 19] - post$m_unknown_parent[, 1, 19])
HPDI(post$m_base[, 1, 19] - post$m_unknown_parent[, 1, 19], prob = 0.9)

# mother dead

mean(post$m_mother_dead[, 1, 1])
HPDI(post$m_mother_dead[, 1, 1], 0.9)

mean(post$m_base[, 1, 1] - post$m_mother_dead[, 1, 1])
HPDI(post$m_base[, 1, 1] - post$m_mother_dead[, 1, 1], 0.9)

mean(post$m_base[, 1, 1])
HPDI(post$m_base[, 1, 1], 0.9)

mean(post$m_mother_dead[, 1, 19])
HPDI(post$m_mother_dead[, 1, 19], 0.9)

mean(post$m_base[, 1, 19] - post$m_mother_dead[, 1, 19])
HPDI(post$m_base[, 1, 19] - post$m_mother_dead[, 1, 19], 0.9)

mean(post$m_base[, 1, 19])
HPDI(post$m_base[, 1, 19], 0.9)

# father model

fit <- readRDS("stanfits/father_survival.rds")

# check fit

png(paste0("output/trace/father_survival.png"), 
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
      file = paste0("output/tables/summary_mother_survival.txt"), 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# extract samples

post <- extract.samples(fit)

# check mother effects

png("output/figures/survival_father_random_effects.png", 
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

post_list <- list(post$m_father_dead,
                  post$m_father_unmarried,
                  post$m_father_married_to_notmother_monogamy,
                  post$m_father_married_to_notmother_polygyny,
                  post$m_father_married_to_mother_polygyny)

type <- c("father deceased", 
          "father unmarried",
          "father married to one \nstep-mother",
          "father polygynously married \n(w/o bio-mother)",
          "father polygynously married \n(with bio-mother)")

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
  
  geom_pointrange(size = 0.6) +
  
  ylim(c(-0.5, 0.2)) +
  
  facet_grid(. ~ type) +
  
  labs(y = "contrast", 
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
                                "father married to one \nstep-mother"  = "cyan4",
                                "father polygynously married \n(w/o bio-mother)" = "sienna4",
                                "father polygynously married \n(with bio-mother)"  = "purple4")) +

  ggtitle("C. Age-specific contrasts in survival for children experiencing different father-states")

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

png("output/figures/survival_big.png",
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
                  rel_widths = c(1.5, 3, 3))

pdf("output/figures/survival.pdf",
    height = 5.5, 
    width = 14.6)

print(tmp2)

dev.off()

# plot separately for supp

pdf("output/figures/survival_age_effects.pdf",
    height = 3.5, 
    width = 10)

print(b)

dev.off()

tmp <- plot_grid(c,
                 d, 
                 rel_widths = c(1.5, 2.5))

pdf("output/figures/survival_birth_effects.pdf",
    height = 3.5, 
    width = 10.5)

print(tmp)

dev.off()

