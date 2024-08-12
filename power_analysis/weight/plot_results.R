# plot results of power analysis for weight outcomes

library(cmdstanr)
library(rethinking)
library(ggplot2)
library(parallel)
library(reshape2)

# load data
load("output/full_res.Rda")

create_plot_with_legend <- function(data, i) {
  
  p <- ggplot(data[data$Variable == i, ]) +
    
    geom_point(aes(x = Age, y = true_values, color = "Simulated Values"), 
               shape = 20, 
               size = 2) +
    
    geom_pointrange(aes(x = Age, y = median_values, ymin = L_values, ymax = H_values, color = "Estimated Values"), 
                    shape = 1, 
                    size = 0.5) +
    
    ylab("Effect Size (log-odds)") +
    
    xlab("Age Category") +
    
    theme_linedraw() +
    
    facet_wrap(. ~ MES, scales = "free") +
    
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    
    scale_color_manual(name = "", 
                       values = c("Simulated Values" = "purple", "Estimated Values" = "orange")) +
    
    theme(legend.position = "right")
  
  return(p)
  
}

for (i in 2:10) {
  
  png(paste0("output/results_", i, ".png"), 
      height = 2000, 
      width = 3000, 
      res = 250)
  
  print(create_plot_with_legend(full_res, i))
  
  dev.off()
  
}
