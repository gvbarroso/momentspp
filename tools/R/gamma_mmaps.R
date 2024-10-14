


suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(scales)
})

setwd("~/Data/bgs_lmr/sim_data/mut_layers/")

emp_map_gamma <- fread("gamma_fitted_emp_mmap.csv")
emp_map_gamma <- emp_map_gamma[3:nrow(emp_map_gamma),]
models <- fread("models.csv")
num_models <- nrow(models)
num_reps <- 10

gamma_fitted_sim <- as.data.frame(matrix(nrow=nrow(models), ncol=3))
names(gamma_fitted_sim) <- c("estimate", "model", "bin_size")

for(i in 1:num_models) {
  
  print(Sys.time())
  cat(paste("Visualizing data from model", i, "...\n"))
  
    
  m1kb <- fread(paste("model_", i, "/rep_", 1, "/maps_1kb.csv.gz", sep=""))
  m10kb <- fread(paste("model_", i, "/rep_", 1, "/maps_10kb.csv.gz", sep=""))
  m100kb <- fread(paste("model_", i, "/rep_", 1, "/maps_100kb.csv.gz",sep=""))
  m1Mb <- fread(paste("model_", i, "/rep_", 1, "/maps_1Mb.csv.gz",sep=""))
    
  g1k <- MASS::fitdistr(m1kb$avg_mut / mean(m1kb$avg_mut), "gamma")
  g10k <- MASS::fitdistr(m10kb$avg_mut / mean(m10kb$avg_mut), "gamma")
  g100k <- MASS::fitdistr(m100kb$avg_mut / mean(m100kb$avg_mut), "gamma")
  g1M <- MASS::fitdistr(m1Mb$avg_mut / mean(m1Mb$avg_mut), "gamma")
  
  gamma_fitted_sim[4 * i + 1,] <- c(g1k$estimate[1], i, 1e+3)
  gamma_fitted_sim[4 * i + 2,] <- c(g10k$estimate[1], i, 1e+4)
  gamma_fitted_sim[4 * i + 3,] <- c(g100k$estimate[1], i, 1e+5)
  gamma_fitted_sim[4 * i + 4,] <- c(g1M$estimate[1], i, 1e+6)
}

tbl <- merge(gamma_fitted_sim, models)

p1 <- ggplot(data=filter(tbl),
             aes(x=bin_size/1e+3, y=estimate, shape=as.factor(mut_spans_small),
                 color=as.factor(shape_rate_u_small))) +
  geom_point(size=2) + theme_bw() + geom_line() +
  facet_grid(+shape_rate_u_large~mut_spans_large) +
  scale_shape_manual(values=c(0, 1, 2), name="Mut. Spans Small") +
  scale_color_discrete(name="Shape_rate Small", type=c("plum3", "seagreen3", "brown1")) +
  scale_x_continuous(breaks=unique(tbl$bin_size/1e+3), trans="log10") +
  scale_y_continuous(breaks=c(2, 10, 100, 1000, 1e+5), trans="log10") +
  labs(title="Gamma Estimates, cols->span muts LARGE, rows->Gamma estimate LARGE",
       x="Map Scale (kb)", y="Shape-Rate Estimate") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom") +
  geom_point(data=emp_map_gamma, col='black', shape=8, size=3)
  
save_plot("Gamma_estimates.png", p1, base_height=8, base_width=12)
