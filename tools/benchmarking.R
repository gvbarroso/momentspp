library(matlib)
library(tidyverse)
library(scales)

setwd("~/Devel/momentspp/benchmarking/")

##############################
#
# timing
#
##############################

pop_scenarios <- 2^(0:3)
MultTimeTable <- as.data.frame(matrix(nrow = length(pop_scenarios), ncol = 5))
names(MultTimeTable) <- c("num_pops", "num_stats", "direct", "eigen", "pseudo-eigen")

for(i in 1:length(pop_scenarios)) {
  stats <- read.table(paste("timing/pops_", pop_scenarios[i], "/final_moments.txt", sep = ""))
  dir <- read.table(paste("timing/pops_", pop_scenarios[i], "/timing_direct.txt", sep = ""))
  eigen <- read.table(paste("timing/pops_", pop_scenarios[i], "/timing_eigen.txt", sep = ""))
  pseudo <- read.table(paste("timing/pops_", pop_scenarios[i], "/timing_pseudo-eigen.txt", sep = ""))
  
  MultTimeTable[i,] <- c(pop_scenarios[i], nrow(stats), dir, eigen, pseudo)
}

molten_time <- pivot_longer(MultTimeTable, cols = c("direct", "eigen", "pseudo-eigen"), names_to = "variable")

bench_plot <- ggplot(data = molten_time, aes(x = num_pops, y = value, shape = variable))
bench_plot <- bench_plot + geom_line(data = molten_time)
bench_plot <- bench_plot + geom_point(aes(shape = variable), size = 4)
bench_plot <- bench_plot + scale_x_continuous(breaks = pop_scenarios, trans="log2") 
bench_plot <- bench_plot + scale_y_continuous(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4), trans="log10")
bench_plot <- bench_plot + scale_shape_manual(values = c(0, 1, 2))
bench_plot <- bench_plot + labs(title = "Benchmark (apply transformation to 10000 generations)", x = "Num Pops.", y = "Time (seconds)") + theme_bw()
bench_plot <- bench_plot + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.position = "bottom")
bench_plot

##############################
#
# 1-pop accuracy
#
##############################

p1_stats <- c("DD", "Dz", "H", "I", "Pi2")
num_stats <- length(p1_stats)

r <- c(1e-8, 1e-7, 1e-6)
mu <- c(1e-8, 1e-7, 1e-6)
Ne <- c(1e+3, 1e+4, 1e+5)

a <- c("low_Ne", "mid_Ne", "high_Ne")
b <- c("low_mu", "mid_mu", "high_mu")
c <- c("low_r", "mid_r", "high_r")

params <- crossing(Ne, mu, r)
num_scenarios <- nrow(params)
  
steadyStateTable <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c), ncol = num_stats))
finalMomentsTable <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c), ncol = num_stats))

names(steadyStateTable) <- p1_stats
names(finalMomentsTable) <- p1_stats

idx = 1
for(x in 1:length(a)) {
  for(y in 1:length(b)) {
    for(z in 1:length(c)) {
      
      s <- t(read.table(paste(paste("accuracy/pops_1", a[x], b[y], c[z], sep = "/"), "/steady_state.txt", sep = ""))$V5)
      f <- t(read.table(paste(paste("accuracy/pops_1", a[x], b[y], c[z], sep = "/"), "/final_moments.txt", sep = ""))$V5)
      
      steadyStateTable[idx,] <- s
      finalMomentsTable[idx,] <- f
      
      idx = idx + 1
    }
  }
}

stats <- bind_cols(bind_rows(steadyStateTable, finalMomentsTable), bind_rows(params, params))
stats$type <- c(rep("steady", num_scenarios), rep("final", num_scenarios))

write.table(stats, "bench_1-pop_moments.csv", sep = ",", row.names = F, quote = F)

dat_stats <- pivot_longer(stats, cols = all_of(p1_stats), names_to = "variable")

p <- ggplot(data = dat_stats, aes(x = variable, y = value, shape = type, color = as.factor(Ne))) + facet_grid(r~mu)
p <- p + geom_point(size = 5) + theme_bw()
p <- p + scale_shape_manual(values = c(0, 1, 2, 3, 4))
p <- p + scale_y_continuous(trans = "log10")
p <- p + labs(title = "1-Epoch Moments x r (rows) and u (cols)", x = "Moment", y = "Value")
p <- p + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), legend.position = "bottom")
p
ggsave("moments_1-pop.pdf", p, device = "pdf", width = 12, height = 12)

##############################
#
# 2-pop accuracy
#
##############################

mu <- 1e-6
r <- c(0, 1e-5)
Ne_0 <- c(1e+3, 1e+4)
Ne_1 <- 1e+4
m_01 <- c(0, 1e-4)
m_10 <- 1e-4
num_epochs <- c(1, 2)

num_gen_epoch_2 <- 5000
Ne_0_bottleneck <- Ne_0 / 10
Ne_1_bottleneck <- Ne_1 / 10

params <- crossing(mu, r, num_epochs, m_01, m_10, Ne_0, Ne_1)
params <- mutate(params, Ne_0_bottleneck = case_when(num_epochs == 1 ~ NA_real_, num_epochs == 2 ~ (Ne_0 / 10)))
params <- mutate(params, Ne_1_bottleneck = case_when(num_epochs == 1 ~ NA_real_, num_epochs == 2 ~ (Ne_1 / 10)))
params <- mutate(params, num_gen_epoch_2 = case_when(num_epochs == 1 ~ NA_real_, num_epochs == 2 ~ 5000))

num_scenarios <- nrow(params)

pop_indices <- c(0, 1)
DD_stats <- NULL
Dz_stats <- NULL
H_stats <- NULL
pi2_stats <- NULL

for(i in 1:length(pop_indices)) {
  for(j in 1:length(pop_indices)) {
    
    DD_stats <- c(DD_stats, paste("DD", pop_indices[i], pop_indices[j], sep = "_"))
    H_stats <- c(H_stats, paste("H", pop_indices[i], pop_indices[j], sep = "_"))
    
    for(k in 1:length(pop_indices)) {
      
      Dz_stats <- c(Dz_stats, paste("Dz", pop_indices[i], pop_indices[j], pop_indices[k], sep = "_"))
      
      for(l in 1:length(pop_indices)) {
        
        pi2_stats <- c(pi2_stats, paste("pi2", pop_indices[i], pop_indices[j], pop_indices[k], pop_indices[l], sep = "_"))
      }
    }
  }
}

p2_stats <- c(DD_stats, Dz_stats, H_stats, "I", pi2_stats)
num_stats <- length(p2_stats)

a <- c("zero_rec", "nonzero_rec")
b <- c("single_epoch", "two_epochs")
c <- c("asymmetric_mig", "symmetric_mig")
d <- c("diff_Ne", "same_Ne")

steadyStateTable <- as.data.frame(matrix(nrow = num_scenarios, ncol = num_stats))
finalMomentsTable <- as.data.frame(matrix(nrow = num_scenarios, ncol = num_stats))

idx = 1
for(x in 1:length(a)) {
  for(y in 1:length(b)) {
    for(z in 1:length(c)) {
      for(w in 1:length(d)) {
        
        steady_state <- read.table(paste(paste("accuracy/pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/steady_state.txt", sep = ""))
        steadyStateTable[idx,] <- t(steady_state$V5)
        
        final_moments <- read.table(paste(paste("accuracy/pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/final_moments.txt", sep = ""))
        finalMomentsTable[idx,] <- t(final_moments$V5)
        
        idx = idx + 1
      }
    }
  }
}

names(steadyStateTable) <- p2_stats
names(finalMomentsTable) <- p2_stats

steadY <- bind_cols(steadyStateTable, params)
finY <- bind_cols(finalMomentsTable, params)

stats <- bind_rows(steadY, finY)
stats$type <- c(rep("steady", num_scenarios), rep("final", num_scenarios))

write.table(stats, "bench_2-pop_moments.csv", sep = ",", row.names = F, quote = F)

dat_stats <- pivot_longer(stats, cols = all_of(p2_stats), names_to = "variable")

# one epoch plot
p1 <- ggplot(data = filter(dat_stats, num_epochs == 1), aes(x = variable, y = value, shape = type, color = as.factor(r))) + facet_grid(Ne_0~m_01)
p1 <- p1 + geom_point(size = 5) + theme_bw()
p1 <- p1 + scale_shape_manual(values = c(0, 1))
p1 <- p1 + labs(title = "1-Epoch Moments x N_0 (rows) and m_01 (cols)", x = "Moment", y = "Value")
p1 <- p1 + scale_y_continuous(trans = "log10")
p1 <- p1 + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p1
ggsave("moments_1-epoch.png", p1, device = "png", width = 12, height = 12)

# two epochs plot
p2 <- ggplot(data = filter(dat_stats, num_epochs == 2), aes(x = variable, y = value, shape = type, color = as.factor(r))) + facet_grid(Ne_0~m_01)
p2 <- p2 + geom_point(size = 5) + theme_bw()
p2 <- p2 + scale_shape_manual(values = c(0, 1))
p2 <- p2 + labs(title = "2-Epoch Moments x N_0 (rows) and m_01 (cols)", x = "Moment", y = "Value")
p2 <- p2 + scale_y_continuous(trans = "log10")
p2 <- p2 + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
p2
ggsave("moments_2-epochs.png", p2, device = "png", width = 12, height = 12)

