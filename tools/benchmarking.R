library(matlib)
library(tidyverse)
library(scales)
library(patchwork) # To display 2 charts together

setwd("Devel/momentspp/benchmarking/")

##############################
#
# timing
#
##############################

MultTimeTable <- as.data.frame(matrix(nrow = 4, ncol = 5))
names(MultTimeTable) <- c("num_pops", "num_stats", "direct", "eigen", "pseudo-eigen")

for(i in 1:4) {
  p <- 2^(i-1)
  stats <- read.table(paste("timing/pops_", p, "/final_moments.txt", sep = ""))
  timing <- read.table(paste("timing/pops_", p, "/timing.txt", sep = ""))
  MultTimeTable[i,] <- c(p, nrow(stats), t(timing))
}

molten_time <- pivot_longer(MultTimeTable, cols = c("direct", "eigen", "pseudo-eigen"), names_to = "variable")

bench_plot <- ggplot(data = molten_bench, aes(x = num_pops, y = value, shape = variable))
bench_plot <- bench_plot + geom_line(data = molten_bench)
bench_plot <- bench_plot + geom_point(aes(shape = variable), size = 4)
bench_plot <- bench_plot + scale_x_continuous(breaks = c(1, 2, 4, 8)) 
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

a <- c("low_Ne", "mid_Ne", "high_Ne")
b <- c("low_mu", "mid_mu", "high_mu")
c <- c("low_r", "mid_r", "high_r")

idx = 1
for(x in 1:length(a)) {
  for(y in 1:length(b)) {
    for(z in 1:length(c)) {
      
      
      idx = idx + 1
    }
  }
}

EigenValuesTableReal$Ne <- c(rep(1e+3, 9), rep(1e+4, 9), rep(1e+5, 9))
EigenValuesTableReal$mu <- c(rep(c(rep(1e-8, 3), rep(1e-7, 3), rep(1e-6, 3)), 3))
EigenValuesTableReal$r <- rep(c(1e-8, 1e-7, 1e-6), 9)

EigenValuesTableReal$theta <- 4 * EigenValuesTableReal$Ne * EigenValuesTableReal$mu
EigenValuesTableReal$rho <- 4 * EigenValuesTableReal$Ne * EigenValuesTableReal$r
EigenValuesTableReal$ratio <- EigenValuesTableReal$theta / EigenValuesTableReal$rho

dat <- pivot_longer(EigenValuesTableReal, cols = starts_with("e"), names_to = "variable")

p <- ggplot(data = dat, aes(x = Ne, y = value, shape = variable, color = log(ratio))) + facet_grid(r~mu)
p <- p + geom_point(size = 3) + theme_bw()
p <- p + scale_shape_manual(values = c(0, 1, 2, 3, 4))
p <- p + scale_x_continuous(trans = "log10", breaks = c(1e+3, 1e+4, 1e+5)) 
p <- p + scale_y_continuous(trans = "log10", breaks = pretty_breaks())
p <- p + labs(title = "Eigenvalues x r (rows) and u (cols)", x = "Ne", y = "Eigenvalue")
p <- p + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = "bottom")
p

ggsave("eigen_vals.pdf", p, device = "pdf", width = 12, height = 12)

##############################
#
# 2-pop accuracy
#
##############################

two_pops_params <- as.data.frame(matrix(nrow = 16, ncol = 10))
names(two_pops_params) <- c("mu", "r", "m_12", "m_21", "N_1_epoch_1", "N_2_epoch_1", "N_1_epoch_2", "N_2_epoch_2", "num_gen_epoch_2", "num_epochs")

two_pops_params$mu <- 1e-08
two_pops_params$r <- c(rep(0, 8), rep(1e-9, 8))
two_pops_params$m_12 <- rep(c(c(1e-5, 1e-5), c(0, 1e-5)), 4)
two_pops_params$m_21 <- 1e-5 
two_pops_params$N_1_epoch_1 <- rep(c(1e+4, 1e+3), 8)
two_pops_params$N_2_epoch_1 <- 1e+4
two_pops_params$num_gen_epoch_2 <- rep(c(rep(5e+3, 4), rep(NA, 4)), 2)
two_pops_params$N_1_epoch_2 <- rep(c(rep(c(1e+3, 1e+2), 2), rep(NA, 4)), 2)
two_pops_params$N_2_epoch_2 <- rep(c(rep(1e+3, 4), rep(NA, 4)), 2)
two_pops_params$num_epochs <- rep(c(rep(2, 4), rep(1, 4)), 2)

write.table(two_pops_params, "bench_params.tsv", sep = ",", row.names = F)

pop_indices <- c(0, 1)
DD_stats <- NULL
Dz_stats <- NULL
H_stats <- NULL
pi2_stats <- NULL

for(i in 1:length(pop_indices)) {
  for(j in 1:length(pop_indices)) {
    
    DD_stats <- c(DD_stats, paste("s_DD", pop_indices[i], pop_indices[j], sep = "_"))
    H_stats <- c(H_stats, paste("s_H", pop_indices[i], pop_indices[j], sep = "_"))
    
    for(k in 1:length(pop_indices)) {
      
      Dz_stats <- c(Dz_stats, paste("s_Dz", pop_indices[i], pop_indices[j], pop_indices[k], sep = "_"))
      
      for(l in 1:length(pop_indices)) {
        
        pi2_stats <- c(pi2_stats, paste("s_pi2", pop_indices[i], pop_indices[j], pop_indices[k], pop_indices[l], sep = "_"))
      }
    }
  }
}

p2_stats <- c(DD_stats, Dz_stats, H_stats, "I", pi2_stats)
num_stats <- length(p2_stats)

a <- c("zero_rec", "nonzero_rec")
b <- c("single_epoch", "two_epochs")
c <- c("symmetric_mig", "asymmetric_mig")
d <- c("same_Ne", "diff_Ne")

# params table (scenarios)
params <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c) * length(d), ncol = 7))
names(params) <- c("N_0", "N_1", "m_01", "m_10", "num_epochs", "r", "mu")

params$N_0 <- rep(c(1e+4, 1e+3), 8)
params$N_1 <- rep(1e+4, 16)
params$m_01 <- rep(c(rep(1e-5, 2), rep(0, 2)), 4)
params$m_10 <- rep(1e-5, 16)
params$num_epochs <- rep(c(rep(2, 4), rep(1, 4)), 2)
params$r <- c(rep(0, 8), rep(1e-8, 8))
params$mu <- 1e-8

EigenValuesTableReal <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c) * length(d), ncol = num_stats))
EigenValuesTableImag <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c) * length(d), ncol = num_stats))
SteadyStateTable <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c) * length(d), ncol = num_stats))
FinalMomentsTable <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c) * length(d), ncol = num_stats))

eigen_list <- list()  # one EigenVectors matrix per scenario

idx = 1
for(x in 1:length(a)) {
  for(y in 1:length(b)) {
    for(z in 1:length(c)) {
      for(w in 1:length(d)) {
        
        reals <- read.table(paste(paste("pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/eigenvalues_real.txt", sep = ""))
        EigenValuesTableReal[idx,] <- t(reals)
        
        imag <- read.table(paste(paste("pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/eigenvalues_imag.txt", sep = ""))
        EigenValuesTableImag[idx,] <- t(imag)
        
        steady_state <- read.table(paste(paste("pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/steady_state.txt", sep = ""))
        SteadyStateTable[idx,] <- t(steady_state$V5)
        
        final_moments <- read.table(paste(paste("pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/final_moments.txt", sep = ""))
        FinalMomentsTable[idx,] <- t(final_moments$V5)
        
        tbl <- as.data.frame(matrix(nrow = num_stats, ncol = num_stats))
        for(i in 1:num_stats) {
          tbl[i, ] <- t(read.table(paste(paste("pops_2", a[x], b[y], c[z], d[w], sep = "/"), "/eigenvector_", i - 1, "_real.txt", sep = "")))
          names(tbl) <- p2_stats
        }
        
        eigen_list[[idx]] <- tbl
        
        idx = idx + 1
      }
    }
  }
}

names(SteadyStateTable) <- p2_stats
names(FinalMomentsTable) <- p2_stats
names(EigenValuesTableReal) <- paste("e", 1:num_stats, sep = "")
names(EigenValuesTableImag) <- paste("e", 1:num_stats, sep = "")

eigen_reals <- bind_cols(EigenValuesTableReal, params)
eigen_imag <- bind_cols(EigenValuesTableImag, params)
steadY <- bind_cols(SteadyStateTable, params)
finY <- bind_cols(FinalMomentsTable, params)

steadY$has_imag <- rowSums(abs(EigenValuesTableImag)) > 0
finY$has_imag <- rowSums(abs(EigenValuesTableImag)) > 0

write.table(steadY, "steadY_tbl.tsv", sep = ",", row.names = F, quote = F)

dat_imag <- pivot_longer(eigen_imag, cols = starts_with("e"), names_to = "variable")
dat_imag_only <- filter(dat_imag, value != 0)
dat_imag_only$value <- abs(dat_imag_only$value) # sends conjugate pairs to the same value in order to plot in log scale

p <- ggplot(data = dat_imag_only, aes(x = r, y = value, color = variable)) + facet_grid(N_0~m_01)
p <- p + geom_point(size = 3, alpha = 0.2) + theme_bw()
p <- p + labs(title = "Eigenvalues x N_0 (rows) and m_01 (cols)", x = "r", y = "Eigenvalue")
p <- p + scale_y_continuous(trans = "log10", breaks = exp(seq(log(1), log(1e-18), length.out=19)))
p <- p + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = "bottom")
p
ggsave("eigen_imag.png", p, device = "png", width = 9, height = 6)

dat_steadY <- pivot_longer(steadY, cols = starts_with("s_"), names_to = "variable")

s <- ggplot(data = dat_steadY, aes(x = variable, y = value, color = variable, shape = has_imag)) + facet_grid(N_0~m_01)
s <- s + geom_point(size = 3, alpha = 0.2) + theme_bw()
s <- s + labs(title = "Stat x N_0 (rows) and m_01 (cols)", x = "stat", y = "eq. value")
s <- s + scale_y_continuous(trans = "log10", breaks = exp(seq(log(1), log(1e-12), length.out=13)))
s <- s + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
s
ggsave("steady.png", s, device = "png")

dat_fin <- pivot_longer(finY, cols = starts_with("s_"), names_to = "variable")

f <- ggplot(data = filter(dat_fin), aes(x = variable, y = value, color = variable, shape = has_imag)) + facet_grid(N_0~m_01)
f <- f + geom_point(size = 3, alpha = 0.5) + theme_bw()
f <- f + labs(title = "Stats x N_0 (rows) and m_01 (cols)", x = "stat", y = "eq. value")
f <- f + scale_y_continuous(trans = "log10", breaks = exp(seq(log(1), log(1e-12), length.out=13)))
f <- f + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
f
ggsave("fin.png", f, device = "png")

f1 <- ggplot(data = filter(dat_fin, num_epochs == 1), aes(x = variable, y = value, color = variable, shape = has_imag)) + facet_grid(N_0~m_01)
f1 <- f1 + geom_point(size = 3, alpha = 0.5) + theme_bw()
f1 <- f1 + labs(title = "Stats [1 epoch] x N_0 (rows) and m_01 (cols)", x = "stat", y = "eq. value")
f1 <- f1 + scale_y_continuous(trans = "log10", breaks = exp(seq(log(1), log(1e-12), length.out=13)))
f1 <- f1 + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
f1
ggsave("fin_1_epoch.png", f1, device = "png")

f2 <- ggplot(data = filter(dat_fin, num_epochs == 2), aes(x = variable, y = value, color = variable, shape = has_imag)) + facet_grid(N_0~m_01)
f2 <- f2 + geom_point(size = 3, alpha = 0.5) + theme_bw()
f2 <- f2 + labs(title = "Stats [2 epochs] x N_0 (rows) and m_01 (cols)", x = "stat", y = "eq. value")
f2 <- f2 + scale_y_continuous(trans = "log10", breaks = exp(seq(log(1), log(1e-12), length.out=13)))
f2 <- f2 + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
f2
ggsave("fin_2_epochs.png", f2, device = "png")

View(filter(dat_fin, num_epochs == 2, m_01 == 1e-5, N_0 == 1e+3))

#######################
#
# misc (mess)
#
#######################

U <- matrix(c(1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 6e-8,
              0, 0, 3e-8, 1, 0,
              0, 0, 0, 0, 1), 5, 5)


R <- matrix(c(1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1), 5, 5)

D <- matrix(c(9.997e-01, 4e-04, 0, 0, 0,
              1e-04, 9.995e-01, 0, 0, 1e-04,
              0, 0, 9.999e-01, 0, 0,
              0, 0, 0, 1, 0,
              1e-04, 0, 0, 0, 9.998e-01), 5, 5)


Z1 <- R %*% D
Z2 <- R %*% U
Z3 <- D %*% R
Z4 <- D %*% U
Z5 <- U %*% R
Z6 <- U %*% D

X1 <- U %*% Z1
X2 <- U %*% Z3
X3 <- R %*% Z4
X4 <- R %*% Z6
X5 <- D %*% Z2
X6 <- D %*% Z5

ev1 <- eigen(X1)
ev2 <- eigen(X2)
ev3 <- eigen(X3)
ev4 <- eigen(X4)
ev5 <- eigen(X5)
ev6 <- eigen(X6)


v <- c(0.1, 0.1, 1e-3, 1, 1e-6)
v <- rep(1, 5)
names(v) <- c("D2", "Dz", "H", "I", "Pi2")

M1 <- X1

for(i in 1:1e+7) {
  #v = X1 %*% v
  #M1 = M1 %*% M1
  W <- W %*% W
}
W
M1
v
X1 %*% ev1$vectors[,1]
M1 %*% ev1$vectors[,1]

M <- matrix(c(9.995001e-01, 9.998000e-05, 0, 0, 9.998000e-05, 
3.999600e-04, 9.994001e-01, 0, 0, 0 ,
0, 0, 9.999000e-01, 2.000000e-04, 0, 
0, 0, 0, 1.000000e+00, 0, 
0, 1.000000e-04, 3.999600e-04, 0, 9.998000e-01), 5, 5)

W <- t(M)
z <- eigen(W)
z


X_inv <- inverse(ev6$vectors)


s <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
N <- c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7)
tbl <- crossing(N, s)
tbl$alpha <- 2 * tbl$N * tbl$s
tbl$p_fix <- tbl$s / (1 - exp(-tbl$alpha))

p_plot <- ggplot(data = tbl, aes(x = alpha, y = p_fix, color = log(N), shape = as.factor(s)))
p_plot <- p_plot + geom_point(size = 3) + theme_bw()
p_plot <- p_plot + scale_shape_manual(values = c(15, 16, 17, 18, 19, 20))
p_plot <- p_plot + scale_x_continuous(trans = "log10", breaks = unique(tbl$alpha)) 
p_plot <- p_plot + scale_y_continuous(trans = "log10", breaks = c(1e-5, 1e-3, 1e-2, 1e-1))
p_plot <- p_plot + labs(title = NULL, x = "Alpha", y = "P(fix)")
p_plot <- p_plot + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = "right")
p_plot
