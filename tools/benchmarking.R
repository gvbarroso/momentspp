library(matlib)
library(tidyverse)
library(scales)
library(patchwork) # To display 2 charts together

setwd("Devel/momentspp/benchmarking/")

two_pops_params <- as.data.frame(matrix(nrow = 16, ncol = 10))
names(two_pops_params) <- c("mu", "r", "m_12", "m_21", "N_1_epoch_1", "N_2_epoch_1", "N_1_epoch_2", "N_2_epoch_2", "num_gen_epoch_2", "num_epochs")

two_pops_params$mu <- 1.5e-08
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

MultTimeTable <- as.data.frame(matrix(nrow = 5, ncol = 5))
names(MultTimeTable) <- c("num_pops", "num_stats", "build_model", "mult_naive", "mult_eigendec")

stats_1 <- read.table("pops_1/moments.txt")
stats_2 <- read.table("pops_2/moments.txt")
stats_4 <- read.table("pops_4/moments.txt")
stats_8 <- read.table("pops_8/moments.txt")
stats_10 <- read.table("pops_10/moments.txt")

MultTimeTable$num_pops <- c(1, 2, 4, 8, 10)
MultTimeTable$num_stats <- c(nrow(stats_1), nrow(stats_2), nrow(stats_4), nrow(stats_8), nrow(stats_10))
MultTimeTable$build_model <- c(0.000308, 0.00212, 0.104, 381, 6.33e+03)
MultTimeTable$mult_naive <- c(0.000587, 0.0256, 4.07, 1.17e+03, 7.07e+03)
MultTimeTable$mult_eigendec <- c(1.7e-05, 0.0001, 0.0174, 63.1, 1.16e+03)

molten_bench <- pivot_longer(MultTimeTable, cols = c("build_model", "mult_naive", "mult_eigendec"), names_to = "variable")

bench_plot <- ggplot(data = molten_bench, aes(x = num_pops, y = value, shape = variable))
bench_plot <- bench_plot + geom_line(data = molten_bench)
bench_plot <- bench_plot + geom_point(aes(shape = variable), size = 4)
bench_plot <- bench_plot + scale_x_continuous(breaks = c(1, 2, 4, 8, 10)) 
bench_plot <- bench_plot + scale_y_continuous(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4), trans="log10")
bench_plot <- bench_plot + scale_shape_manual(values = c(0, 1, 2))
bench_plot <- bench_plot + labs(title = "Benchmark (apply transformation to 1000 generations)", x = "Num Pops.", y = "Time (seconds)") + theme_bw()
bench_plot <- bench_plot + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.position = "bottom")
bench_plot

# single population case
num_stats <- 5 # DD, Dz, H, I, Pi2, == num eigenvalues

a <- c("low_Ne", "mid_Ne", "high_Ne")
b <- c("low_mu", "mid_mu", "high_mu")
c <- c("low_r", "mid_r", "high_r")

EigenValuesTableReal <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c), ncol = num_stats))
EigenVectorsTableReal <- as.data.frame(matrix(nrow = length(a) * length(b) * length(c), ncol = num_stats))

idx = 1
for(x in 1:length(a)) {
  for(y in 1:length(b)) {
    for(z in 1:length(c)) {
      reals <- read.table(paste(paste("pops_1", a[x], b[y], c[z], sep = "/"), "/eigenvalues_real.txt", sep = ""))
      EigenValuesTableReal[idx,] <- t(reals)
      
      vecs <- read.table(paste(paste("pops_1", a[x], b[y], c[z], sep = "/"), "/eigenvector_real.txt", sep = ""))
      EigenVectorsTableReal[idx,] <- t(vecs)
      
      idx = idx + 1
    }
  }
}

EigenValuesTableReal <- select(EigenValuesTableReal, order(colSums(EigenValuesTableReal), decreasing = T))
names(EigenValuesTableReal) <- paste("e", 1:num_stats, sep = "")

EigenValuesTableReal$Ne <- c(rep(1e+6, 9), rep(1e+7, 9), rep(1e+8, 9))
EigenValuesTableReal$mu <- c(rep(c(rep(1e-8, 3), rep(1e-7, 3), rep(1e-6, 3)), 3))
EigenValuesTableReal$r <- rep(c(1e-8, 1e-7, 1e-6), 9)

EigenValuesTableReal$theta <- 4 * EigenValuesTableReal$Ne * EigenValuesTableReal$mu
EigenValuesTableReal$rho <- 4 * EigenValuesTableReal$Ne * EigenValuesTableReal$r
EigenValuesTableReal$ratio <- EigenValuesTableReal$theta / EigenValuesTableReal$rho

dat <- pivot_longer(EigenValuesTableReal, cols = starts_with("e"), names_to = "variable")

p <- ggplot(data = dat, aes(x = Ne, y = value, shape = variable, color = log(ratio))) + facet_grid(r~mu)
p <- p + geom_point(size = 3) + theme_bw()
p <- p + scale_shape_manual(values = c(0, 1, 2, 3, 4))
p <- p + scale_x_continuous(trans = "log10", breaks = c(1e+6, 1e+7, 1e+8)) 
p <- p + scale_y_continuous(trans = "log10", breaks = pretty_breaks())
p <- p + labs(title = "Eigenvalues x r (rows) and u (cols)", x = "Ne", y = "Eigenvalue")
p <- p + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = "bottom")
p

ggsave("eigen_vals.pdf", p, device = "pdf", width = 12, height = 12)




m1 <- lm(e1 ~ (Ne + mu + r), data = EigenValuesTableReal)
m2 <- lm(e2 ~ (Ne + mu + r), data = EigenValuesTableReal)
m3 <- lm(e3 ~ (Ne + mu + r), data = EigenValuesTableReal)
m4 <- lm(e4 ~ (Ne + mu + r), data = EigenValuesTableReal)
m5 <- lm(e5 ~ (Ne + mu + r), data = EigenValuesTableReal)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

filter(EigenValuesTableReal, Ne == 1e+6) %>% e2 / mu

p1_stats <- c("DD", "Dz", "H", "I", "Pi2")



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
