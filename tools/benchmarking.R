library(matlib)
library(reshape2)
library(tidyverse)
library(scales)
library(patchwork) # To display 2 charts together

#######################
#
# timing
#
#######################

setwd("Devel/momentspp/benchmarking/")

MultTimeTable <- as.data.frame(matrix(nrow = 5, ncol = 5))
names(MultTimeTable) <- c("num_pops", "num_stats", "build_model", "mult_naive", "mult_eigendec")

stats_1 <- read.table("pops_1/moments.txt")
stats_2 <- read.table("pops_2/moments.txt")
stats_4 <- read.table("pops_4/moments.txt")
stats_8 <- read.table("pops_8/moments.txt")
stats_10 <- read.table("pops_8/moments.txt")

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

#######################
#
# misc (mess)
#
#######################

U <- matrix(c(1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 4e-4,
              0, 0, 2e-4, 1, 0,
              0, 0, 0, 0, 1), 5, 5)


R <- matrix(c(9.998000e-01, 0, 0, 0, 0,
              0, 9.999000e-01, 0, 0, 0,
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








