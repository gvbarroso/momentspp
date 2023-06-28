
setwd("~/Devel/momentspp/poster/")

library(tidyverse)
library(cowplot)
library(scales)
library(magrittr)

scale.4d <- function(x) sprintf("%.4f", x)

r  <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
s <- c(0, -1e-6, -1e-5, -1e-4, -1e-3)
params <- crossing(-s, r)
params$s <- - params$`-s`
params <- params[,2:3]

commom_stats <- c("DD_0_0", "Dr_0_0_l_0", "Hl_0_0", "Hr_0_0", "pi2_0_0_0_0")
order <- c(2, 5, 10, 20, 100)
model <- 1:30

# constant Ne
vals <- numeric()
for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    tbl <- read.table(paste("constant/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    tbl <- subset(tbl, V1 %in% commom_stats)$V3
    
    vals <- c(vals, tbl)
  }
}

tbl_c <- as.data.frame(vals)
tbl_c$stats <- commom_stats

y <- numeric()
for(i in 1:30) { 
  y <- c(y, rep(i, 25))
}
tbl_c$model <- y

z <- numeric()
for(i in 1:30) { 
  z <- c(z, rep(2, 5), rep(5, 5), rep(10, 5), rep(20, 5), rep(100, 5))
}
tbl_c$order <- z

tbl_c$r <- rep(c(rep(1e-8, 25), rep(1e-7, 25), rep(1e-6, 25), rep(1e-5, 25), rep(1e-4, 25), rep(1e-3, 25)), 5)
tbl_c$s <- 2*5e+4*c(rep(0, 150), rep(-1e-6, 150), rep(-1e-5, 150), rep(-1e-4, 150), rep(-1e-3, 150))
tbl_c$demo <- "constant"

# 10x growth
vals <- numeric()
for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    tbl <- read.table(paste("growth/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    tbl <- subset(tbl, V1 %in% commom_stats)$V3
    
    vals <- c(vals, tbl)
  }
}

tbl_g <- as.data.frame(vals)
tbl_g$stats <- commom_stats

y <- numeric()
for(i in 1:30) { 
  y <- c(y, rep(i, 25))
}
tbl_g$model <- y

z <- numeric()
for(i in 1:30) { 
  z <- c(z, rep(2, 5), rep(5, 5), rep(10, 5), rep(20, 5), rep(100, 5))
}
tbl_g$order <- z

tbl_g$r <- rep(c(rep(1e-8, 25), rep(1e-7, 25), rep(1e-6, 25), rep(1e-5, 25), rep(1e-4, 25), rep(1e-3, 25)), 5)
tbl_g$s <- 2*5e+4*c(rep(0, 150), rep(-1e-6, 150), rep(-1e-5, 150), rep(-1e-4, 150), rep(-1e-3, 150))
tbl_g$demo <- "10x_growth"

tbl <- rbind.data.frame(tbl_c, tbl_g)
write.table(tbl, "sel_models.csv", quote=F, sep=",")

p1 <- ggplot(data=tbl[stats==commom_stats[1],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p1 <- p1 + geom_point(size=3) + theme_bw()
p1 <- p1 + scale_shape_manual(values=(4:0))
p1 <- p1 + scale_color_manual(values=c("black", "red"))
p1 <- p1 + scale_x_log10(breaks = r)
p1 <- p1 + scale_y_log10()
p1 <- p1 + labs(title=NULL, x=NULL, y=expression(D^2), shape="Order (1-2p)")
p1 <- p1 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position="none")


p2 <- ggplot(data=tbl[stats==commom_stats[2],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p2 <- p2 + geom_point(size=3) + theme_bw()
p2 <- p2 + scale_shape_manual(values=(4:0))
p2 <- p2 + scale_color_manual(values=c("black", "red"))
p2 <- p2 + scale_x_log10(breaks = r)
p2 <- p2 + scale_y_log10()
p2 <- p2 + labs(title=NULL, x=NULL, y=expression(Dz), shape="Order (1-2p)")
p2 <- p2 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p3 <- ggplot(data=tbl[stats==commom_stats[3],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p3 <- p3 + geom_point(size=3) + theme_bw()
p3 <- p3 + scale_shape_manual(values=(4:0))
p3 <- p3 + scale_color_manual(values=c("black", "red"))
p3 <- p3 + scale_x_log10(breaks = r)
p3 <- p3 + scale_y_log10()
p3 <- p3 + labs(title=NULL, x=NULL, y=expression(H[l]), shape="Order (1-2p)")
p3 <- p3 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p4 <- ggplot(data=tbl[stats==commom_stats[4],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p4 <- p4 + geom_point(size=3) + theme_bw()
p4 <- p4 + scale_shape_manual(values=(4:0))
p4 <- p4 + scale_color_manual(values=c("black", "red"))
p4 <- p4 + scale_x_log10(breaks = r)
p4 <- p4 + scale_y_log10()
p4 <- p4 + labs(title=NULL, x=NULL, y=expression(H[r]), shape="Order (1-2p)")
p4 <- p4 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p5 <- ggplot(data=tbl[stats==commom_stats[5],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p5 <- p5 + geom_point(size=3) + theme_bw()
p5 <- p5 + scale_shape_manual(values=(4:0))
p5 <- p5 + scale_color_manual(values=c("black", "red"))
p5 <- p5 + scale_x_log10(breaks = r)
p5 <- p5 + scale_y_log10()
p5 <- p5 + labs(title=NULL, x="r", y=expression(pi[2]), shape="Order (1-2p)")
p5 <- p5 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 strip.text.x = element_blank(),
                 legend.position="bottom")


moms <- plot_grid(p1, p2, p3, p4, p5, ncol=1, rel_heights = c(1, 1, 1, 1, 1.55))
save_plot("moms_sel.pdf", moms, base_height=10, base_width=20)

# comparing to moments.TwoLocus
vals <- numeric()
for(i in 1:length(model)) {
    
  tbl <- read.table(paste("constant/m", model[i], "_py.txt", sep=""))$V3
  vals <- c(vals, tbl)
}

tbl_py <- as.data.frame(vals)
tbl_py$stats <- commom_stats

# from 2p(1-p) to p(1-p)
tbl_py[tbl_py$stats=="Hl_0_0", 1] <- tbl_py[tbl_py$stats=="Hl_0_0", 1] / 2
tbl_py[tbl_py$stats=="Hr_0_0", 1] <- tbl_py[tbl_py$stats=="Hr_0_0", 1] / 2

tbl_pp <- tbl_c[tbl_c$order==100,]
tbl_pp$ratio <- tbl_pp$vals / tbl_py$vals

p6 <- ggplot(data=tbl_pp[stats==commom_stats[1],], aes(x=r, y=ratio)) + facet_wrap(~s, nrow=1)
p6 <- p6 + geom_point(size=3) + theme_bw()
p6 <- p6 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p6 <- p6 + labs(title=NULL, x=NULL, y=expression(D^2), shape="Order (1-2p)")
p6 <- p6 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position="none")


p7 <- ggplot(data=tbl_pp[stats==commom_stats[2],], aes(x=r, y=ratio)) + facet_wrap(~s, nrow=1)
p7 <- p7 + geom_point(size=3) + theme_bw()
p7 <- p7 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p7 <- p7 + labs(title=NULL, x=NULL, y=expression(Dz), shape="Order (1-2p)")
p7 <- p7 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p8 <- ggplot(data=tbl_pp[stats==commom_stats[3],], aes(x=r, y=ratio)) + facet_wrap(~s, nrow=1)
p8 <- p8 + geom_point(size=3) + theme_bw()
p8 <- p8 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p8 <- p8 + labs(title=NULL, x=NULL, y=expression(H[l]), shape="Order (1-2p)")
p8 <- p8 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p9 <- ggplot(data=tbl_pp[stats==commom_stats[4],], aes(x=r, y=ratio)) + facet_wrap(~s, nrow=1)
p9 <- p9 + geom_point(size=3) + theme_bw()
p9 <- p9 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p9 <- p9 + labs(title=NULL, x=NULL, y=expression(H[r]), shape="Order (1-2p)")
p9 <- p9 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p10 <- ggplot(data=tbl_pp[stats==commom_stats[5],], aes(x=r, y=ratio)) + facet_wrap(~s, nrow=1)
p10 <- p10 + geom_point(size=3) + theme_bw()
p10 <- p10 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p10 <- p10 + labs(title=NULL, x="r", y=expression(pi[2]), shape="Order (1-2p)")
p10 <- p10 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 strip.text.x = element_blank(),
                 legend.position="bottom")

moms_comp <- plot_grid(p6, p7, p8, p9, p10, ncol=1, rel_heights = c(1, 1, 1, 1, 1.55))
save_plot("moms_++_vs_py.pdf", moms_comp, base_height=10, base_width=20)

