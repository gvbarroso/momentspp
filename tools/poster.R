
setwd("~/Devel/momentspp/poster/selected_1-pop")

library(tidyverse)
library(cowplot)
library(scales)
library(magrittr)

scale.4d <- function(x) sprintf("%.4f", x)

r  <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
s <- c(0, -1e-6, -1e-5, -5e-5, -1e-4, -2e-4, -1e-3)
params <- crossing(-s, r)
params$s <- - params$`-s`
params <- params[,2:3]

commom_stats <- c("DD_0_0", "Dr_0_0_l_0", "Hl_0_0", "Hr_0_0", "pi2_0_0_0_0")
order <- c(2, 5, 10, 20, 100)
model <- 1:42

# constant Ne
vals <- numeric()
for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    x <- read.table(paste("constant/mu_1e-7/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    x <- subset(x, V1 %in% commom_stats)$V3
    
    vals <- c(vals, x)
  }
}

for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    x <- read.table(paste("constant/mu_1e-8/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    x <- subset(x, V1 %in% commom_stats)$V3
    
    vals <- c(vals, x)
  }
}

for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    x <- read.table(paste("constant/mu_1e-9/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    x <- subset(x, V1 %in% commom_stats)$V3
    
    vals <- c(vals, x)
  }
}

tbl_c <- as.data.frame(vals)
tbl_c$stats <- commom_stats
tbl_c$mu <- c(rep(1e-7, 1050), rep(1e-8, 1050), rep(1e-9, 1050))

y <- numeric()
for(i in model) { 
  y <- c(y, rep(i, 25))
}
tbl_c$model <- y

z <- numeric()
for(i in model) { 
  z <- c(z, rep(2, 5), rep(5, 5), rep(10, 5), rep(20, 5), rep(100, 5))
}
tbl_c$order <- z

tbl_c$r <- rep(c(rep(1e-8, 25), rep(1e-7, 25), rep(1e-6, 25), rep(1e-5, 25), rep(1e-4, 25), rep(1e-3, 25)), 7)
tbl_c$s <- 2*5e+4*c(rep(0, 150), rep(-1e-6, 150), rep(-1e-5, 150), rep(-5e-5, 150), rep(-1e-4, 150), rep(-2e-4, 150), rep(-1e-3, 150))
tbl_c$demo <- "constant"

q1 <- ggplot(data=tbl_c[tbl_c$stats==commom_stats[1],], aes(x=r, y=vals, shape=as.factor(order), color=as.factor(mu))) + facet_wrap(~s, nrow=1)
q1 <- q1 + geom_point(size=3) + theme_bw()
q1 <- q1 + scale_shape_manual(values=(4:0))
q1 <- q1 + scale_color_manual(values=c("black", "red", "blue"))
q1 <- q1 + scale_x_log10(breaks = r)
q1 <- q1 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q1 <- q1 + labs(title=NULL, x="r", y=expression(D^2), shape="Order (1-2p)", color="Mu")
q1 <- q1 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 legend.position="bottom")
ggsave("DD_mut_rates.pdf", q1, device="pdf", width=16, height=12)

q2 <- ggplot(data=tbl_c[tbl_c$stats==commom_stats[2],], aes(x=r, y=vals, shape=as.factor(order), color=as.factor(mu))) + facet_wrap(~s, nrow=1)
q2 <- q2 + geom_point(size=3) + theme_bw()
q2 <- q2 + scale_shape_manual(values=(4:0))
q2 <- q2 + scale_color_manual(values=c("black", "red", "blue"))
q2 <- q2 + scale_x_log10(breaks = r)
q2 <- q2 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q2 <- q2 + labs(title=NULL, x="r", y=expression(Dz), shape="Order (1-2p)", color="Mu")
q2 <- q2 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 legend.position="bottom")
ggsave("Dz_mut_rates.pdf", q2, device="pdf", width=16, height=12)

q3 <- ggplot(data=tbl_c[tbl_c$stats==commom_stats[3],], aes(x=r, y=vals, shape=as.factor(order), color=as.factor(mu))) + facet_wrap(~s, nrow=1)
q3 <- q3 + geom_point(size=3) + theme_bw()
q3 <- q3 + scale_shape_manual(values=(4:0))
q3 <- q3 + scale_color_manual(values=c("black", "red", "blue"))
q3 <- q3 + scale_x_log10(breaks = r)
q3 <- q3 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q3 <- q3 + labs(title=NULL, x="r", y=expression(H[l]), shape="Order (1-2p)", color="Mu")
q3 <- q3 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 legend.position="bottom")
ggsave("Hl_mut_rates.pdf", q3, device="pdf", width=16, height=12)

q4 <- ggplot(data=tbl_c[tbl_c$stats==commom_stats[4],], aes(x=r, y=vals, shape=as.factor(order), color=as.factor(mu))) + facet_wrap(~s, nrow=1)
q4 <- q4 + geom_point(size=3) + theme_bw()
q4 <- q4 + scale_shape_manual(values=(4:0))
q4 <- q4 + scale_color_manual(values=c("black", "red", "blue"))
q4 <- q4 + scale_x_log10(breaks = r)
q4 <- q4 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q4 <- q4 + labs(title=NULL, x="r", y=expression(H[r]), shape="Order (1-2p)", color="Mu")
q4 <- q4 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 legend.position="bottom")
ggsave("Hr_mut_rates.pdf", q4, device="pdf", width=16, height=12)

t <- tbl_c[tbl_c$stats==commom_stats[4],]
q4b <- ggplot(data=t[t$mu==1e-7,], aes(x=r, y=vals, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
q4b <- q4b + geom_point(size=3) + theme_bw()
q4b <- q4b + scale_shape_manual(values=(4:0))
q4b <- q4b + scale_x_log10(breaks = r)
q4b <- q4b + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q4b <- q4b + labs(title="Mu=1e-7", x="r", y=expression(H[r]), shape="Order (1-2p)", color="Mu")
q4b <- q4b + theme(axis.title=element_text(size=12),
                   axis.text=element_text(size=10),
                   axis.text.x = element_text(angle = 90),
                   legend.position="bottom")
ggsave("Hr_mu_1e-7.pdf", q4b, device="pdf", width=16, height=12)

q4c <- ggplot(data=t[t$mu==1e-8,], aes(x=r, y=vals, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
q4c <- q4c + geom_point(size=3) + theme_bw()
q4c <- q4c + scale_shape_manual(values=(4:0))
q4c <- q4c + scale_x_log10(breaks = r)
q4c <- q4c + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q4c <- q4c + labs(title="Mu=1e-8", x="r", y=expression(H[r]), shape="Order (1-2p)", color="Mu")
q4c <- q4c + theme(axis.title=element_text(size=12),
                   axis.text=element_text(size=10),
                   axis.text.x = element_text(angle = 90),
                   legend.position="bottom")
ggsave("Hr_mu_1e-8.pdf", q4c, device="pdf", width=16, height=12)

q4d <- ggplot(data=t[t$mu==1e-9,], aes(x=r, y=vals, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
q4d <- q4d + geom_point(size=3) + theme_bw()
q4d <- q4d + scale_shape_manual(values=(4:0))
q4d <- q4d + scale_x_log10(breaks = r)
q4d <- q4d + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q4d <- q4d + labs(title="Mu=1e-9", x="r", y=expression(H[r]), shape="Order (1-2p)", color="Mu")
q4d <- q4d + theme(axis.title=element_text(size=12),
                   axis.text=element_text(size=10),
                   axis.text.x = element_text(angle = 90),
                   legend.position="bottom")
ggsave("Hr_mu_1e-9.pdf", q4d, device="pdf", width=16, height=12)


q5 <- ggplot(data=tbl_c[tbl_c$stats==commom_stats[5],], aes(x=r, y=vals, shape=as.factor(order), color=as.factor(mu))) + facet_wrap(~s, nrow=1)
q5 <- q5 + geom_point(size=3) + theme_bw()
q5 <- q5 + scale_shape_manual(values=(4:0))
q5 <- q5 + scale_color_manual(values=c("black", "red", "blue"))
q5 <- q5 + scale_x_log10(breaks = r)
q5 <- q5 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
q5 <- q5 + labs(title=NULL, x="r", y=expression(pi[2]), shape="Order (1-2p)", color="Mu")
q5 <- q5 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 legend.position="bottom")
ggsave("pi2_mut_rates.pdf", q5, device="pdf", width=16, height=12)

# 10x growth
vals <- numeric()
for(i in 1:length(model)) {
  for(j in 1:length(order)) {
    
    x <- read.table(paste("growth/m", model[i], "_", order[j], "_expectations.txt", sep=""))
    x <- subset(x, V1 %in% commom_stats)$V3
    
    vals <- c(vals, x)
  }
}

tbl_g <- as.data.frame(vals)
tbl_g$stats <- commom_stats

y <- numeric()
for(i in model) { 
  y <- c(y, rep(i, 25))
}
tbl_g$model <- y

z <- numeric()
for(i in model) { 
  z <- c(z, rep(2, 5), rep(5, 5), rep(10, 5), rep(20, 5), rep(100, 5))
}
tbl_g$order <- z

tbl_g$r <- rep(c(rep(1e-8, 25), rep(1e-7, 25), rep(1e-6, 25), rep(1e-5, 25), rep(1e-4, 25), rep(1e-3, 25)), 7)
tbl_g$s <- 2*5e+4*c(rep(0, 150), rep(-1e-6, 150), rep(-1e-5, 150), rep(-5e-5, 150), rep(-1e-4, 150), rep(-2e-4, 150), rep(-1e-3, 150))
tbl_g$demo <- "10x_growth"
tbl_g$mu <- 1e-8

tbl <- rbind.data.frame(tbl_c[tbl_c$mu == 1e-8,], tbl_g)
#write.table(tbl, "sel_models.csv", quote=F, sep=",")

# for normalizing
pi2 <- tbl[tbl$stats==commom_stats[5],]$vals

p1 <- ggplot(data=tbl[tbl$stats==commom_stats[1],], aes(x=r, y=vals / pi2, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p1 <- p1 + geom_point(size=3) + theme_bw()
p1 <- p1 + scale_shape_manual(values=(4:0))
p1 <- p1 + scale_color_manual(values=c("black", "red"))
p1 <- p1 + scale_x_log10(breaks = r)
p1 <- p1 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
p1 <- p1 + labs(title=NULL, x=NULL, y=expression(sigma[d]^2), shape="Order (1-2p)")
p1 <- p1 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position="none")


p2 <- ggplot(data=tbl[tbl$stats==commom_stats[2],], aes(x=r, y=vals / pi2, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p2 <- p2 + geom_point(size=3) + theme_bw()
p2 <- p2 + scale_shape_manual(values=(4:0))
p2 <- p2 + scale_color_manual(values=c("black", "red"))
p2 <- p2 + scale_x_log10(breaks = r)
p2 <- p2 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
p2 <- p2 + labs(title=NULL, x=NULL, y=expression(Dz / pi2), shape="Order (1-2p)")
p2 <- p2 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p3 <- ggplot(data=tbl[tbl$stats==commom_stats[3],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p3 <- p3 + geom_point(size=3) + theme_bw()
p3 <- p3 + scale_shape_manual(values=(4:0))
p3 <- p3 + scale_color_manual(values=c("black", "red"))
p3 <- p3 + scale_x_log10(breaks = r)
p3 <- p3 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
p3 <- p3 + labs(title=NULL, x=NULL, y=expression(H[l]), shape="Order (1-2p)")
p3 <- p3 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p4 <- ggplot(data=tbl[tbl$stats==commom_stats[4],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p4 <- p4 + geom_point(size=3) + theme_bw()
p4 <- p4 + scale_shape_manual(values=(4:0))
p4 <- p4 + scale_color_manual(values=c("black", "red"))
p4 <- p4 + scale_x_log10(breaks = r)
p4 <- p4 + scale_y_log10(labels = function(x) format(x, scientific = TRUE), limits=c(1e-5, 1.02e-3)) 
p4 <- p4 + labs(title=NULL, x=NULL, y=expression(H[r]), shape="Order (1-2p)")
p4 <- p4 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p5 <- ggplot(data=tbl[tbl$stats==commom_stats[5],], aes(x=r, y=vals, shape=as.factor(order), color=demo)) + facet_wrap(~s, nrow=1)
p5 <- p5 + geom_point(size=3) + theme_bw()
p5 <- p5 + scale_shape_manual(values=(4:0))
p5 <- p5 + scale_color_manual(values=c("black", "red"))
p5 <- p5 + scale_x_log10(breaks = r)
p5 <- p5 + scale_y_log10(labels = function(x) format(x, scientific = TRUE))
p5 <- p5 + labs(title=NULL, x="r", y=expression(pi[2]), shape="Order (1-2p)")
p5 <- p5 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 strip.text.x = element_blank(),
                 legend.position="bottom")


moms <- plot_grid(p1, p2, p3, p4, p5, ncol=1, rel_heights = c(1, 1, 1, 1, 1.5))
save_plot("moms_sel.pdf", moms, base_height=10, base_width=30)

# comparing to moments.TwoLocus
vals <- numeric()
for(i in 1:length(model)) {
    
  x <- read.table(paste("constant/mu_1e-8/m", model[i], "_py.txt", sep=""))$V3
  vals <- c(vals, x)
}

tbl_py <- as.data.frame(vals)
tbl_py$stats <- commom_stats

# from 2p(1-p) to p(1-p)
tbl_py[tbl_py$stats=="Hl_0_0", 1] <- tbl_py[tbl_py$stats=="Hl_0_0", 1] / 2
tbl_py[tbl_py$stats=="Hr_0_0", 1] <- tbl_py[tbl_py$stats=="Hr_0_0", 1] / 2

tbl_pp <- tbl[tbl$demo == "constant",] %>% group_by(order) %>% mutate(ratio = vals / tbl_py$vals)
tbl_pp <- tbl_pp[tbl_pp$r < 1e-3,]
tbl_pp <- tbl_pp[tbl_pp$s > -100,]

p6 <- ggplot(data=tbl_pp[tbl_pp$stats==commom_stats[1],], aes(x=r, y=ratio, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
p6 <- p6 + geom_point(size=3) + theme_bw()
p6 <- p6 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p6 <- p6 + scale_shape_manual(values=(4:0))
p6 <- p6 + labs(title=NULL, x=NULL, y=expression(D^2), shape="Order (1-2p)")
p6 <- p6 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position="none")


p7 <- ggplot(data=tbl_pp[tbl_pp$stats==commom_stats[2],], aes(x=r, y=ratio, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
p7 <- p7 + geom_point(size=3) + theme_bw()
p7 <- p7 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p7 <- p7 + scale_shape_manual(values=(4:0))
p7 <- p7 + labs(title=NULL, x=NULL, y=expression(Dz), shape="Order (1-2p)")
p7 <- p7 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p8 <- ggplot(data=tbl_pp[tbl_pp$stats==commom_stats[3],], aes(x=r, y=ratio, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
p8 <- p8 + geom_point(size=3) + theme_bw()
p8 <- p8 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p8 <- p8 + scale_shape_manual(values=(4:0))
p8 <- p8 + labs(title=NULL, x=NULL, y=expression(H[l]), shape="Order (1-2p)")
p8 <- p8 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p9 <- ggplot(data=tbl_pp[tbl_pp$stats==commom_stats[4],], aes(x=r, y=ratio, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
p9 <- p9 + geom_point(size=3) + theme_bw()
p9 <- p9 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p9 <- p9 + scale_shape_manual(values=(4:0))
p9 <- p9 + labs(title=NULL, x=NULL, y=expression(H[r]), shape="Order (1-2p)")
p9 <- p9 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


p10 <- ggplot(data=tbl_pp[tbl_pp$stats==commom_stats[5],], aes(x=r, y=ratio, shape=as.factor(order))) + facet_wrap(~s, nrow=1)
p10 <- p10 + geom_point(size=3) + theme_bw()
p10 <- p10 + scale_x_log10(breaks = r) + scale_y_continuous(labels = scale.4d)
p10 <- p10 + scale_shape_manual(values=(4:0))
p10 <- p10 + labs(title=NULL, x="r", y=expression(pi[2]), shape="Order (1-2p)")
p10 <- p10 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 strip.text.x = element_blank(),
                 legend.position="bottom")

moms_comp <- plot_grid(p6, p7, p8, p9, p10, ncol=1, rel_heights = c(1, 1, 1, 1, 1.3))
save_plot("moms_++_vs_py.pdf", moms_comp, base_height=10, base_width=20)


#

setwd("~/Devel/momentspp/poster/neutral_2-pops")

r  <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
f <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
params <- crossing(f, r)

commom_stats <- read.table("m1_expectations.txt")$V1
model <- 1:42

vals <- numeric()
for(i in 1:length(model)) {
    
  x <- read.table(paste("m", model[i], "_expectations.txt", sep=""))
  vals <- c(vals, x$V3)
}

tbl_n <- as.data.frame(vals)
tbl_n$stats <- commom_stats

y <- numeric()
for(i in model) { 
  y <- c(y, rep(i, length(commom_stats)))
}
tbl_n$model <- y

tbl_n$r <- rep(c(rep(1e-8, length(commom_stats)), rep(1e-7, length(commom_stats)), rep(1e-6, length(commom_stats)), rep(1e-5, length(commom_stats)), rep(1e-4, length(commom_stats)), rep(1e-3, length(commom_stats))), 7)
tbl_n$f <- c(rep(0, 6 * length(commom_stats)), rep(0.1, 6 * length(commom_stats)), rep(0.2, 6 * length(commom_stats)), rep(0.4, 6 * length(commom_stats)), rep(0.6, 6 * length(commom_stats)), rep(0.8, 6 * length(commom_stats)), rep(1, 6 * length(commom_stats)))


r1 <- ggplot(data=tbl_n[tbl_n$stats==commom_stats[1],], aes(x=r, y=vals)) + facet_wrap(~f, nrow=1)
r1 <- r1 + geom_point(size=3) + theme_bw()
r1 <- r1 + scale_x_log10(breaks = r)
r1 <- r1 + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
r1 <- r1 + labs(title=NULL, x=NULL, y=expression(D[11]))
r1 <- r1 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position="none")


r2 <- ggplot(data=tbl_n[tbl_n$stats==commom_stats[2],], aes(x=r, y=vals)) + facet_wrap(~f, nrow=1)
r2 <- r2 + geom_point(size=3) + theme_bw()
r2 <- r2 + scale_x_log10(breaks = r)
r2 <- r2 + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
r2 <- r2 + labs(title=NULL, x=NULL, y=expression(D[1]*D[2]))
r2 <- r2 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")


r3 <- ggplot(data=tbl_n[tbl_n$stats==commom_stats[5],], aes(x=r, y=vals)) + facet_wrap(~f, nrow=1)
r3 <- r3 + geom_point(size=3) + theme_bw()
r3 <- r3 + scale_x_log10(breaks = r)
r3 <- r3 + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
r3 <- r3 + labs(title=NULL, x=NULL, y=expression(D[1]*(1-2*p[1])*(1-2*p[2])))
r3 <- r3 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),   
                 strip.text.x = element_blank(),
                 legend.position="none")

r4 <- ggplot(data=tbl_n[tbl_n$stats==commom_stats[7],], aes(x=r, y=vals)) + facet_wrap(~f, nrow=1)
r4 <- r4 + geom_point(size=3) + theme_bw()
r4 <- r4 + scale_x_log10(breaks = r)
r4 <- r4 + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
r4 <- r4 + labs(title=NULL, x=NULL, y=expression(D[2]*(1-2*p[1])*(1-2*p[1])))
r4 <- r4 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 strip.text.x = element_blank(),
                 legend.position="none")

r5 <- ggplot(data=tbl_n[tbl_n$stats==commom_stats[17],], aes(x=r, y=vals)) + facet_wrap(~f, nrow=1)
r5 <- r5 + geom_point(size=3) + theme_bw()
r5 <- r5 + scale_x_log10(breaks=r)
r5 <- r5 + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
r5 <- r5 + labs(title=NULL, x=NULL, y=expression(pi[1212]))
r5 <- r5 + theme(axis.title=element_text(size=12),
                 axis.text=element_text(size=10),
                 axis.text.x = element_text(angle = 90),
                 strip.text.x = element_blank(),
                 legend.position="none")



moms <- plot_grid(r1, r2, r3, r4, r5, ncol=1, rel_heights = c(1, 1, 1, 1, 1.3))
save_plot("moms_neutral.pdf", moms, base_height=10, base_width=30)


