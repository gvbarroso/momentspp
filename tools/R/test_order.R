#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

# space-separated positional arguments in Rscript --vanilla bgs_pi.R s no
s <- -1e-4#as.numeric(args[1]) 
no <- 30#as.numeric(args[2]) # number of orders tested

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(pracma) # for cubicspline()
  library(bigsnpr) # for seq_log
  library(cowplot)
  library(scales)
  library(RColorBrewer)
  library(shiny)
  library(bslib)
})

setwd("~/Data/momentspp/paper_1/orders/10k/")

params <- fread("params.csv")

r <- unique(params$lookup_r)
u <- unique(params$lookup_u)
o <- 1:no

n_models <- nrow(crossing(s, r, u, o))
hrs <- as.data.frame(matrix(nrow=nrow(params), ncol=6))
hls <- as.data.frame(matrix(nrow=nrow(params), ncol=6))

c <- 1
# loops are nested as in setup_models_1.Rmd and the directory structure
for(i in 1:length(r)) {
  for(j in 1:length(s)) {
    for(k in 1:length(u)) {  
      fnames <- list.files(paste("r_", r[i], "/s_", s[j], "/u_", u[k], sep=""),
                           pattern="*expectations.txt", full.names=TRUE)
      for(l in 1:length(fnames)) {
        order <- as.numeric(strsplit(strsplit(fnames[l], "/")[[1]][4], 
                                 "_")[[1]][1])
        idx <- (c - 1) * length(o) + l

        moms <- fread(fnames[l])
        hls[idx, 1] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
        hls[idx, 2] <- order
        hls[idx, 3] <- c
        hls[idx, 4] <- r[i]
        hls[idx, 5] <- s[j]
        hls[idx, 6] <- u[k]
    
        hrs[idx, 1] <- t(dplyr::select(filter(moms, V1=="Hr_0_0"), V3))
        hrs[idx, 2] <- order
        hrs[idx, 3] <- c
        hrs[idx, 4] <- r[i]
        hrs[idx, 5] <- s[j]
        hrs[idx, 6] <- u[k]
      }
      
      c <- c + 1
    }
  }
}

names(hrs) <-  c("Hr", "Order", "scenario", "r", "s", "u")
names(hls) <-  c("Hl", "Order", "scenario", "r", "s", "u")
        
hrs$N <- 1e+4
hls$N <- 1e+4
hrs$alpha <- 2 * hrs$N * hrs$s
hls$alpha <- 2 * hls$N * hls$s

hls <- hls %>% group_by(scenario) %>% 
  mutate(., Hl_norm=Hl/Hl[length(o)])
hls <- hls %>% group_by(scenario) %>%
  mutate(., Hl_diff=abs(Hl-Hl[length(o)]))

hrs$B <- hrs$Hr / (2 * hrs$N * hrs$u)
hrs <- hrs %>% group_by(scenario) %>% 
  mutate(., Hr_norm=Hr/Hr[length(o)])
hrs <- hrs %>% group_by(scenario) %>%
  mutate(., Hr_diff=abs(Hr-Hr[length(o)]))

hrs <- arrange(hrs, scenario, Order)
hls <- arrange(hls, scenario, Order)

fwrite(hrs, "order_N_1e+4_s_1e-4.csv")

# Hr
rvals <- c(unique(hrs$r)[1], unique(hrs$r)[17])
tbl_hr <- filter(hrs, r %in% c(rvals), u==1e-8)

p1 <- ggplot(data=tbl_hr, aes(x=as.factor(Order), y=Hr)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hr vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hr(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hr_s_", s, ".png", sep=""),
          p1, base_height=8, base_width=16)

p2 <- ggplot(data=tbl_hr, aes(x=as.factor(Order), y=Hr_norm)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hr_norm vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hr_norm(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hr_norm_s_", s, ".png", sep=""),
          p2, base_height=8, base_width=16)

tbl_hr <- filter(tbl_hr, Hr_diff != 0)

p3 <- ggplot(data=tbl_hr, aes(x=as.factor(Order), y=Hr_diff)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hr(O[x]) - Hr(O[max]) vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hr_diff(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hr_diff_s_", s, ".png", sep=""),
          p3, base_height=8, base_width=16)

p4 <- ggplot(data=tbl_hr, aes(x=as.factor(Order), y=B)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="B vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("B(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("B_s_", s, ".png", sep=""),
          p4, base_height=8, base_width=16)

# Hl
rvals <- c(unique(hls$r)[1], unique(hls$r)[17])
tbl_hl <- filter(hls, r %in% c(rvals), u==1e-8)

p5 <- ggplot(data=tbl_hl, aes(x=as.factor(Order), y=Hl)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hl_norm vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hl(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hl_s_", s, ".png", sep=""),
          p5, base_height=8, base_width=16)

p6 <- ggplot(data=tbl_hl, aes(x=as.factor(Order), y=Hl_norm)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hl_norm vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hl_norm(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hl_norm_s_", s, ".png", sep=""), 
          p6, base_height=8, base_width=16)

tbl_hl <- filter(tbl_hl, Hl_diff != 0)

p7 <- ggplot(data=tbl_hl, aes(x=as.factor(Order), y=Hl_diff)) + 
  geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10() +
  labs(title="Hl(O[x]) - Hl(O[max]) vs Order (1-2p), cols->r", 
       x="Order (1-2p)", y=paste("Hl_diff(s=", s, ")", sep="")) +
  scale_x_discrete(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

save_plot(paste("Hl_diff_s_", s, ".png", sep=""),
          p7, base_height=8, base_width=16)
