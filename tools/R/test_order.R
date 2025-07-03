

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

setwd("~/Data/momentspp/paper_1/orders/")

r <- 1e-8
u <- 1e-8
s <- c(-1e-4, -1e-3)
N <- c(1e+4, 1e+5)

params <- crossing(r, u, N, s)
params$model <- 1:nrow(params)
n_models <- nrow(params)

tbl_hrs <- data.table()
tbl_hls <- data.table()

# loops are nested as in setup_models_1.Rmd and the directory structure
for(i in 1:n_models) {
  fnames <- list.files(pattern=paste("model_", i, "_*.*_expectations.txt", sep=""),
                       full.names=TRUE)
  
  hls <- as.data.frame(matrix(nrow=length(fnames), ncol=3))
  hrs <- as.data.frame(matrix(nrow=length(fnames), ncol=3))
  
  for(l in 1:length(fnames)) {
    order <- as.numeric(strsplit(strsplit(strsplit(fnames[l], "/")[[1]][2], "_O_")[[1]][2], "_")[[1]][1])
    moms <- fread(fnames[l])
    hls[l, 1] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
    hls[l, 2] <- order
    hls[l, 3] <- i
    
    hrs[l, 1] <- t(dplyr::select(filter(moms, V1=="Hr_0_0"), V3))
    hrs[l, 2] <- order
    hrs[l, 3] <- i
  }
  
  tbl_hrs <- rbind.data.frame(tbl_hrs, hrs)
  tbl_hls <- rbind.data.frame(tbl_hls, hls)
}

names(tbl_hrs) <-  c("Hr", "Order", "model")
names(tbl_hls) <-  c("Hl", "Order", "model")
        
tbl_hrs <- arrange(tbl_hrs, model, Order)
tbl_hls <- arrange(tbl_hrs, model, Order)

df <- merge(tbl_hrs, params, by="model")
df$alpha <- 2 * df$N * df$s

df$Hr <- df$Hr * 1e+3
df$B <- df$Hr / (2 * df$N * df$u) 
df <- df %>% group_by(model) %>% 
             mutate(., diff=abs(Hr-Hr[length(.)]) / Hr[length(.)])
df$diff

fwrite(df, "tbl_order.csv")

p1 <- ggplot(data=filter(df, N==1e+4, s==-1e-4), 
             aes(x=as.factor(Order), y=diff)) + 
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title=NULL, x="Order (1-2p)", y="Rel. Diff.", sep="") +
  scale_x_discrete(breaks=seq(0, 30, 5)) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

p2 <- ggplot(data=filter(df, N==1e+4, s==-1e-3), 
             aes(x=as.factor(Order), y=diff)) + 
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title=NULL, x="Order (1-2p)", y=NULL, sep="") +
  scale_x_discrete(breaks=seq(10, 100, 10)) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

p3 <- ggplot(data=filter(df, N==1e+5, s==-1e-4),
             aes(x=as.factor(Order), y=diff)) + 
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title=NULL, x="Order (1-2p)", y="Rel. Diff.", sep="") +
  scale_x_discrete(breaks=seq(10, 100, 10)) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

p4 <- ggplot(data=filter(df, N==1e+5, s==-1e-3),
             aes(x=as.factor(Order), y=diff)) + 
  geom_point() + theme_bw() + 
  labs(title=NULL, x="Order (1-2p)", y=NULL, sep="") +
  scale_y_log10(breaks=c(1e-5, 1e-6, 1e-7, 1e-8)) +
  scale_x_discrete(breaks=seq(100, 500, 50)) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        strip.text=element_text(size=14))

p_order <- plot_grid(p1, p2, p3, p4, labels="AUTO", nrow=2)
save_plot("Hr_diff_.png", p_order,base_height=8, base_width=16)
