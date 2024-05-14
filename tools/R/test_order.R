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

setwd("~/Devel/momentspp/debug/")

s <- c(-1e-3, -1e-4, -1e-5)
r <- c(1e-8, 1e-4)
o <- 1:7

params <- crossing(r, s)
params$Model <- 1:nrow(params)
params <- crossing(params, o)

n_models <- nrow(crossing(s, r))
hrs <- as.data.frame(matrix(nrow=nrow(params), ncol=3))
hls <- as.data.frame(matrix(nrow=nrow(params), ncol=3))
for(i in 1:n_models) {
  fnames <- list.files(paste("model_", i, sep=""),
                       pattern="*expectations.txt", full.names=TRUE)
  for(l in 1:length(fnames)) {
    order <- as.numeric(strsplit(strsplit(fnames[l], "/")[[1]][2], 
                                 "_")[[1]][1])
    idx <- (i - 1) * length(o) + l
    cat(paste(idx, "\n"))
    moms <- fread(fnames[l])
    hls[idx, 1] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
    hls[idx, 2] <- order
    hls[idx, 3] <- i
    
    hrs[idx, 1] <- t(dplyr::select(filter(moms, V1=="Hr_0_0"), V3))
    hrs[idx, 2] <- order
    hrs[idx, 3] <- i
  }
}

names(hrs) <-  c("Hr", "Order", "Model")
names(hls) <-  c("Hl", "Order", "Model")
        
params_hl <- cbind.data.frame(dplyr::select(params, -o), dplyr::select(hls, -Model))
params_hr <- cbind.data.frame(dplyr::select(params, -o), dplyr::select(hrs, -Model))

params_hr <- arrange(params_hr, Model, Order)
params_hl <- arrange(params_hl, Model, Order)

params_hr <- params_hr %>% group_by(Model) %>% mutate(., Hr_norm=Hr/Hr[length(o)])
params_hr <- params_hr %>% group_by(Model) %>% mutate(., Hr_diff=abs(Hr-Hr[length(o)]))

params_hr <- filter(params_hr, s != -1e-5, Hr_diff != 0) # tmp
plot_list <- list(length=length(unique(params_hr$s)))
c <- 1

for(sel in unique(params_hr$s)) {
  p <- ggplot(data=filter(params_hr, s==sel),
              aes(x=as.factor(Order), y=Hr_diff)) + 
    geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10()
  if(c==1) {
    p <- p + labs(title="Hr(O[x]) - Hr(O[4Ns]) vs Order (1-2p), cols->r, rows->s", 
         x=NULL, y=paste("Hr(", sel, ")", sep="")) +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          strip.text=element_text(size=14))
  } else if(c==length(unique(params_hr$s))) {
    p <- p + labs(title=NULL, x="Order (1-2p)", y=paste("Hr(", sel, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            strip.text=element_blank())
  } else {
    p <- p + labs(title=NULL, x=NULL, y=paste("Hr(", sel, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            strip.text=element_blank())
  }
  plot_list[[c]] <- p
  c <- c + 1
}

x <- plot_grid(plotlist=plot_list, ncol=1)

plot_list <- list(length=length(unique(params_hr$s)))
c <- 1

for(sel in unique(params_hr$s)) {
  p <- ggplot(data=filter(params_hr, s==sel),
              aes(x=as.factor(Order), y=Hr_norm)) + 
    geom_point() + theme_bw() + facet_wrap(~r) + scale_y_log10()
  if(c==1) {
    p <- p + labs(title="Normalized Hr vs Order (1-2p), cols->r, rows->s", 
                  x=NULL, y=paste("Hr(", sel, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            strip.text=element_text(size=14))
  } else if(c==length(unique(params_hr$s))) {
    p <- p + labs(title=NULL, x="Order (1-2p)", y=paste("Hr(", sel, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            strip.text=element_blank())
  } else {
    p <- p + labs(title=NULL, x=NULL, y=paste("Hr(", sel, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            strip.text=element_blank())
  }
  plot_list[[c]] <- p
  c <- c + 1
}

y <- plot_grid(plotlist=plot_list, ncol=1)

