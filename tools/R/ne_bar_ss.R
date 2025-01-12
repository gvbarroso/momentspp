# Date created: 01/04/2024
# Last modified: 01/04/2024
# Author: Gustavo V. Barroso

#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)


pdf(NULL) # to suppress creation of Rplots.pdf

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(scales)
  library(shiny)
  library(bslib)
})

setwd("~/Data/momentspp/paper_1/4_epoch_model/")

params <- fread("params.csv")

Ns <- unique(params$Ns)
Ts <- unique(params$Ts)
uL <- unique(params$uL)
uR <- unique(params$uR)

ts <- strsplit(Ts, ";")[[1]]
t <- as.integer(strsplit(Ts, ";")[[1]][length(ts)])
  
demo_models <- crossing(Ns, Ts, uL, uR)
num_demo_models <- nrow(demo_models)
sampling_times <- seq(from=t, to=0, by=-500)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=length(sampling_times)))
names(pi0) <- as.character(sampling_times)

fnames <- list.files("demo", pattern="*_e_*.*_time.txt", full.names=TRUE)
hets <- NULL
for(l in 1:length(fnames)) {
  moms <- fread(fnames[l])
  str_sides <- strsplit(fnames[l], "_O_")[[1]]
  order <- as.numeric(strsplit(str_sides[2], "_e")[[1]][1])
  hets <- rbind.data.frame(hets, moms)
}

right <- filter(hets, V1=="Hr_0_0")
right <- right[!duplicated(right$V4)]
pi0[1,] <- t(dplyr::select(filter(right, V1=="Hr_0_0"), V3))

demo_pi0 <- cbind.data.frame(demo_models, pi0)

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(sampling_times),
                      names_to="Generation", values_to="pi0")
m_pi0$Generation <- as.numeric(m_pi0$Generation)
m_pi0$Ne_bar <- m_pi0$pi0 / 2 / m_pi0$uR
fwrite(m_pi0, "demo/m_pi0.csv")

p <- ggplot(data=m_pi0, aes(x=Generation, y=pi0)) + 
  geom_point() + geom_line() + theme_bw() + 
  labs(title=paste("Trajectory of Heterozygosity after size change"), 
       x="Generations ago", y=expression(pi[0])) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=14),
        legend.position="bottom")
  
save_plot("demo/pi0_demo.png", p, base_height=8, base_width=16)

fwrite(dplyr::select(m_pi0, Ne_bar), "demo/ne_bar_pi0.csv")

params_ss <- crossing(params, dplyr::select(m_pi0, c(Generation, Ne_bar, pi0)))
params_ss$Ts <- "0"
#params_ss$Ns <- paste(as.character(params_ss$Ne_bar), ";", sep="")
n_models <- nrow(params_ss)

fwrite(params_ss, "ne_bar_eq/params.csv")

options(digits=16) # precision
for(i in 1:n_models) {
  
  if(i %% 10000 == 0) {
    print(Sys.time())
    cat(paste(i, "\n"))
  }
  
  Ns <- as.numeric(params_ss[i,]$Ne_bar)
  Ts <- as.integer(params_ss[i,]$Ts)
  u <- as.numeric(params_ss[i,]$uL)
  r <- as.numeric(params_ss[i,]$r)
  s <- as.numeric(params_ss[i,]$s)
  
  name <- paste("model_", i, sep="")
  sink(paste("ne_bar_eq/r_", r, "/s_", s,"/", name, ".yaml", sep=""))
  
  cat("time_units: generations\n")
  cat("demes:\n")
  cat("  - name: X\n")
  cat("    epochs:\n")
  for(j in length(Ts):1) {
    cat("      - end_time: ")
    cat(Ts[j])
    cat("\n")
    cat("        start_size: ")
    cat(Ns[j])
    cat("\n")
  }
  cat("metadata:\n")
  cat("  - name: mutation\n")
  cat("    left_factor: 1\n")
  cat("    epochs:\n")
  cat("      - end_time: 0\n")
  cat("        rates: [")
  cat(u)
  cat("]\n")
  cat("  - name: recombination\n")
  cat("    epochs:\n")
  cat("      - end_time: 0\n")
  cat("        rates: [")
  cat(r)
  cat("]\n")
  cat("  - name: selection\n")
  cat("    start_time: .inf\n")
  cat("    epochs:\n")
  cat("      - end_time: 0\n")
  cat("        rates: [")
  cat(format(s, scientific=F))
  cat("]\n")
  
  sink()
}

