# Date created: 01/04/2024
# Last modified: 01/04/2024
# Author: Gustavo V. Barroso

#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

sel <- as.numeric(args[1])
num_iter <- as.numeric(args[2]) # for updating B-vals due to interference

pdf(NULL) # to suppress creation of Rplots.pdf

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

setwd("~/Data/momentspp/paper_1/study")

###############################
#
# Set up
#
###############################

params <- fread("params.csv")

Na <- unique(params$Na)
N1 <- unique(params$N1)
t <- unique(params$t)
uL <- unique(params$uL)

demo_models <- crossing(Na, N1, t, uL)
num_demo_models <- nrow(demo_models)
sampling_times <- seq(from=0, to=t, by=250)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=length(sampling_times)))
names(pi0) <- as.character(sampling_times)

for(i in 1:num_demo_models) {
  
  moms <- read.csv(paste("demo/model_", i, 
                      "_O_2_e_1_expectations.txt", sep=""), sep=" ", header=F)
  
  pi0[i,] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
}

demo_pi0 <- cbind.data.frame(demo_models, pi0)
fwrite(demo_pi0, "demo/demo_pi0.csv")

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(sampling_times),
                      names_to="Generation", values_to="pi0")
m_pi0$Generation <- as.numeric(m_pi0$Generation)

N.labs <- paste("Model ", 1:num_demo_models)
names(N.labs) <- rep(as.character(unique(demo_models$N1)),
                     length(unique(demo_models$uL)))

c <- 1
for(mu in uL) {
  p <- ggplot(data=filter(m_pi0, uL==mu), aes(x=Generation, y=pi0)) + 
    facet_wrap(~N1, labeller=labeller(N1=N.labs[c:(c+1)])) +
    geom_point(size=0.5) + geom_line() + theme_bw() + 
    labs(title=paste("Trajectory of Heterozygosity after size change, u =", mu), 
         x="Generations ago", y=expression(pi[0])) +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          strip.text=element_text(size=14),
          legend.position="bottom")
  
  save_plot(paste("demo/pi0_demo_u_", mu, ".png", sep=""), p, 
            base_height=8, base_width=16)
  
  c <- c + 1
}

d <- ggplot(data=filter(demo_models, uL==uL[1])) + theme_bw() + 
     facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
     geom_segment(aes(x=0, xend=t, y=N1, yend=N1), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t, y=N1, yend=Na), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t+t/10, y=Na, yend=Na), linewidth=1.5) +
     scale_y_log10() +
     scale_x_continuous(breaks=c(0, 1000, t)) +
     labs(title="Demographic Models", 
          x="Generations ago", 
          y="Population Size (N)") +
     theme(axis.title=element_text(size=16), 
           axis.text=element_text(size=12), 
           axis.text.x=element_text(size=12),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           strip.text=element_text(size=14),
           legend.position="bottom")

save_plot("demo/demo_models.png", d, base_height=8, base_width=16)

look_s <- setDT(as.data.frame(unique(params$s)))
look_r <- setDT(as.data.frame(unique(params$r)))

names(look_r) <- "r"
names(look_s) <- "s"

setkey(look_r, r)
setkey(look_s, s)

hl_demo <- fread("hl_time.csv.gz")
hr_demo <- fread("hr_time.csv.gz")

###############################
#
# Single constrained locus
#
###############################

m_het_demo <- pivot_longer(hr_demo, cols=as.character(sampling_times),
                          names_to="Generation", values_to="Hr")
m_het_demo$Generation <- as.numeric(m_het_demo$Generation)
m_het_demo <- setDT(m_het_demo)
setkey(m_het_demo, N1, uL, r, s, Generation)

m_het_demo <- left_join(m_het_demo, dplyr::select(m_pi0,
                        -c(Na, t)), by=c("N1", "Generation", "uL"))
m_het_demo$B <- (m_het_demo$Hr / m_het_demo$pi0)
m_het_demo$uR <- m_het_demo$uL
  
hl_gen <- pivot_longer(hl_demo, cols=as.character(sampling_times),
                       names_to="Generation", values_to="Hl")
hl_gen$Generation <- as.numeric(hl_gen$Generation)
hl_gen <- setDT(hl_gen)
hl_gen <- arrange(hl_gen, Generation)
setkey(hl_gen, N1, uL, r, s, Generation)
hl_gen$uR <- hl_gen$uL

m_het_demo$Hl <- hl_gen$Hl
m_het_demo$piN_pi0 <- m_het_demo$Hl / m_het_demo$pi0
m_het_demo$piN_piS <- m_het_demo$Hl / m_het_demo$Hr

fwrite(m_het_demo, "mpp_stats_demo.csv.gz")

samp_gens <- c(500, 1000, 2500, 5000, 10000, 20000, 100000)
tbl <- as.data.frame(matrix(ncol=3, nrow=length(samp_gens)))
c <- 1
for(g in samp_gens) {
  tmp <- filter(m_het_demo, Generation==g)
  tbl[c,] <- c(unique(tmp$Ne_bar), g)
  c <- c + 1
}
names(tbl) <- c("Ne_bottleneck", "Ne_expansion", "Gen")
fwrite(tbl, "Ne_tbl.csv")

m_demo <- pivot_longer(m_het_demo, cols=c(pi0, Hr, Hl, B, piN_pi0, piN_piS), 
                       names_to="statistic")

svals <-  c(unique(m_demo$s)[1],
            unique(m_demo$s)[10], 
            unique(m_demo$s)[20])

rvals <-  c(unique(m_demo$r)[2],
            unique(m_demo$r)[7], 
            unique(m_demo$r)[12])

for(rec in rvals) {
  for(mom in unique(m_demo$statistic)) {
    
    plot_list <- list(length=length(svals))
    c <- 1
    
    for(sel in svals) {
      p <- ggplot(data=filter(m_demo, statistic==mom, 
                              uL==1e-8, s==sel, r==rec),
                  aes(x=Generation, y=value)) + 
        geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(N1)) +
        scale_x_continuous(breaks=pretty_breaks()) +
        scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none")
      if(c==length(svals)) {
        p <- p + labs(title=NULL, x="Generation", 
                      y=paste(mom, "(s=", sel, ")", sep="")) +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_text(size=12),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="none")
      } else if(c==1) {
        p <- p + labs(title=paste("Temporal dynamics after a size change, u=",
                                  1e-8, "_r=", rec, sep=""),
                      x=NULL, y=paste(mom, "(s=", sel, ")", sep="")) +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_blank(),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="none")
      } else {
        p <- p + labs(title=NULL, x=NULL, y=paste(mom, "(s=", sel, ")", sep="")) +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_blank(),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="none")
      }
      
      plot_list[[c]] <- p
      c <- c + 1
    }
    
    mom_plot <- plot_grid(plotlist=plot_list, ncol=1, align='v')
    save_plot(paste(mom, "_time_u_1e-8_r_", rec, ".png", sep=""),
              mom_plot, base_height=10, base_width=12)
  }
}

r <- ggplot(data=filter(m_demo, Generation==2.5e+5, 
                        s %in% c(unique(m_demo$s)[1],
                                 unique(m_demo$s)[2],
                                 unique(m_demo$s)[3]), 
                                 #unique(m_demo$s)[15],
                                 #unique(m_demo$s)[20]),
                        uL %in% c(unique(m_demo$uL)[1],
                                  unique(m_demo$uL)[5]), 
                        N1==1e+3, statistic=="B"),
       aes(x=r, y=value)) + 
  geom_point() + theme_bw() + geom_line() + 
  facet_grid(as.factor(uL)~as.factor(s)) +
  scale_x_log10(breaks=unique(m_demo$r)) +
  scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none") +
  labs(title="B~r, cols->s, rows->mu", x="r", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("B_r.png", r, base_height=8, base_width=16)

s <- ggplot(data=filter(m_demo, N1==1e+3, Generation==2.5e+5, statistic=="B",
                        uL %in% c(unique(m_demo$uL)[1],
                                  unique(m_demo$uL)[15]),
                        r %in% c(unique(m_demo$r)[1],
                                        unique(m_demo$r)[2], 
                                        unique(m_demo$r)[3],
                                        unique(m_demo$r)[15])),
            aes(x=-s, y=value)) + 
  facet_grid(as.factor(uL)~as.factor(r)) +
  geom_point() + theme_bw() + geom_line() + 
  scale_x_log10(breaks=c(1e-5, 1e-4, 1e-3)) +
  scale_y_log10(breaks=pretty_breaks()) + 
  labs(title=NULL, x="-s", y="B-value") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("B_s.png", s, base_height=10, base_width=16)

s2 <- ggplot(data=filter(m_demo, N1==1e+3, Generation==2.5e+5, 
                        statistic==c("piN_pi0", "piN_piS"),
                        uL %in% c(unique(m_demo$uL)[1],
                                  unique(m_demo$uL)[15]),
                        r %in% c(unique(m_demo$r)[1],
                                 unique(m_demo$r)[5], 
                                 unique(m_demo$r)[10],
                                 unique(m_demo$r)[15])),
            aes(x=-s, y=value, color=statistic)) + 
  facet_grid(as.factor(uL)~as.factor(r)) +
  geom_point() + theme_bw() + geom_line() + 
  scale_x_log10(breaks=c(1e-5, 1e-4, 1e-3)) +
  scale_y_log10(breaks=pretty_breaks()) + 
  labs(title=NULL, x="-s", y="B-value") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("piN_piX_s.png", s2, base_height=10, base_width=16)

d_stats <- filter(m_het_demo, uL==1e-8) %>% 
     group_by(scenario) %>% 
     reframe(N1=N1, lookup_r=lookup_r, lookup_s=lookup_s, Generation=Generation,
             d_pi0=-c(diff(pi0)[1], diff(pi0))/pi0,
             d_Hl=-c(diff(Hl)[1], diff(Hl))/Hl,
             d_Hr=-c(diff(scaled_Hr)[1], diff(scaled_Hr))/scaled_Hr,
             d_B=-c(diff(B)[1], diff(B))/B,
             d_piN_pi0=-c(diff(piN_pi0)[1], diff(piN_pi0))/piN_pi0,
             d_piN_piS=-c(diff(piN_piS)[1], diff(piN_piS))/piN_piS)

m_demo_d <- pivot_longer(d_stats, cols=starts_with("d_"), names_to="statistic")

p <- ggplot(data=filter(m_demo, N1==1e+5, Generation==1e+5, statistic=="B",
                        lookup_r %in% rvals,
                        lookup_s %in% svals),
            aes(x=lookup_u, y=value)) + 
  geom_point() + theme_bw() + geom_line() + 
  facet_grid(as.factor(lookup_r)~as.factor(lookup_s)) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none") +
  labs(title="B~u, cols=s, rows=r", x="u", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
 
for(r in rvals) {
  c <- 1
  plot_list <- list(length=length(unique(m_demo_d$N1)))
  
  for(N in unique(m_demo_d$N1)) {
    p <- ggplot(data=filter(m_demo_d, N1==N, lookup_r==r,
                            lookup_s %in% svals, 
                            statistic %in% c("d_Hl", "d_Hr", "d_pi0")),
                aes(x=Generation, y=value, color=statistic)) + 
      geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(lookup_s)) +
      scale_x_continuous(breaks=pretty_breaks()) +
      scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none")
      if(c==length(unique(m_demo_d$N1))) {
        p <- p + labs(title=NULL, x="Generation", y="Rate of Change") +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_text(size=12),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="bottom")
      } else if(c==1) {
        p <- p + labs(title="Temporal dynamics of statistics after a size change",
                      x=NULL, y="Rate of change") +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_blank(),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="none")
      } else {
        p <- p + labs(title=NULL, x=NULL, y="Rate of change") +
          theme(axis.title=element_text(size=16), 
                axis.text=element_text(size=12), 
                axis.text.x=element_blank(),
                legend.text=element_text(size=16),
                legend.title=element_text(size=16),
                legend.position="none")
      }
    
    plot_list[[c]] <- p
    c <- c + 1
  }
  
  cp <- plot_grid(plotlist=plot_list, ncol=1, align='v')
  save_plot(paste("diffs_stats_time_u_1e-8_r_", r, ".png", sep=""),
            cp, base_height=10, base_width=12)
}

##############################
#
# multiple constrained loci
#
###############################

# reset
m_het_demo$B <- (m_het_demo$Hr / m_het_demo$pi0) ^ 1 
#m_demo <- pivot_longer(m_het_demo, cols=c(
#                       Hl, B, piN_pi0, piN_piS), names_to="statistic")

# interpolates between r values
Bs_demo <- pivot_longer(m_het_demo, cols=B, names_to="statistic") %>%
           filter(., Generation %in% seq(from=0, to=t, by=1e+3)) # speed?
interp <- list(length=length(unique(Bs_demo$lookup_u)) * 
                      length(unique(Bs_demo$lookup_s)) * 
                      length(unique(Bs_demo$N1)) *
                      length(unique(Bs_demo$Generation)))
c <- 1
d <- 1
pb <- txtProgressBar(min=0, max=length(unique(Bs_demo$lookup_u)), style=3)
for(u in unique(Bs_demo$lookup_u)) {
  setTxtProgressBar(pb, d)
  for(N in unique(Bs_demo$N1)) {
    for(s in unique(Bs_demo$lookup_s)) {
      for(gen in unique(Bs_demo$Generation)) {
        tmp <- filter(Bs_demo, lookup_u==u, lookup_s==s, Generation==gen, N1==N)
        rs <- sort(c(seq_log(from=1e-8, to=1e-3, length.out=200),
                     seq(from=1e-4, to=1e-2, length.out=300)))
        Bs_interp <- cubicspline(tmp$lookup_r, tmp$value, rs)
        tmp <- cbind.data.frame(rs, Bs_interp, s, gen, u)
        names(tmp) <- c("r", "B", "s", "Generation", "u")
        tmp$N1 <- N
        
        interp[[c]] <- tmp
        c <- c + 1
      }
    }
  }
  d <- d + 1
}
close(pb)

B_demo_dt <- do.call("rbind", interp)
B_demo_dt <- setDT(B_demo_dt)
setkey(B_demo_dt, u, r, s)
fwrite(B_demo_dt, "B_demo_dt.csv.gz")

# split for faster lookup
dt_u <- data.table(unique(B_demo_dt$u))
dt_r <- data.table(unique(B_demo_dt$r))
dt_s <- data.table(unique(B_demo_dt$s))
names(dt_u) <- "u"
names(dt_r) <- "r"
names(dt_s) <- "s"
setkey(dt_u, u)
setkey(dt_r, r)
setkey(dt_s, s)

# set up chr
num_exons <- 100
ncsl <- rep(0, num_exons)
exon_lengths <- 1000
csl <- rep(exon_lengths, num_exons) 
L <- sum(csl) + sum(ncsl)

maps <- suppressWarnings(c(rbind(csl, ncsl)))
maps <- suppressMessages(setDT(bind_cols(chr="chr1",
    dplyr::lag(cumsum(maps), n=1, default=0),
    dplyr::lag(cumsum(maps), n=0, default=0))))
names(maps) <- c("chr", "start", "end")
maps$midpoint <- maps$start + (maps$end -  maps$start) / 2
maps$s <- c(rbind(rep(sel, length(ncsl)), 0)) 
maps$u <- 1e-8
maps$r <- 1e-8
setkey(maps, start, end)

dt_exons <- filter(maps, s!=0)

jump_length <- 25
neut_pos <- as.numeric(unlist(apply(filter(maps, s==0), 1, 
            function(x) seq(from=x[2], to=x[3], by=jump_length))))
#neut_pos <- L/2 + 1
samp_pos <- sort(c(neut_pos, dt_exons$midpoint)) # sampled sites

focal_mu <- unlist(lapply(samp_pos, function(pos) mean(maps[J(pos), roll=T]$u)))

pos_dt <- suppressMessages(setDT(bind_cols(samp_pos, focal_mu)))
names(pos_dt) <- c("position", "focal_mu")
pos_dt$idx <- 1:nrow(pos_dt) # indexing
setkey(pos_dt, position)

getB <- function(exon_id, samp_id, Ne_f, look_tbl) {
  exon_pos <- dt_exons$midpoint[exon_id]
  idx <- pos_dt[.(exon_pos)]$idx
  focal_s <- dt_exons[exon_id,]$s * B_values[idx] * Ne_f
  focal_mu <- dt_exons[exon_id,]$u * B_values[idx] * Ne_f
  total_r <- abs(pos_dt[samp_id]$cumrec - pos_dt[idx]$cumrec) * Ne_f
  
  closest_r <- dt_r[dt_r[.(total_r), roll="nearest", which=T]]
  
  # finding lower and higher bounds for u (x-axis) and s (y-axis)
  uL <- dt_u[dt_u[.(focal_mu), roll="nearest", which=T]]
  if(uL > focal_mu) {
    uL <- dt_u[dt_u[.(focal_mu), roll="nearest", which=T]-1]
    uH <- dt_u[dt_u[.(focal_mu), roll="nearest", which=T]]
  } else {
    uH <- dt_u[dt_u[.(focal_mu), roll="nearest", which=T] + 1]
  }
  
  sL <- dt_s[dt_s[.(focal_s), roll="nearest", which=T]]
  if(sL > focal_s) {
    sL <- dt_s[dt_s[.(focal_s), roll="nearest", which=T]-1]
    sH <- dt_s[dt_s[.(focal_s), roll="nearest", which=T]]
  } else {
    sH <- dt_s[dt_s[.(focal_s), roll="nearest", which=T] + 1]
  }

  # bilinear interpolation on the square (uL, sL), (uL, sH), (uH, sH), (uL, sL)
  B_ll <- look_tbl[.(uL, closest_r, sL), .(B)]
  B_lh <- look_tbl[.(uL, closest_r, sH), .(B)]
  B_hl <- look_tbl[.(uH, closest_r, sL), .(B)]
  B_hh <- look_tbl[.(uH, closest_r, sH), .(B)]
  
  B <- (B_ll * (uH-focal_mu) * (sH-focal_s) +
        B_hl * (focal_mu-uL)*(sH-focal_s) + 
        B_lh * (uH-focal_mu) * (focal_s-sL) +
        B_hh * (focal_mu-uL) * (focal_s-sL)) / ((uH-uL)*(sH-sL))
  
  return(B ^ exon_lengths)
}

options(datatable.verbose=F)
bmaps_demo <- list(length=length(unique(B_demo_dt$N1)))
mu <- unique(maps$u)
for(n in 1:length(unique(B_demo_dt$N1))) {
  
  N <- unique(B_demo_dt$N1)[n]
  B_maps <- as.data.frame(matrix(ncol=nrow(pos_dt),
      nrow=length(unique(B_demo_dt$Generation))))

  g <- 1
  for(gen in sort(unique(B_demo_dt$Generation), decreasing=T)) {
    
    print(Sys.time())
    cat(paste("Generation", gen, "\n"))
    
    tmp_lookup <- filter(B_demo_dt, Generation==gen, N1==N)
    Ne_bar <- unique(filter(m_het_demo,
                     Generation==gen, N1==N, lookup_u==mu)$pi0) / (2 * mu)
    Ne_ratio <- Ne_bar / 1e+4 # ratio relative to Nanc
    B_values <- rep(1, nrow(pos_dt)) # this generation's B-values 
    gen_g_iters <- as.data.frame(matrix(ncol=nrow(pos_dt), nrow=num_iter))
    
    for(i in 1:num_iter) {
      # approximates cumulative rec by looking at r at sampled sites only
      focal_r <- unlist(lapply(samp_pos,
          function(pos) mean(maps[J(pos), roll=T]$r))) * B_values
      cum_rec <- numeric(length=length(samp_pos))
      cum_rec[1] <- focal_r[1]
      
      for(j in 2:length(samp_pos)) {
        cum_rec[j] <- cum_rec[j-1] + (samp_pos[j] - samp_pos[j-1]) * focal_r[j]
      }
      pos_dt$cumrec <- cum_rec
      
      updated_B_vals <- rep(1, length(B_values))
      pb <- txtProgressBar(min=0, max=nrow(pos_dt), style=3)
      for(j in 1:nrow(pos_dt)) { 
        setTxtProgressBar(pb, j)
        Bs <- sapply(1:nrow(dt_exons), getB, samp_id=j,
                                             Ne_f=Ne_ratio,
                                             look_tbl=tmp_lookup)
        updated_B_vals[j] <- cumprod(Bs)[length(Bs)]
      }
      close(pb)
      
      B_values <- updated_B_vals
      gen_g_iters[i,] <- B_values
    }
    
    B_maps[g, ] <- B_values
    
    names(gen_g_iters) <- pos_dt$position
    gen_g_iters$Iter <- 1:num_iter
    gen_g_iters$s <- sel
    gen_g_iters$N1 <- N
    
    fwrite(gen_g_iters, paste("gen_", gen, 
                              "_iterations_s_", sel, ".csv", sep=""))
    
    m_iters <- pivot_longer(gen_g_iters, 
                            cols=as.character(pos_dt$position),
                            names_to="Pos")
    m_iters$Pos <- as.numeric(m_iters$Pos)
    p <- ggplot(data=m_iters,aes(x=Pos, y=value, color=as.factor(Iter))) + 
      geom_point(aes(alpha=0.5)) + theme_bw() + geom_line() + 
      scale_color_brewer(palette="YlOrRd", name="Iteration") +
      scale_x_continuous(breaks=pretty_breaks()) +
      scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none") +
      labs(title=NULL, x="Pos", y=paste("B (gen ", gen, ")", sep="")) +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            axis.text.x=element_text(size=12),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16),
            legend.position="bottom")
    
    save_plot(paste("iter_Nebar_g_", gen, "_N1_", N, "_s_", sel, ".png", sep=""),
              p, base_width=12, base_height=10)
    
    g <- g + 1
  }

  names(B_maps) <- samp_pos
  B_maps$Generation <- sort(unique(B_demo_dt$Generation), decreasing=T)
  B_maps$N1 <- N
  B_maps$s <- sel
  B_maps$u <- unique(maps$u)
    
  bmaps_demo[[n]] <- B_maps
}

bmaps_demo <- do.call("rbind", bmaps_demo)
m_Bmap_demo <- pivot_longer(bmaps_demo, 
                            cols=as.character(samp_pos), 
                            names_to="Position")
m_Bmap_demo$Position <- as.numeric(m_Bmap_demo$Position)

fwrite(bmaps_demo, paste("bmaps_demo_s_", sel, ".csv", sep=""))
fwrite(m_Bmap_demo, paste("m_Bmap_demo_s_", sel, ".csv", sep=""))

ui <- page_sidebar(
  title="B-value maps over time",
  sidebar=sidebar(
    sliderInput(inputId="Generation",
    label="Generation:",
    min=0,
    max=50000,
    step=10000,
    value=50000),
    checkboxInput("by_demography", "Show Ne", T),
    checkboxInput("by_selection", "Show s", T),
    checkboxInput("by_mutation", "Show u", T),
    checkboxGroupInput(
      "N1", "Ne (color)",
      choices=unique(m_Bmap_demo$N1), 
      selected=unique(m_Bmap_demo$N1[1])
    ),
    checkboxGroupInput(
      "s", "s (shape)",
      choices=unique(m_Bmap_demo$s), 
      selected=unique(m_Bmap_demo$s[1])
    ),
    checkboxGroupInput(
      "u", "u",
      choices=unique(m_Bmap_demo$u), 
      selected=unique(m_Bmap_demo$u[1])
    )
  ),
  plotOutput(outputId="Bmap")
)

server <- function(input, output) {
  output$Bmap <- renderPlot({
    p <- ggplot(data=filter(m_Bmap_demo, 
                            N1==input$N1,
                            s==input$s,
                            u==input$u,
                            Generation==input$Generation),
                aes(x=Position, y=value)) + 
      list(
        geom_point(), geom_line(), theme_bw(),
        scale_y_continuous(breaks=pretty_breaks()),
        if(input$by_demography) aes(color=as.factor(N1)),
        if(input$by_selection) aes(shape=as.factor(s)),
        labs(title=NULL, x="Position", y="B"),
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
      )
    
    p
  }, res = 100)
}

#shinyApp(ui, server)

c <- ggplot(data=m_Bmap_demo, aes(x=Position, y=value, color=Generation)) +
  geom_point() + facet_wrap(~N1) + theme_bw() +
  scale_color_continuous(name="Generation", labels=NULL) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="B-value maps over time", x="Position", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot(paste("Bmap_time_s_", sel, ".png", sep=""), c, 
          base_height=10, base_width=16)

d <- ggplot(data=filter(m_Bmap_demo, Position %in% c(1e+3, 1e+4, 2.5e+4, 5e+4)),
            aes(x=Generation, y=value, color=as.factor(Position))) +
  geom_point() + theme_bw() + facet_wrap(~N1) + geom_line() +
  scale_color_discrete(name="Position") +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="B-values over time", x="Generation", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot(paste("Bpos_time_s_", sel, ".png", sep=""),
          d, base_height=8, base_width=12)

# adds simulation results
sim_exp_strong <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.001.txt")
sim_exp_moderate <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.0001.txt")
sim_exp_weak <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.00001.txt")

sim_exp_strong <- dplyr::select(sim_exp_strong, -V1)
sim_exp_moderate <- dplyr::select(sim_exp_moderate, -V1)
sim_exp_weak <- dplyr::select(sim_exp_weak, -V1)

names(sim_exp_strong) <- as.character(neut_pos)
names(sim_exp_moderate) <- as.character(neut_pos)
names(sim_exp_weak) <- as.character(neut_pos)

sim_exp_tbl <- rbind.data.frame(sim_exp_strong, sim_exp_moderate, sim_exp_weak)
sim_exp_tbl$Generation <- rep(seq(from=0, to=50000, by=1000), 3)
sim_exp_tbl$N1 <- 1e+5
sim_exp_tbl$s <- c(rep(-1e-3, 51), rep(-1e-4, 51), rep(-1e-5, 51))
sim_exp_tbl$u <- 1e-8

# down-sample generations in sim results
sim_exp_tbl <- filter(sim_exp_tbl, Generation %in% unique(bmaps_demo$Generation))

# down-sample sites in m++ results
bmaps_demo_ds <- dplyr::select(bmaps_demo, c(names(sim_exp_tbl)))

sim_exp_tbl$method <- "fwdpy11"
bmaps_demo_ds$method <- "mpp"

df <- rbind.data.frame(filter(bmaps_demo_ds, N1==1e+5), sim_exp_tbl)
m_df <- pivot_longer(df, cols=as.character(neut_pos), names_to="Position")
m_df$Position <- as.numeric(m_df$Position)

sampled_gens <- seq(from=50000, to=0, by=-10000) # downsample gens for speed
plot_list <- list(length=length(sampled_gens))
c <- 1
for(gen in sampled_gens) {
  e <- ggplot(data=filter(m_df, s==sel, Generation==gen), 
              aes(x=Position/1e+3, y=value, color=method)) +
    geom_point() + theme_bw() + 
    scale_y_continuous(breaks=pretty_breaks())
    if(c==length(sampled_gens)) {
      e <- e + labs(title=NULL, x="Position (kb)", 
                    y=paste("B(gen=", gen, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="bottom")
    } else if(c==1) {
      e <- e + labs(title=paste("B-maps over generations"),
                    x=NULL, y=paste("B(gen=", gen, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else {
      e <- e + labs(title=NULL, x=NULL, y=paste("B(gen=", gen, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    }
  
  plot_list[[c]] <- e
  c <- c + 1
}

f <- plot_grid(plotlist=plot_list, ncol=1, align='v',
               rel_heights=c(1.15, 1, 1, 1, 1, 1.65))
save_plot(paste("Bmap_mpp_sims_1_s_", sel, ".png", sep=""),
          f, base_height=16, base_width=14)

g <- ggplot(data=filter(m_df, s==sel, Position %in% c(1e+3, 1e+4, 2.5e+4, 5e+4)),
        aes(x=Generation, y=value, color=as.factor(Position), shape=method)) +
  geom_point(size=2) + theme_bw() + geom_line() +
  scale_color_discrete(name="Position") +
  scale_shape_manual(values=c(0, 1)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="B-value maps over time", x="Generation", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot(paste("Bmap_mpp_sims_2_s_", sel, ".png", sep=""),
          g, base_height=10, base_width=16)

q <- ggplot(data=filter(filter(m_Bmap_demo, Generation==5e+4), s==-1e-3),
       aes(x=Position, y=value, shape=as.factor(N1))) +
  geom_point(size=2) + theme_bw() + geom_line() +
  scale_shape_manual(values=c(0, 1), name="N1") +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Position", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot(paste("Bmap_mpp_demo_s_", sel, ".png", sep=""),
          q, base_height=10, base_width=16)
