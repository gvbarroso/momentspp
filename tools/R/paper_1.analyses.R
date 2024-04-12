# Date created: 01/04/2024
# Last modified: 01/04/2024
# Author: Gustavo V. Barroso

pdf(NULL) # to suppress creation of Rplots.pdf

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(pracma) # for cubicspline()
  library(cowplot)
  library(scales)
  library(RColorBrewer)
})

setwd("~/Data/mpp_data/paper_1/single_neutral")

#########################
#
# Set up
#
#########################

params <- fread("params.csv")

Na <- unique(params$Na)
N1 <- unique(params$N1)
t <- unique(params$t)

demo_models <- crossing(Na, N1, t)
num_demo_models <- nrow(demo_models)
sampling_times <- seq(from=t, to=0, by=-100)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=length(sampling_times)))
names(pi0) <- as.character(sampling_times)

for(i in 1:num_demo_models) {
  
  moms <- read.csv(paste("demo/model_", i, 
                      "_e_1_expectations.txt", sep=""), sep=" ", header=F)
  
  pi0[i,] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
}

demo_pi0 <- cbind.data.frame(demo_models, pi0)
fwrite(demo_pi0, "demo/demo_pi0.csv")

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(sampling_times),
                      names_to="Generation", values_to="pi0")
m_pi0$Generation <- as.numeric(m_pi0$Generation)

N.labs <- paste("Model ", 1:num_demo_models)
names(N.labs) <- as.character(unique(demo_models$N1))

p <- ggplot(data=m_pi0, aes(x=Generation, y=pi0)) + 
  facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
  geom_point(size=0.5) + geom_line() + theme_bw() + 
  labs(title="Trajectory of Heterozygosity after size change", 
       x="Generations ago", y=expression(pi[0])) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=14),
        legend.position="bottom")

save_plot("demo/pi0_demo.png", p, base_height=8, base_width=16)

d <- ggplot(data=demo_models) + theme_bw() + 
     facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
     geom_segment(aes(x=0, xend=t, y=N1, yend=N1), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t, y=N1, yend=Na), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t+5000, y=Na, yend=Na), linewidth=1.5) +
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

# loads tables from setup_models_1.Rmd
lookup_tbl <- fread("lookup_tbl.csv.gz")
setkey(lookup_tbl, lookup_r, lookup_s)

look_s <- setDT(as.data.frame(unique(params$lookup_s)))
look_r <- setDT(as.data.frame(unique(params$lookup_r)))

names(look_r) <- "lookup_r"
names(look_s) <- "lookup_s"

setkey(look_r, lookup_r)
setkey(look_s, lookup_s)

hl_demo <- fread("hl_time.csv.gz")
hr_demo <- fread("hr_time.csv.gz")

plots <- list(length=length(unique(lookup_tbl$N1)))
for(i in 1:length(unique(lookup_tbl$N1))) {
  plots[[i]] <- ggplot(data=filter(lookup_tbl, lookup_r==1e-6,
                                   N1==unique(lookup_tbl$N1)[i]), 
                                   aes(x=-lookup_s, y=hr)) + theme_bw() +
    geom_point() + geom_line() + scale_x_log10() + scale_y_continuous() +
    labs(title=NULL, x="s",  y="Hr") +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          strip.text=element_text(size=14),
          legend.position="bottom")
}

p <- plot_grid(plotlist=plots)
save_plot("hr_s_n.png", p, base_height=8, base_width=16)

###############################
#
# single constrained locus
#
###############################

m_het_demo <- pivot_longer(hr_demo, cols=as.character(sampling_times),
                          names_to="Generation", values_to="Hr")
m_het_demo$Generation <- as.numeric(m_het_demo$Generation)
m_het_demo <- setDT(m_het_demo)
setkey(m_het_demo, N1, lookup_r, lookup_s, Generation)

m_het_demo <- left_join(m_het_demo, dplyr::select(m_pi0, -c(Na, t)),
                       by=c("N1", "Generation"))
m_het_demo$B <- (m_het_demo$Hr / m_het_demo$pi0) ^ 1000 # exon has L sites
m_het_demo$scaled_Hr <- m_het_demo$pi0 * m_het_demo$B
  
hl_gen <- pivot_longer(hl_demo, cols=as.character(sampling_times),
                       names_to="Generation", values_to="Hl")
hl_gen$Generation <- as.numeric(hl_gen$Generation)
hl_gen <- setDT(hl_gen)
hl_gen <- arrange(hl_gen, scenario, Generation)
setkey(hl_gen, N1, lookup_r, lookup_s, Generation)

m_het_demo$Hl <- hl_gen$Hl
m_het_demo$piN_pi0 <- m_het_demo$Hl / m_het_demo$pi0
m_het_demo$piN_piS <- m_het_demo$Hl / m_het_demo$scaled_Hr

fwrite(m_het_demo, "stats_demo.csv")

m_demo <- pivot_longer(m_het_demo, cols=c(pi0, Hr, scaled_Hr, Hl, B, piN_pi0, piN_piS), 
                       names_to="statistic")

svals <- c(-1e-3, -1e-4, -1e-5) # unique(m_tbl$lookup_s), downsampling now
for(mom in unique(m_demo$statistic)) {
  
  plot_list <- list(length=length(svals))
  c <- 1
  
  for(s in svals) {
    p <- ggplot(data=filter(m_demo, statistic==mom, lookup_s==s, lookup_r==1e-10),
                aes(x=Generation, y=value)) + 
      geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(N1)) +
      scale_x_continuous(breaks=pretty_breaks()) +
      scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none")
    if(c==length(svals)) {
      p <- p + labs(title=NULL, x="Generation", y=paste(mom, "(s=", s, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else if(c==1) {
      p <- p + labs(title=paste("Temporal dynamics of", mom, "after a size change"),
                    x=NULL, y=paste(mom, "(s=", s, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else {
      p <- p + labs(title=NULL, x=NULL, y=paste(mom, "(s=", s, ")", sep="")) +
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
  save_plot(paste(mom, "_time.png", sep=""), mom_plot, base_height=10, base_width=12)
}

het_list <- list(length=length(svals))
c <- 1

for(s in svals) {
  p <- ggplot(data=filter(m_demo, statistic==c("scaled_Hr", "pi0"), lookup_s==s, lookup_r==1e-10),
              aes(x=Generation, y=value, color=statistic)) + 
    geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(N1)) +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none")
    labs(title="Temporal dynamics after a size change", x="Generation", y="Value") +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
    
  het_list[[c]] <- p
  c <- c + 1
}

het_plot <- plot_grid(plotlist=het_list, ncol=1, align='v')
save_plot("het_time.png", het_plot, base_height=10, base_width=12)

r <- ggplot(data=filter(m_demo, Generation==50000, N1==1e+5, statistic=="B"),
       aes(x=lookup_r, y=value)) + 
  geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(lookup_s), nrow=1) +
  scale_x_log10(breaks=unique(m_demo$lookup_r)) +
  scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none") +
  labs(title=NULL, x="r", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("B_r.png", r, base_height=7, base_width=10)

d_stats <- m_het_demo %>% 
     group_by(scenario) %>% 
     reframe(N1=N1, lookup_r=lookup_r, lookup_s=lookup_s, Generation=Generation,
             d_pi0=-c(diff(pi0)[1], diff(pi0))/pi0,
             d_Hl=-c(diff(Hl)[1], diff(Hl))/Hl,
             d_Hr=-c(diff(scaled_Hr)[1], diff(scaled_Hr))/scaled_Hr,
             d_B=-c(diff(B)[1], diff(B))/B,
             d_piN_pi0=-c(diff(piN_pi0)[1], diff(piN_pi0))/piN_pi0,
             d_piN_piS=-c(diff(piN_piS)[1], diff(piN_piS))/piN_piS)

m_demo_d <- pivot_longer(d_stats, cols=starts_with("d_"), names_to="statistic")

c <- 1
plot_list <- list(length=length(unique(m_demo_d$N1)))

for(N in unique(m_demo_d$N1)) {
  p <- ggplot(data=filter(m_demo_d, N1==N, lookup_r==1e-6,
                          lookup_s %in% c(-1e-3, -1e-4, -1e-5), 
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
      p <- p + labs(title="Temporal dynamics of statistics after a size change (rows->N1)",
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
save_plot("diffs_stats_time.png", cp, base_height=10, base_width=12)

#####################################
#
# multiple constrained loci
#
#####################################

num_exons <- 1001
ncsl <- rep(50000, num_exons)
exon_lengths <- 1000
csl <- rep(exon_lengths, num_exons) 
L <- sum(csl) + sum(ncsl)

for(N in unique(hr_demo$N1)) {
  
  Bvals_gen <- as.data.frame(matrix(ncol=length(sampling_times), 
                                    nrow=nrow(look_s)))
  names(Bvals_gen) <- as.character(sampling_times)
  
  for(i in 1:nrow(look_s)) {
    
    print(Sys.time())
    cat(paste(i, "\n"))
    
    ss <- rep(as.numeric(look_s[i]), num_exons)
  
    smap <- suppressWarnings(c(rbind(ncsl, csl)))
    smap <- suppressMessages(setDT(bind_cols(chr="chr1",
                                 dplyr::lag(cumsum(smap), n=1, default=0),
                                 dplyr::lag(cumsum(smap), n=0, default=0))))
    names(smap) <- c("chr", "start", "end")
    smap$s <- c(rbind(rep(0, length(csl)), ss))
    setkey(smap, start, end)
  
    dt_neutral <- filter(smap[,2:4], s==0)
    dt_exons <- filter(smap[,2:4], s<0)
    
    rmap <- setDT(bind_cols(chr="chr1", start=0, end=L, r=1e-8)) # flat
    mmap <- setDT(bind_cols(chr="chr1", start=0, end=L, u=1e-8)) # flat
  
    setkey(mmap, start, end)
    setkey(rmap, start, end)
    setkey(smap, start, end)
    setkey(dt_neutral, start, end)
    setkey(dt_exons, start, end)
  
    all_pos <- sort(c(dt_neutral[nrow(dt_neutral)/2]$start + 500, 
                       dt_exons$start + exon_lengths / 2))
    focal_mu <- rep(1e-8, length(all_pos))
    
    pos_dt <- suppressMessages(setDT(bind_cols(all_pos, focal_mu)))
    names(pos_dt) <- c("position", "focal_mu")
    pos_dt$idx <- 1:nrow(pos_dt) # indexing
    setkey(pos_dt, position)
    
    Nanc <- unique(m_het_demo$Na)
    m_hr_N <- filter(m_het_demo, N1==N) 
  
    c <- 1
    for(g in sort(unique(m_hr_N$Generation), decreasing=T)) {
  
      Ne_bar <- unique(filter(m_hr_N, Generation==g)$pi0) / 
                (2*unique(m_hr_N$u))
      
      #num_iter <- 1 NOTE: the following is without interference correction ATM
      B_values <- rep(1, length(all_pos)) # init 
      
      # approximates cumulative rec by looking at r at sampled sites only
      focal_r <- unlist(lapply(all_pos, function(pos) rmap[J(pos), roll=T]$r)) 
      focal_r <- focal_r * B_values
      cum_rec <- numeric(length=length(all_pos))
      cum_rec[1] <- focal_r[1]
      
      for(j in 2:length(all_pos)) {
        cum_rec[j] <- cum_rec[j-1] + (all_pos[j] - all_pos[j-1]) * focal_r[j]
      }
      pos_dt$cumrec <- cum_rec
      
      # pre-computes "effective" rec. dist. between sampled sites and exons
      cr_exons <- pos_dt[position %in% (dt_exons$start + exon_lengths / 2),]
    
      prd <- 1
      if(N <= Nanc) {
        prd <- as.data.frame(2*N*abs(outer(pos_dt$cumrec, cr_exons$cumrec, "-")))
      } else {
        prd <- as.data.frame(2*Ne_bar*abs(outer(pos_dt$cumrec, cr_exons$cumrec, "-")))
      }
      names(prd) <- 1:ncol(prd)
      
      erd <- prd # "effective" rec. distance is inversely proportional to alpha
      for(j in 1:ncol(erd)) { erd[,j] <- erd[,j] / abs(2 * Ne_bar * dt_exons$s[j]) }
      relevant_exons <- apply(erd, 2, function(x) x < 10 * Ne_bar * 1e-3)
      # identifying relevant exons (w.r.t. linked selection) for each sampled site
      exons_per_samp_site <- apply(relevant_exons, 1, function(x) which(x))
      # check if class(exons_per_samp_site) [1] "matrix" "array" 
      
      getB <- function(focal_exon, focal_samp) { # arguments are site indices
        total_r <- prd[focal_samp, focal_exon] / (2 * Ne_bar) 
        if(total_r > 1e-2) { return(1) }
        else {
          exon_pos <- dt_exons$start[focal_exon] + exon_lengths / 2 # midpoint
          idx <- pos_dt[.(exon_pos)]$idx # index within pos_dt / all_pos
          dt_exons[focal_exon,]$s <- dt_exons[focal_exon,]$s * B_values[idx]
          focal_s <- dt_exons[focal_exon,]$s
          
          closest_r <- look_r[look_r[.(total_r), roll="nearest", which=T]]
  
          # both Hr and pi0 are scaled linearly by Ne_bar -> need not include it
          hr <- m_hr_N[.(N, closest_r, focal_s, g)]$Hr 
          pi0 <- m_hr_N[.(N, closest_r, focal_s, g)]$pi0 
          return((hr / pi0) ^ (exon_lengths * B_values[idx]))
        }
      }
      
      tmp <- B_values 
      samp_neut <- filter(pos_dt, position==dt_neutral[nrow(dt_neutral)/2]$start + 500)$idx
      for(k in samp_neut) {
        if(length(exons_per_samp_site[[k]]) > 0) {
          B <- unlist(lapply(exons_per_samp_site[[k]], getB, focal_samp=k))
          tmp[k] <- cumprod(B)[length(B)]
        }
        else {
          tmp[k] <- 1 
        }
      }
  
      B_values <- tmp
      Bvals_gen[i, c] <- B_values[samp_neut]
      
      c <- c + 1
    }
  }
  
  Bvals_gen$s <- look_s$lookup_s
  m_Bs <- pivot_longer(Bvals_gen, cols=as.character(sampling_times))
  
  tbl_gen <- pivot_longer(filter(hl_demo, N1==N, lookup_r==2.001000e-05), 
                          cols=as.character(sampling_times),
                          names_to="Generation", values_to="Hl")
  tbl_gen$Generation <- as.numeric(tbl_gen$Generation)
  tbl_gen <- setDT(tbl_gen)
  setkey(tbl_gen, N1, lookup_r, lookup_s, Generation)
  
  tbl_gen$pi0 <- rep(arrange(filter(m_pi0, N1==N), Generation)$pi0, 7)
  tbl_gen$B <- m_Bs$value
  
  fwrite(tbl_gen, paste("tbl_gen_N1_", N, ".csv", sep=""))
}

tbl_bottleneck <- fread("tbl_gen_N1_1000.csv")
tbl_expansion <- fread("tbl_gen_N1_1e+05.csv")
tbl_demo <- rbind.data.frame(tbl_bottleneck, tbl_expansion)
tbl_demo$Hr <- tbl_demo$B * tbl_demo$pi0
tbl_demo$piN_pi0 <- tbl_demo$Hl / tbl_demo$pi0
tbl_demo$piN_piS <- tbl_demo$Hl / tbl_demo$Hr
  
m_tbl <- pivot_longer(tbl_demo, cols=c(pi0, Hr, Hl, B, piN_pi0, piN_piS),
                      names_to="statistic")
fwrite(m_tbl, "tbl_demo.csv")