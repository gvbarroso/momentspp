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

setwd("~/Data/momentspp/paper_1/single_neutral")

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
  
  moms <- fread(paste("demo/model_", i, 
                      "_e_1_expectations.txt", sep=""))
  
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
  plots[[i]] <- ggplot(data=filter(lookup_tbl, lookup_r==1e-8,
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

x <- plot_grid(plotlist=plots)
save_plot("hr_s_n.png", x, base_height=8, base_width=16)

for(N in unique(hr_demo$N1)) {
  
  Nanc <- unique(m_hr_demo$Na)
  
  m_hr_demo <- pivot_longer(filter(hr_demo, N1==N), 
                            cols=as.character(sampling_times),
                            names_to="Generation", values_to="Hr")
  m_hr_demo$Generation <- as.numeric(m_hr_demo$Generation)
  m_hr_demo <- setDT(m_hr_demo)
  setkey(m_hr_demo, N1, lookup_r, lookup_s, Generation)
  
  m_hr_demo <- left_join(m_hr_demo, dplyr::select(m_pi0, -c(Na, t)),
                         by=c("N1", "Generation"))
  
  num_exons <- 1001
  ncsl <- rep(50000, num_exons)
  exon_lengths <- 1000
  csl <- rep(exon_lengths, num_exons) 
  L <- sum(csl) + sum(ncsl)
  
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
  
    c <- 1
    for(g in sort(unique(m_hr_demo$Generation), decreasing=T)) {
  
      Ne_bar <- unique(filter(m_hr_demo, Generation==g)$pi0) / 
                (2*unique(m_hr_demo$u))
      
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
          hr <- m_hr_demo[.(N, closest_r, focal_s, g)]$Hr 
          pi0 <- m_hr_demo[.(N, closest_r, focal_s, g)]$pi0 
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
  m_Bs <- pivot_longer(Bvals_gen, cols=as.character(sampling_times),
                       names_to="Generation")
  
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
tbl_expansion <- fread("tbl_gen_N1_100000.csv")


pa <- ggplot(data=m_tbl, aes(x=Generation, y=value)) + 
  geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(s)) +
  #scale_fill_manual(values=mycolors) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none") + 
  labs(title="B-value map over time", x="Generation", y="B-value") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("B-vals_time.png", pa, base_height=8, base_width=16)


cat("Performing cubic splines interpolation...")
B_map <- cubicspline(all_pos, B_values, 1:L)
cat("done.\nNow preparing output files.\n")
ones <- rep(1, L)





# focusing on a particular (r, s) pair
closest_s <- look_s[which.min(abs(look_s - -1.5e-4))]
closest_r <- look_r[which.min(abs(look_r - 1e-8))]

plots <- list(length=length(unique(lookup_tbl$N1)))
for(i in 1:length(unique(lookup_tbl$N1))) {
  
  tmp_hl <- filter(m_hl_demo, N1==unique(m_hl_demo$N1)[i], 
                              lookup_r==closest_r, 
                              lookup_s==closest_s)
  tmp_hr <- filter(m_hl_demo, N1==unique(m_hr_demo$N1)[i], 
                              lookup_r==closest_r, 
                              lookup_s==closest_s)
  tmp_pi0 <- filter(m_pi0, N1==unique(m_hr_demo$N1)[i])
  
  
  stats <- tmp_hr
  stats$Hl <- tmp_hl$Hl
  stats$pi0 <- tmp_pi0$pi0
  stats$B <- stats$Hr / stats$pi0
  stats$E <- stats$Hl / stats$pi0
    
  m_stats <- pivot_longer(stats, cols=c(B), names_to="Statistic")
    
  plots[[i]] <- ggplot(data=m_stats, aes(x=Generation,
                                         y=value, 
                                         color=Statistic)) +
                  geom_point() + geom_line() +  theme_bw() + scale_y_log10() +
                  labs(title=NULL, x="Generation", y=NULL) +
                  theme(axis.title=element_text(size=16), 
                        axis.text=element_text(size=12), 
                        axis.text.x=element_text(size=12),
                        legend.text=element_text(size=16),
                        legend.title=element_text(size=16),
                        strip.text=element_text(size=14),
                        legend.position="bottom")
}

x <- plot_grid(plotlist=plots)
save_plot("stats_time.png", x, base_height=8, base_width=16)
