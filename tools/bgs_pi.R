#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(R.utils)
library(tidyverse)
library(data.table)
library(pracma) # cubic splines function
library(cowplot)

num_rounds <- 3 # for updating B-values due to interference selection
jump_length <- as.numeric(args[1]) # sampling dist., for interpolating B-values

# loads results of bgs_maps.R
mmap <- fread("mmap.csv")
rmap <- fread("rmap.csv")
smap <- fread("smap.csv")
lookup_tbl <- fread("lookup_tbl.csv.gz")

dt_neutral <- filter(smap[,2:4], s==0)
dt_exons <- filter(smap, s<0)

setkey(mmap, start, end)
setkey(rmap, start, end)
setkey(smap, start, end)
setkey(dt_neutral, start, end)
setkey(dt_exons, start, end)
setkey(lookup_tbl, lookup_r, lookup_s)

exon_lengths <- unique((dt_exons$end - dt_exons$start))

# finding mu for each sampled pos and the cumulative recombination map
neut_pos <- unlist(apply(dt_neutral, 1, function(x) seq(from=x[1],
                                                        to=x[2], 
                                                        by=jump_length))) + 1
linear_pos <- sort(c(neut_pos, dt_exons$start + exon_lengths / 2))
focal_mu <- unlist(lapply(linear_pos, function(pos) mmap[J(pos), roll=T]$u))
B_values <- rep(1, length(linear_pos)) # init

tbl <- as.data.frame(matrix(ncol=num_rounds, nrow=length(B_values))) # viz
for(i in 1:num_rounds) {
  
  # approximates cumm. rec by looking at r at sampled sites only
  focal_r <- unlist(lapply(linear_pos, function(pos) rmap[J(pos), roll=T]$r)) 
  focal_r <- focal_r * B_values
  cum_rec <- numeric(length=length(linear_pos))
  cum_rec[1] <- focal_r[1]

  for(j in 2:length(linear_pos)) {
    cum_rec[j] <- cum_rec[j-1] + (linear_pos[j] - linear_pos[j-1]) * focal_r[j]
  }

  pos_dt <- setDT(bind_cols(linear_pos, focal_mu, cum_rec))
  pos_dt$idx <- 1:nrow(pos_dt)
  names(pos_dt) <- c("position", "focal_mu", "cumrec", "idx")
  setkey(pos_dt, position)

  # pre-computes "effective" rec. dist. between sampled sites and exons
  cr_exons <- pos_dt[position %in% (dt_exons$start + exon_lengths / 2),]
  prd <- as.data.frame(2 * N * abs(outer(pos_dt$cumrec, cr_exons$cumrec, "-")))
  names(prd) <- 1:ncol(prd)
  erd <- prd # "effective" rec. distance is inversely proportional to alpha
  for(j in 1:ncol(erd)) { erd[,j] <- erd[,j] / abs(2 * N * dt_exons$s[j]) }
  relevant_exons <- apply(erd, 2, function(x) x < 10 * N * 1e-3)
  exons_per_samp_site <- apply(relevant_exons, 1, function(x) which(x))
  
  getB <- function(focal_exon, focal_samp) { 
    total_r <- prd[focal_samp, focal_exon] / (2 * N) 
    if(total_r > 1e-2) { return(1) } 
    else {
      exon_pos <- dt_exons$start[focal_exon] + exon_lengths / 2 # midpoint
      idx <- pos_dt[.(exon_pos)]$idx
      dt_exons[focal_exon,]$s <- dt_exons[focal_exon,]$s * B_values[idx]
      focal_s <- dt_exons[focal_exon,]$s
      
      # if the grid is fine-grained enough, we don't need interpolation
      closest_r <- dt_r[dt_r[.(total_r), roll="nearest", which=T]]
      closest_s <- dt_s[dt_s[.(focal_s), roll="nearest", which=T]]
      hr <- lookup_tbl[.(closest_r, closest_s)]$Hr # lookup_tbl has single u val
      
      # correction for interference on u is applied by scaling exon length
      return((hr / (2 * N * u)) ^ (exon_lengths * B_values[idx]))
    }
  }
  
  # takes time if there are many sampled sites influenced by many exons
  tmp <- B_values
  pb <- txtProgressBar(min=1, max=length(exons_per_samp_site), style=3)
  for(j in 1:length(exons_per_samp_site)) {
    setTxtProgressBar(pb, j)
    if(length(exons_per_samp_site[[j]]) > 0) {
      B <- unlist(lapply(exons_per_samp_site[[j]], getB, focal_samp=j))
      tmp[j] <- cumprod(B)[length(B)]
    }
  }
  close(pb)
  
  B_values <- tmp
  tbl[,i] <- B_values
}

names(tbl) <- c(paste("iter_", rep(1:num_rounds), sep=""))
tbl$pos <- linear_pos
m_tbl <- pivot_longer(tbl, cols=starts_with("iter"), names_to="Iteration")
ggplot(data=m_tbl[5e+4:6e+4,], aes(x=pos, y=value, color=Iteration)) + geom_point(aes(alpha=0.5)) + theme_bw()

# pause for vizualizing the validity of threshold for spotting "relevant" exons
rex <- apply(prd, 2, function(x) x < 4 * N * 1e-2)
ex400 <- apply(rex, 1, function(x) which(x))

Bl400 <- unlist(lapply(ex400[[1]], getB, focal_neutral=1))
Bq400 <- unlist(lapply(ex400[[length(samp_pos)/4]], getB,
                       focal_neutral=length(samp_pos)/4))
Bm400 <- unlist(lapply(ex400[[length(samp_pos)/2]], getB,
                       focal_neutral=length(samp_pos)/2))

Blr <- unlist(lapply(exons_per_samp_neut[[1]], getB, focal_neutral=1))
Bqr <- unlist(lapply(exons_per_samp_neut[[length(samp_pos)/4]], getB,
                     focal_neutral=length(samp_pos)/4))
Bmr <- unlist(lapply(exons_per_samp_neut[[length(samp_pos)/2]], getB,
                     focal_neutral=length(samp_pos)/2))

xl <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bl400))]), Bl400)
xq <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bq400))]), Bq400)
xm <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bm400))]), Bm400)

names(xl) <- c("rec", "B")
names(xq) <- c("rec", "B")
names(xm) <- c("rec", "B")

xl$CummBval <- cumprod(xl$B)
xq$CummBval <- cumprod(xq$B)
xm$CummBval <- cumprod(xm$B)

xl$relevant <- xl$B %in% Blr
xq$relevant <- xq$B %in% Bqr
xm$relevant <- xm$B %in% Bmr

pl <- ggplot(data=xl, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values Left", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

pq <- ggplot(data=xq, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values 1/4 chr", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

pm <- ggplot(data=xm, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values 1/2 chr", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

ql <- ggplot(data=xl, aes(x=rec, y=CummBval, color=relevant)) + geom_point() +
  theme_bw() + scale_y_log10(limits=c(min(xl$CummBval), 1)) +
  labs(title="Cumm. B-value, Left", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

qq <- ggplot(data=xq, aes(x=rec, y=CummBval, color=relevant)) + 
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="Cumm. B-value, 1/4 chr", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

qm <- ggplot(data=xm, aes(x=rec, y=CummBval, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="Cumm. B-value, 1/2 chr", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

cp <- plot_grid(pl, pq, pm, ql, qq, qm, nrow=2)
save_plot(paste("rep_", i, "/Bvals.png", sep=""),
          cp, base_height=10, base_width=12)

# interpolating B-values for every site 1:L can also take some time
B_map <- cubicspline(linear_pos, B_values, 1:L)
ones <- rep(1, L)

# uses closed intervals for mmap & rmap (+1 start coordinates)
pi0s <- invisible(apply(mmap, 1, function(x)
                  ones[(as.numeric(x[2])+1):as.numeric(x[3])] *
                  as.numeric(x[4]) * 2 * N))
pis <- invisible(apply(mmap, 1, function(x)
                 B_map[(as.numeric(x[2])+1):as.numeric(x[3])] *
                 as.numeric(x[4]) * 2 * N))

pi0s <- unlist(pi0s)
pis <- unlist(pis)
hrmap <- setDT(bind_cols(1:L, pis, pi0s, B_map))
names(hrmap) <- c("Pos", "Hr", "pi0", "B")
fwrite(hrmap, "hrmap.csv.gz", compress="gzip") 

# binning 
hrmap$bin <- ((hrmap$Pos - 1) %/% bin_size)
hrmap_1kb <- hrmap %>% group_by(bin) %>%
             summarize(avg_pi=mean(Hr), avg_pi0=mean(pi0), avg_B=mean(B))
hrmap_1kb$bin_10kb <- hrmap_1kb$bin %/% 10
hrmap_1kb$bin_100kb <- hrmap_1kb$bin %/% 100
hrmap_10kb <- hrmap_1kb %>% group_by(bin_10kb) %>% 
              summarise_at(c("avg_pi", "avg_pi0", "avg_B"), mean)
hrmap_100kb <- hrmap_1kb %>% group_by(bin_100kb) %>% 
               summarise_at(c("avg_pi", "avg_pi0", "avg_B"), mean)

fwrite(hrmap_1kb, "hrmap_1kb.csv")
fwrite(hrmap_10kb, "hrmap_10kb.csv")
fwrite(hrmap_100kb, "hrmap_100kb.csv")
