#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(pracma)

jump_length <- as.numeric(args[1]) # sampling dist., for interpolating B-values
num_rounds <- 3 # updating B-values due to interference selection

lookup_tbl <- fread("lookup_tbl.csv")
setkey(lookup_tbl, lookup_r, lookup_s)
mmap <- fread("mmap.csv")
rmap <- fread("rmap.csv")
intervals <- fread("intervals.csv")

dt_neutral <- filter(intervals[,2:4], s==0)
dt_exons <- filter(intervals, s<0)

setkey(intervals, start, end)
setkey(dt_neutral, start, end)
setkey(dt_exons, start, end)

exon_lengths <- unique((dt_exons$end - dt_exons$start))

# finding mu for each sampled pos and the cumulative recombination map
neut_pos <- unlist(apply(dt_neutral, 1, function(x)
                         seq(from=x[1], to=x[2], by=jump_length))) + 1
linear_pos <- sort(c(neut_pos, dt_exons$start + exon_lengths / 2))
focal_mu <- unlist(lapply(linear_pos, function(pos) mmap[J(pos), roll=T]$u))

B_values <- rep(1, length(linear_pos)) # init
for(i in 1:num_rounds) {
  
  focal_r <- unlist(lapply(linear_pos, function(pos) rmap[J(pos), roll=T]$r)) # TODO stop approximating and take weighted avg of rec map between positions
  focal_r <- focal_r * B_values # updates all recs, not just within exons
  cum_rec <- numeric(length=length(linear_pos))
  cum_rec[1] <- focal_r[1]

  for(j in 2:length(linear_pos)) {
    cum_rec[j] <- cum_rec[j-1] + (linear_pos[j] - linear_pos[j-1]) * focal_r[j]
  }

  pos_dt <- setDT(bind_cols(linear_pos, focal_mu, cum_rec))
  names(pos_dt) <- c("position", "focal_mu", "cumrec")
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
      idx <- which(linear_pos==exon_pos)
      dt_exons[focal_exon,]$s <- dt_exons[focal_exon,]$s * B_values[idx]
      focal_s <- dt_exons[focal_exon,]$s
      # if the grid is fine-grained enough, we don't need interpolation
      closest_r <- dt_r[dt_r[.(total_r), roll="nearest", which=T]]
      closest_s <- dt_s[dt_s[.(focal_s), roll="nearest", which=T]]
      hr <- lookup_tbl[.(closest_r, closest_s)]$Hr # lookup_tbl has single u val
      # correction for interference on mu is applied by scaling exon length
      return((hr / (2 * N * u)) ^ (exon_lengths * B_values[idx]))
    }
  }
  
  # gets the B-values per sampled site
  tmp <- B_values
  # this takes time if there are many sampled sites influenced by many exons
  pb <- txtProgressBar(min=1, max=length(exons_per_samp_site), style=3)
  for(j in 1:length(exons_per_samp_site)) {
    setTxtProgressBar(pb, j)
    if(length(exons_per_samp_site[[j]]) > 0) {
      B <- unlist(lapply(exons_per_samp_site[[j]], getB, focal_samp=j))
      #B_values[j] <- cumprod(B)[length(B)]
      tmp[j] <- cumprod(B)[length(B)]
    }
  }
  close(pb)
  B_values <- tmp
  cat(summary(B_values))
}
  
# interpolating B-values for every site 1:L can also take some time
B_map <- cubicspline(linear_pos, B_values, 1:L)
  
# here it's easier to use mmap2 which has closed intervals:
ones <- rep(1, L)
pi0s <- invisible(apply(mmap2, 1, function(x)
                  ones[as.numeric(x[2]):as.numeric(x[3])] *
                  as.numeric(x[4]) * 2 * N))
pis <- invisible(apply(mmap2, 1, function(x)
                 B_map[as.numeric(x[2]):as.numeric(x[3])] *
                 as.numeric(x[4]) * 2 * N))
pi0s <- unlist(pi0s)
pis <- unlist(pis)
hrmap <- setDT(bind_cols(1:L, pis, pi0s, B_map))
names(hrmap) <- c("Pos", "Hr", "pi0", "B")
setkey(hrmap, Pos)
fwrite(hrmap, "hrmap.csv", compress="gzip") 
