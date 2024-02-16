#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(tidyverse)
library(data.table)
library(GenomicRanges)

models <- fread("~/Devel/momentspp/benchmarking/selection/single-pop/landscapes/models.csv")

mdl <- args[1] # which of the models
jump_length <- args[2] # for later interpolating B-values
rec_threshold <- args[3] # for considering whether to include the effect of exons, eg 1e-2 or 1e-3
bin_size <- 1e+4

ncsl <- rgeom(n=models$num_constrained[mdl] + 1, prob=1/models$avg_neutral_lengths[mdl])
csl <- rep(models$exon_lengths[mdl], models$num_constrained[mdl]) 
L <- sum(csl) + sum(ncsl)

ss <- -rgamma(n=length(csl), shape=models$shape_sel[mdl], scale=models$scale_sel[mdl]/models$shape_sel[mdl])
alpha <- 2 * models$N[mdl] * ss

while(any(alpha < -20)) { 
  cat("Selection is too strong in at least one gene! Redrawing...\n")
  ss <- -rgamma(n=length(csl), shape=models$shape_sel[mdl], scale=models$scale_sel[mdl]/models$shape_sel[mdl])
  alpha <- 2 * models$N[mdl] * ss
}

disc_ss <- dt_s[dt_s[J(ss), roll="nearest", which=T]]$lookup_s

intervals <- c(rbind(ncsl, csl))[-2*length(ncsl)]
intervals <- setDT(bind_cols(chr="chr1", dplyr::lag(cumsum(intervals), n=1, default=0), dplyr::lag(cumsum(intervals), n=0, default=0)))
names(intervals) <- c("chr", "start", "end")
intervals$s <- c(rbind(rep(0, length(csl)), disc_ss), 0)
setkey(intervals, start, end)

dt_neutral <- filter(intervals[,2:4], s==0)
setkey(dt_neutral, start, end)
dt_exons <- filter(intervals, s < 0)
dt_exons$u <- 1e-08 
setkey(dt_exons, start, end)

rec_spans <- rgeom(n=ceiling(L/models$avg_rec_spans[mdl]), prob=1/models$avg_rec_spans[mdl])
while(sum(rec_spans) < L) {
  rec_spans <- c(rec_spans, rgeom(n=1, prob=1/models$avg_rec_spans[mdl]))
}
if(sum(rec_spans) > L) {
  rec_spans <- rec_spans[cumsum(rec_spans) < L]
  rec_spans <- c(rec_spans, L - sum(rec_spans))
}
rec_spans <- rec_spans[rec_spans>0]

rmap <- setDT(bind_cols(chr="chr1", dplyr::lag(cumsum(rec_spans), n=1, default=0), dplyr::lag(cumsum(rec_spans), n=0, default=0)))
names(rmap) <- c("chr", "start", "end")
rs <- rgamma(n=nrow(rmap), shape=models$shapes_rec[mdl], scale=models$scales_rec[mdl]/models$shapes_rec[mdl])
disc_rs <- dt_r[dt_r[J(rs), roll="nearest", which=T]]$lookup_r
rmap$r <- disc_rs
setkey(rmap, start, end)

mut_spans <- rgeom(n=ceiling(L/models$avg_mut_spans[mdl]), prob=1/models$avg_mut_spans[mdl])
while(sum(mut_spans) < L) {
  mut_spans <- c(mut_spans, rgeom(n=1, prob=1/models$avg_mut_spans[mdl]))
}
if(sum(mut_spans) > L) {
  mut_spans <- mut_spans[cumsum(mut_spans) < L]
  mut_spans <- c(mut_spans, L - sum(mut_spans))
}
mut_spans <- mut_spans[mut_spans>0]

mu_dist <- round(seq_log(from=1e-10, to=1e-6, length.out=300), digits=11)
dt_u <- data.table(mu_dist)
setkey(dt_u, mu_dist)

mmap <- setDT(bind_cols(chr="chr1", dplyr::lag(cumsum(mut_spans), n=1, default=0), dplyr::lag(cumsum(mut_spans), n=0, default=0)))
names(mmap) <- c("chr", "start", "end")
mus <- rgamma(n=nrow(mmap), shape=models$shapes_mu[mdl], scale=models$scales_mu[mdl]/models$shapes_mu[mdl])
disc_mus <- dt_u[dt_u[J(mus), roll="nearest", which=T]]$mu_dist
mmap$u <- disc_mus
setkey(mmap, start, end)

# finding mu for each sampled pos and the cumulative recombination map over the useful intervals
samp_pos <- unlist(apply(dt_neutral, 1, function(x) seq(from=x[1], to=x[2], by=jump_length))) + 1
linear_pos <- sort(c(samp_pos, dt_exons$start + models$exon_lengths[mdl]/2))
# (the following searches take some time:)
focal_mu <- unlist(lapply(linear_pos, function(pos) mmap[J(pos), roll=T]$u))
focal_r <- unlist(lapply(linear_pos, function(pos) rmap[J(pos), roll=T]$r))
cum_rec <- numeric(length=length(linear_pos))
cum_rec[1] <- focal_r[1]
for(j in 2:length(linear_pos)) { cum_rec[j] <- cum_rec[j-1] + (linear_pos[j] - linear_pos[j-1]) * focal_r[j] }

pos_dt <- setDT(bind_cols(linear_pos, focal_mu, cum_rec))
names(pos_dt) <- c("position", "focal_mu", "cumrec")
setkey(pos_dt, position)

# binning for linear models
bins <- rep(bin_size, floor(L / bin_size))
bins <- bind_cols(chr="chr1", dplyr::lag(cumsum(bins), n=1, default=0), dplyr::lag(cumsum(bins), n=0, default=0))
names(bins) <- c("chr", "start", "end")
bins$start <- bins$start + 1 # so that GRanges understands the intervals -> NOT in BEDGRAPH format
gr_bins <- makeGRangesFromDataFrame(bins)

# binning recombination maps
rmap$start <- rmap$start + 1 # so that GRanges understands the intervals -> NOT in BEDGRAPH format
rec_gr <- makeGRangesFromDataFrame(rmap, keep.extra.columns=T)
hits <- findOverlaps(query=gr_bins, subject=rec_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bins", "rec_bins")
rmap$rec_bins <- 1:nrow(rmap)
rec_bins <- merge(rmap, hit_list) %>% group_by(bins) %>% mutate(avg_rec=weighted.mean(r, end-start))
rec_bins <- rec_bins %>% distinct(bins, .keep_all=T)

# binning mutation maps
mmap$start <- mmap$start + 1 # so that GRanges understands the intervals -> NOT in BEDGRAPH format
mut_gr <- makeGRangesFromDataFrame(mmap, keep.extra.columns=T)
hits <- findOverlaps(query=gr_bins, subject=mut_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bins", "mut_bins")
mmap$mut_bins <- 1:nrow(mmap)
mut_bins <- merge(mmap, hit_list) %>% group_by(bins) %>% mutate(avg_mut=weighted.mean(u, end-start))
mut_bins <- mut_bins %>% distinct(bins, .keep_all=T)

# binning DFE
intervals$start <- intervals$start + 1 # so that GRanges understands the intervals -> NOT in BEDGRAPH format
dfe_gr <- makeGRangesFromDataFrame(intervals, keep.extra.columns=T)
hits <- findOverlaps(query=gr_bins, subject=dfe_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bins", "dfe")
intervals$dfe <- 1:nrow(intervals)
# avg_s encompasses density of constrained sites
dfe_bins <- merge(intervals, hit_list) %>% group_by(bins) %>% mutate(avg_s=weighted.mean(s, end-start)) 
dfe_bins <- dfe_bins %>% distinct(bins, .keep_all=T)

maps_bins <- bind_cols(bins, rec_bins["avg_rec"], mut_bins["avg_mut"], dfe_bins["avg_s"])

# pre-compute recombination distances from sampled neutral sites and exons:
cr_samp <- pos_dt[position %in% samp_pos,]
cr_exons <- pos_dt[position %in% (dt_exons$start + models$exon_lengths[mdl]/2),]
prd <- as.data.frame(abs(outer(cr_samp$cumrec, cr_exons$cumrec, "-"))) # pairwise rec distance
names(prd) <- 1:ncol(prd)
relevant_exons <- apply(prd, 2, function(x) x < rec_threshold)
exons_per_samp_neut <- apply(relevant_exons, 1, function(x) which(x))

getB <- function(focal_exon, x) { # x is the index of focal sampled site
  focal_s <- dt_exons[focal_exon,]$s
  total_r <- prd[x, focal_exon]
  # NOTE: the grid on r is so fine-grained that in practice we don't need interpolation
  closest_r <- dt_r[dt_r[.(total_r), roll="nearest", which=T]]
  hr <- lookup_tbl[.(closest_r, focal_s)]$Hr
  return((hr / (2 * N * u)) ^ models$exon_lengths[mdl])
}

# we finally get the B-values per sampled site!
B_values <- numeric(length=length(exons_per_samp_neut))
for(j in 1:length(exons_per_samp_neut)) {
  B <- unlist(lapply(exons_per_samp_neut[[j]], getB, x=j))
  B_values[j] <- cumprod(B)[length(B)]
}

system.time(B_map <- cubicspline(samp_pos, B_values, 1:L))
pis <- B_map * unlist(lapply(1:L, function(pos) mmap[J(pos), roll=T]$u))
