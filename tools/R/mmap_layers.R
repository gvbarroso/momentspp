suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(MASS)
  library(lmtest)
  library(nlme)
  library(car)
  library(scales)
  library(data.table)
  library(interactions)
  library(GenomicRanges)
})

N <- 1e+5
L <- 3e+7
avg_mut_spans_large <- 1e+5
avg_mut_spans_small <- 10

# large scale, low variance
mut_spans_large <- round(rnorm(n=ceiling(L/avg_mut_spans_large),
    mean=avg_mut_spans_large, sd=avg_mut_spans_large/10))
while(sum(mut_spans_large) < L) {
  mut_spans_large <- c(mut_spans_large, round(rnorm(n=1,
      mean=avg_mut_spans_large, sd=avg_mut_spans_large/10)))
}
if(sum(mut_spans_large) > L) {
  mut_spans_large <- mut_spans_large[cumsum(mut_spans_large) < L]
  mut_spans_large <- c(mut_spans_large, L - sum(mut_spans_large))
}
mut_spans_large <- mut_spans_large[mut_spans_large>0]

mmap_large <- suppressMessages(setDT(bind_cols(chr="chr1",
    dplyr::lag(cumsum(mut_spans_large), n=1, default=1),
    dplyr::lag(cumsum(mut_spans_large), n=0, default=0))))  
names(mmap_large) <- c("chr", "start", "end")
mmap_large$uf <- rgamma(n=length(mut_spans_large), shape=30, rate=30)

MASS::fitdistr(mmap_large$uf, "gamma")

# small scale, high variance
mut_spans_small <- round(rnorm(n=ceiling(L/avg_mut_spans_small),
  mean=avg_mut_spans_small, sd=avg_mut_spans_small/10))
while(sum(mut_spans_small) < L) {
  mut_spans_small <- c(mut_spans_small, round(rnorm(n=1,
      mean=avg_mut_spans_small, sd=avg_mut_spans_small/10)))
}
if(sum(mut_spans_small) > L) {
  mut_spans_small <- mut_spans_small[cumsum(mut_spans_small) < L]
  mut_spans_small <- c(mut_spans_small, L - sum(mut_spans_small))
}
mut_spans_small <- mut_spans_small[mut_spans_small>0]

mmap_small <- suppressMessages(setDT(bind_cols(chr="chr1",
    dplyr::lag(cumsum(mut_spans_small), n=1, default=1),
    dplyr::lag(cumsum(mut_spans_small), n=0, default=0))))  
names(mmap_small) <- c("chr", "start", "end")
mmap_small$uf <- rgamma(n=length(mut_spans_small), shape=3, rate=3)

MASS::fitdistr(mmap_small$uf, "gamma")

gr_small <- makeGRangesFromDataFrame(mmap_small, keep.extra.columns=T)
gr_large <- makeGRangesFromDataFrame(mmap_large, keep.extra.columns=T)

hits <- findOverlaps(query=gr_small, subject=gr_large) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("small_bin", "large_bin")

v <- numeric(length=nrow(mmap_small))
pb <- txtProgressBar(min=1, max=nrow(mmap_large), style=3)
for(i in 1:nrow(mmap_large)) {
  setTxtProgressBar(pb, i)
  tmp <- dplyr::filter(hit_list, large_bin==i)
  v[tmp$small_bin] <- mmap_small[tmp$small_bin, uf] * mmap_large[tmp$large_bin, uf]
}
close(pb)

mmap <- dplyr::select(mmap_small, -uf)
mmap$u <- v * 1e-8

# binning
bin_size <- 1e+3
bins <- rep(bin_size, floor(L / bin_size))
bins <- suppressMessages(bind_cols(chr="chr1",
         dplyr::lag(cumsum(bins), n=1, default=1),
         dplyr::lag(cumsum(bins), n=0, default=0)))
names(bins) <- c("chr", "start", "end")
gr_bins <- makeGRangesFromDataFrame(bins)

mut_gr <- makeGRangesFromDataFrame(mmap, keep.extra.columns=T)
hits <- findOverlaps(query=gr_bins, subject=mut_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bin_1kb", "mut_bins")
mmap$mut_bins <- 1:nrow(mmap)
mmap_1kb <- merge(mmap, hit_list) %>% group_by(bin_1kb) %>%
  mutate(avg_mut=weighted.mean(u, end-start))
mmap_1kb <- mmap_1kb %>% distinct(bin_1kb, .keep_all=T) 
mmap_1kb$start <- bins$start
mmap_1kb$end <- bins$end

mmap_1kb$bin_10kb <- (mmap_1kb$bin_1kb - 1) %/% 10
mmap_1kb$bin_100kb <- (mmap_1kb$bin_1kb - 1) %/% 100
mmap_1kb$bin_1Mb <- (mmap_1kb$bin_1kb - 1) %/% 1000

mmap_10kb <- maps_1kb %>% group_by(bin_10kb) %>% summarise_at("avg_mut", mean)
mmap_100kb <- maps_1kb %>% group_by(bin_100kb) %>% summarise_at("avg_mut", mean)
mmap_1Mb <- maps_1kb %>% group_by(bin_1Mb) %>% summarise_at("avg_mut", mean)

p1 <- ggplot(data=mmap_1kb, aes(x=bin_1kb, y=avg_mut)) + geom_line() + theme_classic() +
  geom_smooth(method="loess", span=0.08, fill='darkred')

p2 <- ggplot(data=mmap_10kb, aes(x=bin_10kb, y=avg_mut)) + geom_line() + theme_classic() +
  geom_smooth(method="loess", span=0.08, fill='darkred')

p3 <- ggplot(data=mmap_100kb, aes(x=bin_100kb, y=avg_mut)) + geom_line() + theme_classic() +
  geom_smooth(method="loess", span=0.08, fill='darkred')

p4 <- ggplot(data=mmap_1Mb, aes(x=bin_1Mb, y=avg_mut)) + geom_line() + theme_classic() 

y1kb <- mmap_1kb$avg_mut / mean(mmap_1kb$avg_mut)
MASS::fitdistr(y1kb, "gamma")

y10kb <- mmap_10kb$avg_mut / mean(mmap_10kb$avg_mut)
MASS::fitdistr(y10kb, "gamma")

y100kb <- mmap_100kb$avg_mut / mean(mmap_100kb$avg_mut)
MASS::fitdistr(y100kb, "gamma")

y1Mb <- mmap_1Mb$avg_mut / mean(mmap_1Mb$avg_mut)
MASS::fitdistr(y1Mb, "gamma")
