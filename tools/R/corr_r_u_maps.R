# Date created: 04/03/2024
# Last modified: 01/10/2024
# Author: Gustavo V. Barroso

#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

m <- as.numeric(args[1]) # index of focal model

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
})

models <- fread("../../models.csv")
cat(paste("Building random genomic landscapes under model", m, "\n"))

bin_size <- 1e+3 # "basal" scale, larger windows are built after
L <- models$L[m]
V <- L - models$num_exons[m] * models$exon_lengths[m]
mult <- models$mult_r_u[m]
  
## elements 
intergenic_lengths <- rgeom(n=models$num_exons[m], prob=models$num_exons[m]/V) 

while(sum(intergenic_lengths) < V) {
  intergenic_lengths <- c(intergenic_lengths, 
                          rgeom(n=1, prob=models$num_exons[m]/V))
}
if(sum(intergenic_lengths) > V) {
  intergenic_lengths <- intergenic_lengths[cumsum(intergenic_lengths) < V]
  intergenic_lengths <- c(intergenic_lengths, V - sum(intergenic_lengths))
}

elements <- suppressWarnings(c(rbind(intergenic_lengths, 
            rep(models$exon_lengths[m], models$num_exons[m]))))
elements <- suppressMessages(setDT(bind_cols(chr="chr1",
		                    dplyr::lag(cumsum(elements), n=1, default=0),
		                    dplyr::lag(cumsum(elements), n=0, default=0))))
names(elements) <- c("chr", "start", "end")
elements <- filter(elements, start < L)
elements$end[nrow(elements)] <- L
elements$selected <- as.integer(elements$end - elements$start == models$exon_lengths[m])

## rec map
rec_spans <- rgeom(n=ceiling(L/models$avg_rec_spans[m]), prob=1/models$avg_rec_spans[m])
while(sum(rec_spans) < L) {
  rec_spans <- c(rec_spans, rgeom(n=1, prob=1/models$avg_rec_spans[m]))
}
if(sum(rec_spans) > L) {
  rec_spans <- rec_spans[cumsum(rec_spans) < L]
  rec_spans <- c(rec_spans, L - sum(rec_spans))
}
rec_spans <- rec_spans[rec_spans>0]

rmap <- suppressMessages(setDT(bind_cols(chr="chr1",
        		  dplyr::lag(cumsum(rec_spans), n=1, default=0),
		          dplyr::lag(cumsum(rec_spans), n=0, default=0))))  
names(rmap) <- c("chr", "start", "end")
rs <- rgamma(n=nrow(rmap),
             shape=models$shapes_rate_r[m],
             rate=models$shapes_rate_r[m])
rmap$r <- rs * models$mean_r[m]
fwrite(rmap, "rec_map.csv.gz")

mut_spans <- rgeom(n=ceiling(L/models$avg_mut_spans[m]), prob=1/models$avg_mut_spans[m])
while(sum(mut_spans) < L) {
  mut_spans <- c(mut_spans, rgeom(n=1, prob=1/models$avg_mut_spans[m]))
}
if(sum(mut_spans) > L) {
  mut_spans <- mut_spans[cumsum(mut_spans) < L]
  mut_spans <- c(mut_spans, L - sum(mut_spans))
}
mut_spans <- mut_spans[mut_spans>0]

mmap <- suppressMessages(setDT(bind_cols(chr="chr1",
         		  dplyr::lag(cumsum(mut_spans), n=1, default=0),
		          dplyr::lag(cumsum(mut_spans), n=0, default=0))))
names(mmap) <- c("chr", "start", "end")
mus <- rgamma(n=nrow(mmap),
            	shape=models$shapes_rate_u[m],
            	rate=models$shapes_rate_u[m])
mmap$u <- mus * models$mean_u[m]

# binning
bins <- rep(bin_size, floor(L / bin_size))
bins <- suppressMessages(bind_cols(chr="chr1",
                  dplyr::lag(cumsum(bins), n=1, default=0),
                  dplyr::lag(cumsum(bins), n=0, default=0)))
names(bins) <- c("chr", "start", "end")
bins$start <- bins$start + 1 # so that GRanges understands the intervals
gr_bins <- makeGRangesFromDataFrame(bins)

rmap2 <- rmap
rmap2$start <- rmap2$start + 1 # so that GRanges understands the intervals
rec_gr <- makeGRangesFromDataFrame(rmap2, keep.extra.columns=T)
hits <- findOverlaps(query=gr_bins, subject=rec_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bins", "rec_bins")
rmap2$rec_bins <- 1:nrow(rmap2)
rec_bins <- merge(rmap2, hit_list) %>%
            group_by(bins) %>%
            mutate(avg_rec=weighted.mean(r, end-start))
rec_bins <- rec_bins %>% distinct(bins, .keep_all=T) 

mmap2 <- mmap 
mmap2$start <- mmap2$start + 1 # so that GRanges understands the intervals
mut_gr <- makeGRangesFromDataFrame(mmap2, keep.extra.columns=T)

# transform mut map with rec map
mmap2$mut_bins <- 1:nrow(mmap2)
hits <- findOverlaps(query=rec_gr, subject=mut_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("rec_bins", "mut_bins")
mmap3 <- merge(mmap2, hit_list) 
mmap3 <- dplyr::select(rmap2, -c("chr", "start", "end")) %>% merge(., mmap3, by="rec_bins")
mmap3$u <- mmap3$u + mmap3$r * mult
mmap3 <- mmap3 %>% distinct(mut_bins, .keep_all=T) %>% dplyr::select(names(mmap))
mmap3$u <- unique(models$mean_u) * (mmap3$u / mean(mmap3$u)) # re-center
mmap2 <- mmap3
mmap2$start <- mmap2$start - 1
fwrite(mmap2, "mut_map.csv.gz")

# binning mut
mut_gr <- makeGRangesFromDataFrame(mmap3, keep.extra.columns=T) # over-write
hits <- findOverlaps(query=gr_bins, subject=mut_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("bins", "mut_bins")
mmap3$mut_bins <- 1:nrow(mmap3)
mut_bins <- merge(mmap3, hit_list, by="mut_bins") %>%
            group_by(bins) %>%
            mutate(avg_mut=weighted.mean(u, end-start))
mut_bins <- mut_bins %>% distinct(bins, .keep_all=T) 

# avg mut rate within each element
elements2 <- elements 
elements2$start <- elements2$start + 1 # so that GRanges understands the intervals
elem_gr <- makeGRangesFromDataFrame(elements2, keep.extra.columns=T)
hits <- findOverlaps(query=elem_gr, subject=mut_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("elems", "mut_bins")
x <- merge(hit_list, mmap3) 
y <- x %>% group_by(elems) %>% mutate(avg_mut=weighted.mean(u, end-start)) 
x <- distinct(y, elems, .keep_all=T)
elements$u <- x$avg_mut

# finally writes elements to file
fwrite(elements, "elements.csv.gz")

maps_1kb <- bind_cols(bins,
                      rec_bins["avg_rec"],
                      mut_bins["avg_mut"])

maps_1kb$start <- maps_1kb$start - 1 # back
maps_1kb$bin <- 1:nrow(maps_1kb)
maps_1kb$bin_10kb <- (maps_1kb$bin - 1) %/% 10
maps_1kb$bin_100kb <- (maps_1kb$bin - 1) %/% 100
maps_1kb$bin_1Mb <- (maps_1kb$bin - 1) %/% 1000
maps_10kb <- maps_1kb %>% group_by(bin_10kb) %>% 
             summarise_at(c("avg_mut", "avg_rec"), mean)
maps_100kb <- maps_1kb %>% group_by(bin_100kb) %>% 
              summarise_at(c("avg_mut", "avg_rec"), mean)
maps_1Mb <- maps_1kb %>% group_by(bin_1Mb) %>% 
            summarise_at(c("avg_mut", "avg_rec"), mean)

bins_1kb <- rep(1e+3, floor(L / 1e+3))
bins_1kb <- suppressMessages(bind_cols(chr="chr1",
  dplyr::lag(cumsum(bins_1kb), n=1, default=0),
  dplyr::lag(cumsum(bins_1kb), n=0, default=0)))
names(bins_1kb) <- c("chr", "start", "end")

bins_10kb <- rep(1e+4, floor(L / 1e+4))
bins_10kb <- suppressMessages(bind_cols(chr="chr1",
  dplyr::lag(cumsum(bins_10kb), n=1, default=0),
  dplyr::lag(cumsum(bins_10kb), n=0, default=0)))
names(bins_10kb) <- c("chr", "start", "end")

bins_100kb <- rep(1e+5, floor(L / 1e+5))
bins_100kb <- suppressMessages(bind_cols(chr="chr1",
  dplyr::lag(cumsum(bins_100kb), n=1, default=0),
  dplyr::lag(cumsum(bins_100kb), n=0, default=0)))
names(bins_100kb) <- c("chr", "start", "end")

bins_1Mb <- rep(1e+6, floor(L / 1e+6))
bins_1Mb <- suppressMessages(bind_cols(chr="chr1",
  dplyr::lag(cumsum(bins_1Mb), n=1, default=0),
  dplyr::lag(cumsum(bins_1Mb), n=0, default=0)))
names(bins_1Mb) <- c("chr", "start", "end")

maps_1kb <- cbind.data.frame(bins_1kb, dplyr::select(maps_1kb, avg_rec, avg_mut))
maps_10kb <- cbind.data.frame(bins_10kb, dplyr::select(maps_10kb, avg_rec, avg_mut))
maps_100kb <- cbind.data.frame(bins_100kb, dplyr::select(maps_100kb, avg_rec, avg_mut))
maps_1Mb <- cbind.data.frame(bins_1Mb, dplyr::select(maps_1Mb, avg_rec, avg_mut))

fwrite(maps_1kb, "maps_1kb.csv.gz")
fwrite(maps_10kb, "maps_10kb.csv.gz")
fwrite(maps_100kb, "maps_100kb.csv.gz")
fwrite(maps_1Mb, "maps_1Mb.csv.gz")