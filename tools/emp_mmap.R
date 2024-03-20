
suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
})

setwd("~/Devel/momentspp/benchmarking/selection/single-pop/SeplyarskiyEtAl2023/")

mmap_chr1 <- fread("VCFs/1_rate_v5.2_TFBS_correction_all.vcf.bgz",
                   showProgress=T, verbose=T, nrows=3*1e+7, 
                   select=c("POS", "REF", "INFO"))

tmp <- data.frame(do.call('rbind', strsplit(mmap_chr1$INFO, ';', fixed=TRUE)))
tmp <- data.frame(do.call('rbind', strsplit(tmp$X2, '=', fixed=TRUE)))
mmap_chr1$mu <- as.numeric(tmp$X2)
rm(tmp)
  
mmap_chr1 <- dplyr::select(mmap_chr1, c(POS, mu))
mmap_chr1 <- mmap_chr1 %>% group_by(POS) %>% mutate(u=sum(mu))
mmap_chr1 <- distinct(mmap_chr1, POS, .keep_all=T)
mmap_chr1 <- dplyr::select(mmap_chr1, c(POS, u))
hist(mmap_chr1$u, nclass=100)

p <- ggplot(data=mmap_chr1[1:1e+6,], aes(x=POS, y=u)) +
  geom_line() + theme_bw() + scale_y_log10() +
  labs("Empirical mutation map", x="Pos", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="top")

save_plot("emp_mmap.png", p, base_height=10, base_width=12)

mmap_chr1$bin_1kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+3
mmap_chr1$bin_10kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+4
mmap_chr1$bin_100kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+5
mmap_chr1$bin_1Mb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+6

mmap_1kb <- mmap_chr1 %>% group_by(bin_1kb) %>% mutate(avg_u=mean(u))
mmap_1kb <- distinct(mmap_1kb, bin_1kb, .keep_all=T)
hist(mmap_1kb$avg_u, nclass=100)

mmap_10kb <- mmap_chr1 %>% group_by(bin_10kb) %>% mutate(avg_u=mean(u))
mmap_10kb <- distinct(mmap_10kb, bin_10kb, .keep_all=T)
hist(mmap_10kb$avg_u, nclass=100)

mmap_100kb <- mmap_chr1 %>% group_by(bin_100kb) %>% mutate(avg_u=mean(u))
mmap_100kb <- distinct(mmap_100kb, bin_100kb, .keep_all=T)
hist(mmap_100kb$avg_u, nclass=20)

mmap_1Mb <- mmap_chr1 %>% group_by(bin_1Mb) %>% mutate(avg_u=mean(u))
mmap_1Mb <- distinct(mmap_1Mb, bin_1Mb, .keep_all=T)
hist(mmap_1Mb$avg_u, nclass=20)

p <- ggplot(data=mmap_1kb, aes(x=bin_1kb, y=avg_u)) +
  geom_line() + theme_bw() + 
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="top")

save_plot("emp_mmap_1kb.png", p, base_height=10, base_width=12)

p <- ggplot(data=mmap_10kb, aes(x=bin_10kb, y=avg_u)) +
  geom_line() + theme_bw() + 
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="top")

save_plot("emp_mmap_10kb.png", p, base_height=10, base_width=12)

p <- ggplot(data=mmap_100kb, aes(x=bin_100kb, y=avg_u)) +
  geom_line(linewidth=4) + theme_bw() + 
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="top")

save_plot("emp_mmap_100kb.png", p, base_height=10, base_width=12)
