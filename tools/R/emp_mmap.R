
suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
})

setwd("~/Data/momentspp/paper_2/")

mmap_chr1 <- fread("de_novo/SeplyarskiyEtAl2023/VCFs/1_rate_v5.2_TFBS_correction_all.vcf.bgz",
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

mmap_chr1$bin_100 <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+2
mmap_chr1$bin_1kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+3
mmap_chr1$bin_10kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+4
mmap_chr1$bin_100kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+5
mmap_chr1$bin_1Mb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+6

mmap_100 <- mmap_chr1 %>% group_by(bin_100) %>% mutate(avg_u=mean(u))
mmap_100 <- distinct(mmap_100, bin_100, .keep_all=T)
hist(mmap_100$avg_u, nclass=100)

res <- fft(mmap_100$u - mean(mmap_100$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

mmap_1kb <- mmap_chr1 %>% group_by(bin_1kb) %>% mutate(avg_u=mean(u))
mmap_1kb <- distinct(mmap_1kb, bin_1kb, .keep_all=T)
hist(mmap_1kb$avg_u, nclass=100)

res <- fft(mmap_1kb$u - mean(mmap_1kb$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

y1kb <- mmap_1kb$avg_u / mean(mmap_1kb$avg_u)
MASS::fitdistr(y1kb, "gamma")

mmap_10kb <- mmap_chr1 %>% group_by(bin_10kb) %>% mutate(avg_u=mean(u))
mmap_10kb <- distinct(mmap_10kb, bin_10kb, .keep_all=T)
hist(mmap_10kb$avg_u, nclass=100)

res <- fft(mmap_10kb$u - mean(mmap_10kb$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

y10kb <- mmap_10kb$avg_u / mean(mmap_10kb$avg_u)
MASS::fitdistr(y10kb, "gamma")

mmap_100kb <- mmap_chr1 %>% group_by(bin_100kb) %>% mutate(avg_u=mean(u))
mmap_100kb <- distinct(mmap_100kb, bin_100kb, .keep_all=T)
hist(mmap_100kb$avg_u, nclass=20)

y100kb <- mmap_100kb$avg_u / mean(mmap_100kb$avg_u)
MASS::fitdistr(y100kb, "gamma")

mmap_1Mb <- mmap_chr1 %>% group_by(bin_1Mb) %>% mutate(avg_u=mean(u))
mmap_1Mb <- distinct(mmap_1Mb, bin_1Mb, .keep_all=T)
hist(mmap_1Mb$avg_u, nclass=20)

y1Mb <- mmap_1Mb$avg_u / mean(mmap_1Mb$avg_u)
MASS::fitdistr(y1Mb, "gamma")

p <- ggplot(data=mmap_1kb, aes(x=bin_1kb, y=avg_u)) +
  geom_line() + theme_bw() + 
  geom_hline(aes(yintercept=mean(mmap_1kb$avg_u), color="red")) +
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="top")

save_plot("emp_mmap_1kb.png", p, base_height=10, base_width=12)

p <- ggplot(data=mmap_10kb, aes(x=bin_10kb, y=avg_u)) +
  geom_line() + theme_bw() + 
  geom_hline(aes(yintercept=mean(mmap_10kb$avg_u), color="red")) +
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="none")

save_plot("emp_mmap_10kb.png", p, base_height=10, base_width=12)

p <- ggplot(data=mmap_100kb, aes(x=bin_100kb, y=avg_u)) +
  geom_line() + theme_bw() + 
  geom_hline(aes(yintercept=mean(mmap_100kb$avg_u), color="red")) +
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="none")

save_plot("emp_mmap_100kb.png", p, base_height=10, base_width=12)

p <- ggplot(data=mmap_1Mb, aes(x=bin_1Mb, y=avg_u)) +
  geom_line() + theme_bw() + 
  geom_hline(aes(yintercept=mean(mmap_1Mb$avg_u), color="red")) +
  labs("Empirical mutation map", x="Bin", y=expression(mu)) + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="none")

save_plot("emp_mmap_1Mb.png", p, base_height=10, base_width=12)

