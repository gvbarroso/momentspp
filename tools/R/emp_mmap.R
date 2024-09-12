
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

gamma_fitted_emp <- as.data.frame(matrix(ncol=2, nrow=6))
names(gamma_fitted_emp) <- c("estimate", "bin_size")

mmap_chr1 <- dplyr::select(mmap_chr1, c(POS, mu))
mmap_chr1 <- mmap_chr1 %>% group_by(POS) %>% mutate(u=sum(mu))
mmap_chr1 <- distinct(mmap_chr1, POS, .keep_all=T)
mmap_chr1 <- dplyr::select(mmap_chr1, c(POS, u))
hist(mmap_chr1$u, nclass=100)

gamma_1bp <- MASS::fitdistr(mmap_chr1$u / mean(mmap_chr1$u), "gamma")
# 1.6297709825   1.6297708331 
gamma_fitted_emp[1,] <- c(gamma_1bp$estimate[1], 1)

mmap_chr1$bin_100 <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+2
mmap_chr1$bin_1kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+3
mmap_chr1$bin_10kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+4
mmap_chr1$bin_100kb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+5
mmap_chr1$bin_1Mb <- (mmap_chr1$POS - mmap_chr1$POS[1]) %/% 1e+6

mmap_100 <- mmap_chr1 %>% group_by(bin_100) %>% mutate(avg_u=mean(u))
mmap_100 <- distinct(mmap_100, bin_100, .keep_all=T)
hist(mmap_100$avg_u, nclass=100)

gamma_100 <- MASS::fitdistr(mmap_100$avg_u / mean(mmap_100$avg_u), "gamma")
# 12.138775750   12.138774581 
gamma_fitted_emp[2,] <- c(gamma_100$estimate[1], 100)

res <- fft(mmap_100$u - mean(mmap_100$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

mmap_1kb <- mmap_chr1 %>% group_by(bin_1kb) %>% mutate(avg_u=mean(u))
mmap_1kb <- distinct(mmap_1kb, bin_1kb, .keep_all=T)
hist(mmap_1kb$avg_u, nclass=100)

gamma_1kb <- MASS::fitdistr(mmap_1kb$avg_u / mean(mmap_1kb$avg_u), "gamma")
# 26.9325749   26.9325712  
gamma_fitted_emp[3,] <- c(gamma_1kb$estimate[1], 1000)

res <- fft(mmap_1kb$u - mean(mmap_1kb$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

mmap_10kb <- mmap_chr1 %>% group_by(bin_10kb) %>% mutate(avg_u=mean(u))
mmap_10kb <- distinct(mmap_10kb, bin_10kb, .keep_all=T)
hist(mmap_10kb$avg_u, nclass=100)

gamma_10kb <- MASS::fitdistr(mmap_10kb$avg_u / mean(mmap_10kb$avg_u), "gamma")
# 42.252085   42.252097 
gamma_fitted_emp[4,] <- c(gamma_10kb$estimate[1], 10000)

res <- fft(mmap_10kb$u - mean(mmap_10kb$u))
ft_tbl <- cbind.data.frame(Mod(res), 1:length(Mod(res)))
names(ft_tbl) <- c("FT", "Index")
ggplot(data=ft_tbl, aes(x=Index, y=FT)) + geom_line() + theme_classic()

mmap_100kb <- mmap_chr1 %>% group_by(bin_100kb) %>% mutate(avg_u=mean(u))
mmap_100kb <- distinct(mmap_100kb, bin_100kb, .keep_all=T)
hist(mmap_100kb$avg_u, nclass=20)

gamma_100kb <- MASS::fitdistr(mmap_100kb$avg_u / mean(mmap_100kb$avg_u), "gamma")
# 57.289101   57.289117
gamma_fitted_emp[5,] <- c(gamma_100kb$estimate[1], 100000)

mmap_1Mb <- mmap_chr1 %>% group_by(bin_1Mb) %>% mutate(avg_u=mean(u))
mmap_1Mb <- distinct(mmap_1Mb, bin_1Mb, .keep_all=T)
hist(mmap_1Mb$avg_u, nclass=20)

gamma_1Mb <- MASS::fitdistr(mmap_1Mb$avg_u / mean(mmap_1Mb$avg_u), "gamma")
# 80.39765   80.39766 
gamma_fitted_emp[6,] <- c(gamma_1Mb$estimate[1], 1000000)

fwrite(gamma_fitted_emp, "gamma_fitted_emp_mmap.csv")

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

