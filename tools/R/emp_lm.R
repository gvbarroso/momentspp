# Date created: 04/03/2024
# Last modified: 08/03/2024
# Author: Gustavo V. Barroso

pdf(NULL) # to suppress creation of Rplots.pdf

setwd("~/Data/momentspp/paper_2/de_novo/")

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(MASS)
  library(lmtest)
  library(nlme)
  library(car)
  library(data.table)
  library(ppcor)
})

print(Sys.time())
cat("Fitting linear models pi ~ u + B...\n")

chr_list <- c(1:22)

tbl_10kb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".10kb.tsv", sep="")
  tbl <- fread(paste("SeplyarskiyEtAl2023/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_10kb <- rbind.data.frame(tbl_10kb, tbl)
}
names(tbl_10kb)[5] <- "diversity"

tbl_100kb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".100kb.tsv", sep="")
  tbl <- fread(paste("SeplyarskiyEtAl2023/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_100kb <- rbind.data.frame(tbl_100kb, tbl)
}
names(tbl_100kb)[5] <- "diversity"

tbl_1Mb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".1Mb.tsv", sep="")
  tbl <- fread(paste("SeplyarskiyEtAl2023/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_1Mb <- rbind.data.frame(tbl_1Mb, tbl)
}
names(tbl_1Mb)[5] <- "diversity"

# standardizing variables to help interpretation of linear coefficients
## 10 kb
x10kb <- na.omit(tbl_10kb)
x10kb$cum_bin <- 1:nrow(x10kb)
x10kb$std_sub <- (x10kb$sub-mean(x10kb$sub)) / sd(x10kb$sub)
x10kb$std_pi <- (x10kb$diversity-mean(x10kb$diversity)) / sd(x10kb$diversity)
x10kb$std_mu <- (x10kb$mu-mean(x10kb$mu)) / sd(x10kb$mu)
x10kb$std_r <- (x10kb$r-mean(x10kb$r)) / sd(x10kb$r)

## 100 kb
x100kb <- na.omit(tbl_100kb)
x100kb$cum_bin <- 1:nrow(x100kb)
x100kb$std_sub <- (x100kb$sub-mean(x100kb$sub)) / sd(x100kb$sub)
x100kb$std_pi <- (x100kb$diversity-mean(x100kb$diversity)) / sd(x100kb$diversity)
x100kb$std_mu <- (x100kb$mu-mean(x100kb$mu)) / sd(x100kb$mu)
x100kb$std_r <- (x100kb$r-mean(x100kb$r)) / sd(x100kb$r)

## 1 Mb
x1Mb <- na.omit(tbl_1Mb)
x1Mb$cum_bin <- 1:nrow(x1Mb)
x1Mb$std_sub <- (x1Mb$sub-mean(x1Mb$sub)) / sd(x1Mb$sub)
x1Mb$std_pi <- (x1Mb$diversity-mean(x1Mb$diversity)) / sd(x1Mb$diversity)
x1Mb$std_mu <- (x1Mb$mu-mean(x1Mb$mu)) / sd(x1Mb$mu)
x1Mb$std_r <- (x1Mb$r-mean(x1Mb$r)) / sd(x1Mb$r)

m_1Mb <- pivot_longer(x1Mb, cols=c("mu", "r", "diversity", "sub"), names_to="var")
ggplot(data=filter(m_1Mb, var=="diversity"),
       aes(x=cum_bin, y=value, color=as.factor(chr))) +
geom_point() + theme_bw() + geom_line() + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  labs(title="Diversity across human chromosomes (1Mb windows)",
       x="Cummulative Bin", y=expression(pi)) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="none")

# diversity
gls.pi.100kb.genome <- gls(std_pi ~ std_mu + std_r, data=x100kb,
                         cor=corAR1(0, ~cum_bin|chr),
                         weights=varPower(0, ~chr),
                         method="ML")
summary(gls.pi.100kb.genome)
vif(gls.pi.100kb.genome)

gls.pi.1Mb.genome <- gls(std_pi ~ std_mu + std_r, data=x1Mb,
                    cor=corAR1(0, ~cum_bin|chr),
                    weights=varPower(0, ~chr),
                    method="ML")
summary(gls.pi.1Mb.genome)
vif(gls.pi.1Mb.genome)

# divergence
gls.sub.100kb.genome <- gls(std_sub ~ std_mu + std_r, data=x100kb,
                           cor=corAR1(0, ~cum_bin|chr),
                           weights=varPower(0, ~chr),
                           method="ML")
summary(gls.sub.100kb.genome)
vif(gls.sub.100kb.genome)

gls.sub.1Mb.genome <- gls(std_sub ~ std_mu + std_r, data=x1Mb,
                         cor=corAR1(0, ~cum_bin|chr),
                         weights=varPower(0, ~chr),
                         method="ML")
summary(gls.sub.1Mb.genome)
vif(gls.sub.1Mb.genome)

