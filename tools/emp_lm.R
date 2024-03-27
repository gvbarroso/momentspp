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
})

print(Sys.time())
cat("Fitting linear models pi ~ u + B...\n")

chr_list <- c(1:4, 17:22)

tbl_10kb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".10kb.tsv", sep="")
  tbl <- fread(paste("primate_data/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_10kb <- rbind.data.frame(tbl_10kb, tbl)
}

tbl_100kb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".100kb.tsv", sep="")
  tbl <- fread(paste("primate_data/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_100kb <- rbind.data.frame(tbl_100kb, tbl)
}

tbl_1Mb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".1Mb.tsv", sep="")
  tbl <- fread(paste("primate_data/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_1Mb <- rbind.data.frame(tbl_1Mb, tbl)
}

# standardizing variables to help interpretation of linear coefficients
## 10 kb
std_10kb <- dplyr::select(tbl_10kb, c(mu, r, pi, sub))
names(std_10kb)[3] <- "diversity"
std_10kb <- as.data.frame(apply(std_10kb[,1:ncol(std_10kb)], 2,
              function(x) (x-mean(x, na.rm=T)) / sd(x, na.rm=T)))
std_10kb$bin <- tbl_10kb$bin
std_10kb$chr <- tbl_10kb$chr

## 100 kb
std_100kb <- dplyr::select(tbl_100kb, c(mu, r, pi, sub))
names(std_100kb)[3] <- "diversity"
std_100kb <- as.data.frame(apply(std_100kb[,1:ncol(std_100kb)], 2,
                function(x) (x-mean(x, na.rm=T)) / sd(x, na.rm=T)))
std_100kb$bin <- tbl_100kb$bin
std_100kb$chr <- tbl_100kb$chr

## 1 Mb
std_1Mb <- dplyr::select(tbl_1Mb, c(mu, r, pi, sub))
names(std_1Mb)[3] <- "diversity"
std_1Mb <- as.data.frame(apply(std_1Mb[,1:ncol(std_1Mb)], 2,
              function(x) (x-mean(x, na.rm=T)) / sd(x, na.rm=T)))
std_1Mb$bin <- tbl_1Mb$bin
std_1Mb$chr <- tbl_1Mb$chr

# linear models
## tables for storing coefficients from linear models
beta_pi_10kb <- as.data.frame(matrix(ncol=2, nrow=length(chr_list)))
names(beta_pi_10kb) <- c("r", "u")
beta_pi_10kb$chr <- chr_list
beta_sub_10kb <- beta_pi_10kb

beta_pi_100kb <- as.data.frame(matrix(ncol=2, nrow=length(chr_list)))
names(beta_pi_100kb) <- c("r", "u")
beta_pi_100kb$chr <- chr_list
beta_sub_100kb <- beta_pi_100kb

beta_pi_1Mb <- as.data.frame(matrix(ncol=2, nrow=length(chr_list)))
names(beta_pi_1Mb) <- c("r", "u")
beta_pi_1Mb$chr <- chr_list
beta_sub_1Mb <- beta_pi_1Mb

pb <- txtProgressBar(min=1, max=length(chr_list) + 1, style=3)
for(i in 1:length(chr_list)) {
  
  setTxtProgressBar(pb, i)
  
  # linear models on pi
  #gls.pi.10kb <- gls(diversity ~ mu + r, data=filter(std_10kb, chr==chr_list[i]),
  #                  cor=corAR1(0, ~bin), method="ML", na.action=na.omit)
  #beta_pi_10kb$u[i] <- gls.pi.10kb$coefficients[2]
  #beta_pi_10kb$r[i] <- gls.pi.10kb$coefficients[3]
  
  gls.pi.100kb <- gls(diversity ~ mu + r, data=na.omit(filter(std_100kb, chr==chr_list[i])),
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit, verbose=T)
  beta_pi_100kb$u[i] <- gls.pi.100kb$coefficients[2]
  beta_pi_100kb$r[i] <- gls.pi.100kb$coefficients[3]
  
  gls.pi.1Mb <- gls(diversity ~ mu + r, data=filter(std_1Mb, chr==chr_list[i]),
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)
  beta_pi_1Mb$u[i] <- gls.pi.1Mb$coefficients[2]
  beta_pi_1Mb$r[i] <- gls.pi.1Mb$coefficients[3]
  
  # linear models on divergence
  #gls.sub.10kb <- gls(sub ~ mu + r, data=filter(std_10kb, chr==chr_list[i]),
  #                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)
  #beta_sub_10kb$u[i] <- gls.sub.10kb$coefficients[2]
  #beta_sub_10kb$r[i] <- gls.sub.10kb$coefficients[3]
  
  gls.sub.100kb <- gls(sub ~ mu + r, data=filter(std_100kb, chr==chr_list[i]),
                       cor=corAR1(0, ~bin), method="ML", na.action=na.omit)
  beta_sub_100kb$u[i] <- gls.sub.100kb$coefficients[2]
  beta_sub_100kb$r[i] <- gls.sub.100kb$coefficients[3]
  
  gls.sub.1Mb <- gls(sub ~ mu + r, data=filter(std_1Mb, chr==chr_list[i]),
                     cor=corAR1(0, ~bin), method="ML", na.action=na.omit)
  beta_sub_1Mb$u[i] <- gls.sub.1Mb$coefficients[2]
  beta_sub_1Mb$r[i] <- gls.sub.1Mb$coefficients[3]
}
close(pb)
