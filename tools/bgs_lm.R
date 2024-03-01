#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(R.utils)
library(tidyverse)
library(MASS)
library(lmtest)
library(nlme)
library(car)
library(data.table)

# loads results of bgs_pi.R
maps_1kb <- fread("maps_1kb.csv")
maps_10kb <- fread("maps_10kb.csv")
maps_100kb <- fread("maps_100kb.csv")

hrmap_1kb <- fread("hrmap_1kb.csv")
hrmap_10kb <- fread("hrmap_10kb.csv")
hrmap_100kb <- fread("hrmap_100kb.csv")

# trimming tails of hrmap to match rate maps
maps_1kb$avg_pi <- hrmap_1kb$avg_pi[1:nrow(maps_1kb)] 
maps_10kb$avg_pi <- hrmap_10kb$avg_pi[1:nrow(maps_10kb)]
maps_100kb$avg_pi <- hrmap_100kb$avg_pi[1:nrow(maps_100kb)] 

# table for storing R^2 values from linear models (rows are bin sizes)
r2_tbl <- as.data.frame(matrix(ncol=5, nrow=3))
names(r2_tbl) <- c("Total", "u", "r", "s", "r:s")

# standardizing variables to help interpretation of linear coefficients
tbl_1kb <- dplyr::select(maps_1kb, c(avg_rec, avg_mut, avg_s, avg_pi))
std_1kb <- as.data.frame(apply(tbl_1kb, 2, function(x) (x-mean(x)) / sd(x)))

tbl_10kb <- dplyr::select(maps_10kb, c(avg_rec, avg_mut, avg_s, avg_pi))
std_10kb <- as.data.frame(apply(tbl_10kb, 2, function(x) (x-mean(x)) / sd(x)))

tbl_100kb <- dplyr::select(maps_100kb, c(avg_rec, avg_mut, avg_s, avg_pi))
std_100kb <- as.data.frame(apply(tbl_100kb, 2, function(x) (x-mean(x)) / sd(x)))

m_1kb_1 <- lm(avg_pi ~ avg_mut + avg_rec + avg_s + avg_rec:avg_s, data=std_1kb)
m_1kb_2 <- lm(avg_pi ~ avg_mut + avg_rec, data=std_1kb)

summary(m_1kb_1)
summary(m_1kb_2)

AIC(m_1kb_1, m_1kb_2)

anova.pi <- Anova(m_1kb_1)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl$Total[1] <- (anova.pi$VarExp[1] +
                    anova.pi$VarExp[2] +
                    anova.pi$VarExp[3] +
                    anova.pi$VarExp[4]) * 100
r2_tbl$u[1] <- anova.pi$VarExp[1] * 100
r2_tbl$r[1] <- anova.pi$VarExp[2] * 100
r2_tbl$s[1] <- anova.pi$VarExp[3] * 100
r2_tbl$`r:s`[1] <- anova.pi$VarExp[4] * 100

m_10kb_1 <- lm(avg_pi~ avg_mut + avg_rec + avg_s + avg_rec:avg_s, data=std_10kb)
m_10kb_2 <- lm(avg_pi~ avg_mut + avg_rec, data=std_10kb)

summary(m_10kb_1)
summary(m_10kb_2)

AIC(m_10kb_1, m_10kb_2) 

anova.pi <- Anova(m_10kb_1)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl$Total[2] <- (anova.pi$VarExp[1] + 
                    anova.pi$VarExp[2] +
                    anova.pi$VarExp[3] +
                    anova.pi$VarExp[4]) * 100
r2_tbl$u[2] <- anova.pi$VarExp[1] * 100
r2_tbl$r[2] <- anova.pi$VarExp[2] * 100
r2_tbl$s[2] <- anova.pi$VarExp[3] * 100
r2_tbl$`r:s`[2] <- anova.pi$VarExp[4] * 100

m_100kb_1 <- lm(avg_pi~ avg_mut + avg_rec +avg_s + avg_rec:avg_s,data=std_100kb)
m_100kb_2 <- lm(avg_pi~ avg_mut + avg_rec, data=std_100kb)

summary(m_100kb_1)
summary(m_100kb_2)

AIC(m_100kb_1, m_100kb_2) 

anova.pi <- Anova(m_100kb_1)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl$Total[3] <- (anova.pi$VarExp[1] + 
                    anova.pi$VarExp[2] +
                    anova.pi$VarExp[3] +
                    anova.pi$VarExp[4]) * 100
r2_tbl$u[3] <- anova.pi$VarExp[1] * 100
r2_tbl$r[3] <- anova.pi$VarExp[2] * 100
r2_tbl$s[3] <- anova.pi$VarExp[3] * 100
r2_tbl$`r:s`[3] <- anova.pi$VarExp[4] * 100

fwrite(r2_tbl, "r2_tbl.csv")
