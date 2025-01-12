# Date created: 04/03/2024
# Last modified: 18/12/2024
# Author: Gustavo V. Barroso

pdf(NULL) # to suppress creation of Rplots.pdf

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(MASS)
  library(lmtest)
  library(nlme)
  library(car)
  library(data.table)
  library(interactions)
})

print(Sys.time())
cat("Fitting linear models pi ~ u + B...\n")

# loads results of bgs_pi.R
maps_1kb <- fread("tbl_1kb.csv.gz")
maps_10kb <- fread("tbl_10kb.csv.gz")
maps_100kb <- fread("tbl_100kb.csv.gz")
maps_1Mb <- fread("tbl_1Mb.csv.gz")

# tables for storing R^2 values from linear models (rows are bin sizes)
## expected pi
r2_tbl_exp <- as.data.frame(matrix(ncol=3, nrow=4))
names(r2_tbl_exp) <- c("u", "B", "B:u")
r2_tbl_exp$scale <- c(1e+3, 1e+4, 1e+5, 1e+6)

## simulated pi
r2_tbl_sim <- as.data.frame(matrix(ncol=3, nrow=4))
names(r2_tbl_sim) <- c("u", "B", "B:u")
r2_tbl_sim$scale <- c(1e+3, 1e+4, 1e+5, 1e+6)

# table for storing total R^2 values straight from the lm function (rows are bin sizes)
r2_raw <- as.data.frame(matrix(ncol=2, nrow=4))
names(r2_raw) <- c("r2_typeI_exppi", "r2_typeI_simpi")
r2_raw$scale <- c(1e+3, 1e+4, 1e+5, 1e+6)

# tables for storing coefficients from linear models (rows are bin sizes)
## expected pi
coeff_tbl_exp <- as.data.frame(matrix(ncol=3, nrow=4))
names(coeff_tbl_exp) <- c("u", "B", "B:u")
coeff_tbl_exp$scale <- c(1e+3, 1e+4, 1e+5, 1e+6)

## simulated pi
coeff_tbl_sim <- as.data.frame(matrix(ncol=3, nrow=4))
names(coeff_tbl_sim) <- c("u", "B", "B:u")
coeff_tbl_sim$scale <- c(1e+3, 1e+4, 1e+5, 1e+6)

# standardizing variables to help interpretation of linear coefficients
tbl_1kb <- dplyr::select(maps_1kb, c(avg_mut, B, avg_pi, sim_pi))
std_1kb <- as.data.frame(apply(tbl_1kb, 2, function(x) (x-mean(x)) / sd(x)))
std_1kb$bin <- 1:nrow(std_1kb)

tbl_10kb <- dplyr::select(maps_10kb, c(avg_mut, B, avg_pi, sim_pi))
std_10kb <- as.data.frame(apply(tbl_10kb, 2, function(x) (x-mean(x)) / sd(x)))
std_10kb$bin <- 1:nrow(std_10kb)

tbl_100kb <- dplyr::select(maps_100kb, c(avg_mut, B, avg_pi, sim_pi))
std_100kb <- as.data.frame(apply(tbl_100kb, 2, function(x) (x-mean(x)) / sd(x)))
std_100kb$bin <- 1:nrow(std_100kb)

tbl_1Mb <- dplyr::select(maps_1Mb, c(avg_mut, B, avg_pi, sim_pi))
std_1Mb <- as.data.frame(apply(tbl_1Mb, 2, function(x) (x-mean(x)) / sd(x)))
std_1Mb$bin <- 1:nrow(std_1Mb)

# 1 kb
## expected pi
m_1kb_exp <- lm(avg_pi ~ (B + avg_mut) ^ 2, std_1kb)
summary(m_1kb_exp)
r2_raw$r2_typeI_exppi[1] <- summary(m_1kb_exp)$r.squared
vif(m_1kb_exp, type = 'predictor')

coeff_tbl_exp$u[1] <- m_1kb_exp$coefficients[3]
coeff_tbl_exp$B[1] <- m_1kb_exp$coefficients[2]
coeff_tbl_exp$`B:u`[1] <- m_1kb_exp$coefficients[4]

anova.pi <- Anova(m_1kb_exp)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_exp$`B:u`[1] <- anova.pi$VarExp[3] * 100
r2_tbl_exp$u[1] <- anova.pi$VarExp[2] * 100
r2_tbl_exp$B[1] <- anova.pi$VarExp[1] * 100

## simulated pi with mutational variance
m_1kb_sim <- lm(sim_pi ~ (B + avg_mut) ^ 2, std_1kb)
summary(m_1kb_sim)
r2_raw$r2_typeI_simpi[1] <- summary(m_1kb_sim)$r.squared
vif(m_1kb_sim, type = 'predictor')

coeff_tbl_sim$u[1] <- m_1kb_sim$coefficients[3]
coeff_tbl_sim$B[1] <- m_1kb_sim$coefficients[2]
coeff_tbl_sim$`B:u`[1] <- m_1kb_sim$coefficients[4]

anova.pi <- Anova(m_1kb_sim)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_sim$`B:u`[1] <- anova.pi$VarExp[3] * 100
r2_tbl_sim$u[1] <- anova.pi$VarExp[2] * 100
r2_tbl_sim$B[1] <- anova.pi$VarExp[1] * 100

# 10 kb
## expected pi
m_10kb_exp <- lm(avg_pi ~ (B + avg_mut) ^ 2, std_10kb)
summary(m_10kb_exp)
r2_raw$r2_typeI_exppi[2] <- summary(m_10kb_exp)$r.squared
vif(m_10kb_exp, type = 'predictor')

coeff_tbl_exp$u[2] <- m_10kb_exp$coefficients[3]
coeff_tbl_exp$B[2] <- m_10kb_exp$coefficients[2]
coeff_tbl_exp$`B:u`[2] <- m_10kb_exp$coefficients[4]

anova.pi <- Anova(m_10kb_exp)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_exp$`B:u`[2] <- anova.pi$VarExp[3] * 100
r2_tbl_exp$u[2] <- anova.pi$VarExp[2] * 100
r2_tbl_exp$B[2] <- anova.pi$VarExp[1] * 100

## simulated pi with mutational variance
m_10kb_sim <- lm(sim_pi ~ (B + avg_mut) ^ 2, std_10kb)
summary(m_10kb_sim)
r2_raw$r2_typeI_simpi[2] <- summary(m_10kb_sim)$r.squared
vif(m_10kb_sim, type = 'predictor')

coeff_tbl_sim$u[2] <- m_10kb_sim$coefficients[3]
coeff_tbl_sim$B[2] <- m_10kb_sim$coefficients[2]
coeff_tbl_sim$`B:u`[2] <- m_10kb_sim$coefficients[4]

anova.pi <- Anova(m_10kb_sim)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_sim$`B:u`[2] <- anova.pi$VarExp[3] * 100
r2_tbl_sim$u[2] <- anova.pi$VarExp[2] * 100
r2_tbl_sim$B[2] <- anova.pi$VarExp[1] * 100

# 100 kb
## expected pi
m_100kb_exp <- lm(avg_pi ~ (B + avg_mut) ^ 2, std_100kb)
summary(m_100kb_exp)
r2_raw$r2_typeI_exppi[3] <- summary(m_100kb_exp)$r.squared
vif(m_100kb_exp, type = 'predictor')

coeff_tbl_exp$u[3] <- m_100kb_exp$coefficients[3]
coeff_tbl_exp$B[3] <- m_100kb_exp$coefficients[2]
coeff_tbl_exp$`B:u`[3] <- m_100kb_exp$coefficients[4]

anova.pi <- Anova(m_100kb_exp)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_exp$`B:u`[3] <- anova.pi$VarExp[3] * 100
r2_tbl_exp$u[3] <- anova.pi$VarExp[2] * 100
r2_tbl_exp$B[3] <- anova.pi$VarExp[1] * 100

## simulated pi with mutational variance
m_100kb_sim <- lm(sim_pi ~ (B + avg_mut) ^ 2, std_100kb)
summary(m_100kb_sim)
r2_raw$r2_typeI_simpi[3] <- summary(m_100kb_sim)$r.squared
vif(m_100kb_sim, type = 'predictor')

coeff_tbl_sim$u[3] <- m_100kb_sim$coefficients[3]
coeff_tbl_sim$B[3] <- m_100kb_sim$coefficients[2]
coeff_tbl_sim$`B:u`[3] <- m_100kb_sim$coefficients[4]

anova.pi <- Anova(m_100kb_sim)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_sim$`B:u`[3] <- anova.pi$VarExp[3] * 100
r2_tbl_sim$u[3] <- anova.pi$VarExp[2] * 100
r2_tbl_sim$B[3] <- anova.pi$VarExp[1] * 100

# 1 Mb
## expected pi
m_1Mb_exp <- lm(avg_pi ~ (B + avg_mut) ^ 2, data=std_1Mb)
summary(m_1Mb_exp)
r2_raw$r2_typeI_exppi[4] <- summary(m_1Mb_exp)$r.squared
vif(m_1Mb_exp, type = 'predictor')

coeff_tbl_exp$u[4] <- m_1Mb_exp$coefficients[3]
coeff_tbl_exp$B[4] <- m_1Mb_exp$coefficients[2]
coeff_tbl_exp$`B:u`[4] <- m_1Mb_exp$coefficients[4]

anova.pi <- Anova(m_1Mb_exp)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_exp$`B:u`[4] <- anova.pi$VarExp[3] * 100
r2_tbl_exp$u[4] <- anova.pi$VarExp[2] * 100
r2_tbl_exp$B[4] <- anova.pi$VarExp[1] * 100

## simulated pi with mutational variance
m_1Mb_sim <- lm(sim_pi ~ (B + avg_mut) ^ 2, std_1Mb)
summary(m_1Mb_sim)
r2_raw$r2_typeI_simpi[4] <- summary(m_1Mb_sim)$r.squared
vif(m_1Mb_sim, type = 'predictor')

coeff_tbl_sim$u[4] <- m_1Mb_sim$coefficients[3]
coeff_tbl_sim$B[4] <- m_1Mb_sim$coefficients[2]
coeff_tbl_sim$`B:u`[4] <- m_1Mb_sim$coefficients[4]

anova.pi <- Anova(m_1Mb_sim)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl_sim$`B:u`[4] <- anova.pi$VarExp[3] * 100
r2_tbl_sim$u[4] <- anova.pi$VarExp[2] * 100
r2_tbl_sim$B[4] <- anova.pi$VarExp[1] * 100

fwrite(coeff_tbl_sim, "coeff_tbl_sim.csv.gz")
fwrite(coeff_tbl_exp, "coeff_tbl_exp.csv.gz")

fwrite(r2_tbl_sim, "r2_tbl_sim.csv.gz")
fwrite(r2_tbl_exp, "r2_tbl_exp.csv.gz")

fwrite(r2_raw, paste("r2_raw_tbl.csv.gz", sep=""))

cat("bgs_lm.R is done!\n\n")