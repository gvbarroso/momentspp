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

chr_list <- c(17:18, 21:22)

tbl_100kb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".100kb.tsv", sep="")
  tbl <- fread(paste("primate_data/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_100kb <- rbind.data.frame(tbl, tbl_100kb)
}

tbl_1Mb <- data.frame()
for(i in chr_list) {
  name <- paste("chr", i, ".1Mb.tsv", sep="")
  tbl <- fread(paste("primate_data/", name, sep=""))
  tbl$bin <- 1:nrow(tbl)
  tbl$chr <- i
  tbl_1Mb <- rbind.data.frame(tbl, tbl_1Mb)
}

# table for storing R^2 values from linear models (rows are bin sizes)
r2_tbl <- as.data.frame(matrix(ncol=3, nrow=3))
names(r2_tbl) <- c("Total", "u", "B")
r2_tbl$scale <- c(1e+4, 1e+5, 1e+6)

# NOTE linear models technically cannot be fit when the mutation map is flat 
# because in this case pi==B (and the residual sum of squares is zero)

# standardizing variables to help interpretation of linear coefficients
std_100kb <- dplyr::select(tbl_100kb, c(mu, r, pi))
names(std_100kb)[3] <- "diversity"
std_100kb <- as.data.frame(apply(std_100kb[,1:ncol(std_100kb)], 2,
                function(x) (x-mean(x, na.rm=T)) / sd(x, na.rm=T)))
std_100kb$bin <- tbl_100kb$bin
std_100kb$chr <- tbl_100kb$chr

# standardizing variables to help interpretation of linear coefficients
std_1Mb <- dplyr::select(tbl_1Mb, c(mu, r, pi))
names(std_1Mb)[3] <- "diversity"
std_1Mb <- as.data.frame(apply(std_1Mb[,1:ncol(std_1Mb)], 2,
                         function(x) (x-mean(x, na.rm=T)) / sd(x, na.rm=T)))
std_1Mb$bin <- tbl_1Mb$bin
std_1Mb$chr <- tbl_1Mb$chr

m.pi.1Mb.1 <- lm(diversity ~ (mu + r + chr)^2, data=std_1Mb)
summary(m.pi.1Mb.1)

m.pi.100kb.1 <- lm(diversity ~ (mu + r + chr)^2, data=std_100kb)
summary(m.pi.100kb.1)

m.pi.100kb.2 <- lm(diversity ~ (r + chr)^2, data=std_100kb)
summary(m.pi.100kb.2)

plot(resid(m.pi.100kb.1)~fitted(m.pi.100kb.1))
hist(resid(m.pi.100kb.1), nclass = 30)
dwtest(m.pi.100kb.1)
hmctest(m.pi.100kb.1, nsim=10000) 

g.pi.100kb.1 <- gls(diversity ~ (mu + r + mu:r), data=std_100kb, 
                    cor=corAR1(0, ~bin|chr), method="ML", na.action=na.omit)

summary(g.pi.100kb.1)
vif(g.pi.100kb.1)

g.pi.100kb.2 <- gls(diversity ~ (mu + r + mu:r), data=std_100kb, 
                  weights=varPower(0, ~mu|chr), cor=corAR1(0, ~bin|chr), 
                  method="ML", na.action=na.omit)

summary(g.pi.100kb.2)
vif(g.pi.100kb.2)

AIC(g.pi.100kb.1, g.pi.100kb.2)
AIC(m.pi.100kb, g.pi.100kb.1)

anova.pi <- Anova(m_1kb)
apiss <- anova.pi$"Sum Sq"
anova.pi$VarExp <- apiss / sum(apiss)

r2_tbl$Total[1] <- (anova.pi$VarExp[1] + anova.pi$VarExp[2]) * 100
r2_tbl$u[1] <- anova.pi$VarExp[1] * 100
r2_tbl$B[1] <- anova.pi$VarExp[2] * 100


m.pi.100kb <- lm(diversity ~ (mu + r)^2, data=filter(std_100kb, chr==17))
summary(m.pi.100kb)

g.pi.100kb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_100kb, chr==17), 
                   cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.100kb.3)
vif(g.pi.100kb.3)

m.pi.100kb <- lm(diversity ~ (mu + r)^2, data=filter(std_100kb, chr==18))
summary(m.pi.100kb)

g.pi.100kb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_100kb, chr==18), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.100kb.3)
vif(g.pi.100kb.3)

m.pi.100kb <- lm(diversity ~ (mu + r)^2, data=filter(std_100kb, chr==21))
summary(m.pi.100kb)

g.pi.100kb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_100kb, chr==21), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.100kb.3)
vif(g.pi.100kb.3)

m.pi.100kb <- lm(diversity ~ (mu + r)^2, data=filter(std_100kb, chr==22))
summary(m.pi.100kb)

g.pi.100kb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_100kb, chr==22), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.100kb.3)
vif(g.pi.100kb.3)


#
m.pi.1Mb <- lm(diversity ~ (mu + r)^2, data=filter(std_1Mb, chr==17))
summary(m.pi.1Mb)

g.pi.1Mb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_1Mb, chr==17), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.1Mb.3)
vif(g.pi.1Mb.3)

m.pi.1Mb <- lm(diversity ~ (mu + r)^2, data=filter(std_1Mb, chr==18))
summary(m.pi.1Mb)

g.pi.1Mb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_1Mb, chr==18), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.1Mb.3)
vif(g.pi.1Mb.3)

m.pi.1Mb <- lm(diversity ~ (mu + r)^2, data=filter(std_1Mb, chr==21))
summary(m.pi.1Mb)

g.pi.1Mb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_1Mb, chr==21), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.1Mb.3)
vif(g.pi.1Mb.3)

m.pi.1Mb <- lm(diversity ~ (mu + r)^2, data=filter(std_1Mb, chr==22))
summary(m.pi.1Mb)

g.pi.1Mb.3 <- gls(diversity ~ (mu + r + mu:r), data=filter(std_1Mb, chr==22), 
                    cor=corAR1(0, ~bin), method="ML", na.action=na.omit)

summary(g.pi.1Mb.3)
vif(g.pi.1Mb.3)
