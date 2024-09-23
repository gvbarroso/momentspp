
library(tidyverse)
library(MASS)
library(ppcor)
library(car)
library(scales)


# when z = x * y, what can linear regression tell us about these (pairwise) relationships?
s <- 1e+6 # sample size

# parameters of our toy models
means_u <- 100
means_B <- c(100, 200, 500)
sds_u <- c(1, 2, 5)
sds_B <- c(1, 2, 5)
  
params <- crossing(means_u, means_B, sds_u, sds_B)
n_models <- nrow(params)
params$model <- 1:nrow(params)

# table to store estimates
estimates <- as.data.frame(matrix(nrow=n_models, ncol=14))
names(estimates) <- c("cor_u", "cor_B", 
                      "cor_u_partial", "cor_B_partial",
                      "coeff_u_m1", "coeff_B_m1",
                      "r2_u_m1", "r2_B_m1",
                      "coeff_u_m2", "coeff_B_m2",
                      "r2_u_m2", "r2_B_m2",
                      "coeff_u_m3", "coeff_B_m3")
estimates$model <- 1:nrow(estimates)

pb <- txtProgressBar(min=0, max=n_models, style=3)
for(i in 1:n_models) {

  setTxtProgressBar(pb, i)
  # build toy model
  u <- rnorm(s, params$means_u[i], params$sds_u[i])
  B <- rnorm(s, params$means_B[i], params$sds_B[i])
  p <- u * B + rnorm(s, 0, 1) # adds noise to avoid perfect fit
  
  t1 <- cor.test(p, u)
  t2 <- cor.test(p, B)
  
  estimates$cor_u[i] <- t1$estimate
  estimates$cor_B[i] <- t2$estimate
  
  t3 <- pcor.test(p, u, B)
  t4 <- pcor.test(p, B, u)
  
  estimates$cor_u_partial[i] <- t3$estimate
  estimates$cor_B_partial[i] <- t4$estimate
  
  tbl <- cbind.data.frame(p, u, B)
  std <- as.data.frame(apply(tbl, 2, function(x) (x-mean(x)) / sd(x)))
  
  # "naive' linear model
  m1 <- lm(p ~ u + B, tbl)
  summary(m1)
  
  estimates$coeff_u_m1[i] <- m1$coefficients[2]
  estimates$coeff_B_m1[i] <- m1$coefficients[3]
  
  anova_m1 <- Anova(m1)
  apiss <- anova_m1$"Sum Sq"
  anova_m1$VarExp <- apiss / sum(apiss)
  anova_m1
  
  estimates$r2_u_m1[i] <- anova_m1$VarExp[1]
  estimates$r2_B_m1[i] <- anova_m1$VarExp[2]
  
  # linear model on log scale
  m2 <- lm(log(p) ~ log(u) + log(B), tbl)
  summary(m2)
  
  estimates$coeff_u_m2[i] <- m2$coefficients[2]
  estimates$coeff_B_m2[i] <- m2$coefficients[3]
  
  anova_m2 <- Anova(m2)
  apiss <- anova_m2$"Sum Sq"
  anova_m2$VarExp <- apiss / sum(apiss)
  anova_m2
  
  estimates$r2_u_m2[i] <- anova_m2$VarExp[1]
  estimates$r2_B_m2[i] <- anova_m2$VarExp[2]
  
  # linear model with standardized variables
  m3 <- lm(p ~ u + B, std)
  summary(m3)
  
  estimates$coeff_u_m3[i] <- m3$coefficients[2]
  estimates$coeff_B_m3[i] <- m3$coefficients[3]
}
close(pb)

tbl <- merge(estimates, params, by="model")

tbl_cm1 <- pivot_longer(tbl, cols=c("coeff_u_m1", "coeff_B_m1"))

p1 <- ggplot(data=tbl_cm1, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Coeff. M1 (linear)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p1

tbl_r2m1 <- pivot_longer(tbl, cols=c("r2_u_m1", "r2_B_m1"))

p2 <- ggplot(data=tbl_r2m1, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Var. Explained M1 (linear)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p2

tbl_cm2 <- pivot_longer(tbl, cols=c("coeff_u_m2", "coeff_B_m2"))

p3 <- ggplot(data=tbl_cm2, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Coeff. M2 (log)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p3

tbl_r2m2 <- pivot_longer(tbl, cols=c("r2_u_m2", "r2_B_m2"))

p4 <- ggplot(data=tbl_r2m2, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Var. Explained M2 (log)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p4

tbl_cm3 <- pivot_longer(tbl, cols=c("coeff_u_m3", "coeff_B_m3"))

p5 <- ggplot(data=tbl_cm3, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Coeff. M3 (z-scores)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p5

tbl_cor <- pivot_longer(tbl, cols=c("cor_u", "cor_B"))

p6 <- ggplot(data=tbl_cor, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="Cor (p, x)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p6

tbl_cor_p <- pivot_longer(tbl, cols=c("cor_u_partial", "cor_B_partial"))

p7 <- ggplot(data=tbl_cor_p, aes(x=means_B, y=value, shape=name)) +
  geom_point(size=3) + theme_bw() + geom_line() +
  facet_grid(+sds_u~sds_B) +
  scale_shape_manual(values=c(0, 1), name=NULL, labels=c("B", "u")) +
  scale_y_log10(breaks=pretty_breaks()) +
  labs(title="rows->SD(u), cols->SD(B)", x="Mean (B)", y="P. Cor (p, x, z)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p7
