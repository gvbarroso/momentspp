#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(R.utils)
library(tidyverse)
library(cowplot)
library(MASS)
library(lmtest)
library(nlme)
library(car)
library(data.table)
library(scales)

models <- fread("models.csv.gz")
num_models <- 36 # nrow(models)
num_reps <- 10 # as.numeric(args[1])

# tables for storing R^2 values from linear models
r2_tbl_1kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))
r2_tbl_10kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))
r2_tbl_100kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))

names(r2_tbl_1kb) <- c("Total", "u", "r", "s", "r:s")
names(r2_tbl_10kb) <- c("Total", "u", "r", "s", "r:s")
names(r2_tbl_100kb) <- c("Total", "u", "r", "s", "r:s")

for(i in 1:num_models) {
  # standardizing variables to help interpretation of linear coefficients
  maps_1kb$std_s <- (maps_1kb$avg_s - mean(maps_1kb$avg_s)) / 
                    sd(maps_1kb$avg_s)
  maps_1kb$std_mu <- (maps_1kb$avg_mut - mean(maps_1kb$avg_mut)) /
                     sd(maps_1kb$avg_mut)
  maps_1kb$std_r <- (maps_1kb$avg_rec - mean(maps_1kb$avg_rec)) / 
                    sd(maps_1kb$avg_rec)
  maps_1kb$std_pi <- (maps_1kb$avg_pi - mean(maps_1kb$avg_pi)) / 
                     sd(maps_1kb$avg_pi)
  
  m_1kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s, data=maps_1kb)
  m_1kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_1kb)
  
  summary(m_1kb_1)
  summary(m_1kb_2)
  
  comp <- AIC(m_1kb_1, m_1kb_2) # model comparison
  
  anova.pi <- Anova(m_1kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_1kb$Total[i] <- (anova.pi$VarExp[1] +
                            anova.pi$VarExp[2] +
                            anova.pi$VarExp[3] +
                            anova.pi$VarExp[4]) * 100
  r2_tbl_1kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_1kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_1kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_1kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  
  # builds maps at larger genomic scales
  maps_1kb$bin_10kb <- (maps_1kb$bin - 1) %/% 10
  maps_1kb$bin_100kb <- (maps_1kb$bin - 1) %/% 100
  maps_10kb <- maps_1kb %>%
    group_by(bin_10kb) %>%
    summarise_at(c("std_pi", "std_s", "std_mu", "std_r",
                   "avg_pi", "avg_s", "avg_mut", "avg_rec"), mean)
  
  maps_100kb <- maps_1kb %>%
    group_by(bin_100kb) %>%
    summarise_at(c("std_pi", "std_s", "std_mu", "std_r",
                   "avg_pi", "avg_s", "avg_mut", "avg_rec"), mean)
  
  # 10 kb
  m_10kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s, data=maps_10kb)
  m_10kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_10kb)
  
  summary(m_10kb_1)
  summary(m_10kb_2)
  
  AIC(m_10kb_1, m_10kb_2) # model comparison
  
  anova.pi <- Anova(m_10kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_10kb$Total[i] <- (anova.pi$VarExp[1] + 
                             anova.pi$VarExp[2] +
                             anova.pi$VarExp[3] +
                             anova.pi$VarExp[4]) * 100
  r2_tbl_10kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_10kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_10kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_10kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  
  # 100 kb
  m_100kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s,data=maps_100kb)
  m_100kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_100kb)
  
  summary(m_100kb_1)
  summary(m_100kb_2)
  
  AIC(m_100kb_1, m_100kb_2) # model comparison
  
  anova.pi <- Anova(m_100kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_100kb$Total[i] <- (anova.pi$VarExp[1] + 
                              anova.pi$VarExp[2] +
                              anova.pi$VarExp[3] +
                              anova.pi$VarExp[4]) * 100
  r2_tbl_100kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_100kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_100kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_100kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  r2_mdl_1kb <- fread(paste("model_", i, "/r2_1kb.csv", sep=""))
  r2_mdl_10kb <- fread(paste("model_", i, "/r2_10kb.csv", sep=""))
  r2_mdl_100kb <- fread(paste("model_", i, "/r2_100kb.csv", sep=""))
  
  tbl <- rbind.data.frame(r2_mdl_1kb, r2_mdl_10kb, r2_mdl_100kb)
  tbl$model <- i
  tbl$rep <- rep(1:10, 3)
  tbl$scale <- c(rep(1, 10), rep(100, 10), rep(1000, 10)) # in kb
  
  r2_tbl <- rbind.data.frame(r2_tbl, tbl)
}

r2s <- r2_tbl %>% group_by(model, scale) %>% 
                  summarize_at(vars(Total, u, r, s, `r:s`), 
                                    list(avg=mean, sd=sd))

for(i in 1:num_models) {
  molten_avgs <- pivot_longer(r2s, cols=ends_with("avg"), names_to="variable")
  p <- ggplot(data=filter(molten_avgs, model==i), 
              aes(x=scale, y=value, shape=variable)) +
       geom_line() + theme_bw() + geom_point(size=4) +
       scale_shape_manual(values=c(0, 8, 2, 3, 1), name=NULL) +
       scale_x_continuous(breaks=unique(molten_avgs$scale), trans="log10") +
       scale_y_continuous(breaks=pretty_breaks()) +
    labs(title=paste("Model", i, "(Mean over replicates)"),
         x="Map Scale (kb)", y="Variance Explained (%)") +
       theme(axis.title=element_text(size=16), 
             axis.text=element_text(size=12), 
             axis.text.x=element_text(size=12),
             legend.text=element_text(size=16),
             legend.position="bottom")
  
  molten_sds <- pivot_longer(r2s, cols=ends_with("sd"), names_to="variable")
  q <- ggplot(data=filter(molten_sds, model==i), 
              aes(x=scale, y=value, shape=variable)) +
    geom_line() + theme_bw() + geom_point(size=4) +
    scale_shape_manual(values=c(0, 8, 2, 3, 1), name=NULL) +
    scale_x_continuous(breaks=unique(molten_sds$scale), trans="log10") +
    scale_y_continuous(breaks=pretty_breaks()) +
    labs(title=paste("Model", i, "(SD over replicates)"),
         x="Map Scale (kb)", y=NULL) +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16), 
          legend.position="bottom")
  
  save_plot(paste("model_", i, "_r2.png", sep=""), plot_grid(p, q, nrow=1), 
            base_height=10, base_width=15)
}
  
fwrite(r2_tbl_1kb, "r2_1kb.csv")
fwrite(r2_tbl_10kb, "r2_10kb.csv")
fwrite(r2_tbl_100kb, "r2_100kb.csv")

# the plots above show us that by and large, sd over replicates is negligible
# and u and s are the only meaningful variables affecting pi in these models
# thus we simplify:
r2s_f <- dplyr::select(r2s, -starts_with("r")) %>% dplyr::select(., -ends_with("sd"))
models_f <- dplyr::select(models, c(scale_sel, shapes_mu, shapes_rec, scales_rec, model))
r2_models <- merge(models_f, r2s_f, by="model")
m_r2_models <- pivot_longer(r2_models, cols=ends_with("avg"), names_to="variable")

p <- ggplot(data=filter(m_r2_models, shapes_rec==1), # focus shape_rec==1
            aes(x=scale, y=value, color=as.factor(scale_sel), shape=variable)) +
  geom_line() + geom_point(size=4) + theme_bw() + 
  facet_grid(+scales_rec~shapes_mu) +
  scale_shape_manual(values=c(0, 1, 2), name=NULL) +
  scale_color_discrete(name="Mean(s)", type=c("plum3", "seagreen3")) +
  scale_x_continuous(breaks=unique(m_r2_models$scale), trans="log10") +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=paste("Models with r ~ exp(), mean", expression(R^2),
                   "over replicates"),
       x="Map Scale (kb)", y="Variance Explained (%)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("r2_plot_clean.png", p, base_height=8, base_width=10)

q <- ggplot(data=filter(m_r2_models, shapes_rec==1, variable=="u_avg"),
            aes(x=shapes_mu, y=value, color=as.factor(scale_sel),
                shape=as.factor(scales_rec))) +
  geom_line() + geom_point(size=4) + theme_bw() + 
  facet_wrap(~scale) +
  scale_shape_manual(values=c(0, 1, 2), name="Mean(r)") +
  scale_color_discrete(name="Mean(s)", type=c("plum3", "seagreen3")) +
  scale_x_continuous(breaks=unique(m_r2_models$shapes_mu)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=paste("Models with r ~ exp(), mean", expression(R^2),
                   "over replicates, cols. = genomic scale (kb)"),
       x="Shape of mutation rate distribution",
       y="Variance Explained by Mutation (%)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
q
save_plot("r2_u_only.png", q, base_height=8, base_width=16)

# visualization of shapes of mutation rate distributions
mut_spans <- rgeom(n=100, prob=1/models$avg_mut_spans[1])
mut_spans <- mut_spans[mut_spans>0]

avg_mu <- unique(models$scales_mu)
shapes <- unique(models$shapes_mu)
mut_rates_3 <- rgamma(n=length(mut_spans), shapes[1], rate=shapes[1]) * avg_mu
mut_rates_10 <- rgamma(n=length(mut_spans), shapes[2], rate=shapes[2]) * avg_mu
mut_rates_30 <- rgamma(n=length(mut_spans), shapes[3], rate=shapes[3]) * avg_mu


mmap <- setDT(bind_cols(dplyr::lag(cumsum(mut_spans), n=1, default=0),
                        dplyr::lag(cumsum(mut_spans), n=0, default=0)))
names(mmap) <- c("start", "end")
mmap$`3` <- mut_rates_3
mmap$`10` <- mut_rates_10
mmap$`30` <- mut_rates_30

molten_map <- pivot_longer(mmap, cols=c("3", "10", "30"), names_to="shape")
p <- ggplot(data=molten_map, aes(x=start/1e+3, y=value,
                                 color=as.factor(as.numeric(shape)))) +
     geom_step() + theme_bw() +
     scale_color_discrete(name="Shape",
                          type=c("plum3", "seagreen3", "salmon3")) +
     scale_x_continuous(breaks=pretty_breaks()) +
     scale_y_continuous(breaks=pretty_breaks(), trans="log10") +
     labs(title="Mutation landscapes for different shape parameters",
          x="Position (kb)", y=expression(mu)) +
     theme(axis.title=element_text(size=16), 
           axis.text=element_text(size=12), 
           axis.text.x=element_text(size=16),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           legend.position="bottom")
p

q <- ggplot(data=molten_map, aes(x=value, fill=as.factor(as.numeric(shape)))) + 
     geom_density(alpha=0.5) + theme_bw() +
     scale_fill_discrete(name="Shape",
                          type=c("plum3", "seagreen3", "salmon3")) +
     scale_x_continuous(trans="log10") +
     scale_y_continuous(breaks=pretty_breaks()) +
     labs(title="Density of mutation rate distributions",
          x="Rate", y="Density") +
     theme(axis.title=element_text(size=16), 
           axis.text=element_text(size=12), 
           axis.text.x=element_text(size=12),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           legend.position="bottom")
q

save_plot("mut_rates_viz.png", plot_grid(p, q, nrow=1),
          base_height=12, base_width=24)

# meta linear model
r2_models <- merge(select(models, c(model, scale_sel, shapes_mu, shapes_rec, scales_rec)),
                   select(filter(r2s, scale==1), c(model, u_avg, s_avg, r_avg, `r:s_avg`)),
                   by="model")
std_r2 <- as.data.frame(apply(r2_models[,2:ncol(r2_models)], 2, 
                              function(x) (x-mean(x)) / sd(x)))

m_mu <- lm(u_avg ~ scale_sel + shapes_mu + shapes_rec + scales_rec, data=std_r2)
m_s <- lm(s_avg ~ scale_sel + shapes_mu + shapes_rec + scales_rec, data=std_r2)
m_r <- lm(r_avg ~ scale_sel + shapes_mu + shapes_rec + scales_rec, data=std_r2)

summary(m_mu)
summary(m_s)
summary(m_r)