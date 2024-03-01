#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(R.utils)
library(tidyverse)
library(data.table)
library(GenomicRanges)

# loads the main tables
models <- fread("models.csv.gz")
num_models <- nrow(models)
num_reps <- 10

# if scripts bgs_maps.R, bgs_pi.R and bgs_lm.R have been run within rep folders:
for(i in 1:num_models) {
  for(j in 1:num_reps) {
    
    smap <- fread(paste("/model_", i, "/rep_", j, "/hrmap.csv.gz", sep=""))
    rmap <- fread(paste("/model_", i, "/rep_", j, "/rmap.csv.gz", sep=""))
    mmap <- fread(paste("/model_", i, "/rep_", j, "/mmap.csv.gz", sep=""))
    hrmap <- fread(paste("/model_", i, "/rep_", j, "/hrmap.csv.gz", sep=""))
    
    m1kb <- fread(paste("/model_", i, "/rep_", j, "/map_1kb.csv.gz", sep=""))
    m10kb <- fread(paste("/model_", i, "/rep_", j, "/map_10kb.csv.gz", sep=""))
    m100kb <- fread(paste("/model_", i, "/rep_", j, "/map_100kb.csv.gz",sep=""))
    
    r2_tbl <- fread(paste("/model_", i, "/rep_", j, "/r2_tbl.csv", sep=""))
    
    # plots single-nucleotide diversity
    seg <- filter(smap, s<0, start >= 0, end <= 1e+6)
    wsize <- seg$end[nrow(seg)] - seg$start[1]
    barcode <- ggplot(data=seg) +
      geom_segment(aes(x=start, xend=end, y=1, yend=1, size=3, color=s)) +
      theme_void() +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            legend.position="none")
    
    molten_div <- pivot_longer(hrmap, cols=c("Hr", "pi0"), names_to="Var")
    ppi <- ggplot(data=molten_div[1:seg$end[nrow(seg)],], 
                  aes(x=Pos, y=value, color=Var)) +
      geom_line(linewidth=1.05) + theme_bw() + scale_y_log10() +
      scale_color_discrete(type=c("plum3", "seagreen3"),
                           name=NULL, 
                           labels=c(expression(pi), expression(pi[0]))) +
      labs(title=paste("Diversity of a ", wsize, 
                       "-bp window of model ",i, sep=""),
           x=NULL, y="Pairwise Diversity") + 
      theme(axis.title=element_text(size=16),
            axis.text=element_text(size=12),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            legend.position="top")
    
    ppi2 <- plot_grid(ppi, barcode, ncol=1, align='v', axis='l',
                      rel_heights=c(1, 0.1))
    save_plot(paste("rep_", i, "/pis.png", sep=""), 
              ppi2, base_height=10, base_width=12) 
    
    # plots binned maps 
    pi_1kb <- ggplot(data=maps_1kb, aes(x=(bin -1) * 1e+4, y=avg_pi)) +
      geom_line(data=maps_1kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scale.4d) +
      labs(title="1 kb", x=NULL, y=expression(pi)) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=12),
            axis.title.y=element_text(size=16),
            plot.title=element_text(size=20),
            legend.position="none") 
    
    u_1kb <- ggplot(data=maps_1kb, aes(x=(bin -1) * 1e+4, y=avg_mut)) + 
      geom_line(data=maps_1kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y=expression(mu)) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=12),
            axis.title.y=element_text(size=16),
            plot.title=element_text(size=20), 
            legend.position="none")
    
    s_1kb <- ggplot(data=maps_1kb, aes(x=(bin -1) * 1e+4, y=avg_s)) + 
      geom_line(data=maps_1kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y="s") +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=12),
            axis.title.y=element_text(size=16),
            plot.title=element_text(size=20),
            legend.position="none")
    
    r_1kb <- ggplot(data=maps_1kb, aes(x=(bin -1) * 1e+4, y=avg_rec)) + 
      geom_line(data=maps_1kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x="Position", y="r") +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            axis.title.y=element_text(size=16),
            plot.title=element_text(size=20),
            legend.position="none")
    
    lands_1kb <- plot_grid(pi_1kb, u_1kb, s_1kb, r_1kb, align='v', ncol=1)
    
    pi_10kb <- ggplot(data=maps_10kb, aes(x=(bin_10kb -1) * 1e+5, y=avg_pi)) +
      geom_line(data=maps_10kb) + theme_bw() + 
      scale_y_continuous(breaks=pretty_breaks(), labels=scale.4d) +
      labs(title="10 kb", x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position = "none") 
    
    u_10kb <- ggplot(data=maps_10kb, aes(x=(bin_10kb -1) * 1e+5, y=avg_mut)) + 
      geom_line(data=maps_10kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none")
    
    s_10kb <- ggplot(data=maps_10kb, aes(x=(bin_10kb -1) * 1e+5, y=avg_s)) +
      geom_line(data=maps_10kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none")
    
    r_10kb <- ggplot(data=maps_10kb, aes(x=(bin_10kb -1) * 1e+5, y=avg_rec)) +
      geom_line(data=maps_10kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x="Position", y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_text(size=12),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20), legend.position="none")
    
    lands_10kb <- plot_grid(pi_10kb, u_10kb, s_10kb, r_10kb, align='v', ncol=1)
    
    pi_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_pi))+
      geom_line(data=maps_100kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scale.4d) +
      labs(title="100 kb", x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none") 
    
    u_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_mut))+
      geom_line(data=maps_100kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none")
    
    s_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_s)) +
      geom_line(data=maps_100kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x=NULL, y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none")
    
    r_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_rec)) +
      geom_line(data=maps_100kb) + theme_bw() +
      scale_y_continuous(breaks=pretty_breaks(), labels=scientific) +
      labs(title=NULL, x="Position", y=NULL) +
      theme(axis.title.x=element_text(size=16),
            axis.text.x=element_text(size=12),
            axis.text.y=element_blank(),
            plot.title=element_text(size=20),
            legend.position="none")
    
    lands_100kb <- plot_grid(pi_100kb,u_100kb,s_100kb,r_100kb, align='v',ncol=1)
    
    lands_scales <- plot_grid(lands_1kb, lands_10kb, lands_100kb, nrow=1)
    save_plot(paste("rep_", i, "/maps.png", sep=""),
              lands_scales, base_height=12, base_width=16)
    
r2s <- r2_tbl %>% group_by(model, scale) %>% 
       summarize_at(vars(Total, u, r, s, `r:s`), 
               list(avg=mean, sd=sd))


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


# the plots above show us that by and large, sd over replicates is negligible
# and u and s are the only meaningful variables affecting pi in these models
# thus we simplify:
r2s_f <- dplyr::select(r2s, -starts_with("r")) %>% 
         dplyr::select(., -ends_with("sd"))
models_f <- dplyr::select(models, c(scale_sel,
                                    shapes_mu,
                                    shapes_rec, 
                                    scales_rec, 
                                    model))
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



# check the correlation between B-values and pi
for(i in 1:num_models) {
  for(j in 1:num_models) {
    hrmap <- fread(paste("model_", i, "/rep_", j, "/hrmap.csv", sep=""))
    hrmap$B <- hrmap$Hr / hrmap$pi0
    hrmap$bin <- (hrmap$Pos - 1) %/% 1e+5
    hrbins <- hrmap %>% group_by(bin) %>% summarize_at(vars(B, Hr, pi0), mean)
    cor.test(hrbins$B, hrbins$Hr)
    plot(hrbins$Hr, hrbins$B)
  }
}


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