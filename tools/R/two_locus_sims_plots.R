

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(scales)
  library(bslib)
  library(purrr)
})

setwd("~/Data/momentspp/D_iagnostics/selected/")

sim_params <- "Ne=1e+4_s=-1e-4_u=1e-8_r=1e-4_L=1_G=1e+9"
type <- "Sel" # Neu || Sel

stats_files <- list.files(pattern="stats_test_*", full.names=TRUE)
stats <- list(length=length(stats_files))

pb <- txtProgressBar(min=0, max=length(stats_files), style=3)
for(i in 1:length(stats_files)) {
  setTxtProgressBar(pb, i)
  stats[[i]] <- fread(stats_files[i]) %>% dplyr::select(., -V2) %>% dplyr::filter(., V1!="random_seed")
}
close(pb)

stats <- bind_rows(stats)
names(stats) <- c("Statistic", "Value")
avg_stats <- stats %>% group_by(Statistic) %>% summarize(mean(Value))

# D profile
d_files <- list.files(pattern="tbl_D_test_*", full.names=TRUE)
Ds <- list(length=length(d_files))

pb <- txtProgressBar(min=0, max=length(d_files), style=3)
for(i in 1:length(d_files)) {
  setTxtProgressBar(pb, i)
  Ds[[i]] <- fread(d_files[i]) 
}
close(pb)

Ds <- bind_rows(Ds)

get_D_prime <- function(D, p, q) {
  d_max = 1
  if (D < 0)
    d_max = min(p * q, (1-p) * (1-q))
  
  else if (D > 0)
    d_max = min(p * (1-q), (1-p) * q)
  
  else
    return(NA)
  
  return(D / d_max)
}

Ds <- Ds %>% mutate(D_prime=pmap_dbl(., get_D_prime))

Ds$Dr <- Ds$D * (1 - 2 * Ds$q)
Ds$Dl <- Ds$D * (1 - 2 * Ds$p)
Ds$Dz <- Ds$D * (1 - 2 * Ds$q) * (1 - 2 * Ds$p)

fwrite(Ds, paste("Ds_", sim_params, ".csv.gz", sep=""))
Ds <- fread(paste("Ds_", sim_params, ".csv.gz", sep=""))

summary(Ds)

Ds_q <- Ds %>% group_by(q) %>% summarize_at(vars(D, D_prime, Dr, Dl, Dz, p), list(avg=mean))
Ds_p <- Ds %>% group_by(p) %>% summarize_at(vars(D, D_prime, Dr, Dl, Dz, q), list(avg=mean))

#Dx <- Ds %>% group_by(p, q) %>% summarize_at(vars(D, Dr, Dl, Dz), list(avg=mean))
#D_p_1 <- ggplot(data=Dx, aes(x=p, y=D_avg, color=q)) + 
#  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
#  scale_colour_gradientn(colours = terrain.colors(10)) +
#  geom_hline(yintercept=0, color="red", linetype="dashed") + 
#  scale_x_continuous(breaks=pretty_breaks()) +
#  scale_y_continuous(breaks=pretty_breaks()) +
#  labs(title=sim_params, x="p", y="D") +
#  theme(axis.title=element_text(size=16), 
#        axis.text=element_text(size=12), 
#        axis.text.x=element_text(size=12), 
#        axis.text.y=element_text(size=12),
#        legend.position="bottom")
#ggsave("../Dx_p_S.pdf", D_p_1, dpi=500, width=10, height=7)

# nominal distributions
D_p <- ggplot(data=Ds_p, aes(x=p, y=D_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x=NULL, y="D") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="none")
#ggsave(paste("../D_p_nominal_", type,".pdf", sep=""), D_p, dpi=500, width=10, height=7)

D_q <- ggplot(data=Ds_q, aes(x=q, y=D_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="none")
#ggsave(paste("../D_q_nominal_", type,".pdf", sep=""), D_q, dpi=500, width=10, height=7)

Dr_q <- ggplot(data=Ds_q, aes(x=q, y=Dr_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10), breaks=c(0.03, 0.037, 0.045)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x=NULL, y=NULL) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="none")
#ggsave(paste("../Dr_q_nominal_", type,".pdf", sep=""), Dr_q, dpi=500, width=10, height=7)

Dr_p <- ggplot(data=Ds_p, aes(x=p, y=Dr_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10), breaks=c(0.05, 0.1, 0.15)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x=NULL, y="D(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="none")
#ggsave(paste("../Dr_p_nominal_", type,".pdf", sep=""), Dr_p, dpi=500, width=10, height=7)

Dl_q <- ggplot(data=Ds_q, aes(x=q, y=Dl_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10), breaks=c(0.03, 0.037, 0.045)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="q", y=NULL) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
#ggsave(paste("../Dl_q_nominal_", type,".pdf", sep=""), Dl_q, dpi=500, width=10, height=7)

Dl_p <- ggplot(data=Ds_p, aes(x=p, y=Dl_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10), breaks=c(0.05, 0.1, 0.15)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="p", y="D(1-2p)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
#ggsave(paste("../Dl_p_nominal_", type,".pdf", sep=""), Dl_p, dpi=500, width=10, height=7)

pgd <- plot_grid(D_p, D_q, Dr_p, Dr_q, Dl_p, Dl_q,
                 rel_heights=c(1, 1, 1.35), 
                 nrow=3, labels="AUTO", label_size=18, align='v')
save_plot("../D_sims.pdf", pgd, base_height=15, base_width=16)

Dz_q <- ggplot(data=Ds_q, aes(x=q, y=Dz_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="q", y="D(1-2p)(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dz_q_nominal_", type,".pdf", sep=""), Dz_q, dpi=500, width=10, height=7)

Dz_p <- ggplot(data=Ds_p, aes(x=p, y=Dz_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="p", y="D(1-2p)(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dz_p_nominal_", type,".pdf", sep=""), Dz_p, dpi=500, width=10, height=7)

# weighted distributions
Ds_samp <- slice_sample(Ds, n=2e+5)
D_qs <- ggplot(data=Ds_samp, aes(x=q, y=D, color=p)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=sim_params, x="q", y="D") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../D_q_sample_", type,".pdf", sep=""), D_qs, dpi=500, width=10, height=7)

D_ps <- ggplot(data=Ds_samp, aes(x=p, y=D, color=q)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=sim_params, x="p", y="D") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../D_p_sample_", type,".pdf", sep=""), D_ps, dpi=500, width=10, height=7)

Dr_q <- ggplot(data=Ds_samp, aes(x=q, y=Dr, color=p)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=sim_params, x="q", y="D(1-2q)") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../D_q_dist_", type, ".pdf", sep=""), D_q, dpi=500, width=10, height=7)

##########################
#
# c-space time series
#
##########################

c_files <- list.files(pattern="hap_trajectories*", full.names=TRUE)
c_files <- c_files[1:5e+3]
cs <- list(length=length(c_files))

pb <- txtProgressBar(min=0, max=length(c_files), style=3)
for(i in 1:length(c_files)) {
  setTxtProgressBar(pb, i)
  cs[[i]] <- fread(c_files[i], header=T, colClasses = 'character', fill=Inf) 
}
close(pb)

x <- bind_rows(cs)
cs <- x

cs$I <- as.integer(cs$fAB != "removed")
cs$Series <- with(cs, ave(I, cumsum(I == 0), FUN = cumsum))
y <- filter(cs, fAB!="removed")

#z <- y[,1:4]
#x <- lapply(z, strsplit, split= "= ")
#c1 <- unlist(x$fAB)
#fAB = as.numeric(c1[seq_len(length(c1)) %% 2 != 1])
#c2 <- unlist(x$fAb)
#fAb = as.numeric(c2[seq_len(length(c2)) %% 2 != 1])
#c3 <- unlist(x$faB)
#faB = as.numeric(c3[seq_len(length(c3)) %% 2 != 1])
#c4 <- unlist(x$fab)
#fab = as.numeric(c4[seq_len(length(c4)) %% 2 != 1])

y$fAB <- as.numeric(y$fAB)
y$fAb <- as.numeric(y$fAb)
y$faB <- as.numeric(y$faB)
y$fab <- as.numeric(y$fab)

hap_df <- dplyr::select(y, c(fAB, fAb, faB, fab, Series))
names(hap_df)[5] <- "Generation"
hap_df$p <- hap_df$fAB + hap_df$fAb
hap_df$q <- hap_df$fAB + hap_df$faB
hap_df$D <- hap_df$fAB * hap_df$fab - hap_df$faB * hap_df$fAb
hap_df$Dr <- (hap_df$fAB * hap_df$fab - hap_df$faB * hap_df$fAb) * (1 - 2 * hap_df$q)

hap_df_avg <- na.omit(hap_df) %>% group_by(Generation) %>% 
  summarize_at(vars(fAB, fAb, faB, fab, p, q, D, Dr), list(avg=mean))

m_hap_series <- pivot_longer(hap_df_avg, cols=starts_with("f"), names_to="hap")

p0 <- ggplot(data=m_hap_series, aes(x=Generation, y=value, color=hap)) + 
  theme_bw() + geom_line() +
  scale_color_discrete(name="Hap", type=c("plum3", "seagreen3", "brown1", "orange")) +
  geom_hline(yintercept=0, color="black", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="Hap. Freq.") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
save_plot("../haps_time.pdf", p0, base_height=8, base_width=10)

p1 <- ggplot(data=hap_df_avg, aes(x=Generation, y=D_avg)) + 
  theme_bw() + geom_line() +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="D") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
save_plot("../D_time.pdf", p1, base_height=8, base_width=10)

p2 <- ggplot(data=hap_df_avg, aes(x=Generation, y=Dr_avg)) + 
  theme_bw() + geom_line() +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="D(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
save_plot("../Dr_time.pdf", p2, base_height=8, base_width=10)

m_D_time <- pivot_longer(hap_df_avg, cols=starts_with("D"), names_to="Stat")
q0 <- ggplot(data=m_D_time, aes(x=Generation, y=value, color=Stat)) + 
  theme_bw() + geom_line() +
  scale_color_discrete(name="Stat", type=c("plum3", "seagreen3")) +
  geom_hline(yintercept=0, color="black", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="Value") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="none")
save_plot("../D_stats_time.pdf", q0, base_height=10, base_width=12)

p3 <- ggplot(data=hap_df_avg, aes(x=Generation, y=p_avg)) + 
  theme_bw() + geom_line() +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="Mean p") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
save_plot("../p_time.pdf", p3, base_height=10, base_width=12)

p4 <- ggplot(data=hap_df_avg, aes(x=Generation, y=q_avg)) + 
  theme_bw() + geom_line() + 
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="mean q") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
save_plot("../q_time.pdf", p4, base_height=10, base_width=12)

m_hap_df <- pivot_longer(na.omit(hap_df), cols=c(fAB, fAb, faB, fab), names_to="Hap")

p5 <- ggplot(data=m_hap_df[1:3e+5,], aes(x=Generation, y=value, color=Hap)) + 
  theme_bw() + geom_point(shape=21, alpha=0.75) + 
  scale_color_discrete(name="Hap", type=c("plum3", "seagreen3", "brown1", "orange")) +
  geom_hline(yintercept=0, color="grey20", linetype="dashed") + 
  labs(title=NULL, x="Generation", y="Hap. Freq.") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_log10() +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")

p6 <- ggplot(data=hap_df[1:1e+6,], aes(x=Generation, y=D, color=q)) + 
  theme_bw() + geom_point(shape=21, alpha=1) + 
  geom_smooth(color="seagreen3") + 
  geom_hline(yintercept=0, color="grey20", linetype="dashed") + 
  labs(title=NULL, x="Generation", y="D") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous() +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
