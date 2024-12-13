

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(scales)
  library(bslib)
  library(purrr)
})

setwd("~/Data/momentspp/D_iagnostics/neutral/")

sim_params <- "Ne=1e+4_s=0_u=1e-6_r= 1e-4_L=10_G=1e+7"
type <- "Neu" # Neu || Sel

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
#ggsave("../Dx_p_S.png", D_p_1, dpi=500, width=10, height=7)

# nominal distributions
D_p <- ggplot(data=Ds_p, aes(x=p, y=D_avg, color=q_avg)) + 
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
ggsave(paste("../D_p_nominal_", type,".png", sep=""), D_p, dpi=500, width=10, height=7)

D_q <- ggplot(data=Ds_q, aes(x=q, y=D_avg, color=p_avg)) + 
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
ggsave(paste("../D_q_nominal_", type,".png", sep=""), D_q, dpi=500, width=10, height=7)

Dr_q <- ggplot(data=Ds_q, aes(x=q, y=Dr_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="q", y="D(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dr_q_nominal_", type,".png", sep=""), Dr_q, dpi=500, width=10, height=7)

Dr_p <- ggplot(data=Ds_p, aes(x=p, y=Dr_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="p", y="D(1-2q)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dr_p_nominal_", type,".png", sep=""), Dr_p, dpi=500, width=10, height=7)

Dl_q <- ggplot(data=Ds_q, aes(x=q, y=Dl_avg, color=p_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="q", y="D(1-2p)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dl_q_nominal_", type,".png", sep=""), Dl_q, dpi=500, width=10, height=7)

Dl_p <- ggplot(data=Ds_p, aes(x=p, y=Dl_avg, color=q_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="black") + 
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=sim_params, x="p", y="D(1-2p)") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
ggsave(paste("../Dl_p_nominal_", type,".png", sep=""), Dl_p, dpi=500, width=10, height=7)

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
ggsave(paste("../Dz_q_nominal_", type,".png", sep=""), Dz_q, dpi=500, width=10, height=7)

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
ggsave(paste("../Dz_p_nominal_", type,".png", sep=""), Dz_p, dpi=500, width=10, height=7)

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
ggsave(paste("../D_q_sample_", type,".png", sep=""), D_qs, dpi=500, width=10, height=7)

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
ggsave(paste("../D_p_sample_", type,".png", sep=""), D_ps, dpi=500, width=10, height=7)

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
ggsave("../D_q_dist_N.png", D_q, dpi=500, width=10, height=7)

##########################
#
# c-space time series
#
##########################

c_files <- list.files(pattern="hap_trajectories*", full.names=TRUE)
#c_files <- c_files[1:1e+4]
cs <- list(length=length(c_files))

pb <- txtProgressBar(min=0, max=length(c_files), style=3)
for(i in 1:length(c_files)) {
  setTxtProgressBar(pb, i)
  cs[[i]] <- fread(c_files[i], header=T, colClasses = 'character') 
}
close(pb)

x <- bind_rows(cs)
cs <- x

cs$I <- as.integer(cs$fAB != "removed")
cs$Series <- with(cs, ave(I, cumsum(I == 0), FUN = cumsum))
y <- filter(cs, fAB!="removed")
z <- y[,1:4]
q <- lapply(z, strsplit, split= "= ")

c1 <- unlist(q$fAB)
fAB = as.numeric(c1[seq_len(length(c1)) %% 2 != 1])
c2 <- unlist(q$fAb)
fAb = as.numeric(c2[seq_len(length(c2)) %% 2 != 1])
c3 <- unlist(q$faB)
faB = as.numeric(c3[seq_len(length(c3)) %% 2 != 1])
c4 <- unlist(q$fab)
fab = as.numeric(c4[seq_len(length(c4)) %% 2 != 1])

hap_df <- cbind.data.frame(fAB, fAb, faB, fab, y$Series)
names(hap_df)[5] <- "Generation"
hap_df$p <- hap_df$fAB + hap_df$fAb
hap_df$q <- hap_df$fAB + hap_df$faB
hap_df$D <- hap_df$fAB * hap_df$fab - hap_df$faB * hap_df$fAb

hap_df_avg <- hap_df %>% group_by(Generation) %>% 
  summarize_at(vars(fAB, fAb, faB, fab, p, q, D), list(avg=mean))

p1 <- ggplot(data=hap_df_avg, aes(x=Generation, y=D_avg)) + 
  theme_bw() + geom_line() + #geom_smooth(color="black") + 
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  #geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="D") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")

p2 <- ggplot(data=hap_df_avg, aes(x=Generation, y=p_avg)) + 
  theme_bw() + geom_line() + #geom_smooth(color="black") + 
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  #geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="mean p") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")

p3 <- ggplot(data=hap_df_avg, aes(x=Generation, y=q_avg)) + 
  theme_bw() + geom_line() + #geom_smooth(color="black") + 
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  #geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="Generations after 2nd mutation", y="mean q") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")

p4 <- ggplot(data=hap_df, aes(x=q, y=D)) + 
  theme_bw() + geom_line() + #geom_smooth(color="black") + 
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  #geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="q", y="D") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")

hap_df_avg_2 <- hap_df %>% group_by(q) %>% 
  summarize_at(vars(fAB, fAb, faB, fab, p, D), list(avg=mean))

p5 <- ggplot(data=hap_df_avg_2, aes(x=q, y=D_avg)) + 
  theme_bw() + geom_point(shape=21) + geom_smooth(color="blue") + 
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") + 
  labs(title=NULL, x="q", y="D") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        legend.position="bottom")
