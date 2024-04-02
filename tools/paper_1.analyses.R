# Date created: 01/04/2024
# Last modified: 01/04/2024
# Author: Gustavo V. Barroso

pdf(NULL) # to suppress creation of Rplots.pdf

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(pracma) # for cubicspline()
  library(cowplot)
  library(scales)
  library(RColorBrewer)
})

setwd("~/Data/momentspp/paper_1/")

params <- fread("pilot/params.csv")

Na <- unique(params$Na)
N1 <- unique(params$N1)
tchange <- unique(params$t)

demo_models <- crossing(Nanc, N1, tchange)
num_demo_models <- nrow(demo_models)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=251))
names(pi0) <- as.character(seq(from=tchange, to=0, by=-10))

for(i in 1:num_demo_models) {
  
  moms <- fread(paste("pilot/demo_neutral/model_", i, 
                      "_e_1_expectations.txt", sep=""))
  
  pi0[i,] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
}

demo_pi0 <- cbind.data.frame(demo_models, pi0)
fwrite(demo_pi0, "pilot/demo_neutral/demo_pi0.csv")

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(seq(from=tchange, to=0, by=-10)),
                      names_to="Generation", values_to="pi0")
m_pi0$Generation <- as.numeric(m_pi0$Generation)

N.labs <- paste("Model ", 1:num_demo_models)
names(N.labs) <- as.character(unique(demo_models$N1))

p <- ggplot(data=m_pi0, aes(x=Generation, y=pi0)) + 
  facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
  geom_point(size=0.5) + geom_line() + theme_bw() + 
  labs(title="Trajectory of Heterozygosity after size change", 
       x="Generations ago", y=expression(pi[0])) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=14),
        legend.position="bottom")

save_plot("pilot/demo_neutral/pi0_demo.png", p, base_height=8, base_width=16)

d <- ggplot(data=demo_models) + theme_bw() + 
     facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
     geom_segment(aes(x=0, xend=tchange, y=N1, yend=N1), linewidth=1.5) +
     geom_segment(aes(x=tchange, xend=tchange, y=N1, yend=Nanc), linewidth=1.5) +
     geom_segment(aes(x=tchange, xend=tchange+100, y=Nanc, yend=Nanc), linewidth=1.5) +
     scale_y_log10() +
     scale_x_continuous(breaks=c(0, 1000, tchange)) +
     labs(title="Demographic Models", 
          x="Generations ago", 
          y="Population Size (N)") +
     theme(axis.title=element_text(size=16), 
           axis.text=element_text(size=12), 
           axis.text.x=element_text(size=12),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           strip.text=element_text(size=14),
           legend.position="bottom")

save_plot("pilot/demo_neutral/demo_models.png",
          d, base_height=8, base_width=16)

# loads tables from setup_models_1.Rmd
lookup_tbl <- fread("pilot/lookup_tbl.csv.gz")
look_s <- unique(params$lookup_s)
look_r <- unique(params$lookup_r)

hl_demo <- fread("pilot/hl_time.csv.gz")
hr_demo <- fread("pilot/hr_time.csv.gz")

plots <- list(length=length(unique(lookup_tbl$N1)))
for(i in 1:length(unique(lookup_tbl$N1))) {
  plots[[i]] <- ggplot(data=filter(lookup_tbl, lookup_r==1e-9,
                                    N1==unique(lookup_tbl$N1)[i]), 
                                    aes(x=-lookup_s, y=hr)) + theme_bw() +
    geom_point() + geom_line() + scale_x_log10() + scale_y_continuous() +
    labs(title=NULL, x="s",  y="Hr") +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          strip.text=element_text(size=14),
          legend.position="bottom")
}

x <- plot_grid(plotlist=plots)
save_plot("pilot/hr_s_n.png", x, base_height=8, base_width=16)

#m_hl_demo <- pivot_longer(hl_demo, 
#                          cols=as.character(seq(from=tchange, to=0, by=-10)),
#                          names_to="Generation", values_to="Hl")
#m_hl_demo$Generation <- as.numeric(m_hl_demo$Generation)
#m_hl_demo <- setDT(m_hl_demo)
#setkey(m_hl_demo, N1, lookup_r, lookup_s)

m_hr_demo <- pivot_longer(hr_demo, 
                          cols=as.character(seq(from=tchange, to=0, by=-10)),
                          names_to="Generation", values_to="Hr")
m_hr_demo$Generation <- as.numeric(m_hr_demo$Generation)
m_hr_demo <- setDT(m_hr_demo)
setkey(m_hr_demo, N1, lookup_r, lookup_s)

b <- left_join(m_hr_demo, m_pi0, by=c("N1", "Generation"))
b$B <- b$Hr / b$pi0


# TODO simulate chr with genes and rec map,
# then plot the B landscape over time!


# focusing on a particular (r, s) pair
closest_s <- look_s[which.min(abs(look_s - -1.5e-4))]
closest_r <- look_r[which.min(abs(look_r - 1e-8))]

plots <- list(length=length(unique(lookup_tbl$N1)))
for(i in 1:length(unique(lookup_tbl$N1))) {
  
  tmp_hl <- filter(m_hl_demo, N1==unique(m_hl_demo$N1)[i], 
                              lookup_r==closest_r, 
                              lookup_s==closest_s)
  tmp_hr <- filter(m_hl_demo, N1==unique(m_hr_demo$N1)[i], 
                              lookup_r==closest_r, 
                              lookup_s==closest_s)
  tmp_pi0 <- filter(m_pi0, N1==unique(m_hr_demo$N1)[i])
  
  
  stats <- tmp_hr
  stats$Hl <- tmp_hl$Hl
  stats$pi0 <- tmp_pi0$pi0
  stats$B <- stats$Hr / stats$pi0
  stats$E <- stats$Hl / stats$pi0
    
  m_stats <- pivot_longer(stats, cols=c(B), names_to="Statistic")
    
  plots[[i]] <- ggplot(data=m_stats, aes(x=Generation,
                                         y=value, 
                                         color=Statistic)) +
                  geom_point() + geom_line() +  theme_bw() + scale_y_log10() +
                  labs(title=NULL, x="Generation", y=NULL) +
                  theme(axis.title=element_text(size=16), 
                        axis.text=element_text(size=12), 
                        axis.text.x=element_text(size=12),
                        legend.text=element_text(size=16),
                        legend.title=element_text(size=16),
                        strip.text=element_text(size=14),
                        legend.position="bottom")
}

x <- plot_grid(plotlist=plots)
save_plot("pilot/stats_time.png", x, base_height=8, base_width=16)
