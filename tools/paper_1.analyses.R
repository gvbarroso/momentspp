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

Nanc <- unique(params$Na)
N1 <- unique(params$N1)
tsplit <- unique(params$t)

demo_models <- crossing(Nanc, N1, tsplit)
num_demo_models <- nrow(demo_models)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=251))
names(pi0) <- as.character(seq(from=tsplit, to=0, by=-10))

for(i in 1:num_demo_models) {
  
  moms <- fread(paste("pilot/demo_neutral/model_", i, 
                      "_e_1_expectations.txt", sep=""))
  
  pi0[i,] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
}

demo_pi0 <- cbind.data.frame(demo_models, pi0)
fwrite(demo_pi0, "pilot/demo_neutral/demo_pi0.csv")

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(seq(from=tsplit, to=0, by=-10)),
                      names_to="gen")
m_pi0$gen <- as.numeric(m_pi0$gen)

p <- ggplot(data=m_pi0, aes(x=gen, y=value, color=as.factor(N1))) + 
  geom_point(size=0.5) + geom_line() + theme_bw() + 
  labs(title="Trajectory of diversity after split 2500 generations ago", 
       x="gen", y=expression(pi[0])) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("pilot/demo_neutral/pi0_demo.png", p, base_height=8, base_width=16)

d1 <- ggplot(data=demo_models) + theme_bw() + facet_wrap(~N1) +
      geom_segment(aes(x=0, xend=tsplit, y=N1, yend=N1), linewidth=1.5) +
      geom_segment(aes(x=tsplit, xend=tsplit, y=N1, yend=Nanc), linewidth=1.5) +
      geom_segment(aes(x=tsplit, xend=2e+4, y=Nanc, yend=Nanc), linewidth=1.5) +
      scale_y_log10() +
      labs(title="Demographic Models", x="Generation", y="N") +
      theme(axis.title=element_text(size=16), 
            axis.text=element_text(size=12), 
            axis.text.x=element_text(size=12),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16),
            legend.position="bottom")

save_plot("pilot/demo_neutral/demo_models.png",
          d1, base_height=8, base_width=16)
    