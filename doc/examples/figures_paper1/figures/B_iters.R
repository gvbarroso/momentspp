library(RColorBrewer)
library(tidyverse)
library(data.table)
library(scales)
library(cowplot)


B_iters <- fread("../data/B_iters.txt", header=F) # obtained with bgshr
names(B_iters) <- as.character(1:ncol(B_iters))
B_iters$iteration <- 1:nrow(B_iters)

m_B_iters <- pivot_longer(B_iters, cols=1:ncol(B_iters)-1, names_to="Pos")
m_B_iters$Pos <- as.numeric(m_B_iters$Pos)

p <- ggplot(data=m_B_iters, aes(x=Pos, y=value, color=as.factor(iteration))) +
  geom_line() + geom_point() + theme_bw() + 
  scale_color_discrete(name="Iteration", type=brewer.pal(10, "Spectral")) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Position (kb)", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p
ggsave("B_iters.pdf", width=10, height=8)
