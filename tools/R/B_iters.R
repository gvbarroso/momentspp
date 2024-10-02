library(RColorBrewer)

setwd("~/Data/mpp_figures/figures/")

B_iters <- fread("B_iters.txt", header=F)  # from Python
names(B_iters) <- as.character(1:ncol(B_iters))
B_iters$iteration <- 1:nrow(B_iters)

m_B_iters <- pivot_longer(B_iters, cols=1:ncol(B_iters)-1, names_to="Pos")
m_B_iters$Pos <- as.numeric(m_B_iters$Pos)

p <- ggplot(data=m_B_iters, aes(x=Pos, y=value, color=as.factor(iteration))) +
  geom_line() + geom_point() + theme_bw() + 
  scale_color_discrete(name="Iteration") +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Position (kb)", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
