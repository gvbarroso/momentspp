suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(MASS)
  library(lmtest)
  library(nlme)
  library(car)
  library(scales)
  library(data.table)
  library(interactions)
})

N=1e+5
L=1e+7

Bs <- 0.8 * rnorm(n=L, mean=1, sd=0.1)
us <- 1e-8 * rgamma(n=L, shape=10, rate=10)
Tmrcas <- rgeom(n=L, 1/2/N)

pi_1 <- us * Bs * 2 * N # expected pi
pi_2 <- us * Bs * Tmrcas # with genealogical variance (drift)
pi_3 <- rbinom(n=L, size=1, prob=pi_2) # with mutational variance on top of drift

pis <- cbind.data.frame(pi_1, pi_2, pi_3)
pis$bin <- (1:nrow(pis) - 1) %/% 1e+3

pis_1kb <- pis %>% group_by(bin) %>% summarise_at(c("pi_1", "pi_2", "pi_3"), mean)
m_pis <- pivot_longer(pis_1kb, cols=starts_with("pi_"), names_to="diversity")

p <- ggplot(data=filter(m_pis, diversity != "pi_4"), aes(x=bin, y=value, color=diversity)) +
  theme_bw() + geom_point() + #geom_line() + 
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="1 kb", x="Pos (kb)", y=expression(pi)) +
  theme(axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=16),
        plot.title=element_text(size=20),
        legend.position="bottom") 
p

# increasing diploid sample size S
S <- 20
pb <- txtProgressBar(min=1, max=S, style=3)
for(i in 1:(S-1)) {
  setTxtProgressBar(pb, i)
  Tmrcas <- rgeom(n=L, 1/2/N)
  x <- us * Bs * Tmrcas 
  y <- rbinom(n=L, size=1, prob=x)
  
  pi_2 <- pi_2 + x
  pi_3 <- pi_3 + y
}
close(pb)

pi_2 <- pi_2 / S
pi_3 <- pi_3 / S

pis <- cbind.data.frame(pi_1, pi_2, pi_3)
pis$bin <- (1:nrow(pis) - 1) %/% 1e+3

pis_1kb <- pis %>% group_by(bin) %>% summarise_at(c("pi_1", "pi_2", "pi_3"), mean)
m_pis <- pivot_longer(pis_1kb, cols=starts_with("pi_"), names_to="diversity")
  
p <- ggplot(data=filter(m_pis, diversity != "pi_4"), aes(x=bin, y=value, color=diversity)) +
  theme_bw() + geom_point() + #geom_line() + 
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="1 kb", x="Pos (kb)", y=expression(pi)) +
  theme(axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=16),
        plot.title=element_text(size=20),
        legend.position="bottom") 
p

