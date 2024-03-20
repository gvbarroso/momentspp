setwd("~/Devel/momentspp/benchmarking/selection/single-pop/test_data/model_B/")

j500 <- fread("j_500/hrmap.csv.gz")
j1000 <- fread("j_1000/hrmap.csv.gz")
j2000 <- fread("j_2000/hrmap.csv.gz")

cor.test(j500$B, j2000$B)

Bjs <- cbind.data.frame(j500$B, j1000$B, j2000$B)
names(Bjs) <- c("500", "1000", "2000")
Bjs$pos <- 1:nrow(Bjs)

molten_div <- pivot_longer(Bjs, cols=c("500", "1000", "2000"), names_to="J")
ppi <- ggplot(data=molten_div[1:1e+6,], 
              aes(x=pos, y=value, color=J)) +
  geom_line(linewidth=1.05) + theme_bw() + scale_y_log10() +
  labs(title="B-values for different sampling schemes", x="Position", y="B") + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        legend.position="top")

# j500
## 7198 sampled sites
## run-time B's: 85 min
## run-time cubic splines: 1min

# j1000 + ds=5
## 3985 sampled sites
## run-time B's: 50 min
## run-time cubic splines: 30s

# j2000 + ds=5
## 2962 sampled sites
## run-time B's: 37 min
## run-time cubic splines: 30s

