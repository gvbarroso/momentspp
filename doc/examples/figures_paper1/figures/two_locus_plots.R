
library(tidyverse)
library(scales)
library(bigsnpr) # for seq_log
library(data.table)
library(cowplot)

lookup_tbl <- fread("../lookup_tables/lookup_tbl_two-locus.csv.gz")
lookup_tbl$Ns[which(lookup_tbl$Ns == "1000;10000")] <- "10-fold bottleneck"
lookup_tbl$Ns[which(lookup_tbl$Ns == "1e+05;10000")] <- "10-fold expansion"
scaleFUN <- function(x) sprintf("%.5f", x)

df <- dplyr::filter(lookup_tbl, s==c(-1e-3, -1e-4, -1e-5))
p1 <- ggplot(data=filter(df, s==-1e-3),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=c(0, 30000, 60000, 90000, 120000)) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

p2 <- ggplot(data=filter(df, s==-1e-4),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

p3 <- ggplot(data=filter(df, s==-1e-5),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x="Generations ago", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="bottom")


cp1 <- plot_grid(p1, p2, p3, ncol=1, rel_heights=c(1, 1, 1.35), labels="AUTO", align='v')
cp1
save_plot("two-locus_time_B.pdf", cp1, base_height=15, base_width=12)

q1 <- ggplot(data=filter(lookup_tbl, s==-1e-3), aes(x=Generation, y=Hl)) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=c(0, 30000, 60000, 90000, 120000)) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="Hl") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

q2 <- ggplot(data=filter(lookup_tbl, s==-1e-4), aes(x=Generation, y=Hl)) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=c(0, 30000, 60000, 90000, 120000)) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="Hl") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

q3 <- ggplot(data=filter(lookup_tbl, s==-1e-5), aes(x=Generation, y=Hl)) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=c(0, 30000, 60000, 90000, 120000)) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x="Generations ago", y="Hl") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

cq1 <- plot_grid(q1, q2, q3, ncol=1, rel_heights=c(1, 1, 1.35), labels="AUTO", align='v')
cq1
save_plot("Hl_time.pdf", cq1, base_height=15, base_width=12)

# looking for an extended period of time in the expansion scenario
lookup_tbl <- fread("../lookup_tables/lookup_tbl_two-locus_extended.csv.gz")
lookup_tbl$Ns[which(lookup_tbl$Ns == "1e+05;10000")] <- "10-fold expansion"

p4 <- ggplot(data=filter(lookup_tbl, s==-1e-3),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) + facet_grid(~s~Ns) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

p5 <- ggplot(data=filter(lookup_tbl, s==-1e-4),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) +facet_grid(~s~Ns) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x=NULL, y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="none")

p6 <- ggplot(data=filter(lookup_tbl, s==-1e-5),
             aes(x=Generation, y=B^1000, color=as.factor(r))) + 
  theme_bw() + geom_line(linewidth=1.25) +  facet_grid(~s~Ns) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks(), labels=scaleFUN) +
  labs(title=NULL, x="Generations ago", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12),
        strip.text.x=element_text(size=13),
        strip.text.y=element_text(size=13),
        legend.position="bottom") +
  guides(color=guide_legend(title="r"))

cp2 <- plot_grid(p4, p5, p6, ncol=1, rel_heights=c(1, 1, 1.35),
                 labels="AUTO", align='v')
cp2
save_plot("two-locus_time_exp_extended.pdf", cp2, base_height=15, base_width=8)
