# first we check the correlation between B-values and pi
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

# plots
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

pi_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_pi)) +
  geom_line(data=maps_100kb) + theme_bw() +
  scale_y_continuous(breaks=pretty_breaks(), labels=scale.4d) +
  labs(title="100 kb", x=NULL, y=NULL) +
  theme(axis.title.x=element_text(size=16),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=20),
        legend.position="none") 

u_100kb <- ggplot(data=maps_100kb, aes(x=(bin_100kb -1) * 1e+6, y=avg_mut)) +
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

lands_100kb <- plot_grid(pi_100kb,u_100kb,s_100kb, r_100kb, align='v', ncol=1)

lands_scales <- plot_grid(lands_1kb, lands_10kb, lands_100kb, nrow=1)
save_plot(paste("rep_", i, "/maps.png", sep=""),
          lands_scales, base_height=12, base_width=16)

# using max rec distance (r==1e-2) to visually check validity of rec threshold
rex <- apply(prd, 2, function(x) x < 4 * N * 1e-2)
ex400 <- apply(rex, 1, function(x) which(x))

Bl400 <- unlist(lapply(ex400[[1]], getB, focal_neutral=1))
Bq400 <- unlist(lapply(ex400[[length(samp_pos)/4]], getB,
                       focal_neutral=length(samp_pos)/4))
Bm400 <- unlist(lapply(ex400[[length(samp_pos)/2]], getB,
                       focal_neutral=length(samp_pos)/2))

Blr <- unlist(lapply(exons_per_samp_neut[[1]], getB, focal_neutral=1))
Bqr <- unlist(lapply(exons_per_samp_neut[[length(samp_pos)/4]], getB,
                     focal_neutral=length(samp_pos)/4))
Bmr <- unlist(lapply(exons_per_samp_neut[[length(samp_pos)/2]], getB,
                     focal_neutral=length(samp_pos)/2))

xl <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bl400))]), Bl400)
xq <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bq400))]), Bq400)
xm <- cbind.data.frame(as.numeric(prd[1, as.numeric(names(Bm400))]), Bm400)

names(xl) <- c("rec", "B")
names(xq) <- c("rec", "B")
names(xm) <- c("rec", "B")

xl$CummBval <- cumprod(xl$B)
xq$CummBval <- cumprod(xq$B)
xm$CummBval <- cumprod(xm$B)

xl$relevant <- xl$B %in% Blr
xq$relevant <- xq$B %in% Bqr
xm$relevant <- xm$B %in% Bmr

pl <- ggplot(data=xl, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values Left", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

pq <- ggplot(data=xq, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values 1/4 chr", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

pm <- ggplot(data=xm, aes(x=rec, y=B, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="B-values 1/2 chr", x=NULL, y=NULL) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="none")

ql <- ggplot(data=xl, aes(x=rec, y=CummBval, color=relevant)) + geom_point() +
  theme_bw() + scale_y_log10(limits=c(min(xl$CummBval), 1)) +
  labs(title="Cumm. B-value, Left", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

qq <- ggplot(data=xq, aes(x=rec, y=CummBval, color=relevant)) + 
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="Cumm. B-value, 1/4 chr", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

qm <- ggplot(data=xm, aes(x=rec, y=CummBval, color=relevant)) +
  geom_point() + theme_bw() + scale_y_log10() +
  labs(title="Cumm. B-value, 1/2 chr", x="Cummulative Rho", y="B") +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=90, size=12, vjust=0.5, hjust=1.0),
        legend.position="bottom")

cp <- plot_grid(pl, pq, pm, ql, qq, qm, nrow=2)
save_plot(paste("rep_", i, "/Bvals.png", sep=""),
          cp, base_height=10, base_width=12)

# plotting diversity per site for an arbitrary segment of the chr
seg <- filter(dt_exons, start >= 0, end <= 1e+6)
wsize <- seg$end[nrow(seg)] - seg$start[1]
barcode <- ggplot(data=seg) +
  geom_segment(aes(x=start, xend=end, y=1, yend=1, size=3, color=s))+
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