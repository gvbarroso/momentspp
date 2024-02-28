#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

library(tidyverse)
library(cowplot)
library(MASS)
library(lmtest)
library(nlme)
library(car)
library(data.table)
library(GenomicRanges)
library(pracma)
library(bigsnpr)
library(scales)

scale.4d <- function(x) sprintf("%.4f", x) # for plotting

# loads the main tables
models <- fread("models.csv")
lookup_tbl <- fread("lookup_tbl.csv")
setkey(lookup_tbl, lookup_r, lookup_s)

# for discretizing Gamma distributions and make look-up faster
dt_r <- data.table(unique(lookup_tbl$lookup_r))
dt_s <- data.table(unique(lookup_tbl$lookup_s))
names(dt_r) <- "lookup_r"
names(dt_s) <- "lookup_s"
setkey(dt_r, lookup_r)
setkey(dt_s, lookup_s)

m <- as.numeric(args[1]) # model index

# "global" parameters
num_reps <- 10
jump_length <- 500 # sampling distance, for interpolating B-values
bin_size <- 1e+3 # "basal" scale, larger windows are built after
N <- unique(lookup_tbl$N)
u <- unique(lookup_tbl$u)

# tables for storing R^2 values from linear models
r2_tbl_1kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))
r2_tbl_10kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))
r2_tbl_100kb <- as.data.frame(matrix(ncol=5, nrow=num_reps))

names(r2_tbl_1kb) <- c("Total", "u", "r", "s", "r:s")
names(r2_tbl_10kb) <- c("Total", "u", "r", "s", "r:s")
names(r2_tbl_100kb) <- c("Total", "u", "r", "s", "r:s")

for(i in 1:num_reps) {
  
  cat(paste("Running analyses for model ", m, 
            ", replicate ", i, " [", Sys.time(), "]\n", sep=""))
  
  ncsl <- rgeom(n=models$num_constrained[m] + 1, 
                prob=1/models$avg_neutral_lengths[m])
  csl <- rep(models$exon_lengths[m], models$num_constrained[m]) 
  L <- sum(csl) + sum(ncsl)

  ss <- -rgamma(n=length(csl), shape=models$shape_sel[m],
                scale=models$scale_sel[m]/models$shape_sel[m])
  alpha <- 2 * models$N[m] * ss

  while(any(alpha < -20)) { 
    cat("Selection is too strong in at least one gene! Redrawing...\n")
    ss <- -rgamma(n=length(csl), shape=models$shape_sel[m],
                  scale=models$scale_sel[m]/models$shape_sel[m])
    alpha <- 2 * models$N[m] * ss
  }

  disc_ss <- dt_s[dt_s[J(ss), roll="nearest", which=T]]$lookup_s # discretize

  intervals <- c(rbind(ncsl, csl))[-2*length(ncsl)]
  intervals <- setDT(bind_cols(chr="chr1",
			                         dplyr::lag(cumsum(intervals), n=1, default=0),
			                         dplyr::lag(cumsum(intervals), n=0, default=0)))
  names(intervals) <- c("chr", "start", "end")
  intervals$s <- c(rbind(rep(0, length(csl)), disc_ss), 0)
  setkey(intervals, start, end)
  fwrite(intervals, paste("rep_", i, "/intervals.csv", sep=""), sep=",")

  dt_neutral <- filter(intervals[,2:4], s==0)
  setkey(dt_neutral, start, end)
  dt_exons <- filter(intervals, s<0)
  dt_exons$u <- u # all exons have same mut rate
  setkey(dt_exons, start, end)
  fwrite(dt_exons, paste("rep_", i, "/exons.csv", sep=""), sep=",")

  rec_spans <- rgeom(n=ceiling(L/models$avg_rec_spans[m]),
	 	                 prob=1/models$avg_rec_spans[m])
  while(sum(rec_spans) < L) {
    rec_spans <- c(rec_spans, rgeom(n=1, prob=1/models$avg_rec_spans[m]))
  }
  if(sum(rec_spans) > L) {
    rec_spans <- rec_spans[cumsum(rec_spans) < L]
    rec_spans <- c(rec_spans, L - sum(rec_spans))
  }
  rec_spans <- rec_spans[rec_spans>0]

  rmap <- setDT(bind_cols(chr="chr1",
	        		  dplyr::lag(cumsum(rec_spans), n=1, default=0),
			          dplyr::lag(cumsum(rec_spans), n=0, default=0)))  
  names(rmap) <- c("chr", "start", "end")
  rs <- rgamma(n=nrow(rmap),
	            shape=models$shapes_rec[m],
	            scale=models$scales_rec[m]/models$shapes_rec[m])
  disc_rs <- dt_r[dt_r[J(rs), roll="nearest", which=T]]$lookup_r # discretize
  rmap$r <- disc_rs
  setkey(rmap, start, end)
  fwrite(rmap, paste("rep_", i, "/rmap.csv", sep=""), sep=",")

  mut_spans <- rgeom(n=ceiling(L/models$avg_mut_spans[m]),
	            	     prob=1/models$avg_mut_spans[m])
  while(sum(mut_spans) < L) {
    mut_spans <- c(mut_spans, rgeom(n=1, prob=1/models$avg_mut_spans[m]))
  }
  if(sum(mut_spans) > L) {
    mut_spans <- mut_spans[cumsum(mut_spans) < L]
    mut_spans <- c(mut_spans, L - sum(mut_spans))
  }
  mut_spans <- mut_spans[mut_spans>0]

  mu_dist <- round(seq_log(from=1e-10, to=1e-6, length.out=300), digits=11)
  dt_u <- data.table(mu_dist)
  setkey(dt_u, mu_dist)

  mmap <- setDT(bind_cols(chr="chr1",
	         		  dplyr::lag(cumsum(mut_spans), n=1, default=0),
			          dplyr::lag(cumsum(mut_spans), n=0, default=0)))
  names(mmap) <- c("chr", "start", "end")
  mus <- rgamma(n=nrow(mmap),
	            	shape=models$shapes_mu[m],
	            	scale=models$scales_mu[m]/models$shapes_mu[m])
  disc_mus <- dt_u[dt_u[J(mus), roll="nearest", which=T]]$mu_dist # discretize
  mmap$u <- disc_mus
  setkey(mmap, start, end)
  fwrite(mmap, paste("rep_", i, "/mmap.csv", sep=""), sep=",")

  # finding mu for each sampled pos and the cumulative recombination map
  samp_pos <- unlist(apply(dt_neutral, 1, function(x)
	                  		   seq(from=x[1], to=x[2], by=jump_length))) + 1
  linear_pos <- sort(c(samp_pos, dt_exons$start + models$exon_lengths[m]/2))
  focal_mu <- unlist(lapply(linear_pos, function(pos) mmap[J(pos), roll=T]$u))
  focal_r <- unlist(lapply(linear_pos, function(pos) rmap[J(pos), roll=T]$r))
  cum_rec <- numeric(length=length(linear_pos))
  cum_rec[1] <- focal_r[1]

  for(j in 2:length(linear_pos)) {
    cum_rec[j] <- cum_rec[j-1] + (linear_pos[j] - linear_pos[j-1]) * focal_r[j]
  }

  pos_dt <- setDT(bind_cols(linear_pos, focal_mu, cum_rec))
  names(pos_dt) <- c("position", "focal_mu", "cumrec")
  setkey(pos_dt, position)

  # defines bins for linear models
  bins <- rep(bin_size, floor(L / bin_size))
  bins <- bind_cols(chr="chr1",
		                dplyr::lag(cumsum(bins), n=1, default=0),
	     	            dplyr::lag(cumsum(bins), n=0, default=0))
  names(bins) <- c("chr", "start", "end")
  bins$start <- bins$start + 1 # so that GRanges understands the intervals
  gr_bins <- makeGRangesFromDataFrame(bins)
  
  # binning recombination maps
  rmap2 <- rmap
  rmap2$start <- rmap2$start + 1 # so that GRanges understands the intervals
  rec_gr <- makeGRangesFromDataFrame(rmap2, keep.extra.columns=T)
  hits <- findOverlaps(query=gr_bins, subject=rec_gr) 
  hit_list <- as.data.frame(hits)
  names(hit_list) <- c("bins", "rec_bins")
  rmap2$rec_bins <- 1:nrow(rmap2)
  rec_bins <- merge(rmap2, hit_list) %>%
	                  group_by(bins) %>%
	                  mutate(avg_rec=weighted.mean(r, end-start))
  rec_bins <- rec_bins %>% distinct(bins, .keep_all=T)
  
  # binning mutation maps
  mmap2 <- mmap
  mmap2$start <- mmap2$start + 1 # so that GRanges understands the intervals
  mut_gr <- makeGRangesFromDataFrame(mmap2, keep.extra.columns=T)
  hits <- findOverlaps(query=gr_bins, subject=mut_gr) 
  hit_list <- as.data.frame(hits)
  names(hit_list) <- c("bins", "mut_bins")
  mmap2$mut_bins <- 1:nrow(mmap2)
  mut_bins <- merge(mmap2, hit_list) %>%
	            group_by(bins) %>%
	            mutate(avg_mut=weighted.mean(u, end-start))
  mut_bins <- mut_bins %>% distinct(bins, .keep_all=T)
  
  # binning DFE
  intervals2 <- intervals
  intervals2$start <- intervals2$start + 1 # so that GRanges understands
  dfe_gr <- makeGRangesFromDataFrame(intervals2, keep.extra.columns=T)
  hits <- findOverlaps(query=gr_bins, subject=dfe_gr) 
  hit_list <- as.data.frame(hits)
  names(hit_list) <- c("bins", "dfe")
  intervals2$dfe <- 1:nrow(intervals2)
  # NOTE: avg_s encompasses local density of constrained sites
  dfe_bins <- merge(intervals2, hit_list) %>%
	                  group_by(bins) %>%
	                  mutate(avg_s=weighted.mean(s, end-start)) 
  dfe_bins <- dfe_bins %>% distinct(bins, .keep_all=T)
  
  maps_1kb <- bind_cols(bins,
                        rec_bins["avg_rec"],
                        mut_bins["avg_mut"],
                        dfe_bins["avg_s"])
  
  # pre-computes "effective" rec. dist. between sampled neutral sites and exons
  cr_samp <- pos_dt[position %in% samp_pos,]
  cr_exons <- pos_dt[position %in% (dt_exons$start + models$exon_lengths[m]/2),]
  prd <- as.data.frame(2 * N * abs(outer(cr_samp$cumrec, cr_exons$cumrec, "-")))
  names(prd) <- 1:ncol(prd)
  erd <- prd # "effective" rec. distance is inversely proportional to alpha
  for(j in 1:ncol(erd)) { erd[,j] <- erd[,j] / abs(2 * N * dt_exons$s[j]) }
  relevant_exons <- apply(erd, 2, function(x) x < 10 * N * 1e-3)
  exons_per_samp_neut <- apply(relevant_exons, 1, function(x) which(x))
  
  getB <- function(focal_exon, focal_neutral) { 
    focal_s <- dt_exons[focal_exon,]$s
    total_r <- prd[focal_neutral, focal_exon] / (2 * N) # TODO check if total_r remains reasonable for mpp (< 1e-2) since erd divides by alpha
    # the grid on r is fine-grained enough that we don't need interpolation
    closest_r <- dt_r[dt_r[.(total_r), roll="nearest", which=T]]
    hr <- lookup_tbl[.(closest_r, focal_s)]$Hr # lookup table has single u value
    return((hr / (2 * N * u)) ^ models$exon_lengths[m])
  }
  
  # we finally get the B-values per sampled site! 
  # this takes time if there are many sampled sites influenced by many exons
  B_values <- numeric(length=length(exons_per_samp_neut))
  for(j in 1:length(exons_per_samp_neut)) {
    if(length(exons_per_samp_neut[[j]]) > 0) {
       B <- unlist(lapply(exons_per_samp_neut[[j]], getB, focal_neutral=j))
       B_values[j] <- cumprod(B)[length(B)]
    }
    else { # if no exons within relevant distance
      B_values[j] <- 1
    }
  }
  
  # using max rec distance (r==1e-2) to check validity of threshold used above
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
  
  # getting B-values per site can also take some time:
  B_map <- cubicspline(samp_pos, B_values, 1:L)
 
  # here it's easier to use mmap2 which has closed intervals:
  ones <- rep(1, L)
  pi0s <- invisible(apply(mmap2, 1, function(x)
			             ones[as.numeric(x[2]):as.numeric(x[3])] *
			                  as.numeric(x[4]) * 2 * N))
  pis <- invisible(apply(mmap2, 1, function(x)
			             B_map[as.numeric(x[2]):as.numeric(x[3])] *
			                   as.numeric(x[4]) * 2 * N))
  pi0s <- unlist(pi0s)
  pis <- unlist(pis)
  hrmap <- setDT(bind_cols(1:L, pis, pi0s, B_map))
  names(hrmap) <- c("Pos", "Hr", "pi0", "B")
  setkey(hrmap, Pos)
  fwrite(hrmap, paste("rep_", i, "/hrmap.csv", sep=""), compress="gzip")
  
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
  
  hrmap$bin <- ((hrmap$Pos - 1) %/% bin_size)
  hrmap <- hrmap[1:maps_1kb$end[nrow(maps_1kb)],] # trims the tail (convenience)
  avgs <- hrmap %>% group_by(bin) %>%
          summarize(avg_pi=mean(Hr), avg_pi0=mean(pi0))
  maps_1kb$avg_pi <- avgs$avg_pi
  maps_1kb$bin <- 1:nrow(maps_1kb)
  
  # standardizing variables to help interpretation of linear coefficients
  maps_1kb$std_s <- (maps_1kb$avg_s - mean(maps_1kb$avg_s)) / 
                    sd(maps_1kb$avg_s)
  maps_1kb$std_mu <- (maps_1kb$avg_mut - mean(maps_1kb$avg_mut)) /
                     sd(maps_1kb$avg_mut)
  maps_1kb$std_r <- (maps_1kb$avg_rec - mean(maps_1kb$avg_rec)) / 
                    sd(maps_1kb$avg_rec)
  maps_1kb$std_pi <- (maps_1kb$avg_pi - mean(maps_1kb$avg_pi)) / 
                     sd(maps_1kb$avg_pi)
  
  m_1kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s, data=maps_1kb)
  m_1kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_1kb)
  
  summary(m_1kb_1)
  summary(m_1kb_2)
  
  comp <- AIC(m_1kb_1, m_1kb_2) # model comparison
  
  anova.pi <- Anova(m_1kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_1kb$Total[i] <- (anova.pi$VarExp[1] +
			                    anova.pi$VarExp[2] +
			                    anova.pi$VarExp[3] +
			                    anova.pi$VarExp[4]) * 100
  r2_tbl_1kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_1kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_1kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_1kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  
  # builds maps at larger genomic scales
  maps_1kb$bin_10kb <- (maps_1kb$bin - 1) %/% 10
  maps_1kb$bin_100kb <- (maps_1kb$bin - 1) %/% 100
  maps_10kb <- maps_1kb %>%
	       group_by(bin_10kb) %>%
	       summarise_at(c("std_pi", "std_s", "std_mu", "std_r",
	        	            "avg_pi", "avg_s", "avg_mut", "avg_rec"), mean)

  maps_100kb <- maps_1kb %>%
	        group_by(bin_100kb) %>%
	       	summarise_at(c("std_pi", "std_s", "std_mu", "std_r",
			                   "avg_pi", "avg_s", "avg_mut", "avg_rec"), mean)
  
  # 10 kb
  m_10kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s, data=maps_10kb)
  m_10kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_10kb)
  
  summary(m_10kb_1)
  summary(m_10kb_2)
  
  AIC(m_10kb_1, m_10kb_2) # model comparison
  
  anova.pi <- Anova(m_10kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_10kb$Total[i] <- (anova.pi$VarExp[1] + 
			                     anova.pi$VarExp[2] +
			                     anova.pi$VarExp[3] +
			                     anova.pi$VarExp[4]) * 100
  r2_tbl_10kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_10kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_10kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_10kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  
  # 100 kb
  m_100kb_1 <- lm(std_pi ~ std_mu + std_r + std_s + std_r:std_s,data=maps_100kb)
  m_100kb_2 <- lm(std_pi ~ std_mu + std_r, data=maps_100kb)
  
  summary(m_100kb_1)
  summary(m_100kb_2)
  
  AIC(m_100kb_1, m_100kb_2) # model comparison
  
  anova.pi <- Anova(m_100kb_1)
  apiss <- anova.pi$"Sum Sq"
  anova.pi$VarExp <- apiss / sum(apiss)
  
  r2_tbl_100kb$Total[i] <- (anova.pi$VarExp[1] + 
			                      anova.pi$VarExp[2] +
			                      anova.pi$VarExp[3] +
			                      anova.pi$VarExp[4]) * 100
  r2_tbl_100kb$u[i] <- anova.pi$VarExp[1] * 100
  r2_tbl_100kb$r[i] <- anova.pi$VarExp[2] * 100
  r2_tbl_100kb$s[i] <- anova.pi$VarExp[3] * 100
  r2_tbl_100kb$`r:s`[i] <- anova.pi$VarExp[4] * 100
  
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
  
  # 10 kb
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
  
  # 100 kb
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
}

fwrite(r2_tbl_1kb, "r2_1kb.csv")
fwrite(r2_tbl_10kb, "r2_10kb.csv")
fwrite(r2_tbl_100kb, "r2_100kb.csv")

