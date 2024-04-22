# Date created: 01/04/2024
# Last modified: 01/04/2024
# Author: Gustavo V. Barroso

pdf(NULL) # to suppress creation of Rplots.pdf

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(pracma) # for cubicspline()
  library(bigsnpr) # for seq_log
  library(cowplot)
  library(scales)
  library(RColorBrewer)
  library(shiny)
  library(bslib)
})

setwd("~/Data/momentspp/paper_1/single_neutral")

###############################
#
# Set up
#
###############################

params <- fread("params.csv")

Na <- unique(params$Na)
N1 <- unique(params$N1)
t <- unique(params$t)
u <- unique(params$u)

demo_models <- crossing(u, Na, N1, t)
num_demo_models <- nrow(demo_models)
sampling_times <- seq(from=t, to=0, by=-1000)
pi0 <- as.data.frame(matrix(nrow=num_demo_models, ncol=length(sampling_times)))
names(pi0) <- as.character(sampling_times)

for(i in 1:num_demo_models) {
  
  moms <- read.csv(paste("demo/model_", i, 
                      "_e_1_expectations.txt", sep=""), sep=" ", header=F)
  
  pi0[i,] <- t(dplyr::select(filter(moms, V1=="Hl_0_0"), V3))
}

demo_pi0 <- cbind.data.frame(demo_models, pi0)
fwrite(demo_pi0, "demo/demo_pi0.csv")

m_pi0 <- pivot_longer(demo_pi0, 
                      cols=as.character(sampling_times),
                      names_to="Generation", values_to="pi0")
m_pi0$Generation <- as.numeric(m_pi0$Generation)

N.labs <- paste("Model ", 1:num_demo_models)
names(N.labs) <- rep(as.character(unique(demo_models$N1)),
                     length(unique(demo_models$u)))

plot_list <- list(length=length(u))
c <- 1
for(mu in u) {
  p <- ggplot(data=filter(m_pi0, u==mu), aes(x=Generation, y=pi0)) + 
    facet_wrap(~N1, labeller=labeller(N1=N.labs[1:2])) +
    geom_point(size=0.5) + geom_line() + theme_bw() + 
    labs(title=paste("Trajectory of Heterozygosity after size change, u =", mu), 
         x="Generations ago", y=expression(pi[0])) +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          strip.text=element_text(size=14),
          legend.position="bottom")
  
  plot_list[[c]] <- p
  c <- c + 1
}

pi0_demo <- plot_grid(plotlist=plot_list, ncol=length(u))
save_plot("demo/pi0_demo.png", pi0_demo, base_height=8, base_width=16)

d <- ggplot(data=demo_models) + theme_bw() + 
     facet_wrap(~N1, labeller=labeller(N1=N.labs)) +
     geom_segment(aes(x=0, xend=t, y=N1, yend=N1), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t, y=N1, yend=Na), linewidth=1.5) +
     geom_segment(aes(x=t, xend=t+t/10, y=Na, yend=Na), linewidth=1.5) +
     scale_y_log10() +
     scale_x_continuous(breaks=c(0, 1000, t)) +
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

save_plot("demo/demo_models.png", d, base_height=8, base_width=16)

look_s <- setDT(as.data.frame(unique(params$lookup_s)))
look_r <- setDT(as.data.frame(unique(params$lookup_r)))

names(look_r) <- "lookup_r"
names(look_s) <- "lookup_s"

setkey(look_r, lookup_r)
setkey(look_s, lookup_s)

hl_demo <- fread("hl_time.csv.gz")
hr_demo <- fread("hr_time.csv.gz")

###############################
#
# Single constrained locus
#
###############################

m_het_demo <- pivot_longer(hr_demo, cols=as.character(sampling_times),
                          names_to="Generation", values_to="Hr")
m_het_demo$Generation <- as.numeric(m_het_demo$Generation)
m_het_demo <- setDT(m_het_demo)
setkey(m_het_demo, N1, lookup_r, lookup_s, Generation)

m_het_demo <- left_join(m_het_demo, dplyr::select(m_pi0,
                        -c(Na, t)), by=c("N1", "Generation", "u"))
m_het_demo$B <- (m_het_demo$Hr / m_het_demo$pi0) ^ 1000 # exons have L sites
m_het_demo$scaled_Hr <- m_het_demo$pi0 * m_het_demo$B
  
hl_gen <- pivot_longer(hl_demo, cols=as.character(sampling_times),
                       names_to="Generation", values_to="Hl")
hl_gen$Generation <- as.numeric(hl_gen$Generation)
hl_gen <- setDT(hl_gen)
hl_gen <- arrange(hl_gen, scenario, Generation)
setkey(hl_gen, N1, lookup_r, lookup_s, Generation)

m_het_demo$Hl <- hl_gen$Hl
m_het_demo$piN_pi0 <- m_het_demo$Hl / m_het_demo$pi0
m_het_demo$piN_piS <- m_het_demo$Hl / m_het_demo$scaled_Hr

fwrite(m_het_demo, "stats_demo.csv")

m_demo <- pivot_longer(m_het_demo, cols=c(pi0, Hr, scaled_Hr, 
                       Hl, B, piN_pi0, piN_piS), names_to="statistic")

svals <- unique(m_demo$lookup_s)
for(mom in unique(m_demo$statistic)) {
  
  plot_list <- list(length=length(svals))
  c <- 1
  
  for(s in svals) {
    p <- ggplot(data=filter(m_demo, statistic==mom, u==1e-8,
                            lookup_s==s, lookup_r==1e-8),
                aes(x=Generation, y=value)) + 
      geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(N1)) +
      scale_x_continuous(breaks=pretty_breaks()) +
      scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none")
    if(c==length(svals)) {
      p <- p + labs(title=NULL, x="Generation", 
                    y=paste(mom, "(s=", s, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else if(c==1) {
      p <- p + labs(title=paste("Temporal dynamics after a size change"),
                    x=NULL, y=paste(mom, "(s=", s, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else {
      p <- p + labs(title=NULL, x=NULL, y=paste(mom, "(s=", s, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    }
    
    plot_list[[c]] <- p
    c <- c + 1
  }
  
  mom_plot <- plot_grid(plotlist=plot_list, ncol=1, align='v')
  save_plot(paste(mom, "_time_u_", 1e-8, ".png", sep=""),
            mom_plot, base_height=10, base_width=12)
}

het_list <- list(length=length(svals))
c <- 1

for(s in svals) {
  p <- ggplot(data=filter(m_demo, statistic==c("scaled_Hr", "pi0"),
                          u==1e-8, lookup_s==s, lookup_r==1e-8),
              aes(x=Generation, y=value, color=statistic)) + 
    geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(N1)) +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_log10(breaks=pretty_breaks()) + guides(alpha="none")
    labs(title="Temporal dynamics after a size change",
         x="Generation", y="Value") +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
    
  het_list[[c]] <- p
  c <- c + 1
}

het_plot <- plot_grid(plotlist=het_list, ncol=1, align='v')
save_plot(paste("het_time_u_", 1e-8, ".png", sep=""),
          het_plot, base_height=10, base_width=12)

r <- ggplot(data=filter(m_demo, Generation==50000,
                        u==1e-8, N1==1e+5, statistic=="B"),
       aes(x=lookup_r, y=value)) + 
  geom_point() + theme_bw() + geom_line() + 
  facet_wrap(~as.factor(lookup_s), nrow=1) +
  scale_x_log10(breaks=unique(m_demo$lookup_r)) +
  scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none") +
  labs(title=NULL, x="r", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot(paste("B_r_u_", 1e-8, ".png", sep=""), r,
          base_height=7, base_width=10)

d_stats <- filter(m_het_demo, u==1e-8) %>% 
     group_by(scenario) %>% 
     reframe(N1=N1, lookup_r=lookup_r, lookup_s=lookup_s, Generation=Generation,
             d_pi0=-c(diff(pi0)[1], diff(pi0))/pi0,
             d_Hl=-c(diff(Hl)[1], diff(Hl))/Hl,
             d_Hr=-c(diff(scaled_Hr)[1], diff(scaled_Hr))/scaled_Hr,
             d_B=-c(diff(B)[1], diff(B))/B,
             d_piN_pi0=-c(diff(piN_pi0)[1], diff(piN_pi0))/piN_pi0,
             d_piN_piS=-c(diff(piN_piS)[1], diff(piN_piS))/piN_piS)

m_demo_d <- pivot_longer(d_stats, cols=starts_with("d_"), names_to="statistic")

c <- 1
plot_list <- list(length=length(unique(m_demo_d$N1)))

for(N in unique(m_demo_d$N1)) {
  p <- ggplot(data=filter(m_demo_d, N1==N, lookup_r==1e-8,
                          lookup_s %in% c(-1e-3, -1e-4, -1e-5), 
                          statistic %in% c("d_Hl", "d_Hr", "d_pi0")),
              aes(x=Generation, y=value, color=statistic)) + 
    geom_point() + theme_bw() + geom_line() + facet_wrap(~as.factor(lookup_s)) +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_continuous(breaks=pretty_breaks()) + guides(alpha="none")
    if(c==length(unique(m_demo_d$N1))) {
      p <- p + labs(title=NULL, x="Generation", y="Rate of Change") +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="bottom")
    } else if(c==1) {
      p <- p + labs(title="Temporal dynamics of statistics after a size change",
                    x=NULL, y="Rate of change") +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else {
      p <- p + labs(title=NULL, x=NULL, y="Rate of change") +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    }
  
  plot_list[[c]] <- p
  c <- c + 1
}

cp <- plot_grid(plotlist=plot_list, ncol=1, align='v')
save_plot(paste("diffs_stats_time_u_", 1e-8, ".png", sep=""),
          cp, base_height=10, base_width=12)

###############################
#
# multiple constrained loci part 1
#
###############################

# reset
m_het_demo$B <- (m_het_demo$Hr / m_het_demo$pi0) ^ 1 
m_demo <- pivot_longer(m_het_demo, cols=c(pi0, Hr, scaled_Hr, 
                       Hl, B, piN_pi0, piN_piS), names_to="statistic")

# assembles lookup tables with continuous recombination from 1e-8 to 1e-2
Bs_bottleneck <- filter(m_demo, N1==1e+3, statistic=="B")
interp_bottleneck <- data.frame()
  
for(mu in unique(Bs_bottleneck$u)) {
  for(s in unique(Bs_bottleneck$lookup_s)) {
    print(Sys.time())
    cat(paste(s, "\n"))
    for(gen in unique(Bs_bottleneck$Generation)) {
      tmp <- filter(Bs_bottleneck, lookup_s==s, Generation==gen, u==mu)
      Bs_interp <- cubicspline(tmp$lookup_r, tmp$value,
                               seq_log(from=1e-8, to=1e-2, length.out=1000))
      tmp <- cbind.data.frame(seq_log(from=1e-8, to=1e-2, length.out=1000),
                              Bs_interp, s, gen, mu)
      names(tmp) <- c("r", "B", "s", "Generation", "u")
      
      interp_bottleneck <- rbind.data.frame(interp_bottleneck, tmp)
    }
  }
}

dt_bottleneck <- setDT(interp_bottleneck)
setkey(dt_bottleneck, r, s)

Bs_expansion <- filter(m_demo, N1==1e+5, statistic=="B")
interp_expansion <- data.frame()

for(mu in unique(Bs_expansion$u)) {
  for(s in unique(Bs_expansion$lookup_s)) {
    print(Sys.time())
    cat(paste(s, "\n"))
    for(gen in unique(Bs_expansion$Generation)) {
      tmp <- filter(Bs_expansion, lookup_s==s, Generation==gen, u==mu)
      Bs_interp <- cubicspline(tmp$lookup_r, tmp$value,
                               seq_log(from=1e-8, to=1e-2, length.out=1000))
      tmp <- cbind.data.frame(seq_log(from=1e-8, to=1e-2, length.out=1000),
                              Bs_interp, s, gen, mu)
      names(tmp) <- c("r", "B", "s", "Generation", "u")
      
      interp_expansion <- rbind.data.frame(interp_expansion, tmp)
    }
  }
}

dt_expansion <- setDT(interp_expansion)
setkey(dt_expansion, r, s)

# thinning for speed
bt_df <- filter(dt_bottleneck, Generation %in% seq(from=0, to=50000, by=1000))
exp_df <- filter(dt_expansion, Generation %in% seq(from=0, to=50000, by=1000))

# set up chr
num_exons <- 100
ncsl <- rep(0, num_exons)
exon_lengths <- 1000
csl <- rep(exon_lengths, num_exons) 
L <- sum(csl) + sum(ncsl)

maps <- suppressWarnings(c(rbind(csl, ncsl)))
maps <- suppressMessages(setDT(bind_cols(chr="chr1",
    dplyr::lag(cumsum(maps), n=1, default=0),
    dplyr::lag(cumsum(maps), n=0, default=0))))
names(maps) <- c("chr", "start", "end")
maps$s <- c(rbind(rep(1, length(ncsl)), 0)) # indicator, whether is selected
maps$r <- 1e-8 * (maps$end - maps$start)
maps$cum_r <- cumsum(maps$r)
setkey(maps, start, end)

dt_neut <- filter(maps, s==0)
dt_exons <- filter(maps, s==1)

# DFE proportions
dt_exons$p1 <- runif(nrow(dt_exons), min=0, max=1)
p2 <- numeric(length=nrow(dt_exons))
for(i in 1:nrow(dt_exons)) {
  p2[i] <- runif(1, min=0, max=1-dt_exons$p1[i])
}
dt_exons$p2 <- p2
dt_exons$p3 <- 1 - dt_exons$p2 - dt_exons$p1

dt_exons$p1 <- 1/3
dt_exons$p2 <- 1/3
dt_exons$p3 <- 1/3
dt_exons$avg_s <- dt_exons$p1 * -1e-3 + dt_exons$p2 * -1e-4 + dt_exons$p3 * -1e-5

getB <- function(exon_id, samp_id, look_tbl) { # id, id, data.table for generation
  focal_r <- abs(dt_neut[samp_id]$cum_r - dt_exons[exon_id]$cum_r) 
  if(focal_r > 1e-2) { 
    return(1)
  } else {
    focal_B <- look_tbl[look_tbl[.(focal_r), roll="nearest", which=T]]$B
    B <- focal_B ^ exon_lengths # (exon_lengths * dt_exons$p1[exon_id])) NOTE
    
    return(B)
  }
}

bmaps_demo <- data.frame()

for(mu in unique(bt_df$u)) {
  
  # one table per selection strength
  B_maps_bt_strong <- as.data.frame(matrix(ncol=nrow(dt_neut),
      nrow=length(unique(bt_df$Generation))))
  
  B_maps_exp_strong <- as.data.frame(matrix(ncol=nrow(dt_neut),
      nrow=length(unique(exp_df$Generation))))
  
  B_maps_bt_moderate <- B_maps_bt_strong
  B_maps_bt_weak <- B_maps_bt_strong
  
  B_maps_exp_moderate <- B_maps_exp_strong
  B_maps_exp_weak <- B_maps_exp_strong
  
  c <- 1
  for(g in sort(unique(bt_df$Generation), decreasing=F)) {
    
    print(Sys.time())
    cat(paste(g, "\n"))
    
    bt_tmp <- filter(bt_df, Generation==g, u==mu)
    exp_tmp <- filter(exp_df, Generation==g, u==mu)
    
    bt_strong <- filter(bt_tmp, s==-1e-3)
    bt_moderate <- filter(bt_tmp, s==-1e-4)
    bt_weak <- filter(bt_tmp, s==-1e-5)
    
    exp_strong <- filter(exp_tmp, s==-1e-3)
    exp_moderate <- filter(exp_tmp, s==-1e-4)
    exp_weak <- filter(exp_tmp, s==-1e-5)
  
    for(i in 1:nrow(dt_neut)) { # for each sampled position
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=bt_strong)
      B_bt_strong <- cumprod(b)[length(b)]
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=bt_moderate)
      B_bt_moderate <- cumprod(b)[length(b)]
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=bt_weak)
      B_bt_weak <- cumprod(b)[length(b)]
      
      B_maps_bt_strong[c, i] <- B_bt_strong
      B_maps_bt_moderate[c, i] <- B_bt_moderate
      B_maps_bt_weak[c, i] <- B_bt_weak
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=exp_strong)
      B_exp_strong <- cumprod(b)[length(b)]
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=exp_moderate)
      B_exp_moderate <- cumprod(b)[length(b)]
      
      b <- sapply(1:nrow(dt_exons), getB, samp_id=i, look_tbl=exp_weak)
      B_exp_weak <- cumprod(b)[length(b)]
      
      B_maps_exp_strong[c, i] <- B_exp_strong
      B_maps_exp_moderate[c, i] <- B_exp_moderate
      B_maps_exp_weak[c, i] <- B_exp_weak
    }
    
    c <- c + 1
  }

  names(B_maps_bt_strong) <- dt_neut$end
  B_maps_bt_strong$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_bt_strong$N1 <- 1e+3
  B_maps_bt_strong$s <- -1e-3
  
  names(B_maps_bt_moderate) <- dt_neut$end
  B_maps_bt_moderate$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_bt_moderate$N1 <- 1e+3
  B_maps_bt_moderate$s <- -1e-4
  
  names(B_maps_bt_weak) <- dt_neut$end
  B_maps_bt_weak$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_bt_weak$N1 <- 1e+3
  B_maps_bt_weak$s <- -1e-5
  
  B_maps_bottleneck <- rbind.data.frame(B_maps_bt_strong,
                                        B_maps_bt_moderate,
                                        B_maps_bt_weak)

  names(B_maps_exp_strong) <- dt_neut$end
  B_maps_exp_strong$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_exp_strong$N1 <- 1e+5
  B_maps_exp_strong$s <- -1e-3
  
  names(B_maps_exp_moderate) <- dt_neut$end
  B_maps_exp_moderate$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_exp_moderate$N1 <- 1e+5
  B_maps_exp_moderate$s <- -1e-4
  
  names(B_maps_exp_weak) <- dt_neut$end
  B_maps_exp_weak$Generation <- sort(unique(bt_df$Generation), decreasing=F)
  B_maps_exp_weak$N1 <- 1e+5
  B_maps_exp_weak$s <- -1e-5
  
  B_maps_expansion <- rbind.data.frame(B_maps_exp_strong,
                                       B_maps_exp_moderate,
                                       B_maps_exp_weak)
  bmaps_demo_u <- rbind.data.frame(B_maps_bottleneck, B_maps_expansion)
  bmaps_demo_u$u <- mu
  
  bmaps_demo <- rbind.data.frame(bmaps_demo, bmaps_demo_u)
}

m_Bmap_demo <- pivot_longer(bmaps_demo, cols=as.character(dt_neut$end), 
                            names_to="Position")
m_Bmap_demo$Position <- as.numeric(m_Bmap_demo$Position)

fwrite(bmaps_demo, "bmaps_demo.csv")
fwrite(m_Bmap_demo, "m_Bmap_demo.csv")

m_Bmap_demo <- fread("m_Bmap_demo.csv")
bmaps_demo <- fread("bmaps_demo.csv")

ui <- page_sidebar(
  title="B-value maps over time",
  sidebar=sidebar(
    sliderInput(inputId="Generation",
    label="Generation:",
    min=0,
    max=50000,
    step=1000,
    value=50000),
    checkboxInput("by_demography", "Show Ne", T),
    checkboxInput("by_selection", "Show s", T),
    checkboxInput("by_mutation", "Show u", T),
    checkboxGroupInput(
      "N1", "Ne (color)",
      choices=unique(m_Bmap_demo$N1), 
      selected=unique(m_Bmap_demo$N1[1])
    ),
    checkboxGroupInput(
      "s", "s (shape)",
      choices=unique(m_Bmap_demo$s), 
      selected=unique(m_Bmap_demo$s[1])
    ),
    checkboxGroupInput(
      "u", "u",
      choices=unique(m_Bmap_demo$u), 
      selected=unique(m_Bmap_demo$u[1])
    )
  ),
  plotOutput(outputId="Bmap")
)

server <- function(input, output) {
  output$Bmap <- renderPlot({
    p <- ggplot(data=filter(m_Bmap_demo, 
                            N1==input$N1,
                            s==input$s,
                            u==input$u,
                            Generation==input$Generation),
                aes(x=Position, y=value)) + 
      list(
        geom_point(), geom_line(), theme_bw(),
        scale_y_continuous(breaks=pretty_breaks()),
        if(input$by_demography) aes(color=as.factor(N1)),
        if(input$by_selection) aes(shape=as.factor(s)),
        labs(title=NULL, x="Position", y="B"),
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
      )
    
    p
  }, res = 100)
}

shinyApp(ui, server)

for(sel in svals) {
  c <- ggplot(data=filter(m_Bmap_demo, s==sel),
              aes(x=Position, y=value, color=Generation)) +
    geom_point() + facet_wrap(~N1) + theme_bw() +
    scale_color_continuous(name="Generation") +
    scale_y_continuous(breaks=pretty_breaks()) +
    labs(title="B-value maps over time", x="Position", y="B") +
    theme(axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
  
  save_plot(paste("Bmap_time_s_", sel, ".png", sep=""), c, 
            base_height=10, base_width=16)
}

d <- ggplot(data=filter(m_Bmap_demo, u==1e-8,
                        Position %in% c(1e+3, 1e+4, 2.5e+4, 5e+4)),
            aes(x=Generation, y=value, color=as.factor(Position))) +
  geom_point() + theme_bw() + facet_grid(N1~s) + geom_line() +
  scale_color_discrete(name="Position") +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="B-values over time", x="Generation", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("Bpos_time_u_1e-8.png", d, base_height=8, base_width=12)

# adds simulation results
sim_exp_strong <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.001.txt")
sim_exp_moderate <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.0001.txt")
sim_exp_weak <- fread("fwdpy11_Bvalues.expansion.Q4.s_-0.00001.txt")

sim_exp_strong <- dplyr::select(sim_exp_strong, -V1)
sim_exp_moderate <- dplyr::select(sim_exp_moderate, -V1)
sim_exp_weak <- dplyr::select(sim_exp_weak, -V1)

names(sim_exp_strong) <- as.character(dt_neut$end)
names(sim_exp_moderate) <- as.character(dt_neut$end)
names(sim_exp_weak) <- as.character(dt_neut$end)

sim_exp_tbl <- rbind.data.frame(sim_exp_strong, sim_exp_moderate, sim_exp_weak)
sim_exp_tbl$Generation <- rep(seq(from=0, to=50000, by=1000), 3)
sim_exp_tbl$N1 <- 1e+5
sim_exp_tbl$s <- c(rep(-1e-3, 51), rep(-1e-4, 51), rep(-1e-5, 51))
sim_exp_tbl$u <- 1e-8
#sim_exp_tbl$Q <- 4

m_sim_exp <- pivot_longer(sim_exp_tbl, cols=1:100, names_to="Position")
m_sim_exp$Position <- as.numeric(m_sim_exp$Position)
m_sim_exp$method <- "fwdpy11"

y <- ggplot(data=filter(m_sim_exp, Position %in% c(1e+3, 1e+4, 2.5e+4, 5e+4)),
            aes(x=Generation, y=value, color=as.factor(Position))) +
  geom_point(size=2) + theme_bw() + facet_wrap(~s) + 
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="Simulations w/ expansion, u=1e-8", x="Generation", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("fwd.png", y, base_height=8, base_width=12)

m_Bmap_demo$method <- "mpp"
#bmaps_demo <- relocate(bmaps_demo, c(Generation, N1), .after=last_col())
#names(sim_exp_tbl) <- names(bmaps_demo)
df <- rbind.data.frame(filter(bmaps_demo, N1==1e+5, u==1e-8), sim_exp_tbl)

m_df <- pivot_longer(df, cols=as.character(dt_neut$end), names_to="Position")
m_df$Position <- as.numeric(m_df$Position)

sampled_gens <- seq(from=50000, to=0, by=-10000)
plot_list <- list(length=length(sampled_gens))
c <- 1
for(g in sampled_gens) {
  e <- ggplot(data=filter(m_df, Generation==g), 
              aes(x=Position/1e+3, y=value, color=method)) +
    geom_point() + facet_wrap(~s) + theme_bw() + 
    scale_y_continuous(breaks=pretty_breaks())
    if(c==length(sampled_gens)) {
      e <- e + labs(title=NULL, x="Position (kb)", 
                    y=paste("B(gen=", g, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_text(size=12),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="bottom")
    } else if(c==1) {
      e <- e + labs(title=paste("B-maps over generations"),
                    x=NULL, y=paste("B(gen=", g, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    } else {
      e <- e + labs(title=NULL, x=NULL, y=paste("B(gen=", g, ")", sep="")) +
        theme(axis.title=element_text(size=16), 
              axis.text=element_text(size=12), 
              axis.text.x=element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              legend.position="none")
    }
  
  plot_list[[c]] <- e
  c <- c + 1
}

f <- plot_grid(plotlist=plot_list, ncol=1, align='v',
               rel_heights=c(1.15, 1, 1, 1, 1, 1.65))
save_plot("Bmap_mpp_sims_1.png", f, base_height=16, base_width=14)

g <- ggplot(data=filter(m_df, Position %in% c(1e+3, 1e+4, 2.5e+4, 5e+4)),
            aes(x=Generation, y=value, color=as.factor(Position), shape=method)) +
  geom_point(size=2) + theme_bw() + facet_wrap(~s) + geom_line() +
  scale_color_discrete(name="Position") +
  scale_shape_manual(values=c(0, 1)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title="B-value maps over time", x="Generation", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")

save_plot("Bmap_mpp_sims_2.png", g, base_height=10, base_width=16)

x <- filter(m_Bmap_demo, Generation==5e+4)
ggplot(data=filter(x, s==-1e-3), aes(x=Position, y=value, shape=as.factor(N1))) +
  geom_point(size=2) + theme_bw() + facet_wrap(~s) + geom_line() +
  scale_shape_manual(values=c(0, 1)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Pos", y="B") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
#####################################
#
# multiple constrained loci part 2
#
#####################################

num_exons <- 1001
ncsl <- rep(50000, num_exons)
exon_lengths <- 1000
csl <- rep(exon_lengths, num_exons) 
L <- sum(csl) + sum(ncsl)

for(N in unique(hr_demo$N1)) {
  
  Bvals_gen <- as.data.frame(matrix(ncol=length(sampling_times), 
                                    nrow=nrow(look_s)))
  names(Bvals_gen) <- as.character(sampling_times)
  
  for(i in 1:nrow(look_s)) {
    
    print(Sys.time())
    cat(paste(i, "\n"))
    
    ss <- rep(as.numeric(look_s[i]), num_exons)
  
    smap <- suppressWarnings(c(rbind(ncsl, csl)))
    smap <- suppressMessages(setDT(bind_cols(chr="chr1",
                                 dplyr::lag(cumsum(smap), n=1, default=0),
                                 dplyr::lag(cumsum(smap), n=0, default=0))))
    names(smap) <- c("chr", "start", "end")
    smap$s <- c(rbind(rep(0, length(csl)), ss))
    setkey(smap, start, end)
  
    dt_neutral <- filter(smap[,2:4], s==0)
    dt_exons <- filter(smap[,2:4], s<0)
    
    rmap <- setDT(bind_cols(chr="chr1", start=0, end=L, r=1e-8)) # flat
    mmap <- setDT(bind_cols(chr="chr1", start=0, end=L, u=1e-8)) # flat
  
    setkey(mmap, start, end)
    setkey(rmap, start, end)
    setkey(smap, start, end)
    setkey(dt_neutral, start, end)
    setkey(dt_exons, start, end)
  
    all_pos <- sort(c(dt_neutral[nrow(dt_neutral)/2]$start + 500, 
                       dt_exons$start + exon_lengths / 2))
    focal_mu <- rep(1e-8, length(all_pos))
    
    pos_dt <- suppressMessages(setDT(bind_cols(all_pos, focal_mu)))
    names(pos_dt) <- c("position", "focal_mu")
    pos_dt$idx <- 1:nrow(pos_dt) # indexing
    setkey(pos_dt, position)
    
    Nanc <- unique(m_het_demo$Na)
    m_hr_N <- filter(m_het_demo, N1==N) 
  
    c <- 1
    for(g in sort(unique(m_hr_N$Generation), decreasing=T)) {
  
      Ne_bar <- unique(filter(m_hr_N, Generation==g)$pi0) / 
                (2*unique(m_hr_N$u))
      
      #num_iter <- 1 NOTE: the following is without interference correction ATM
      B_values <- rep(1, length(all_pos)) # init 
      
      # approximates cumulative rec by looking at r at sampled sites only
      focal_r <- unlist(lapply(all_pos, function(pos) rmap[J(pos), roll=T]$r)) 
      focal_r <- focal_r * B_values
      cum_rec <- numeric(length=length(all_pos))
      cum_rec[1] <- focal_r[1]
      
      for(j in 2:length(all_pos)) {
        cum_rec[j] <- cum_rec[j-1] + (all_pos[j] - all_pos[j-1]) * focal_r[j]
      }
      pos_dt$cumrec <- cum_rec
      
      # pre-computes "effective" rec. dist. between sampled sites and exons
      cr_exons <- pos_dt[position %in% (dt_exons$start + exon_lengths / 2),]
    
      prd <- 1
      if(N <= Nanc) {
        prd <- as.data.frame(2*N*abs(outer(pos_dt$cumrec,
                                           cr_exons$cumrec, "-")))
      } else {
        prd <- as.data.frame(2*Ne_bar*abs(outer(pos_dt$cumrec,
                                                cr_exons$cumrec, "-")))
      }
      names(prd) <- 1:ncol(prd)
      
      erd <- prd # "effective" rec. distance is inversely proportional to alpha
      for(j in 1:ncol(erd)) { erd[,j] <- erd[,j] / 
        abs(2 * Ne_bar * dt_exons$s[j]) }
      relevant_exons <- apply(erd, 2, function(x) x < 10 * Ne_bar * 1e-3)
      # identifying relevant exons (w.r.t. linked selection) for each sampled site
      exons_per_samp_site <- apply(relevant_exons, 1, function(x) which(x))
      # check if class(exons_per_samp_site) [1] "matrix" "array" 
      
      getB <- function(focal_exon, focal_samp) { # arguments are site indices
        total_r <- prd[focal_samp, focal_exon] / (2 * Ne_bar) 
        if(total_r > 1e-2) { return(1) }
        else {
          exon_pos <- dt_exons$start[focal_exon] + exon_lengths / 2 # midpoint
          idx <- pos_dt[.(exon_pos)]$idx # index within pos_dt / all_pos
          dt_exons[focal_exon,]$s <- dt_exons[focal_exon,]$s * B_values[idx]
          focal_s <- dt_exons[focal_exon,]$s
          
          closest_r <- look_r[look_r[.(total_r), roll="nearest", which=T]]
  
          # both Hr and pi0 are scaled linearly by Ne_bar -> need not include it
          hr <- m_hr_N[.(N, closest_r, focal_s, g)]$Hr 
          pi0 <- m_hr_N[.(N, closest_r, focal_s, g)]$pi0 
          return((hr / pi0) ^ (exon_lengths * B_values[idx]))
        }
      }
      
      tmp <- B_values 
      samp_neut <- filter(pos_dt, 
          position==dt_neutral[nrow(dt_neutral)/2]$start + 500)$idx
      
      for(k in samp_neut) {
        if(length(exons_per_samp_site[[k]]) > 0) {
          B <- unlist(lapply(exons_per_samp_site[[k]], getB, focal_samp=k))
          tmp[k] <- cumprod(B)[length(B)]
        }
        else {
          tmp[k] <- 1 
        }
      }
  
      B_values <- tmp
      Bvals_gen[i, c] <- B_values[samp_neut]
      
      c <- c + 1
    }
  }
  
  Bvals_gen$s <- look_s$lookup_s
  m_Bs <- pivot_longer(Bvals_gen, cols=as.character(sampling_times))
  
  tbl_gen <- pivot_longer(filter(hl_demo, N1==N, lookup_r==2.001000e-05), 
                          cols=as.character(sampling_times),
                          names_to="Generation", values_to="Hl")
  tbl_gen$Generation <- as.numeric(tbl_gen$Generation)
  tbl_gen <- setDT(tbl_gen)
  setkey(tbl_gen, N1, lookup_r, lookup_s, Generation)
  
  tbl_gen$pi0 <- rep(arrange(filter(m_pi0, N1==N), Generation)$pi0, 7)
  tbl_gen$B <- m_Bs$value
  
  fwrite(tbl_gen, paste("tbl_gen_N1_", N, ".csv", sep=""))
}

tbl_bottleneck <- fread("tbl_gen_N1_1000.csv")
tbl_expansion <- fread("tbl_gen_N1_1e+05.csv")
tbl_demo <- rbind.data.frame(tbl_bottleneck, tbl_expansion)
tbl_demo$Hr <- tbl_demo$B * tbl_demo$pi0
tbl_demo$piN_pi0 <- tbl_demo$Hl / tbl_demo$pi0
tbl_demo$piN_piS <- tbl_demo$Hl / tbl_demo$Hr
  
m_tbl <- pivot_longer(tbl_demo, cols=c(pi0, Hr, Hl, B, piN_pi0, piN_piS),
                      names_to="statistic")
fwrite(m_tbl, "tbl_demo.csv")