maps$cum_r_scaled_bt <- maps$cum_r
maps$cum_r_scaled_exp <- maps$cum_r

## @@
focal_r <- abs(dt_neut[samp_id]$cum_r_scaled - dt_exons[exon_id]$cum_r_scaled) 

## @@
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
    
    Ne_bar_bt <- unique(filter(m_het_demo, Generation==g, u==mu, N1==1e+3)$pi0)/
      (2 * mu)
    
    Ne_bar_exp <- unique(filter(m_het_demo, Generation==g, u==mu, N1==1e+5)$pi0)/
      (2 * mu)
    
    Ne_ratio_bt <- Ne_bar_bt / 1e+4
    Ne_ratio_exp <- Ne_bar_exp / 1e+4
    
    dt_neut$cum_r_scaled <- 
      
      bt_strong <- filter(bt_tmp, s==-1e-3)
    bt_moderate <- filter(bt_tmp, s==-1e-4)
    bt_weak <- filter(bt_tmp, s==-1e-5)
    
    exp_strong <- filter(exp_tmp, s==-1e-3)
    exp_moderate <- filter(exp_tmp, s==-1e-4)
    exp_weak <- filter(exp_tmp, s==-1e-5)
    
    B_values <- rep(1, nrow(dt_neut)) # init 
    
    for(x in 1:num_iter) {
      #approximates cumulative rec by looking at r at sampled sites only
      focal_r <- unlist(lapply(samp_pos, function(pos) rmap[J(pos), roll=T]$r)) 
      focal_r <- focal_r * B_values
      cum_rec <- numeric(length=length(samp_pos))
      cum_rec[1] <- focal_r[1]
      
      for(j in 2:length(samp_pos)) {
        cum_rec[j] <- cum_rec[j-1] + (samp_pos[j] - samp_pos[j-1]) * focal_r[j]
      }
      pos_dt$cumrec <- cum_rec
      
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