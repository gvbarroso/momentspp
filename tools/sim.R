
library(tidyverse)
library(cowplot)

# parameters
u <- 1e-8
r <- 1e-7
s <- -1e-4

# vector of relative frequencies of haplotypes c11, c12, c13, c14
x1 <- c(5000, 3000, 1500, 500) # init

# vector of relative frequencies of haplotypes c21, c22, c23, c24
x2 <- c(5000, 3000, 1500, 500) # init

# we work in absolute counts vectors to ensure they remain integers throughout evolution

N1 <- sum(x1)
N2 <- sum(x2)

# operators in c-space
mutate_det <- function(x, u) {
  x[1] <- x[1] - 2 * u * x[1]
  x[2] <- x[2] - u * x[2] + u * x[1]
  x[3] <- x[3] - u * x[3] + u * x[1]
  x[4] <- x[4] + u * x[2] + u * x[3]
  
  return(x)
}

mutate_random <- function(x, u) {
  n_ab <- x[1]
  n_Ab <- x[2]
  n_aB <- x[3]
  
  n_ab_mut <- rbinom(1, n_ab, 2 * u)
  n_Ab_mut <- rbinom(1, n_Ab, u)
  n_aB_mut <- rbinom(1, n_aB, u)
  
  to_Ab <- rbinom(1, n_ab_mut, 0.5)
  to_aB <- n_ab_mut - to_Ab
  
  x[1] <- x[1] - n_ab_mut
  x[2] <- x[2] - n_Ab_mut  + to_Ab
  x[3] <- x[3] - n_aB_mut + to_aB
  x[4] <- x[4] + n_Ab_mut + n_aB_mut
  
  return(x)
}

recombine_det <- function(x, r) {
  x[1] <- x[1] - r * (x[1] * x[4] - x[2] * x[3])
  x[2] <- x[2] + r * (x[1] * x[4] - x[2] * x[3])
  x[3] <- x[3] + r * (x[1] * x[4] - x[2] * x[3])
  x[4] <- x[4] - r * (x[1] * x[4] - x[2] * x[3])
  
  return(x)
}

recombine_random <- function(x, r) {
  
  n_ab <- x[1]
  n_Ab <- x[2]
  n_aB <- x[3]
  n_AB <- x[3]
  
  N <- sum(x)
  p <- (n_Ab + n_AB) / N
  q <- (n_aB + n_AB) / N
  
  n_ab_rec <- rbinom(1, n_ab, r)
  n_Ab_rec <- rbinom(1, n_Ab, r)
  n_aB_rec <- rbinom(1, n_aB, r)
  n_AB_rec <- rbinom(1, n_AB, r)
  
  ab_to_Ab <- 0.5 * rbinom(1, n_ab_rec, p)
  ab_to_aB <- 0.5 * rbinom(1, n_ab_rec, q)
  ab_stay_ab <- n_ab_rec - ab_to_Ab - ab_to_aB
  
  x[1] <- x[1] - n_ab_rec
  x[2] <- x[2] - n_Ab_mut  + to_Ab
  x[3] <- x[3] - n_aB_mut + to_aB
  x[4] <- x[4] + n_Ab_mut + n_aB_mut
  
  return(x)
}

select_det <- function(x, s) { 
  x[1] <- x[1] - s * x[1] * (x[2] + x[4])
  x[2] <- x[2] + s * (1 - x[2] - x[4])
  x[3] <- x[3] - s * x[1] * (x[2] + x[4])
  x[4] <- x[4] + s * (1 - x[2] - x[4])
  
  return(x)
}

drift <- function(x, N) { 
  return(rmultinom(n=1, size=2*N, prob=x) / (2 * N))
}

evolve <- function(x, N, u, r, s) {
  return(select_det(recombine_det(mutate_random(drift(x, N), u, N), r), s))
}

# 1-population
x1 <- c(0.6, 0.25, 0.1, 0.05)

# num gens spent in each epoch
epochs_gen <- c(1e+3, 5e+3, 4e+3)

# pop sizes per epoch
epochs_pop <- c(1e+4, 2e+4, 1e+5)

burn_in <- function(x, N, u, r, s, tol) {
  
  prev <- x # init
  
  repeat {
    new <- evolve(prev, N, u, r, s)
    diff <- sum((new - prev) ^ 2)
    
    cat(diff)
    cat("\n")
    
    if(sum < tol) { 
      break
    }
    
    else {
      prev <- new
    }
  }
  
  return(new)
}
