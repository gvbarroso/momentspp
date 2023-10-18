
library(tidyverse)
library(cowplot)

# operators in c-space

## deterministic (working on expectations)
### neglects genealogical, mutational, recombinational variance
mutate_det <- function(x, u) {
  
  x[1] <- x[1] - 2 * u * x[1]
  x[2] <- x[2] - u * x[2] + u * x[1]
  x[3] <- x[3] - u * x[3] + u * x[1]
  x[4] <- x[4] + u * x[2] + u * x[3]
  
  return(x)
}

recombine_det <- function(x, r) {
  
  x[1] <- x[1] - r * (x[1] * x[4] - x[2] * x[3])
  x[2] <- x[2] + r * (x[1] * x[4] - x[2] * x[3])
  x[3] <- x[3] + r * (x[1] * x[4] - x[2] * x[3])
  x[4] <- x[4] - r * (x[1] * x[4] - x[2] * x[3])
  
  return(x)
}

select_det <- function(x, s) { 
 
  y <- x / sum(x)
  
  y[1] <- y[1] - s * y[1] * (y[2] + y[4])
  y[2] <- y[2] + s * (1 - y[2] - y[4])
  y[3] <- y[3] - s * y[1] * (y[2] + y[4])
  y[4] <- y[4] + s * (1 - y[2] - y[4])
  
  return(x)
}

## random
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
  
  ab_to_Ab <- round(0.5 * rbinom(1, n_ab_rec, p), 0)
  ab_to_aB <- round(0.5 * rbinom(1, n_ab_rec, q), 0)
  ab_stay_ab <- n_ab - ab_to_Ab - ab_to_aB

  Ab_to_ab <- round(0.5 * rbinom(1, n_Ab_rec, 1-p), 0)
  Ab_to_AB <- round(0.5 * rbinom(1, n_ab_rec, q), 0)
  Ab_stay_Ab <- n_Ab - Ab_to_ab - Ab_to_AB

  aB_to_AB <- round(0.5 * rbinom(1, n_aB_rec, p), 0)
  aB_to_ab <- round(0.5 * rbinom(1, n_ab_rec, 1-q), 0)
  aB_stay_aB <- n_aB - aB_to_AB - aB_to_ab

  AB_to_aB <- round(0.5 * rbinom(1, n_aB_rec, 1-p), 0)
  AB_to_Ab <- round(0.5 * rbinom(1, n_ab_rec, 1-q), 0)
  AB_stay_AB <- n_AB - AB_to_aB - AB_to_Ab

  x[1] <- ab_stay_ab + Ab_to_ab + aB_to_ab
  x[2] <- Ab_stay_Ab + ab_to_Ab + AB_to_Ab
  x[3] <- aB_stay_aB + ab_to_aB + AB_to_aB
  x[4] <- AB_stay_AB + Ab_to_AB + aB_to_AB
  
  return(x)
}

drift <- function(x, N) { 
  return(rmultinom(n=1, size=N, prob=x))
}

## joint operators
evolve_det <- function(x, N, u, r, s) {
  return(select_det(recombine_det(mutate_det(drift(x, N), u), r), s))
}

evolve_random <- function(x, N, u, r, s) {
  return(select_det(recombine_random(mutate_random(drift(x, N), u), r), s))
}

burn_in <- function(x, N, u, r, s, tol) {
  
  prev <- x # init
  
  repeat {
    new <- evolve_random(prev, N, u, r, s)
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

# parameters

## we work in absolute counts vectors to ensure they remain integers throughout evolution

## vector of relative frequencies of haplotypes c11, c12, c13, c14
x1 <- c(5000, 3000, 1500, 500) # init

## vector of relative frequencies of haplotypes c21, c22, c23, c24
x2 <- c(5000, 3000, 1500, 500) # init

## working with haplotype counts, therefore, no factor of 2 in pop. sizes
N1 <- sum(x1)
N2 <- sum(x2)

## 1-population
u <- 1e-7
r <- 1e-6
s <- -1e-4

## num gens spent in each epoch
epochs_gen <- c(1e+3, 5e+3, 4e+3)

## pop sizes per epoch
epochs_pop <- c(1e+4, 2e+4, 1e+5)

v <- burn_in(x1, N, u, r, s, 1e-2)
