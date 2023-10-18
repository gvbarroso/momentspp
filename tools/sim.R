
library(tidyverse)
library(cowplot)

# operators in c-space

## deterministic operators (working on expectations)
## neglect genealogical, mutational etc variance and operate on rel. freqs. x
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
 
  y[1] <- y[1] - s * y[1] * (y[2] + y[4])
  y[2] <- y[2] + s * (1 - y[2] - y[4])
  y[3] <- y[3] - s * y[1] * (y[2] + y[4])
  y[4] <- y[4] + s * (1 - y[2] - y[4])
  
  return(x)
}

## random operators (use prob. dist.)
## operate on abs. freqs x because rbinom() & rmultinom() require integers
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
  
  ab_to_Ab <- 0.5 * rbinom(1, n_ab_rec, p)
  ab_to_aB <- 0.5 * rbinom(1, n_ab_rec, q)
  # integers
  ab_to_Ab <- round(jitter(ab_to_Ab, amount = 1e-6))
  ab_to_aB <- round(jitter(ab_to_aB, amount = 1e-6))
  
  ab_stay_ab <- n_ab - ab_to_Ab - ab_to_aB

  Ab_to_ab <- 0.5 * rbinom(1, n_Ab_rec, 1-p)
  Ab_to_AB <- 0.5 * rbinom(1, n_Ab_rec, q)
  # integers
  Ab_to_ab <- round(jitter(Ab_to_ab, amount = 1e-6))
  Ab_to_AB <- round(jitter(Ab_to_AB, amount = 1e-6))
  
  Ab_stay_Ab <- n_Ab - Ab_to_ab - Ab_to_AB

  aB_to_AB <- 0.5 * rbinom(1, n_aB_rec, p)
  aB_to_ab <- 0.5 * rbinom(1, n_aB_rec, 1-q)
  # integers
  aB_to_AB <- round(jitter(aB_to_AB, amount = 1e-6))
  aB_to_ab <- round(jitter(aB_to_ab, amount = 1e-6))
  
  aB_stay_aB <- n_aB - aB_to_AB - aB_to_ab

  AB_to_aB <- 0.5 * rbinom(1, n_AB_rec, 1-p)
  AB_to_Ab <- 0.5 * rbinom(1, n_AB_rec, 1-q)
  # integers
  AB_to_aB <- round(jitter(AB_to_aB, amount = 1e-6))
  AB_to_Ab <- round(jitter(AB_to_Ab, amount = 1e-6))
  
  AB_stay_AB <- n_AB - AB_to_aB - AB_to_Ab

  x[1] <- ab_stay_ab + Ab_to_ab + aB_to_ab
  x[2] <- Ab_stay_Ab + ab_to_Ab + AB_to_Ab
  x[3] <- aB_stay_aB + ab_to_aB + AB_to_aB
  x[4] <- AB_stay_AB + Ab_to_AB + aB_to_AB
  
  return(x)
}

select_random <- function(x, s) { 
  
  n_ab <- x[1]
  n_Ab <- x[2]
  n_aB <- x[3]
  n_AB <- x[4]
  
  Ab_dead <- rbinom(1, n_Ab, abs(s))
  Ab_replacement <- rmultinom(1, Ab_dead, prob=x)
    
  AB_dead <- rbinom(1, n_AB, abs(s))
  AB_replacement <- rmultinom(1, AB_dead, prob=x)
  
  x[1] <- n_ab + Ab_replacement[1] + AB_replacement[1]
  x[2] <- n_Ab - Ab_dead + Ab_replacement[2] + AB_replacement[2]
  x[3] <- n_aB + Ab_replacement[3] + AB_replacement[3]
  x[4] <- n_AB - AB_dead + AB_replacement[4] + AB_replacement[4]
  
  return(x)
}

drift_rel <- function(x, N) { 
  return (rmultinom(n=1, size=N, prob=x))
}

drift_abs <- function(x) {
  N <- sum(x)
  return (rmultinom(n=1, size=N, prob=x))
}

## joint operators
evolve_det <- function(x, N, u, r, s) {
  return (select_det(recombine_det(mutate_det(drift_rel(x, N), u), r), s))
}

evolve_random <- function(x, u, r, s) {
  return (select_random(recombine_random(mutate_random(drift_abs(x), u), r), s))
}

evolve_neutral <- function(x, u, r, s) {
  return (recombine_random(mutate_random(drift_abs(x), u), r))
}

burn_in <- function(x, u, r, s, tol) {
  
  N <- sum(x)
  prev <- x # init
  
  repeat {
    new <- evolve_neutral(prev, N, u, r)
    diff <- sum((new / N - prev / N) ^ 2)
    
    cat(diff)
    cat("\n")
    
    if(diff < tol) { 
      break
    }
    
    else {
      prev <- new
    }
  }
  
  return(new)
}

# parameters

## vector of relative frequencies of haplotypes c11, c12, c13, c14
x1 <- c(5000, 3000, 1500, 2500) # init

## vector of relative frequencies of haplotypes c21, c22, c23, c24
x2 <- c(5000, 3000, 1500, 500) # init

## working with haplotype counts, therefore, no factor of 2 in pop. sizes
N1 <- sum(x1)
N2 <- sum(x2)

## 1-population
u <- 1e-6
r <- 1e-5
s <- -1e-4

## num gens spent in each epoch
epochs_gen <- c(1e+3, 5e+3, 4e+3)

## pop sizes per epoch
epochs_pop <- c(1e+4, 2e+4, 1e+5)

v <- burn_in(x1, u, r, s, tol=1e-7)
v

for(i in 1:1e+5) {
  x1 <- recombine_random(x1, 1e-4)
  cat(x1)
  cat("\n")
}

for(i in 1:1e+5) {
  x1 <- drift_abs(x1)
  cat(x1)
  cat("\n")
}

for(i in 1:1e+5) {
  x1 <- mutate_random(drift_abs(x), u)
  cat(x1)
  cat("\n")
}
