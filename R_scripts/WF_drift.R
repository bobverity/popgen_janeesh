# WF_drift.R
#
# Author: Bob Verity
# Date: 2023-03-06
#
# Purpose:
# Simulate drift under the Wright-Fisher model with mutation. Plot the allele
# frequency over time.
#
# ------------------------------------------------------------------

# define parameters
N <- 100
mu <- 0.001
t_max <- 1e3
p_init <- 0.5

# vector for storing allele frequencies in each generation
p <- rep(NA, t_max)
p[1] <- p_init

# simulate allele frequencies in all generations
for (i in 2:t_max) {
  p_mut <- p[i - 1]*(1 - mu) + (1 - p[i - 1])*mu
  p[i] <- rbinom(1, 2*N, p_mut) / (2*N) 
}

plot(p, type = 'l', ylim = c(0, 1))
