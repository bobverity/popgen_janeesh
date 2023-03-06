# WF_equilibrium.R
#
# Author: Bob Verity
# Date: 2023-03-06
#
# Purpose:
# Simulate drift under the Wright-Fisher model with mutation. Explore the
# equilibrium allele frequency distribution.
#
# ------------------------------------------------------------------

# define parameters
N <- 1e2
mu <- 0.001
t_max <- 1e4  # this should be at least 5/mu to ensure approximately at equilibrium
p_init <- 0.5
reps <- 1e4

# vector for storing final allele frequency over reps
p <- rep(p_init, reps)

# simulate drift and mutation
for (i in 2:t_max) {
  
  # write progress to console
  if ((i %% 1e3) == 0) {
    message(i)
  }
  
  # update allele frequencies
  p_mut <- p*(1 - mu) + (1 - p)*mu
  p <- rbinom(reps, 2*N, p_mut) / (2*N)
}

# plot simulated equilibrium allele frequency distribution
hist(p, probability = TRUE, breaks = seq(0, 1, 0.01))

# overlay analytical solution
x <- seq(0, 1, 0.01)
theta <- 4*N*mu
fx <- dbeta(x, shape1 = theta, shape2 = theta)

lines(x, fx, col = 2, lwd = 2)



