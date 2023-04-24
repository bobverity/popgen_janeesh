# theta_posterior.R
#
# Author: Bob Verity
# Date: 2023-04-24
#
# Purpose:
# Simulate data from a beta-binomial distribution (assuming equilibrium under a
# WF model). Infer theta (scaled mutation rate) using Bayes' theorem.
#
# ------------------------------------------------------------------

#set.seed(6)

################################################################
# SIMULATE DATA

# draw from beta-binomial distribution
rbbinom <- function(n, theta, p_A) {
  p <- rbeta(1, shape1 = theta*p_A, shape2 = theta*(1 - p_A))
  n_A <- rbinom(1, size = n, prob = p)
  n_B <- n - n_A
  c(n_A, n_B)
}

# simulate data from a known theta
theta_true <- 1
n <- 1000
p_A_true <- 0.5

sim_data <- rbbinom(n = n, theta = theta_true, p_A = p_A_true)
sim_data

################################################################
# INFER THETA

# log-likelihood of observed data
ll <- function(n_A, n_B, theta, p_A, p_B = 1 - p_A) {
  lgamma(theta) + lgamma(n_A + theta*p_A) + lgamma(n_B + theta*p_B) - (lgamma(n_A + n_B + theta) + lgamma(theta*p_A) + lgamma(theta*p_B))
}

# plot the likelihood distribution of theta
theta_vec <- seq(0, 20, l = 1001)
ll_x <- ll(n_A = sim_data[1], n_B = sim_data[2], theta = theta_vec, p_A = p_A_true)
f_x <- exp(ll_x)

# define a prior and normalise
prior_theta <- dlnorm(theta_vec, meanlog = 2, sdlog = 2)
prior_theta <- dexp(theta_vec, rate = 1)
prior_theta <- prior_theta / sum(prior_theta)

# calculate posterior and normalise
post_theta <- f_x * prior_theta
post_theta <- post_theta / sum(post_theta, na.rm = TRUE)

# plot prior, likelihood and posterior. Vertical line indicates true value used
# to generate data
plot(theta_vec, prior_theta, type = 'l', ylim = c(0, 0.02), lty = 2)
lines(theta_vec, f_x / sum(f_x, na.rm = TRUE), col = 3)
lines(theta_vec, post_theta)
abline(v = theta_true, col = 2)
