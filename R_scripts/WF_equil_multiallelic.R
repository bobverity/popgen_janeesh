# MULTIALLELIC

# code to model how the AF change over time in a multiallelic case with 3
# alleles without mutations occurring and with mutations 


#-------------------------------------------------------------------------------

## NO MUTATION

# pop size 
N <- 1e2
# length of simulation 
t_max <- 1e3
# initial AF of A 
p_0 <- c(0.4, 0.3, 0.3)
k <- length(p_0)

# store all the AF for each time point 
p <- matrix(NA, nrow = t_max, ncol = k)
p[1,] <- p_0

for (i in 2:t_max) {
  # prob - AF A, AF not A 
  n <- as.vector(rmultinom(n = 1, size = 2*N, prob = p[i-1,]))
  p[i,] <- n / (2*N)
}

plot(0, type = "n", xlim = c(0, t_max), ylim = c(0, 1))
for (i in 1:k) {
  lines(p[,i], col = i)
}


#-------------------------------------------------------------------------------

## MUTATION

# define parameters
N <- 1e2 # pop size 
t_max <- 1e3  # this should be at least 5/mu to ensure approximately at equilibrium
p_0 <- c(1, 0, 0)
k <- length(p_0)
mu <- 0.001
z <- c(1/3, 1/3, 1/3)

# store all the AF for each time point 
p <- matrix(NA, nrow = t_max, ncol = k)
p[1,] <- p_0

# simulate drift and mutation
for (i in 2:t_max) {
  
  # update allele frequencies
  p_new <- p[i-1,]*(1 - mu) + mu*z
  n <- rmultinom(n = 1, size = 2*N, prob = p_new)
  p[i,] <- n / (2*N)
}

plot(0, type = "n", xlim = c(0, t_max), ylim = c(0, 1))
for (i in 1:k) {
  lines(p[,i], col = i)
}
