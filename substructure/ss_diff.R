#################### substructure differentials solution ########################

# substituting complex parameter values - trying to simplfy everything 

# parameter values selected 
# within deme population size 
N <- 30000
# migration rate 
m <- 1e-4
m <- 1/(4*N)
# number of demes 
K <- 2

# alpha and beta from the quadratic equation 
alpha <- (4*N*m*K + K - 1) / (2*N*(K - 1))
beta <- sqrt(alpha^2 - 4*m/(N*(K - 1)))

# eigenvalues calculated 
l1 <- 0.5*(-alpha + beta)
l2 <- 0.5*(-alpha - beta)

# substituting complex matrix values 
a <- -2*m - 1/(2*N)
b <- 2*m/(K - 1)

# eta2 of the eigenvectors in terms of gamma and episilon 
gamma <- (l1 - a)/b
eps <- (l2 - a)/b

# constants of the solution  
c1 <- eps / (eps - gamma)
c2 <- gamma / (gamma - eps)

# ------------------------------------------------------------------------------

# solved equations and plotting 

# number of generations 
t <- seq(0, 1e5, l = 1001)
t <- seq(0, 1e6, l = 1001)

# solved equation for S(t) - prob uncoal and same deme 
St <- c1*exp(t*l1) + c2*exp(t*l2)
# solved equation for D(t) - prob of uncoal and different demes 
Dt <- c1*exp(t*l1)*gamma  + c2*exp(t*l2)*eps

plot(t, St, type = 'l', ylim = c(0, 1))
lines(t, Dt, col = 2)

# probability of being uncoal ever (same or diff demes)
Ut <- St + Dt
plot(t, Ut, type = 'l', ylim = c(0, 1))
lines(t, St, col = 2)
lines(t, Dt, col = 3)
# probability of coal event if pop is panmictic
ft <- exp(-t/(2*K*N))
lines(t, ft, col = 2, lty = 2)

# ------------------------------------------------------------------------------

# coal rates and Ne

# proportion of uncoal due to diff demes 
Zt <- Dt / (Dt + St)
plot(t, Zt, type = 'l', ylim = c(0, 1))

# proportion uncoal due to same demes x coal rate = effective coal rate 
coal_eff <- (St / (St + Dt))/(2*N)
plot(t, coal_eff, type = 'l')

# coal rate for pop size N
abline(h = 1/(2*N), col = 2)
# coal rate if was panmictic (N x K)
abline(h = 1/(2*N*K), col = 2)

# effective pop size 
Neff <- N*(1 + Dt / St)
plot(Neff, type = 'l', ylim = c(0, max(Neff)))
plot(Neff, type = 'l')
# population size within deme 
abline(h = N, col = 2)
# population size whole system - panmictic 
abline(h = N*K, col = 2)
