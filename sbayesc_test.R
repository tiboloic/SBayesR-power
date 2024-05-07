# LT 15/04/2024
# 
# Test of numerical integration for SBayesC power
# compare with Monte Carlo simulations

# first scenario, looking at all SNPs (alpha_v = 0)

# parameters
n <- 30000
h2 <- 0.5
pi <- 0.01
m <- 1e6

# number of samples for MC integration
N = 1e3

# draw beta
betas <- rnorm(N, 0, sqrt(h2 / m / pi))
vs <- betas ^ 2

# plot distribution of v
hist(vs, prob = TRUE, breaks = 100)
plot(function(x) dnorm(sqrt(x), 0, sqrt(h2 / m / pi)) / sqrt(x),
     xlim = c(min(vs), max(vs)), add = TRUE)

# test naive MC integration drawing from chi2
P_v <- sapply(vs, function(v) mean(rchisq(N, 1, n * v / (1 - h2)) > 0.9))
mean(P_v)

# test inner integral
fP_mc <- function(P, n, h2, m, pi) {
  # draw v
  betas <- rnorm(N, 0, sqrt(h2 / m / pi))
  vs <- betas ^ 2
  z <- P_2_z(P, n = n, h2 = h2, m = m, pi = pi)
  mean(dchisq(z, 1, n * vs / (1 - h2)))
}

fP_int <- function(P, n, h2, m, pi){
  z <- P_2_z(P, n = n, h2 = h2, m = m, pi = pi)
  fP(0, z, n, h2, m, pi) 
}

fP_mc(0.3, 1e5, 0.5, m = m, pi = pi)
fP_int(0.3, 1e5, 0.5, m = m, pi = pi)  
# do match

# Monte Carlo integration of power
pow_mc <- function(P, n, h2, m, pi, N = 1000) {
  # draw P
  Ps = runif(N, min = P, max = 1)
  # for each P calculate f(P)
  (1 - P) * mean(sapply(Ps, function(P) fP_mc(P, n, h2, m, pi)))
}

pow_mc(0.1, n,h2, m, pi)
pow(0.1, 0, n = n, h2 = h2, m = m, pi = pi)              

# do not match


