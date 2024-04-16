# LT 15/04/2024
# 
# Test of numerical integration for SBayesC power
# compare with Monte Carlo simulations

# first scenario, looking at all SNPs (alpha_v = 0)

# parameters
n <- 10000
h2 <- 0.5
pi <- 0.01
m <- 1e6

# number of samples for MC integration
N = 1e3

# draw beta
betas <- rnorm(N, 0, h2 / m / pi)
vs <- betas ^ 2

# draw from non central chi-square
P_v <- sapply(vs, function(v) mean(rchisq(N, 1, n * v / (1 - h2)) > 0.9))
mean(P_v)

pow(0.9, 0, n = n, h2 = h2, pi =pi, m=m)              

# do not match
# test inner integral
