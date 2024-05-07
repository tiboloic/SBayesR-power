# MC integration tests

# integral dnorm, 0 +inf
# rewrite y = exp(-x)

us <- runif(N)
mean(1/us*dnorm(-log(us)))

# double integral
as <- runif(N)
bs <- runif(N)
mean(dchisq(-log(as), 1, -log(bs) * n / (1 - h2)) / as * dnorm(sqrt(-log(bs)), 0, sqrt(h2 / m / pi)) / sqrt(-log(bs)) / bs)



# simpler try
f <- function(x) {
  dnorm(x[1], 0, 1) + dnorm(x[2], 0, 1)
}

cubintegrate(f, lower=c(0,0), upper=c(Inf, Inf), method = "pcubature")
cubintegrate(dnorm, lower=0, upper=Inf, method = "pcubature")

g <- function(x) {
  dchisq(x[1], 1, x[2] * n / (1 - h2)) * dnorm(sqrt(x[2]), 0, sqrt(h2 / m / pi)) / sqrt(x[2])
}
cubintegrate(g, lower=c(0,0), upper=c(1, 1), method = "pcubature")
g_2 <-  function(x) {
  dchisq((log(x[1]) + log(1-x[1]) -log(A))/B, 1, x[2] *  n / (1 - h2)) * 
    dnorm(sqrt(x[2]), 0, sqrt(h2 / m / pi)) / sqrt(x[2]) / B / x[1] / (1-x[1])
}
cubintegrate(g_2, lower=c(0,0), upper=c(1, +Inf), method = "hcubature")
hcubature(g_2, c(0,0), c(1, +Inf), tol = 1e-10, absError= 1e-10)

# analytical example
# a compound random variable:exponential distributed with a rate from a [0,1] uniform

# analytical solution
exp_unif_anal <- function(x) {
  1 + exp(-x)/x -1/x
}

# MC integration by compound drawing
mean(rexp(10000, runif(10000))<2)

# MC integration by drawing from mixing distribution

# MC integration by double integral
as= runif(1000)
bs = runif(1000)
mean(as * exp(-as * bs))

# new try on dnorm(sqrt(x))/sqrt(x)
# 

h <- function(x, sig2) {
  dnorm(sqrt(x), 0, sqrt(sig2))/sqrt(x)
}
cubintegrate(h, lower=0, upper=Inf, method = "hcubature", sig2 = 0.007)

h <- function(y) exp(log_f_reparam(30, y = y, n = n, h2 = h2, m = m, pi = pi))
cubintegrate(h, lower=-30, upper=30, method = "hcubature")

integrate(
  function(y) exp(log_f_reparam(30, y = y, n = n, h2 = h2, m = m, pi = pi)),
  -15, 15, abs.tol = 0)$value

# reparam double integration with cubature
j <- function(x, o, h2, q, pi) exp(log_f_reparam(x[1], y = x[2], n = o, h2 = h2, m = q, pi = pi))
cubintegrate(j, lower=c(0,-30), upper=c(Inf,30), method = "hcubature",
             o = n, h2 = h2, q = m, pi = pi)

pow_cub <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  cubintegrate(j, lower = c(P_2_z(P0, n, h2, m, pi),-30), upper = c(Inf,30),
               method = "hcubature", absTol = 1e-12, maxEval = 1e7,#.Machine$double.eps,
               o = n, h2 = h2, q = m, pi = pi)
}

k <- function(x, o, h2, q, pi) {
  P <- x[1]
  v <- x[2]
  
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  uP <- (log(P) - log((1-P) * A)) / B
  ncp <- n * v / (1 - h2)
  exp(ldchi(uP, ncp) + lnor(v, sqrt(h2 / m / pi)) -log(P) -log(1-P) -log(B))
}

pow_cub_2 <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  cubintegrate(k, lower = c(P0, 0), upper = c(1, +Inf),
               method = "hcubature", o = n, h2 = h2, q = m, pi = pi)
}

# compare pow and pow_cub
pow(0.9, 30000, 0.5, 1e5, 0.001)
pow_cub(0.9, n=30000, h2 = 0.5, m = 1e5, pi=0.001)
# same. pow_cub is way faster

k <-  function(x, n=30000, h2 = 0.5, m = 1e5, pi=0.001) {
  P <- x[1]
  y <- x[2]
  
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  uP <- log(P / (1-P) / A) / B
  ncp <- - n * log(y) / (1 - h2)
  exp(ldchi(uP, ncp) + lnor(-log(y), sqrt(h2 / m / pi)) -log(B * P * (1-P) * y))
}
  
}