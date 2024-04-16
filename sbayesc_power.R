# LT 12/04/2024
# 
# SBayesC power calculations



Bs <- function(n, h2, gamma) {
  # eq. 9  
  lambda <- (1 - h2) / (gamma * h2)
  
  # eq. 10
  C <- n + lambda
  
  # eq. 21
  B <- 0.5 * n * (1 - h2) / C
  
  B
}

## functions to handle re parametrisation

P_2_z <- function(P, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  (log(P) - log(1-P) - log(A)) / B
}

z_2_P <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  1 - 1 / (1 + A * exp(B * z))
}

# numerically robust non central chi-square with 1 degree of freedom
ldchi <- function(x, lambda) {
  log(0.5) -0.5 * (sqrt(lambda) - sqrt(x))^2 - 0.25 * log(x/lambda) +
    log(besselI(sqrt(x * lambda), -0.5, expon.scaled = TRUE))
}

# log pdf of dnorm(sqrt(x), sigma) * sqrt(x)
ldnor <- function(x, sig) {
  0.5 * log(x) - log(sig) - 0.5 * log(2 * 3.14159265359) - 0.5 * x / sig ^ 2
}

# f(P|v) * f(v)
log_f <- function(z, y, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  ldchi(z, exp(y)) + 
    ldnor((1 - h2) / n * exp(y), sqrt(h2 / m / pi))
}

F_V <- function(x, sig, pi) {
  ifelse(x > 0,
         (1 - pi) + pi * (pnorm(sqrt(x), 0, sig) - pnorm(-sqrt(x), 0, sig)),
         0)
}

# generic function to get the quantile function of an arbitrary cdf
inverse = function (f, lower = 0, upper = 1) {
  function (y, ...) uniroot(function (x) f(x, ...) - y,
                       lower = lower, upper = upper,
                       maxiter = 100000)$root
}

# TODO write a safer alternative for q_V, dealing with the lower bound for v (1 - pi)
q_V <- function(p, sig, pi)
  ifelse (p < 1 - pi, 0, inverse(F_V)(p, sig = sig, pi = pi))



# f(P) (v is integrated out)
# re parametrised in z instead of P
# p = proportion of genetic variance threshold 0 <= p <= 1  
fP <- function(p, z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  
  # calculate alpha_v see eq. 48
  alpha_v <- q_V(p, sig = sqrt( (1 - h2) / m / pi), pi = pi)
  alpha_y <- log(n * alpha_v / (1 - h2))
  sapply(z, function(z) {
    integrate(
      function(y) exp(log_f(z, y = y, n = n, h2 = h2, m = m, pi = pi)),
     max(alpha_y, -20), 20, abs.tol = 0)$value
  })
}

# v  = 1 - % of variance explained:
# v = 0 => all SNPS, whatever the % variance they explain
# PIP: threshold PIP above which we consider variant is causal
pow <- function(PIP, v, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  integrate(fP, P_2_z(PIP, n=n, h2 = h2, m = m, pi = pi), Inf, 
            p = v, n = n, h2 = h2, m = m, pi = pi, 
            abs.tol = 0, subdivisions = 10000)$value
}

