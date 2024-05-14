# LT 13/05/2024
# 
# 

Bs <- function(n, h2, gamma) {
  # eq. 9  
  lambda <- (1 - h2) / (gamma * h2)
  
  # eq. 10
  C <- n + lambda
  
  # eq. 21
  B <- 0.5 * n * (1 - h2) / C
  
  B
}

z_2_P <- function(z, n, h2 = 0.5, gammas, pis) {
  lambda <- (1 - h2) / h2 / gammas
  C <- n + lambda
  A <- pi / pi[1] * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  1 - 1 / (1 + sum(A[-1] * exp(B[-1] * z)))
}

P_2_z <- function(P, n, h2 = 0.5, gammas, pis) {
  Pmin <- z_2_P(0, n, h2, gammas, pis)
  P <- max(P, Pmin)
  uniroot(function (x, ...) z_2_P(x, ...) - P,
          lower = 0, upper = 1e6,
          n = n, h2 = h2, gammas = gammas, pis = pis,
          tol = .Machine$double.eps,
          maxiter = 100000)$root
}

# try with optimize. Not better
P_2_z_2 <- function(P, n, h2 = 0.5, gammas, pis) {
  Pmin <- z_2_P(0, n, h2, gammas, pis)
  P <- max(P, Pmin)
  optimize(function (x, ...) (z_2_P(-log(x), ...) - P)^2,
          lower = 0, upper = 1,
          n = n, h2 = h2, gammas = gammas, pis = pis)
}

# numerically robust non central chi-square with 1 degree of freedom
ldchi <- function(x, ncp) {
  if (ncp == 0) {
    dchisq(x, 1, log = TRUE)
  } else {
    log(0.5) -0.5 * (sqrt(ncp) - sqrt(x))^2 - 0.25 * log(x/ncp) +
      log(besselI(sqrt(x * ncp), -0.5, expon.scaled = TRUE))
  }
}

# log pdf of dnorm(sqrt(x), sigma) * sqrt(x)
ldnor <- function(x, sig) {
  0.5 * log(x) - log(sig) - 0.5 * log(2 * 3.14159265359) - 0.5 * x / sig ^ 2
}

# f(P|v) * f(v)
log_f_i <- function(i, z, y, n, h2 = 0.5, gammas, pis)
{
  ldchi(z, exp(y)) + 
    ldnor((1 - h2) / n * exp(y), sqrt(h2 * gammas[i]))
}

f <- function(z, y, n , h2, gammas, pis) {
  sum(sapply(2:5, function(i)
    pis[i] * exp(log_f_i(i, z, y, n, h2, gammas, pis))))
}

f_2d <- function(x, n , h2, gammas, pis) {
  z <- x[1]
  y <- x[2]
  sum(sapply(2:5, function(i)
    pis[i] * exp(log_f_i(i, z, y, n, h2, gammas, pis))))
}

pow_cub <- function(P0, n, h2, gammas, pis) {
  hcubature(f_2d, c(P_2_z(P0, n, h2, gammas, pis), -15), c(+Inf, 15),
            n = n, h2 = h2, gammas = gammas, pis = pis)
}

fP <- function(z, n, h2 = 0.5, gammas, pis) {
  sapply(z, function(z) {
    integrate(
,
      -20, 20, abs.tol = 0)$value
  })
}

pow <- function(P0, n, h2 = 0.5, gammas, pis) {
  integrate(fP, P_2_z(P0, n, h2, gammas, pis), Inf,
            n = n, h2 = h2, gammas = gammas, pis = pis, abs.tol=0)$value
}
