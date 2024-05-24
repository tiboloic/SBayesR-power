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

z_2_P_R <- function(z, n, h2 = 0.5, gammas, pis) {
  lambda <- (1 - h2) / h2 / gammas
  C <- n + lambda
  A <- pis / pis[1] * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  1 - 1 / (1 + sum(A[-1] * exp(B[-1] * z)))
  #cat(lambda, '\n')
  #cat (A, '\n')
  #cat(B, '\n')
  #cat(C, '\n')
  }

P_2_z_R <- function(P, n, h2 = 0.5, gammas, pis) {
  Pmin <- z_2_P_R(0, n, h2, gammas, pis)
  P <- max(P, Pmin)
  uniroot(function (x, ...) z_2_P_R(x, ...) - P,
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

ldchi <- function(x, ncp) {
  if (length(ncp) > 1) {
    res <- numeric(length(ncp))
    if (any(ncp == 0)) {
      if (length(x) > 1) sub_x <- x[ncp == 0]
      else sub_x <- x
      res[ncp == 0] <- dchisq(sub_x, 1, log = TRUE)
    }
    if (any(ncp > 0)) {
      if(length(x) > 1) x <- x[ncp > 0]
      res[ncp > 0] <- log(0.5) -0.5 * (sqrt(ncp[ncp>0]) - sqrt(x))^2 - 0.25 *
        log(x/ncp[ncp>0]) + 
        log(besselI(sqrt(x * ncp[ncp>0]), -0.5, expon.scaled = TRUE))
    }
    res
  } else {
    if(ncp == 0)
      dchisq(x, 1, log = TRUE)
    else
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
  m <- sapply(2:length(gammas), function(i)
    pis[i] * exp(log_f_i(i, z, y, n, h2, gammas, pis))) / sum(pis[-1])
  
  if (is.null(dim(m)))
     sum(m)
  else
    rowSums(m)
}

f_2d <- function(x, n , h2, gammas, pis) {
  z <- x[1]
  y <- x[2]
  sum(sapply(2:5, function(i)
    pis[i] * exp(log_f_i(i, z, y, n, h2, gammas, pis))))  / sum(pis[-1])
}

pow_cub <- function(P0, n, h2, gammas, pis) {
  hcubature(f_2d, c(P_2_z_R(P0, n, h2, gammas, pis), -15), c(+Inf, 15),
            n = n, h2 = h2, gammas = gammas, pis = pis, tol = 1e-10)
}

fP_R <- function(zs, n, h2 = 0.5, gammas, pis) {
  sapply(zs, function(z) {
    integrate(function(y, ...) f(y = y, ...),
              -20, 20,
              z = z, n = n, h2 = h2, gammas = gammas, pis = pis,
              abs.tol = 0)$value  
  })
}

pow_R <- function(P0, n, h2 = 0.5, gammas, pis) {
  integrate(fP_R, P_2_z_R(P0, n, h2, gammas, pis), Inf,
            n = n, h2 = h2, gammas = gammas, pis = pis, abs.tol=0)$value
}

#. Monte Carlo integration
pow_R_mc <- function(P0, n, h2 = 0.5, gammas, pis, N = 10000) {
  
  if (length(gammas) == 2) {
    deltas <- 2
  } else {
    deltas <- sample(2:length(gammas), N, replace = TRUE, prob = pis[-1])
  }

  sds <- sqrt(gammas * h2)
  vs <- rnorm(N, 0, sds[deltas]) ^ 2
  ncp <- n * vs / (1 - h2)
  Ps <- sapply(rchisq(N, 1, ncp), z_2_P_R,
               n = n, h2 = h2, gammas = gammas, pis = pis)
  mean(Ps > P0)
}

# check on SBayesC : gamma = c(0, 1/m/pi), pis = c(1-pi, pi)
pow(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_R(0.2, n=30000, h2 = 0.5, gammas=c(0, 1/1e6/0.001), pis=c(1-0.001, 0.001))
pow_R_mc(0.2, n=30000, h2 = 0.5, gammas=c(0, 1/1e6/0.001), pis=c(1-0.001, 0.001))

# check on a 5 components mixture
pow_R(0.9, n=3000, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_mc(0.2, n=30000, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
