# LT 13/05/2024
# 
# 
library(cubature)

z_2_P_R <- function(z, n, h2 = 0.5, gammas, pis) {
  lambda <- (1 - h2) / h2 / gammas
  C <- n + lambda
  A <- pis / pis[1] * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  1 - 1 / (1 + sum(A[-1] * exp(B[-1] * z)))
}

get_max_z <- function(P, n, h2 = 0.5, gammas, pis) {
  lambda <- (1 - h2) / h2 / gammas
  C <- n + lambda
  A <- pis / pis[1] * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  # based on avoiding averflow exp(z_max * B) = +Inf
  500 / max(B[-1])
}

P_2_z_R <- function(P, n, h2 = 0.5, gammas, pis) {
  Pmin <- z_2_P_R(0, n, h2, gammas, pis)
  P <- max(P, Pmin)
  uniroot(function (x, ...) z_2_P_R(x, ...) - P,
          lower = 0, upper = get_max_z(P, n, h2, gammas, pis),
          n = n, h2 = h2, gammas = gammas, pis = pis,
          tol = .Machine$double.eps,
          maxiter = 100000)$root
}

# try with optimize. Not better
P_2_z_2 <- function(P, n, h2 = 0.5, gammas, pis) {
  Pmin <- z_2_P(0, n, h2, gammas, pis)
  P <- max(P, Pmin)
  optimize(function (x, ...) (z_2_P_R(x, ...) - P)^2,
          n = n, h2 = h2, gammas = gammas, pis = pis,
          lower = 0, upper = 1e3, tol = .Machine$double.eps)
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

# vectorized for efficiency
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
  
  # special case when only 2 mixtures (SBayesC)
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
pow_R_mc <- function(P0, n, h2 = 0.5, gammas, pis, N = 1e6) {
  
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
  power <- mean(Ps > P0)
  cat("power = ", power, " with standard error se = ", sd(Ps > P0)/sqrt(N), "\n")
  return(power)
}

# new implementation based on the cdf of the non central chi square

f_R_e <- function(v, z0, n, h2 = 0.5, gammas, pis) {
  ncp <- n * v / (1 - h2)
  chi_int <-  pchisq(z0, 1, ncp, lower.tail = FALSE)
  m <- sapply(2:length(gammas), function(i)
    pis[i] / sum(pis[-1]) * chi_int *
      dnorm(sqrt(v), 0, sqrt(gammas[i] * h2)) / sqrt(v))
  
  # special case when only 2 mixtures (SBayesC)
  if (is.null(dim(m)))
    sum(m)
  else
    rowSums(m)
}

pow_R_e <-  function(P0, n, h2 = 0.5, gammas, pis) {
  z0 <- P_2_z_R(P0, n, h2, gammas, pis)
  sds <- sqrt(gammas * h2)
  v_max <- 10 * max(sds)
  
  #integrate(f_R_e, lower = 0, upper = v_max,
  #          z0 = z0, n = n, h2 = h2, gammas = gammas, pis = pis, abs.tol=0)$value
  hcubature(f_R_e, 0, v_max,
            z0 = z0, n = n, h2 = h2, gammas = gammas, pis = pis, absError = 0)$integral
}

f_prop_var_R <- function(v, z0, n, h2 = 0.5, gammas, pis) {
  ncp <- n * v / (1 - h2)
  chi_int <-  pchisq(z0, 1, ncp, lower.tail = FALSE)
  m <- sapply(2:length(gammas), function(i)
    pis[i] / sum(pis[-1]) * chi_int *
      dnorm(sqrt(v), 0, sqrt(gammas[i] * h2)) * sqrt(v))
  
  # special case when only 2 mixtures (SBayesC)
  if (is.null(dim(m)))
    sum(m)
  else
    rowSums(m)
}

# P(v) to calculate expectation of v
f_exp_v <- function(v, h2 = 0.5, gammas, pis) {
  m <- sapply(2:length(gammas), function(i)
    pis[i] / sum(pis[-1]) * dnorm(sqrt(v), 0, sqrt(gammas[i] * h2)) * sqrt(v))
  
  # special case when only 2 mixtures (SBayesC)
  if (is.null(dim(m)))
    sum(m)
  else
    rowSums(m)
}

exp_v <- function(h2 = 0.5, gammas, pis) {
  sds <- sqrt(gammas * h2)
  v_max <- 10 * max(sds)
  
  integrate(f_exp_v, lower = 0, upper = v_max,
            h2 = h2, gammas = gammas, pis = pis, abs.tol=0)$value
  
}

prop_var_R <- function(P0, n, h2 = 0.5, gammas, pis) {
  z0 <- P_2_z_R(P0, n, h2, gammas, pis)
  sds <- sqrt(gammas * h2)
  v_max <- 10 * max(sds)
  
  integrate(f_prop_var_R, lower = 0, upper = v_max,
            z0 = z0, n = n, h2 = h2, gammas = gammas, pis = pis, abs.tol=0)$value /
    exp_v(h2, gammas, pis)
}

exp_v_mc <- function(h2 = 0.5, gammas, pis, N = 100000) {
  
  if (length(gammas) == 2) {
    deltas <- 2
  } else {
    deltas <- sample(2:length(gammas), N, replace = TRUE, prob = pis[-1])
  }
  
  sds <- sqrt(gammas * h2)
  vs <- rnorm(N, 0, sds[deltas]) ^ 2
  Ev <- mean(vs)
  cat("expectation of v = ", Ev, " with standard error se = ", sd(vs)/sqrt(N), "\n")
  return(Ev)
}

calc_power <- pow_R_e
calc_power_mc <- pow_R_mc 

# check on SBayesC : gamma = c(0, 1/m/pi), pis = c(1-pi, pi)
#pow(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_R(0.2, n=30000, h2 = 0.5, gammas=c(0, 1/1e6/0.001), pis=c(1-0.001, 0.001))
pow_R_mc(0.2, n=30000, h2 = 0.5, gammas=c(0, 1/1e6/0.001), pis=c(1-0.001, 0.001))
pow_R_e(0.2, n=30000, h2 = 0.5, gammas=c(0, 1/1e6/0.001), pis=c(1-0.001, 0.001))

# check on a 5 components mixture
pow_R(0.9, n=30000, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_mc(0.9, n=30000, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_e(0.9, n=30000, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))

# try for large n
pow_R(0.9, n=1e6, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_mc(0.9, n=1e6, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_e(0.9, n=1e6, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))

# very large n, small h2
pow_R_mc(0.9, n=1e7, h2 = 0.05, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
pow_R_e(0.9, n=1e7, h2 = 0.05, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))

# very large n, small h2
pow_R_mc(0.9, n=1e7, h2 = 0.05, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.9, 0.04, 0.03, 0.02, 0.01))
pow_R_e(0.9, n=1e7, h2 = 0.05, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.9, 0.04, 0.03, 0.02, 0.01))


# test expectation of v
exp_v(h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
exp_v_mc(h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))

exp_v(h2 = 0.3, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.9, 0.04, 0.03, 0.02, 0.01))
exp_v_mc(h2 = 0.3, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.9, 0.04, 0.03, 0.02, 0.01))
exp_v_mc_2(h2 = 0.3, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.9, 0.04, 0.03, 0.02, 0.01))

prop_var_R(0.9, 1e6, h2 = 0.5, gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2), pis =c (0.99, 0.004, 0.003, 0.002, 0.001))
