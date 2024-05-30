# LT 9/04/2024
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
  
u <- function(PIP, n, h2 = 0.5,
              gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
              pi = c(0.99, 0.004, 0.003, 0.002, 0.001)) {
  # eq. 9  
  lambda <- (1 - h2) / (gamma * h2)
  
  # eq. 10
  C <- n + lambda
  
  # eq. 20
  A <- pi / pi[1] * sqrt(lambda) / sqrt(C)
  
  # eq. 21
  B <- 0.5 * n * (1 - h2) / C
  
  # eq. 36
  z <- (4 * log(PIP) - 4 * log(1 - PIP) + sum(log(pi[-1]) - 4 * log(1 - pi[1]) -
                                                sum(log(A[-1])))) / sum(B[-1])
  return(z)
}

ncp <- function(v, h2 = 0.5, n) {
  # non-central parameter of the \Chi^2
  ncp <-  n * v / (1 - h2)
  ncp
}

# integrand 
f <- function(v, i, PIP = 0.9, n = 100, h2 = 0.5,
              gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
              pi = c(0.99, 0.004, 0.003, 0.002, 0.001))
{
  lam <- ncp(v, h2, n)
  cat("non central parameter: ", lam, "\n")
  
  uP <- u(PIP, n, h2, gamma, pi)
  cat("u(P): ", uP, "\n")
  
  0.5 * exp(-0.5 * (sqrt(lam) - sqrt(uP))^2) * (lam/uP) ^ 0.25 *
    besselI(sqrt(uP * lam), -0.5, expon.scaled = TRUE) / sqrt(v) *
    dnorm(sqrt(v), 0, sqrt(gamma[i] * h2))
}

# plot terms to see
# density is concentrate around 0+

# numerically robust non central chi-square with 1 degree of freedom
# # vectorized in x and ncp
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

lnor <- function(x, sig) {
  - 0.5 * log(x) - log(sig) - 0.5 * log(2 * 3.14159265359) - 0.5 * x / sig^2
}

# use log scale to avoid underflow and overflow
log_f <- function(v, i, PIP = 0.9, n = 1e8, h2 = 0.5,
                  gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                  pi = c(0.99, 0.004, 0.003, 0.002, 0.001))
{
  lam <- ncp(v, h2, n)
  
  uP <- u(PIP, n, h2, gamma, pi)
  
  ldchi(uP, lam) + lnor(v, sqrt(gamma[i] * h2)) 

}

# still numerically unstable
#test <- integrate(f, 0, Inf, i = 4)
test2 <- integrate( function(x) exp(log_f(x, i = 4)), 0, Inf)  
  
# try change of variable 
log_f_reparam <- function(y, i = 5, PIP = 0.9, n = 1e8, h2 = 0.5,
                          gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                          pi = c(0.99, 0.004, 0.003, 0.002, 0.001))
{
  uP <- u(PIP, n, h2, gamma, pi)
  
  ldchi(uP, exp(y)) + lnor(sqrt((1 - h2) / n) * exp(0.5 * y), sqrt(gamma[i] * h2)) +
    0.5 * log((1 - h2)) - 0.5 * log(n) + 0.5 * y
} 

# cannot integrate on rreal line
test2 <- integrate( function(x) exp(log_f_reparam(x)), -Inf, Inf) 

# but is now numerically stable on a reasonable interval
test3 <- integrate( function(x) exp(log_f_reparam(x)), -15, 15)

# f(P)
P_0 <- function(PIP = 0.9, n = 1e8, h2 = 0.5,
                gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                pi = c(0.99, 0.004, 0.003, 0.002, 0.001))
{
  sum(sapply(2:5, function(i) {
    int <- integrate(function(x) exp(log_f_reparam(x,
                                          i = i,
                                          PIP = PIP,
                                          n = n,
                                          h2 = h2,
                                          gamma = gamma,
                                          pi = pi)), -15, 15, abs.tol=0)
    if (int$message != "OK")
      stop(paste0("Inner integral failed with message: ", int$message))
    
    4 * pi[i] / sum(Bs(n, h2, gamma)) / PIP / (1 - PIP) * int$value
    #4 * pi[i] / sum(Bs(n, h2, gamma))  * int$value
  }))
}

# plot pdf of P
plot(function(x) sapply(x, function(p) P_0(p)), xlim=c(0.01,0.99))

# try to integrate P

test <- integrate(function(x) sapply(x, function(p) P_0(PIP=p)), 0.5, 1)


g <- function(PIP, v = 0.01, n = 1e8, h2 = 0.5,
              gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
              pi = c(0.9, 0.04, 0.03, 0.02, 0.01)) {
  sum(sapply(2:length(gamma), function(i) {
    4 * pi[i] / sum(Bs(n, h2, gamma)) / PIP / (1-PIP) *
      exp(log_f(v=v, i=i, PIP=PIP, n=n, h2=h2, gamma=gamma, pi=pi))
    }))
}

###################
###################
###################

# simple BayesC
h <- function(PIP, n = 3000, v = 0.01, h2 = 0.5, gamma=c(0,1), pi = c(0.99, 0.01)) {
  lam <- ncp(v, h2, n)
  uP <- u(PIP, n, h2, gamma, pi)
  B <- Bs(n, h2, gamma)
  4 / B[2] * pi[2] / PIP / (1-PIP) * dchisq(uP, 1, lam) * dnorm(sqrt(v), 0, sqrt(gamma[2] * h2)) 
}
plot(function(x) sapply(x, function(p) h(p)), xlim=c(0.01,0.99))
plot(function(x) sapply(x, function(p) g(p, n = 3000, v = 0.01, h2 = 0.5, gamma=c(0,1), pi = c(0.99, 0.01))), xlim=c(0.01,0.99))


# new implementation as a function of z
# use log scale to avoid underflow and overflow
log_fz <- function(v, i, z = 1, n = 1e8, h2 = 0.5,
                  gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                  pi = c(0.99, 0.004, 0.003, 0.002, 0.001))
{
  lam <- ncp(v, h2, n)
  
  ldchi(z, lam) + lnor(v, sqrt(gamma[i] * h2)) 
  
}

plot(function(x) sapply(x, function(z) exp(log_fz(z = z, v = 0.01, i = 5,  n = 10, h2 = 0.5,
                                                  gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                                                  pi = c(0.99, 0.004, 0.003, 0.002, 0.001)))), xlim=c(0,3))

# test using JZ's dPIP for BayesC
# 

P_2_z <- function(P, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  Pmin <- 1 - 1 / (1 + A)
  P <- max(P, Pmin)
  (log(P) - log(1-P) - log(A)) / B
}

z_2_P <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  1 - 1 / (1 + A * exp(B * z))
}

dPIP <- function(P, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  ncp = n * v / (1-h2)
  z = (log(P) - log(1-P) - log(A)) / B
  dchisq(z, 1, ncp ) / B / P / (1 - P)
}
integrate(dPIP, 0, 1, n = 3000, v = 0.01)

dPIP_reparam <- function(z, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  ncp = n * v / (1 - h2)
  dchisq(z, 1, ncp)
}
integrate(dPIP_reparam, 0, +Inf, n = 3000, v = 0.01)
# both integrate to 1

integrate(dPIP_reparam, P_2_z(0.5, n=3000, h2, m, pi), +Inf, n = 3000, v = 0.01)
integrate(dPIP, 0.5, 1, n = 3000, v = 0.01)

integrate(dPIP_reparam, P_2_z(0.2, n=3000, h2, m, pi), +Inf, n = 3000, v = 0.01)
integrate(dPIP, 0.2, 1, n = 3000, v = 0.01)
# same

# compare dPIP and dPIP reparam
integrate(dPIP_reparam, P_2_z(0.1, n = 3000), P_2_z(0.9, n=3000), n = 3000, v = 0.01)
integrate(dPIP, 0.1, 0.9, n = 3000, v = 0.01)
# same !

# now test ldchi
ldchi(1, 10)
dchisq(1, 1, 10, log=TRUE)

# test lnor
lnor(10, .1)
dnorm(sqrt(10),0, .1, log=TRUE)- 0.5 * log(10)
dchisq(10/0.1^2, 1,log=TRUE) - 2 * log(0.1)
# all the same

# naive implementation
dPIP_reparam <- function(z, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  #lambda <- (1 - h2) / h2 * m * pi
  #C <- n + lambda
  #A <- pi / (1 - pi) * sqrt(lambda / C)
  #B <- 0.5 * n * (1 - h2) / C
  ncp = n * v / (1 - h2)
  dchisq(z, 1, ncp) * dnorm(sqrt(v), 0, sqrt(h2 / m / pi)) / sqrt(v)
}

dPIP_stable <- function(z, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  ncp = n * v / (1 - h2)
  exp(ldchi(z, ncp) + lnor(v,sqrt(h2 / m / pi)))
}
# plot as a function of z
plot(function(x) sapply(x, function(z) dPIP_reparam(z, n = 3000, v = 1e-4)), xlim = c(0.1,5))

# plot as a function of v
plot(function(x) sapply(x, function(v) dPIP_reparam(P_2_z(0.9, n=3000), n = 3000, v)), xlim = c(0.000001,0.002))
plot(function(x) sapply(x, function(v) dPIP_stable(P_2_z(0.9, n=3000), n = 3000, v)), xlim = c(0.000001,0.002), col = 2, add=T)
# numerical errors ? seems strange
plot(function(x) sapply(x, function(v) dPIP_reparam(20, n = 3000, v)), xlim = c(0.000001,0.002))
plot(function(x) sapply(x, function(v) dPIP_stable(20, n = 3000, v)), xlim = c(0.000001,0.002), col=2, add=T)


naive_inner <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01)
  sapply(z, function(z) 
    integrate(function(x) dPIP_reparam(z=z, n=n, v = x, m=m, pi=pi), 0, +Inf, abs.tol = 0, subdivisions=10000)$value
  )


plot(function(x) sapply(x, function(z) naive_inner(z, n = 30000)), xlim = c(0.1,1))
plot(function(x) sapply(x, function(z) naive_inner_stable(z, n = 30000)), xlim = c(0.1,1))

# mc integration


# try integrate naive_inner ?
naive_outer <- function(n, h2 = 0.5, m = 1e6, pi = 0.01) {
  integrate(function(x) naive_inner(z=x, n=n, m=m, pi=pi), 0, +Inf, abs.tol = 0)$value
}

# BayesC power calculation
# z = u(P)
# y = ncp
log_f_reparam <- function(z, y, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  ldchi(z, exp(y)) + 
    lnor((1 - h2) / n * exp(y), sqrt(h2 / m / pi)) +
    log(1 - h2) - log(n) + y
} 
f_reparam <- function(z, y, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  exp(ldchi(z, exp(y)) + 
    lnor((1 - h2) / n * exp(y), sqrt(h2 / m / pi)) +
    log(1 - h2) - log(n) + y)
} 

fP <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  sapply(z, function(z) {
  integrate(
    function(y) exp(log_f_reparam(z, y = y, n = n, h2 = h2, m = m, pi = pi)),
    -15, 15, abs.tol = 0)$value
  })
}
plot(function(x) sapply(x, function(p) fP(z, n = 3000)), xlim = c(0.01,1))
plot(function(x) sapply(x, function(p) fP(P_2_z(x,n=3000), n = 3000)), xlim = c(0.01,1))

pow <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  integrate(fP, P_2_z(P0, n, h2, m, pi), Inf, n = n, h2 = h2, m = m, pi = pi, abs.tol=0)$value
}

j <- function(x, o, h2, q, pi) exp(log_f_reparam(x[1], y = x[2], n = o, h2 = h2, m = q, pi = pi))
pow_cub <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  hcubature(j, lower = c(P_2_z(P0, n, h2, m, pi),-30), upper = c(Inf,30),
            o = n, h2 = h2, q = m, pi = pi,
            tol = sqrt(.Machine$double.eps))
}

pow_cub2 <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  cubintegrate(j, lower = c(P_2_z(P0, n, h2, m, pi),-30), upper = c(Inf,30),
               method = "hcubature", absTol = 1e-12, maxEval = 1e7,#.Machine$double.eps,
               o = n, h2 = h2, q = m, pi = pi)
}

# example where it fails:
mean(sapply(runif(1000, 0.00502515, 1), function(p) fP(P_2_z(p , n=30000, h2, m, pi=0.01), n=30000, h2, m, pi=0.001)))
mean(sapply(runif(1000, 0, 1), function(u) fP(-log(u), n=30000, h2, m, pi=0.001)))
pow(0.1, n=30000, h2 = 0.1, m = 1e5, pi=0.001)
pow_cub(0.1, n=30000, h2 = 0.1, m = 1e5, pi=0.001)

# check fP with MC integration
f_mc_1 <- function(P, n, h2, m, pi) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  zs <- rnorm(N, 0, sqrt(h2 / m / pi))^2
  2 * mean(dchisq(log(P/(1-P)/A) / B, 1, n * zs / (1 - h2))) / P / (1 - P) / B
}
fP_1 <- function(P, n, h2, m, pi) {
  fP(P_2_z(P , n, h2, m, pi), n, h2, m, pi)
}

f_mc_1(0.1, 30000, h2, m, pi)
fP_1(0.1, 30000, h2, m, pi)
# not the same !


# compare 
# compare.MC integration and numerical integration on fP reparametrized in z
f_mc_2 <- function(z, n, h2, m, pi) {
  vs <- rnorm(N, 0, sqrt(h2 / m / pi))^2
  mean(dchisq(z, 1, n * vs / (1 - h2)))
}

fP(1, n = 3000, h2, m, pi)
f_mc_2(1, n = 3000, h2, m, pi)
fP(1, n = 3000, h2, m=1e5, pi)
f_mc_2(1, n = 3000, h2, m=1e5, pi)
# they match

# new attempt to get power by MC integraation
# following JZ's approach
f_mc_3 <- function(P0, n, h2, m, pi, N = 1000) {
  vs <- rnorm(N, 0, sqrt(h2 / m / pi)) ^2
  
  mean(sapply(vs, function(v) {
    ncp <- n * v / (1 - h2)
    mean(rchisq(N, 1, ncp) > P_2_z(P0, n, h2, m, pi))}))
}
f_mc_4 <- function(P0, n, h2, m, pi, N = 1000) {
  vs <- rnorm(N, 0, sqrt(h2 / m / pi)) ^2
  
  mean(sapply(vs, function(v) {
    ncp <- n * v / (1 - h2)
    mean(z_2_P(rchisq(N, 1, ncp), n, h2, m , pi) > P0)}))
}
pow(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_cub(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_3(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_4(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)

# power calculation fail for very large SNP effect for a few SNPS
pow(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
pow_cub(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
f_mc_3(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
f_mc_4(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)

# but OK for "decent" parameters
pow(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
pow_cub(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
f_mc_3(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
f_mc_4(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)

# new try with different integration order
P_v <- function(y, P0, n, h2, m, pi) {
  # vectorize for integration
  sapply(y, function(this_y)
  integrate(f_reparam, P_2_z(P0, n, h2, m, pi), Inf, 
            y = this_y, n = n, h2 = h2, m = m, pi = pi,
            abs.tol = 0)$value)
}
P_v_cub <- function(y, P0, n_, h2, m, pi) {
  # vectorize for integration
  sapply(y, function(this_y)
    hcubature(f_reparam, P_2_z(P0, n_, h2, m, pi), Inf,
              y = this_y, n = n_, h2 = h2, m = m, pi = pi,
              tol = sqrt(.Machine$double.eps))$integral)
}
pow_inv <- function(P0, n, h2, m, pi) {
  integrate(P_v,
            -20, +20, P0 = P0, n = n, h2 = h2, m = m, pi = pi,
            abs.tol = 0, subdivisions = 10000)
}
pow_inv_cub <- function(P0, n, h2, m, pi) {
  hcubature(P_v_cub, -20, +20,
            P0 = P0, n_ = n, h2 = h2, m = m, pi = pi,
            tol = sqrt(.Machine$double.eps))
}
pow(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_inv(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_inv_cub(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_3(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_4(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)

# power calculation fail for very large SNP effect for a few SNPS
pow(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
pow_inv(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
pow_inv_cub(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
pow_cub(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
f_mc_3(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)
f_mc_4(0.2, n=30000, h2 = 0.7, m = 5e3, pi=0.01)

# but OK for "decent" parameters
pow(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
pow_inv(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
pow_cub(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
pow_inv_cub(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
f_mc_3(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)
f_mc_4(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)

# another try with double chi-square
g <- function(z, v, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  exp(ldchi(z, n * v/(1 - h2)) + 
        dchisq(v * m * pi / h2, 1, log = TRUE) - log(m * pi / h2))
} 
g_cub <- function(x, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  z <- x[1]; v <- x[2];
  exp(ldchi(z, n * v / (1 - h2)) + 
        dchisq(v * m * pi / h2, 1, log = TRUE) + log(m * pi / h2))
} 
pow_g_cub <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  hcubature(g_cub, lower = c(P_2_z(P0, n, h2, m, pi), 0), upper = c(Inf,Inf),
          n = n, h2 = h2, m = m, pi = pi,
          tol = sqrt(.Machine$double.eps))
}
h_cub <- function(x, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  z <- x[1]; y <- x[2];
  exp(dchisq(z, 1, n * y * h2 / (1 - h2) / m / pi, log = TRUE) + 
        dchisq(y, 1, log = TRUE))
}
h_cub_ap <- function(x, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  z <- x[1]; y <- x[2];
  exp(ldchi_ap(z, n * y * h2 / (1 - h2) / m / pi) + 
        dchisq(y, 1, log = TRUE))
}
pow_h_cub <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  hcubature(h_cub, lower = c(P_2_z(P0, n, h2, m, pi), 0), upper = c(Inf,Inf),
            n = n, h2 = h2, m = m, pi = pi,
            tol = sqrt(.Machine$double.eps))
}
pow_h_cub_ap <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  hcubature(h_cub_ap, lower = c(P_2_z(P0, n, h2, m, pi), 0), upper = c(Inf,Inf),
            n = n, h2 = h2, m = m, pi = pi,
            tol = sqrt(.Machine$double.eps))
}

pow_g_cub(0.2, n=30000, h2 = 0.7, m = 1e5, pi=0.01)

pow(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_inv(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_inv_cub(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_g_cub(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
pow_h_cub(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_3(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)
f_mc_4(0.2, n=30000, h2 = 0.5, m = 1e6, pi=0.001)

# experiment using pracma, with a change of variables to (0,1)

i_prac <- function(u, t, n, h2, m, pi) {
  exp(ldchi(-log(t), -n * log(u) / (1 - h2)) - log(t) - 0.5 * log(-log(u)) - log(u) +
        dnorm(sqrt(-log(u)), sqrt(h2 / m / pi), log = TRUE))
}
i_prac <- function(u, t, n, h2, m, pi) {
  exp(ldchi(-log(t), -n * log(u) / (1 - h2)) - log(t) - log(u) +
        lnor(-log(u), sqrt(h2 / m / pi)))
}
# test with normal integration, compare to fP
integrate(i_prac, 0, 1, t = exp(-100), n=30000, h2 = 0.5, m = 1e6, pi=0.001, abs.tol = 0)
cuhre(i_prac, nComp = 1, 0, 1, t = exp(-100), n=30000, h2 = 0.5, m = 1e6, pi=0.001, absTol = 0)

fP(100, 30000, 0.5, 1e6, 0.001)
# there is inconsistency between integration methods

ldchi_ap <- function(x, ncp) {
  dnorm(x, 1 + ncp, sqrt(2 + 4 * ncp), log = TRUE)
}

f_reparam_ap <- function(z, y, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
  exp(ldchi_ap(z, exp(y)) + 
        lnor((1 - h2) / n * exp(y), sqrt(h2 / m / pi)) +
        log(1 - h2) - log(n) + y)
} 

fP_ap <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  sapply(z, function(z) {
    integrate(
      function(y) f_reparam(z, y = y, n = n, h2 = h2, m = m, pi = pi),
      -15, 15, abs.tol = 0)$value
  })
}

pow_ap <- function(P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  integrate(fP_ap, P_2_z(P0, n, h2, m, pi), Inf, n = n, h2 = h2, m = m, pi = pi, abs.tol=0)$value
}

# Compare fP and fP_ap
fP(10, 1e6)
fP_ap(10, 1e6)

# compare power
pow(0.2, n=1e6, h2 = 0.5, m = 1e6, pi=0.001)     
pow_ap(0.2, n=1e6, h2 = 0.5, m = 1e6, pi=0.001)
# works well but is not fastre

f_ap2 <- function(z, v, n = 300, h2 = 0.5, m = 1e6, pi = 0.01) {
  ncp <- n * v / (1 - h2)
  sigma_chi <- 2 + 4 * ncp
  1 / 2 / 3.14159265359 / sqrt(v) / sqrt(sigma_chi * h2 / pi / m) *
    exp(-v / 2 / h2 * pi * m - (z - 1 - ncp)^2 / 2 / sigma_chi)
}

f_ap_test <- function(z, v, n = 300, h2 = 0.5, m = 1e6, pi = 0.01) {
  ncp <- n * v / (1 - h2)
  sigma_chi <- 2 + 4 * ncp
  dnorm(z, 1 + ncp, sqrt(2 + 4 * ncp)) / sqrt(v) *
    dnorm(sqrt(v), 0, sqrt(h2 / pi / m))
}

fP_ap2 <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  sapply(z, function(z) {
    integrate(
      function(v) f_ap2(z, v = v, n = n, h2 = h2, m = m, pi = pi),
      0, 100, abs.tol = 0)$value
  })
}

# same with change of variable
f_ap3 <- function(z, y, n = 300, h2 = 0.5, m = 1e6, pi = 0.01) {
  sigma_chi <- 2 + 4 * exp(y)
  sqrtv <- sqrt((1 - h2) * exp(y) / n)
  exp(dnorm(z, 1 + exp(y), sqrt(sigma_chi), log = TRUE) +
        dnorm(sqrtv, 0, sqrt(h2 / m / pi)) + log(sqrtv))
}



fP(10, 1e6)
fP_ap(10, 1e6)
fP_ap2(10, 1e6)

# Eureka ?

# inner function nc chi square
h_e <- function(z, l, n = 300, h2 = 0.5, m = 1e6, pi = 0.01)
{
#  cat("z: ", z, "l:", l, "\n")
  exp(ldchi(z, l))
}
p_e <- function(y, z0, n, h2, m, pi) {
  integrate(h_e, lower = z0, upper = 10000, 
            l = y, n = n, h2 = h2, m = m, pi = pi,
            abs.tol = 0)$value
}

q_e <- function(v, z0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  ncp <- n * v / (1 - h2)
  dnorm(sqrt(v), 0, sqrt(h2/ m / pi)) / sqrt(v) *
    pchisq(z0, 1, ncp, lower.tail = FALSE)
}

q_e_2 <- function(v, P0, n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  exp(lnor(v, sqrt(h2/ m / pi)) +
    pchisq(P_2_z(P0, n, h2, m, pi), 1, n * v / (1 - h2), lower.tail = FALSE, log.p = TRUE))
}
pow_e <- function(P0, n, h2, m, pi) {

  z0 <- P_2_z(P0, n, h2, m, pi)
  #ncp_max <- 10000
  #v_max <- ncp_max * (1 - h2) / n
  #v_max <- min(v_max, 10 * sqrt(h2 / m / pi))
  # limit integration to 10 standard deviations
  v_max <- 10 * sqrt(h2 / m / pi)
  
  cat("z0: ", z0, "ncp_max: ", ncp_max, "v_max: ", v_max, "\n")
  integrate(q_e, lower = 0, upper = v_max,
            z0 = z0,  n = n, h2 = h2, m = m, pi = pi,
            subdivisions = 10000, abs.tol = 0)$value
}

pow_e_cub <- function(P0, n, h2, m, pi) {
  
  z0 <- P_2_z(P0, n, h2, m, pi)
  #ncp_max <- 10000
  #v_max <- ncp_max * (1 - h2) / n
  #v_max <- min(v_max, 10 * sqrt(h2 / m / pi))
  # limit integration to 10 standard deviations
  v_max <- 10 * sqrt(h2 / m / pi)
  
  cat("z0: ", z0, "ncp_max: ", "unused", "v_max: ", v_max, "\n")
  hcubature(q_e, lower = 0, upper = Inf,
            z0 = z0,  n = n, h2 = h2, m = m, pi = pi,
            vectorInterface = TRUE, absError = 0)
}
# per sSNP variance explained
