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

# re implement robust non central chi-square with 1 degree of freedom
ldchi <- function(x, lambda) {
  log(0.5) -0.5 * (sqrt(lambda) - sqrt(x))^2 - 0.25 * log(x/lambda) +
    log(besselI(sqrt(x * lambda), -0.5, expon.scaled = TRUE))
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

# cut integral in 2
P_0_1 <- function(PIP = 0.9, n = 1e8, h2 = 0.5,
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
                                                   pi = pi)), -10, 10, abs.tol=0)
    if (int$message != "OK")
      stop(paste0("Inner integral failed with message: ", int$message))
    
    4 * pi[i] / sum(Bs(n, h2, gamma)) / PIP * int$value
  }))
}

g <- function(PIP, v = 0.01, n = 1e8, h2 = 0.5,
              gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
              pi = c(0.9, 0.04, 0.03, 0.02, 0.01)) {
  sum(sapply(2:length(gamma), function(i) {
    4 * pi[i] / sum(Bs(n, h2, gamma)) / PIP / (1-PIP) *
      exp(log_f(v=v, i=i, PIP=PIP, n=n, h2=h2, gamma=gamma, pi=pi))
    }))
}

# simple BayesC
h <- function(PIP, n = 3000, v = 0.01, h2 = 0.5, gamma=c(0,1), pi = c(0.99, 0.01)) {
  lam <- ncp(v, h2, n)
  uP <- u(PIP, n, h2, gamma, pi)
  B <- Bs(n, h2, gamma)
  4 / B[2] * pi[2] / PIP / (1-PIP) * dchisq(uP, 1, lam) * dnorm(sqrt(v), 0, sqrt(gamma[2] * h2)) 
}
plot(function(x) sapply(x, function(p) h(p)), xlim=c(0.01,0.99))
plot(function(x) sapply(x, function(p) g(p, n = 3000, v = 0.01, h2 = 0.5, gamma=c(0,1), pi = c(0.99, 0.01))), xlim=c(0.01,0.99))


# new implementation as a fonction of z
# 
fz <- function() {
  
}

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

P_2_z <- function(P, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  (log(P) - log(1-P) - log(A)) / B
}

z_2_P <- function(z, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
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

# compare dPIP and dPIP reparam
integrate(dPIP_reparam, P_2_z(0.1, n = 3000, v = 0.01), +Inf, n = 3000, v = 0.01)
integrate(dPIP, 0.1, 1, n = 3000, v = 0.01)
# same !

# now test ldchi
ldchi(1, 10)
dchisq(1, 1, 10, log=TRUE)

# test lnor
lnor(10,1)
dnorm(sqrt(10),0,1, log=TRUE)- 0.5 * log(10)

# naive implementation
dPIP_reparam <- function(z, n, v, h2 = 0.5, m = 1e6, pi = 0.01) {
  lambda <- (1 - h2) / h2 * m * pi
  C <- n + lambda
  A <- pi / (1 - pi) * sqrt(lambda / C)
  B <- 0.5 * n * (1 - h2) / C
  ncp = n * v / (1 - h2)
  dchisq(z, 1, ncp) * dnorm(sqrt(v), 0, sqrt(h2 / m / pi)) / sqrt(v)
}
plot(function(x) sapply(x, function(z) dPIP_reparam(z, n = 3000, v = 1e-4)), xlim = c(0.1,5))

naive_inner <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01)
  sapply(z, function(z) 
    integrate(function(x) dPIP_reparam(z=z, n=n, v = x, m=m, pi=pi), 0, +Inf, abs.tol = 0, subdivisions=1000)$value
  )

plot(function(x) sapply(x, function(z) naive_inner(z, n = 30000)), xlim = c(0.1,10))

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

fP <- function(z, n, h2 = 0.5, m = 1e6, pi = 0.01) {
  sapply(z, function(z) {
  integrate(
    function(y) exp(log_f_reparam(z, y = y, n = n, h2 = h2, m = m, pi = pi)),
    -15, 15, abs.tol = 0)$value
  })
}
plot(function(x) sapply(x, function(z) fP(z, n = 3000)), xlim = c(0.01,10))

pow <- function(n = 3000, h2 = 0.5, m = 1e6, pi = 0.01) {
  integrate(fP, 0.000, Inf, n = n, h2 = h2, m = m, pi = pi, abs.tol=0)$value
}
