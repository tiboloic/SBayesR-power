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

