# MC integration tests

# integral dnorm, 0 +inf
# rewrite y = exp(-x)

us <- runif(N)
mean(1/us*dnorm(-log(us)))

# double integral
as <- runif(N)
bs <- runif(N)
mean(dchisq(-log(as), 1, -log(bs) * n / (1 - h2)) / as * dnorm(sqrt(-log(bs)), 0, sqrt(h2 / m / pi)) / sqrt(-log(bs)) / bs)
     