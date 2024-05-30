# LT 23/05/2024
# 
#  Figure 7:power and proportion of variance
#  

# height
h2 <- 0.68
pis <- c(9262822.000000, 60758.699219, 2834.479980, 190.014999, 0.135000)
gammas <- c(0, 1e-5, 1e-4, 1e-3, 1e-2)

# sample sizes
ns <- 10 ^ seq(4, 7, 0.1)

ht_power <- sapply(ns, function(n) pow_R_e(0.9, n = n, h2 = h2,
                               gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                               pis = pis))
ht_power_mc <- sapply(ns, function(n) pow_R_mc(0.9, n = n, h2 = h2,
                                           gammas = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                                           pis = pis))
plot(ht_power ~ ns, log = 'xy', type = "l", col = 2)
lines(ht_power_mc ~ ns, col = 3)
