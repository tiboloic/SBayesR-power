# LT 23/05/2024
# 
#  Figure 7:power and proportion of variance
#  

# sample sizes
ns <- 10 ^ seq(4, 7, 0.1)

# height
h2 <- 0.68
ms <- c(9262822.000000, 60758.699219, 2834.479980, 190.014999, 0.135000)
m_tot <- sum(ms)
pis <- ms / m_tot
gammas <- c(0, 1e-5, 1e-4, 1e-3, 1e-2)

ht_power <- m_tot * (1 - pis[1]) * sapply(ns, function(n) pow_R_e(0.9, n = n, h2 = h2,
                               gammas = gammas, pis = pis))
ht_power_mc <- m_tot * (1 - pis[1]) * sapply(ns, function(n) pow_R_mc(0.9, n = n, h2 = h2,
                                           gammas = gammas, pis = pis))
plot(ht_power ~ ns, log = 'xy', type = "l", col = 2)
lines(ht_power_mc ~ ns, col = 3)

# proportion of variance explained
ht_var_exp <- sapply(ns, function(n) prop_var_R(0.9, n = n, h2 = h2,
                                                   gammas = gammas, pis = pis))
plot(ht_var_exp ~ ns, log = 'x', type = "l", col = 2)

# SimLD
h2 <- 0.5
ms <- c(1028069, 0, 10145, 4, 0)
m_tot <- sum(ms)
pis <- ms / m_tot

simld_power <- sapply(ns, function(n) pow_R_e(0.9, n = n, h2 = h2,
                                           gammas = gammas, pis = pis))
# problem with 0 pis

# SCZ
h2 <- 0.2633
ms <- c(7276536.7375, 77483.0705, 2495.0030, 2.6315, 0.0575)
m_tot <- sum(ms)
pis <- ms / m_tot

scz_power <- m_tot * (1 - pis[1]) * sapply(ns, function(n) pow_R_e(0.9, n = n, h2 = h2,
                                              gammas = gammas, pis = pis))
plot(scz_power ~ ns, log = 'xy', type = "l", col = 2)

