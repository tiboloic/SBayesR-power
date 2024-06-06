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


# final figure including the dots indicating the SNP and credible set outcome
# 
pred <- function(trait, h2, m, kappa = 0,
                 gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2)) {
  m_tot <- sum(m)
  pis <- m / m_tot
  ns <- 10 ^ seq(4, 7, 0.1)
  power <-  sapply(ns, function(n) pow_R_e(0.9, n = n, h2 = h2,
                                           gammas = gammas, pis = pis))
  numFM <- m_tot * (1 - pis[1]) * power
  PVE <- sapply(ns, function(n) prop_var_R(0.9, n = n, h2 = h2,
                                        gammas = gammas, pis = pis))
  data.frame(Trait = trait, Kappa = kappa, N = ns, Power = power,
             NumFM = numFM, PVE = PVE)
}

pred2 <- function(trait, h2, m, kappa = 0,
                 gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2)) {
  m_tot <- sum(m)
  pis <- m / m_tot
  ns <- 10 ^ seq(4, 7, 0.1)
  power <-  sapply(ns, function(n) pow_R_mc(0.9, n = n, h2 = h2,
                                           gammas = gammas, pis = pis))
  numFM <- m_tot * (1 - pis[1]) * power
  PVE <- sapply(ns, function(n) prop_var_R(0.9, n = n, h2 = h2,
                                           gammas = gammas, pis = pis))
  data.frame(Trait = trait, Kappa = kappa, N = ns, Power = power,
             NumFM = numFM, PVE = PVE)
}
## prediction based on genetic architecture estimates
height.mc <- pred(trait = "Height", h2 = 0.68, m = c(9262822.000000, 60758.699219, 2834.479980, 190.014999, 0.135000))
scz.mc <- pred(trait = "SCZ", h2 = 0.2633, m = c(7276536.7375, 77483.0705, 2495.0030, 2.6315, 0.0575))
bip.mc <- pred(trait = "BIP", h2 = 0.1158, m = c(7301780.265, 49350.143, 5383.282, 4.282, 0.028))
bip.gk <- pred(trait = "BIP", h2 = 0.1158, m = c(7301780.265, 49350.143, 5383.282, 4.282, 0.028))
mdd.mc <- pred(trait = "MDD", h2 = 0.0752, m = c(7267886.1945, 86987.6690, 1642.5530, 1.5810, 0.0025))  
cd.mc <- pred(trait = "CD", h2 = 0.290922, m = c(9273676.000000, 50740.875000, 1912.599976, 269.864990, 6.315000))  
cd.mc2 <- pred2(trait = "CD", h2 = 0.290922, m = c(9273676.000000, 50740.875000, 1912.599976, 269.864990, 6.315000))  

simLD.mc <- pred(trait = "SimLD", h2 = 0.5, m = c(1028069, 0, 10145, 4, 0))
simLE.mc <- pred(trait = "SimLE", h2 = 0.5, m = c(53243.781250, 0, 10186.235352, 2.985000, 0))

res.mc = rbind(height.mc, cd.mc, scz.mc, bip.mc, mdd.mc, simLD.mc, simLE.mc)
res.mc$Trait = factor(res.mc$Trait, levels=c("SimLD", "SimLE","Height","CD","SCZ","BIP","MDD"))

sim <- res.mc
sim$Kappa <- as.factor(sim$Kappa)
sim <- subset(sim, Trait != "SimLE")
sim$Trait = factor(sim$Trait, levels=c( "SimLD","Height","CD","SCZ","BIP","MDD"))

#### LD-based approach fine-mapping results
## observed results for these traits with different GWAS datasets
## PEP = 0.7
res.obs = data.frame()
res.obs = rbind(res.obs, c("Height", 346385, "SNP+CS", 224+1703, 0.0320, 0.3517))
res.obs = rbind(res.obs, c("CD", 52965, "SNP+CS", 19+72, 0.0016, 0.1931))
res.obs = rbind(res.obs, c("SCZ", 225255, "SNP+CS", 13+200, 0.0024, 0.0417))
res.obs = rbind(res.obs, c("BIP", 269864, "SNP+CS", 5+122, 0.0019, 0.0350))
res.obs = rbind(res.obs, c("MDD", 1480703, "SNP+CS", 28+596, 0.0068, 0.0573))
res.obs = rbind(res.obs, c("SimLD", 100000, "SNP+CS", 92+912, 0.0928, 0.3168))
names(res.obs) = c("Trait", "N", "Outcome", "NumFM","Power", "PVE")
res.obs$N = as.numeric(res.obs$N)
res.obs$NumFM = as.numeric(res.obs$NumFM)
res.obs$Power = as.numeric(res.obs$Power)
res.obs$PVE = as.numeric(res.obs$PVE)
res.obs$Trait = factor(res.obs$Trait, levels=c( "SimLD","Height","CD","SCZ","BIP","MDD"))

p1 = ggplot(subset(sim, Kappa==0), aes(N, NumFM)) + geom_line(aes(col=Trait, linetype=Trait)) + 
  geom_point(data=subset(res.obs, Outcome=="SNP+CS"), aes(col=Trait)) +
  scale_y_log10(breaks = 10^seq(0,5,1), labels = 10^seq(3,8,1)/1000) + 
  scale_x_log10(breaks = 10^seq(3,8,1), labels = 10^seq(3,8,1)/1000) +
  scale_shape_manual(values=c(3, 16)) +
  ylab("Number of identified causal variants (thousand)") + xlab("GWAS sample size (thousand)") +
  annotation_logticks(sides="bl") +
  coord_cartesian(ylim = c(1, 10^5))
ggsave("nfm_traits_02June24.pdf", p1, width=8, height=6)

p2 = ggplot(subset(sim, Kappa==0), aes(N, Power)) + geom_line(aes(col=Trait, linetype=Trait)) + 
  geom_point(data=subset(res.obs, Outcome=="SNP+CS"), aes(col=Trait)) +
  scale_x_log10(breaks = 10^seq(3,8,1), labels = 10^seq(3,8,1)/1000) +
  annotation_logticks(sides="b") +
  ylab("Power") + xlab("GWAS sample size (thousand)")
ggsave("power_traits_02June24.pdf", p2, width=8, height=6)

p3 = ggplot(subset(sim, Kappa==0), aes(N, PVE)) + geom_line(aes(col=Trait, linetype=Trait)) + 
  geom_point(data=subset(res.obs, Outcome=="SNP+CS"), aes(col=Trait)) +
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  scale_x_log10(breaks = 10^seq(3,8,1), labels = 10^seq(3,8,1)/1000) + 
  annotation_logticks(sides="b") +
  ylab("Proportion of SNP-based heritability") + xlab("GWAS sample size (thousand)")
ggsave("pgv_traits_02June24.pdf", p3, width=8, height=6)