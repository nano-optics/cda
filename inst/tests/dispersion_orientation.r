library(cda)

library(cda)
library(ggplot2)

wvl <- seq(400,900, by=1)
gold <- epsAu(wvl)


clust <- makeRodChain(N=2, pitch=200e-3)

n <- 1.5
k0 <- 2 * pi/gold$wavelength
kn <- k0 * n
invalpha <- make.invalpha(clust, gold, polarizability.fun = polarizability.ellipsoid, 
                          n = n, kuwata = TRUE)

Angles <- rep(seq(0, pi/2, length=12), 3)
Axes <- rep(letters[24:26], each=12)

clust <- makeRodChain(N=1)
clust <- makeRodChain(N=1, a=0.03, b=0.03, c=0.05)
clust <- makeRodChain(N=1, a=0.03, b=0.05, c=0.03)

results <- dispersion_spectrum(clust, Angles, Axes, gold, 
                               polarisation="linear", progress=TRUE)

test <- melt(results, meas="value")

ggplot(subset(test, type == "extinction"), 
       aes(wavelength, value, colour=angles, group=angles)) +
  facet_grid(axes ~ polarisation, scales="free") +
  geom_path()
