## solving CD equations for two polarisations and a range of angles

library(cda)
library(ggplot2)

wvl <- seq(400,900, by=1)
gold <- epsAu(wvl)

clust <- makeSpheresCluster(N=5, radius=5e-3, R0=12e-3, pitch=15e-3,
                            delta=pi/2, right=TRUE)

## clust <- makeRodChain(N=2)

Angles <- as.matrix(expand.grid(theta=pi/2, psi=seq(0, pi/2, length=10), phi=0))
Angles <- as.matrix(expand.grid(theta=pi/2, psi=seq(0, pi, length=20), phi=seq(0, pi, length=20)))
## results <- dispersion_spectrum(clust, Angles, gold, polarisation="linear", progress=FALSE)
results <- dispersion_spectrum(clust, Angles, gold, polarisation="circular", progress=FALSE)
test <- melt(results, meas="value")

## ggplot(test, aes(wavelength, value, colour=angles.psi, group=angles.psi)) +
##   facet_grid(type ~ polarisation, scales="free") +
##   geom_path()

ggplot(subset(test, polarisation == "CD"), aes(wavelength, value,
                      colour=angles.phi, group=interaction(angles.psi, angles.phi))) +
  ## facet_grid(type ~ polarisation, scales="free") +
  geom_line()
