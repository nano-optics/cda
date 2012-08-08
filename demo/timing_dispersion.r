## testing the direct inversion of the interaction matrix vs solving the linear system

library(cda)
library(ggplot2)

wvl <- seq(400,900, by=1)
gold <- epsAu(wvl)

clust1 <- makeRodChain(N=2)
clust2 <- makeRodChain(N=20)
Angles1 <- as.matrix(expand.grid(theta=0, psi=seq(0, pi, length=2), phi=pi/2))
Angles2 <- as.matrix(expand.grid(theta=0, psi=seq(0, pi, length=500), phi=pi/2))
results <- dispersion_spectrum(clust1, Angles1, gold, progress=FALSE)
test <- melt(results, meas="value")

ggplot(test, aes(wavelength, value, colour=angles.psi, group=angles.psi)) +
  facet_grid(type ~ polarisation, scales="free") +
  geom_path()

## small number of angles, small system
system.time(results1 <- dispersion_spectrum(clust1, Angles1, gold, progress=FALSE, invert=TRUE))
system.time(results2 <- dispersion_spectrum(clust1, Angles1, gold, progress=FALSE, invert=FALSE))


## large number of angles, small system
system.time(results1 <- dispersion_spectrum(clust1, Angles2, gold, progress=FALSE, invert=TRUE))
system.time(results2 <- dispersion_spectrum(clust1, Angles2, gold, progress=FALSE, invert=FALSE))

## small number of angles, large system
system.time(results1 <- dispersion_spectrum(clust2, Angles1, gold, progress=FALSE, invert=TRUE))
system.time(results2 <- dispersion_spectrum(clust2, Angles1, gold, progress=FALSE, invert=FALSE))

