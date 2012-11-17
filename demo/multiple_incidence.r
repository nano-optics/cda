## solving CD equations for two polarisations and a range of angles

library(cda)
library(ggplot2)

wvl <- seq(400,900, by=1)
gold <- epsAu(wvl)

clust <- makeSpheresCluster(N=5, radius=5e-3, R0=12e-3, pitch=15e-3,
                            delta=pi/2, right=TRUE)

if(interactive() && require(rgl)){ # display RGL window
open3d()
rgl.ellipsoids(clust$r, clust$sizes, clust$angles, col="gold")
axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')

# rgl.snapshot( "helix-rgl.png", fmt="png", top=TRUE )
}

angles <- rep(seq(0, pi/2, length=10), 3)
axis <- rep(letters[24:26],each=10)
results <- dispersion_spectrum(clust, angles, axis = axis, gold, 
                               polarisation="circular", progress=FALSE)
test <- melt(results, meas="value")

ggplot(subset(test, polarisation == "CD"), aes(wavelength, value,
                      colour=angles)) +
 facet_grid(axis ~ polarisation, scales="free") +
  geom_line()
