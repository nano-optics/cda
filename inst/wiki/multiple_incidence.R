
## @knitr load
library(cda)
library(rgl)
library(ggplot2)


## @knitr setup

rgl_annotate = function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
}
theme_set(theme_minimal())


## @knitr cluster

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# two clusters
cl <- cluster_chain(N=2, pitch=100)

cl2 <- cluster_helix(5, R0=200, pitch=150, 
                          delta=pi/2, delta0=0, right=TRUE,
                          a=50, b=50/2, c=50/2,
                          angles="helix")
hel <- helix(N = 5, R0 = 200, pitch = 150, delta = pi/2, 
             delta0 = 0, right = TRUE)

# visualise
rgl.ellipsoids(cl2$r, cl2$sizes, cl2$angles, col="gold")
lines3d(hel$smooth, lwd=1, col="red")
  shift <- cbind(rep(1000, nrow(cl$r)),0,0)
rgl.ellipsoids(cl$r+cbind(rep(-500, 2),0,0), cl$sizes, cl$angles, col="gold")


## @knitr cd

Angles <- rep(seq(0, pi/2, length=12), 3)
Axes <- rep(c('x','y','z'), each=12)

results <- dispersion_spectrum(cl, Angles, Axes, gold, 
                               polarisation="linear")

test <- melt(results, meas="value")

ggplot(subset(test, type == "extinction"), 
       aes(wavelength, value, colour=angles, group=angles)) +
  facet_grid(axis ~ polarisation, scales="free") +
  geom_line() +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="incident angle")


## @knitr comparison

variables <- expand.grid(Angles = seq(0, 2*pi, length=36),
                         Axes = c('x','y','z'))

results <- dispersion_spectrum(cl2, variables$Angles, variables$Axes, gold, 
                               polarisation="circular")
average <- circular_dichroism_spectrum(cl2, material = gold)

ggplot(subset(results, polarisation == "CD"), aes(wavelength, value)) +
  facet_grid(axis ~ polarisation, scales="free") +
  geom_line(aes(colour=angles, group=angles)) +
  geom_line(data=subset(average, type == "CD" & variable =="extinction"), linetype=2, size=1.2) +
  labs(y=expression(sigma[CD]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="incident angle")


