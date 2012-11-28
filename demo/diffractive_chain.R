
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

cl <- cluster_chain(N=10, pitch=500)
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)


## @knitr comparison

chain <- function(N, pitch = 500, ...){
  cl <- cluster_chain(N, pitch, ...)
  linear_extinction_spectrum(cluster = cl, material = gold)
}
  
params <- data.frame(N=c(1, 10, 50))
comparison <- mdply(params, chain)

p <- 
  ggplot(data=comparison)+
labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm),
       colour = expression(N), linetype=expression(polarisation))+
  geom_line(aes(wavelength, value, linetype=variable,
                colour=factor(N),
                group=interaction(N,variable)))

p


