
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

# define a cluster of particles
cl <- list(r = rbind(c(0, 0, 0),
                      c(0, 0, 200)),
            angles = rbind(c(0, 0, 0),
                           c(pi/4, 0, 0)),
            sizes = rbind(c(40, 20, 20),
                          c(40, 20, 20)))

# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()



## @knitr linear
# calculate extinction spectrum at fixed incidence

linear <- linear_extinction_spectrum(cl, gold)
ggplot(linear, aes(wavelength, value, color=variable)) + geom_path()



## @knitr oa
circular <- circular_dichroism_spectrum(cl, gold)

ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~variable, scales="free") + geom_path()



