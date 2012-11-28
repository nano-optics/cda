
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

# define a dimer with dihedral angle
cl <- cluster_dimer(d=100, 
                    dihedral=45*pi/180, alpha1=10*pi/180, alpha2=0,
                    a=35, b=12)

# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()



## @knitr comparison
  
dimer <- function(dihedral=45, ...){
  cl <- cluster_dimer(dihedral = dihedral * pi/180, ...)
  circular_dichroism_spectrum(cl, material = gold)
  
}
params <- data.frame(dihedral=c(45, 30, 10, 0, -10, -30,-45),
                     d=100, a=35, b=12, 
                     alpha1=10*pi/180, alpha2=0)
comparison <- mdply(params, dimer)

p <- 
  ggplot(data=comparison) + 
  facet_grid(type~variable, scales="free") +
  geom_line(aes(wavelength, value, 
                colour=factor(dihedral))) +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="dihedral angle") +
         scale_colour_brewer(type="div", palette=3)

p


