library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 900))

cl <- makeDimerDihedral(d=100e-3, 
              dihedral=0*pi/180, alpha1=10*pi/180, alpha2=0,
              a=35e-3, b=12e-3, 
              right=TRUE)

# visualise
#rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")

linear <- linear_extinction_spectrum(cl, gold)
p1 <- ggplot(linear, aes(wavelength, value, color=variable)) + geom_path()

circular <- circular_dichroism_spectrum(cl, gold,N=300, averaging="grid")

p2 <- ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~variable, scales="free") + geom_path()

p2