library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 900))

cl <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12, 
              right=TRUE)

# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")

linear <- linear_extinction_spectrum(cl, gold)
p1 <- ggplot(linear, aes(wavelength, value, color=variable)) + geom_path()

circular <- circular_dichroism_spectrum(cl, gold)

p2 <- ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~variable, scales="free") + geom_path()

p2