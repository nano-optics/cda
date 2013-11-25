require(cda)

cl <- cluster_dimer(dihedral = 45 * pi/180, d=100, a=35, b=12, 
                    alpha1=10*pi/180, alpha2=0)

test <- dispersion_spectrum(cl, epsAu(300:900), angles=seq(0,pi/2,length=10), axis="y")

require(ggplot2)
ggplot(test, aes(wavelength, value, colour=angles,group=angles)) +
  geom_line() + facet_grid(type~polarisation)