library(cda)
library(ggplot2)
library(reshape2)
library(plyr)


# dielectric function
wvl <- seq(400, 700)
gold <- epsAu(wvl)

dimer <- function(dihedral=45, ...){
  cl <- cluster_dimer(dihedral = dihedral * pi/180, ...)
  circular_dichroism_spectrum(cl, material = gold)
  
}
params <- data.frame(dihedral=c(0),
                     d=200, a=50, b=30, 
                     alpha1=0, alpha2=0)
comparison <- mdply(params, dimer)

p <- 
  ggplot(data=subset(comparison, type!="CD")) + 
  facet_grid(type~., scales="free") +
  geom_line(aes(wavelength, value, 
                colour=variable)) +
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="dihedral angle") 

p
