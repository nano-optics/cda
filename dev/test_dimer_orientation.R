library(cda)
library(ggplot2)
library(reshape2)
library(plyr)


# dielectric function
wvl <- seq(300, 600)
gold <- epsAg(wvl)

dimer <- function(dihedral=45, ...){
  cl <- cluster_dimer(dihedral = dihedral * pi/180, ...)
  dispersion_spectrum(cl, material= gold, medium = 1.33, angles=c(pi/2,pi/2),
                      axes = c("x","z"), polarisation = c("linear"))
  
}
params <- data.frame(dihedral=c(0),
                     d=40, a=30, b=30, 
                     alpha1=0, alpha2=0)
comparison <- mdply(params, dimer)

p <- 
  ggplot(data=subset(comparison, type=="extinction")) + 
  facet_grid(polarisation~., scales="free") +
  geom_line(aes(wavelength, value,group=axes, 
                colour=axes)) +
  geom_line(aes(wavelength, value,group=axes, 
                colour=axes),data=ref,lty=2) +
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="pol") 

p
