library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 900))

cl <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12, 
              right=TRUE)

params <- expand.grid(N=c(500, 1000, 5000),
                       averaging=c("grid", "GL", "QMC"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl, material=gold)


p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=N, group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))+
         theme_minimal()

p
