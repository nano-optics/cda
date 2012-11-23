library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 600, by=100))

cl <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12, 
              right=TRUE)

beta <- inverse_polarizability(cl, gold)
head(beta)

