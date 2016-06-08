## ---- setup,echo=FALSE ---------------
library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
theme_set(theme_bw() + theme(strip.background=element_blank()))

## ---- demo ---------------

# dielectric function
material <- epsAu(seq(500, 800))

# cluster geometry
cl <- cluster_dimer(d = 50, a=30, b=10, c=10, dihedral = pi/4)

# simulation
results <- spectrum_oa(cl, material)


## ---- plot,echo=TRUE,fig.width=8 ---------------
p <- ggplot(results, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line() + 
  labs(x = "wavelength /nm", y = expression(sigma/nm^2), 
       colour = "variable")

p
