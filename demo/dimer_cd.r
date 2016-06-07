library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


# dielectric function
gold <- epsAu(seq(400, 700))


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
cl <- cluster_dimer(d = 50, a=30, b=10, c=10, dihedral = pi/4)
gold <- epsAu(seq(500, 800))
circular <- spectrum_oa(cl, gold, Nsca=10)

p <- ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line() + 
  labs(x = "wavelength /nm", y = expression(sigma/nm^2), 
       colour = "variable")

p
