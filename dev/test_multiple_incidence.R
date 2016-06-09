
## ----load,message=FALSE--------------------------------------------------
library(cd)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="multiple-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

cl <- cluster_dimer(d=100)
cl2 <- cluster_dimer(d=100)
## ----cd,echo=TRUE,tidy=FALSE,fig.path="multiple-",fig.width=8------------

Incidence <- rep(seq(0, pi/2, length=12), 3)
Axes <- rep(c('x','y','z'), each=12)

results <- spectrum_dispersion_ellipsoids(cl, gold, Incidence=Incidence, Axes = Axes, 
                               polarisation="linear")

test <- melt(results, meas="value")

ggplot(subset(test, type == "extinction"), 
       aes(wavelength, value, colour=Incidence, group=Incidence)) +
  facet_grid(Axes ~ polarisation, scales="free") +
  geom_line() +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="incident angle")


## ----comparison,echo=TRUE,tidy=FALSE,fig.path="multiple-",fig.width=8----

params <- expand.grid(Incidence = seq(0, 2*pi, length=36),
                         Axes = c('x','y','z'))

results <- spectrum_dispersion_ellipsoids(cl2, gold, Incidence=params$Incidence, Axes = params$Axes, 
                                          polarisation="circular")

average <- spectrum_oa_ellipsoids(cl2, material = gold)

ggplot(subset(results, polarisation == "CD"), aes(wavelength, value)) +
  facet_grid(Axes ~ polarisation, scales="free") +
  geom_line(aes(colour=Incidence, group=Incidence)) +
  geom_line(data=subset(average, type == "CD" & variable =="extinction"), linetype=2, size=1.2) +
  labs(y=expression(sigma[CD]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="incident angle")


