## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)
library(dielectric)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="basic-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a cluster of particles
cl <- cluster_dimer()
# visualise
visualise(cl)



## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

linear <- cda::spectrum_dispersion(cl, gold)
ggplot(linear, aes(wavelength, value, linetype=type)) +
  facet_wrap(~polarisation) + geom_path()


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
circular <- cda::spectrum_oa(cl, gold)

ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line()



