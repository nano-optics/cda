## ----load,message=FALSE--------------------------------------------------
library(cd)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="basic-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a cluster of particles
cl <- list(positions = cbind(c(0, 0, 0),
                      c(0, 0, -50)),
            angles = cbind(c(0, 0, 0),
                           c(pi/4, 0, 0)),
            sizes = cbind(c(30, 10, 10),
                          c(30, 10, 10)))

# visualise
visualise_cluster(cl)


## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

linear <- spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 100)
linear2 <- spectrum_dispersion_ellipsoids(cl, gold, method="ls")
d1 <- subset(linear, type=="cross-section")
d2 <- subset(linear2, type=="cross-section")

ggplot(d1, aes(wavelength, value, colour=variable)) +
  facet_wrap(~polarisation) + geom_line() + geom_line(data=d2, lty=2)


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
circular <- spectrum_oa_ellipsoids(cl, gold, method="oos")

ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line()



