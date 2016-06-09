## testing largeish array with OS formulation

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
cl <- cluster_chain(1000, pitch = 500, a = 50, b = 30, c = 30)
  

# visualise
# visualise_cluster(cl)


## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

library(microbenchmark)

# wvl <- c(400, 900)
# gold <- epsAu(wvl)
# microbenchmark(linear = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 3),
# linear2 = spectrum_dispersion_ellipsoids(cl, gold, method="ls"),
# times = 10L)

wvl <- c(400, 900)
gold <- epsAu(wvl)
microbenchmark(linear = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 5),
               linear2 = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 10),
               times = 1L)


# 
# 
# wvl <- seq(550, 750, length=80)
# gold <- epsAu(wvl)
# linear = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 6)
# linear2 = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 8)
# # linear2 = spectrum_dispersion_ellipsoids(cl, gold, method="ls")
# d1 <- subset(linear, type=="cross-section")
# d2 <- subset(linear2, type=="cross-section")
# 
# ggplot(d1, aes(wavelength, value, colour=variable)) +
#   facet_wrap(~polarisation) + geom_line() + 
#   geom_line(data=d2, lty=2) +
#   theme()
