## ----load,message=FALSE--------------------------------------------------
library(cd)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="basic-"----

# dielectric function
wvl <- seq(400, 900)
#wvl <- seq(400, 402)
gold <- epsAu(wvl)

# define a cluster of particles

cl <- list(positions = cbind(c(0, 0, 0),
                      c(0, 0, -45)),
            angles = cbind(c(0, 0, 0),
                           c(pi/4, 0, 0)),
            sizes = cbind(c(30, 10, 10),
                          c(30, 10, 10)))
# cl <- cluster_helix(N=5, a=10, b=10, c=30, R0=100, pitch = 300, delta = pi/4)
# visualise
 # visualise_cluster(cl)


## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------
# 
# linear <- spectrum_dispersion(cl, gold)
# ggplot(linear, aes(wavelength, value, linetype=type)) +
#   facet_wrap(~polarisation) + geom_path()

Incidence <- rep(seq(0,pi,length=20), 3)
Axes <- rep(c("x", "y", "z"), each=20)
## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
c1 <- spectrum_dispersion_ellipsoids(cl, gold, medium=1.33, 
                                     polarisation = "circular",
                                     # polarisation = "linear",
                                     Incidence=Incidence, Axes=Axes,
                                     method = 'ls')


ggplot(c1, aes(wavelength, value, colour=Incidence,group=Incidence)) + 
  facet_grid(variable~Axes+polarisation, scales="free") +
  geom_line()

