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
cl <- cluster_helix(N=5, a=10, b=10, c=30, R0=100, pitch = 300, delta = pi/4)
# visualise
 # visualise_cluster(cl)


## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------
# 
# linear <- spectrum_dispersion(cl, gold)
# ggplot(linear, aes(wavelength, value, linetype=type)) +
#   facet_wrap(~polarisation) + geom_path()


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
c1 <- spectrum_oa_ellipsoids(cl, gold, medium=1.33, quadrature = "gl", Nq=36, method = 'ls')
c2 <- spectrum_oa_ellipsoids(cl, gold, medium=1.33, quadrature = "gl", Nq=36, 
                                   tol=1e-4, method = 'oos', maxiter=50)

ggplot(c1, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line()+ geom_line(data=c2, linetype=2)


