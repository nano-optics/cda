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
gold <- epsAg(wvl)

cluster_array <- function(N, pitch=500, d=100, a=30, b=10, c=b){
  
  xyz1 <- expand.grid(x=seq_len(N)* pitch, y=seq_len(N)* pitch, z=0)
  xyz2 <- expand.grid(x=seq_len(N)* pitch - d, y= seq_len(N)* pitch, z=0)
  positions <- cbind(t(as.matrix(xyz1)), t(as.matrix(xyz2)))
  
  sizes <- cd:::equal_sizes(a=a, b=b, c=c, N=2*N^2)
  angles <- cd:::equal_angles(0,0,0, N=2*N^2)
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = angles),
            class="cluster")
  
}

# define a cluster of particles
cl <- cluster_array(15, pitch = 500, d=160, a = 60, b = 60, c = 15)
  

# visualise
# visualise_cluster(cl)


## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------



wvl <- seq(500, 1000, length=80)
silver <- epsAg(wvl)
# linear = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 6, medium=1.515, progress = TRUE)
# linear2 = spectrum_dispersion_ellipsoids(cl, gold, method="oos", maxiter = 8, medium=1.515, progress = TRUE)
linear2 = spectrum_dispersion_ellipsoids(cl, silver, method="ls", medium=1.515, progress = TRUE)
d1 <- subset(linear, type=="cross-section")
d2 <- subset(linear2, type=="cross-section")

ggplot(d1, aes(wavelength, value, colour=variable)) +
  facet_wrap(~polarisation) + 
  # geom_line() + 
  geom_line(data=d2, lty=2) +
  theme()
