library(cd)
library(rgl)

wavelength <- seq(300, 800)
test <- polarizability_dye(wavelength, a=1, b=0.5, c=0.2)

test2 <- polarizability_ellipsoid(wavelength, epsAg(wavelength)$epsilon, 
                                  a = 0.8, b=0.8, c=0.8)

matplot(wavelength, Im(test2), t="l", col=1)
matlines(wavelength, Im(test), t="l", col=2)

cl <- cluster_dimer(d=10, a=1, b=1)
res <- spectrum_shell(cluster = cl, core=FALSE, 
                      wavelength = seq(300,600,length=100), 
                      method = "ls", cg=FALSE, iterative=FALSE,
                      averaging = "QMC", Nquad = 50,
                      medium=1.33, progress = FALSE)
# 
# ggplot(res,aes(wavelength, value, colour=variable))+geom_line()

# rgl.spheres(0,0,0,15)
# rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="red")

simulation <- function(d=10,  
                       dipole = c('uniaxial','isotropic','planar'), 
                       medium=1.33, 
                       wavelength=seq(300,600,length=100),...){
  dipole <- match.arg(dipole)
  if(dipole=="uniaxial"){
    a <- 1; b <- 0
  } else if(dipole=="isotropic"){
    a <- 1; b <- 1
  } else if(dipole=="planar"){
    a <- 0.1; b <- 1
  } 
  
  cl <- cluster_dimer(d=d, a=a, b=b, ...)
  
  res <- spectrum_shell(cluster = cl, core=FALSE, wavelength = wavelength, 
                         method = "ls", cg=FALSE, iterative=FALSE,
                         averaging = "QMC", Nquad = 50,
                         medium=medium, progress = FALSE)
  
  res
}

params <- expand.grid(d=c(2, 3,5, 10, 100), 
                      dipole = c('uniaxial','isotropic','planar'), 
                      stringsAsFactors = FALSE)

library(ggplot2)
library(plyr)
all <- mdply(params, simulation, .progress = "text")


ggplot(subset(all,variable=="absorption"), 
       aes(wavelength, value, colour=factor(d), 
                group=d)) + 
  geom_line() + 
  facet_wrap(~dipole,scales="free",ncol=3) +
  theme()

# ggsave("coreshell.pdf",width=10,height=4)
