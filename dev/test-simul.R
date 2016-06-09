library(cd)
library(rgl)

wavelength <- seq(300, 800)
test <- polarizability_dye(wavelength, a=1, b=0.5, c=0.2)

test2 <- polarizability_ellipsoid(wavelength, epsAg(wavelength)$epsilon, 
                                  a = 0.8, b=0.8, c=0.8)

matplot(wavelength, Im(test2), t="l", col=1)
matlines(wavelength, Im(test), t="l", col=2)



wavelength <- seq(450, 600, length=50)
N <- 50
cl <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=1, c=1, 
                    position = "quasi-random", orientation = "radial")

cl1 <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=0, c=0, 
                    position = "quasi-random", orientation = "radial")

cl2 <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=0, c=0, 
                     position = "quasi-random", orientation = "flat")


# # cl1 <- cl
# rgl.spheres(0,0,0,15)
# rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="red")

simulation <- function(N=10, orientation="radial", 
                       dipole = c('uniaxial','isotropic','planar'), 
                       core=FALSE, 
                       medium=1.33, R0=20, ...){
  dipole <- match.arg(dipole)
  if(dipole=="uniaxial"){
    a <- 1; b <- 0
  } else if(dipole=="isotropic"){
    a <- 1; b <- 1
  } else if(dipole=="planar"){
    a <- 0.1; b <- 1
  } 
  
  cl <- cluster_shell(N = N, d=1, R0 = R0,  a=a, b=b, c=b, 
                      position = "quasi-random",
                      orientation = orientation, ...)
  
  res <- spectrum_shell(cluster = cl, core=core, wavelength = wavelength, 
                         method = "ls", cg=FALSE, iterative=FALSE,
                         averaging = "QMC", Nquad = 50,
                         medium=medium, progress = FALSE)
  
  res
}
N=c(1, 10, 30, 50, 70)

params <- expand.grid(N=c(10, 50, 100, 150, 200), 
                      orientation=c("radial", "flat"),
                      dipole = c('uniaxial','isotropic','planar'), 
                      # core=c(TRUE, FALSE),
                      stringsAsFactors = FALSE)

library(ggplot2)
library(plyr)
all <- mdply(params, simulation, .progress = "text",core=TRUE, seed=234)


ggplot(subset(all, variable == "absorption"), 
       aes(wavelength, value/N, colour=factor(N), 
                group=interaction(N,dipole))) + 
  geom_line() + 
  facet_wrap(~dipole+orientation,scales="free",ncol=3) +
  theme()

# ggsave("coreshell.pdf",width=10,height=4)
