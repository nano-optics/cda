
## ----load,message=FALSE, echo=1:6----------------------------------------
require(cda)
require(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, tidy=FALSE, fig.path='averaging-'--------------------------
gold <- epsAu(seq(400, 800, length=200))
gold2 <- epsAu(seq(400, 800, length=100))

cl <- cluster_helix(N=10, R0=50, pitch=100, 
                     delta=pi/3, delta0=0, right=TRUE,
                     a=20, b=10, c=10,
                     angles='helix')
library(rgl)
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()

params <- expand.grid(Nquad=c(2, 50, 100),
                      nmax = 50,
                      tol=c(1e-10),
                      born = c(TRUE, FALSE),
                      stringsAsFactors=FALSE)

timer <- function(...){
  
  system.time(circular_dichroism_spectrum(cluster=cl, material=gold, 
                                          averaging="QMC", ...))
}
test <- mdply(params, timer)

ggplot(test, aes(Nquad, elapsed, linetype = born)) +
  geom_line() + facet_wrap(~nmax)

## ----comparison, tidy=FALSE, fig.path='averaging-'-----------------------
params <- expand.grid(Nquad=50,
                      nmax = c(2, 5, 10, 20),
                      tol=c(1e-5, 1e-6),
                      born = c(TRUE, FALSE),
                      stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, 
                    cluster=cl, material=gold, averaging="QMC", cg=TRUE)

cheap <- circular_dichroism_spectrum(cluster=cl, 
                                     material=gold, averaging="QMC",Nquad=50, cg=FALSE)

p <- 
  ggplot(subset(comparison, type == "cross section" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_wrap(~tol)+
  geom_path(aes(colour=factor(nmax), linetype=born))+
  geom_path(data=subset(cheap, type == "cross section" & variable == "extinction"), linetype=2)+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))
# 
p

p2 <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_wrap(~tol)+
  geom_path(aes(colour=factor(nmax), linetype=born))+
  geom_path(data=subset(cheap, type == "CD" & variable == "extinction"), linetype=2)+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p2

