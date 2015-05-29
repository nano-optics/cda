
## ----load,message=FALSE, echo=1:6----------------------------------------
require(cda)
require(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, tidy=FALSE, fig.path='averaging-'--------------------------
gold <- epsAu(seq(400, 800, length=200))

cl <- cluster_helix(N=10, R0=500, pitch=1000, 
                     delta=pi/7, delta0=0, right=TRUE,
                     a=100, b=50, c=50,
                     angles='helix')


## ----comparison, tidy=FALSE, fig.path='averaging-'-----------------------
params <- expand.grid(Nquad=100,
                      nmax = c(5, 10, 15),
                      tol=c(1e-4, 1e-5),
                      stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, 
                    cluster=cl, material=gold, averaging="QMC", cg=TRUE)

cheap <- circular_dichroism_spectrum(cluster=cl, 
                                     material=gold, averaging="QMC",Nquad=100, cg=FALSE)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_wrap(~tol)+
  geom_path(aes(colour=factor(nmax), group=nmax))+
  geom_path(data=subset(cheap, type == "CD" & variable == "extinction"), linetype=2)+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p

