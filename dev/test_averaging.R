
## ----load,message=FALSE, echo=1:6----------------------------------------
require(cd)
require(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, tidy=FALSE, fig.path='averaging-'--------------------------
gold <- epsAu(seq(600, 800))

cl <- cluster_dimer(d=100)

# achiral cluster (plane of symmetry)
cl2 <- cluster_dimer(d=100, dihedral=0, alpha2 = pi/4)


## ----comparison, tidy=FALSE, fig.path='averaging-'-----------------------
params <- expand.grid(Nq=c(5, 10, 50, 100),
                       quadrature=c("random", "gl", "qmc"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, spectrum_oa, cluster=cl, material=gold)
cheap <- spectrum_oa(cluster=cl, material=gold, quadrature="cheap")
converged <- spectrum_oa(cluster=cl, material=gold, quadrature="qmc", Nq=500)

p <- 
  ggplot(subset(comparison, type == "dichroism" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_grid(quadrature~., scales="free")+
  geom_line(aes(colour=factor(Nq), group=Nq))+
  geom_line(data=subset(cheap, type == "dichroism" & variable == "extinction"), linetype=2)+
  geom_line(data=subset(converged, type == "dichroism" & variable == "extinction"))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


## ----achiral, tidy=FALSE, fig.path='averaging-'--------------------------
params <- expand.grid(Nq=c(50, 100, 500),
                      quadrature=c("random", "gl", "qmc"),
                      stringsAsFactors=FALSE)

comparison <- mdply(params, spectrum_oa, cluster=cl2, material=gold)
cheap <- spectrum_oa(cluster=cl2, material=gold, quadrature="cheap")
converged <- spectrum_oa(cluster=cl2, material=gold, quadrature="gl", Nq=500)

p <- 
  ggplot(subset(comparison, type == "dichroism" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_grid(quadrature~., scales="free")+
  geom_line(aes(colour=factor(Nq), group=Nq))+
  geom_line(data=subset(cheap, type == "dichroism" & variable == "extinction"), linetype=2)+
  geom_line(data=subset(converged, type == "dichroism" & variable == "extinction"))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


