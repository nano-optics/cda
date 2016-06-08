
## ----load,message=FALSE, echo=1:6----------------------------------------
require(cda)
require(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, tidy=FALSE, fig.path='averaging-'--------------------------
gold <- epsAu(seq(600, 800))

cl <- cluster_dimer(d=100, 
              dihedral=10*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)

# achiral cluster (plane of symmetry)
cl2 <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)


## ----comparison, tidy=FALSE, fig.path='averaging-'-----------------------
params <- expand.grid(Nq=c(5, 25, 50, 75, 100, 500),
                      quadrature=c("gl", "qmc", "random"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, spectrum_oa, cluster=cl, material=gold)
cheap <- spectrum_oa(cluster=cl, material=gold, quadrature="cheap")
converged <- spectrum_oa(cluster=cl, material=gold, quadrature="qmc", Nq=1000)

p <- 
  ggplot(subset(comparison, type == "dichroism" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_grid(quadrature~., scales="free")+
  geom_path(aes(colour=factor(Nq), group=Nq))+
  geom_path(data=subset(cheap, type == "dichroism" & variable == "extinction"), linetype=2)+
  geom_path(data=subset(converged, type == "dichroism" & variable == "extinction"), linetype=3)+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p

