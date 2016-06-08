library(cda)
library(rgl)
library(ggplot2)
library(dielectric)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))

# dielectric function
wvl <- seq(450, 750, length=200)
gold <- epsAu(wvl)

# two clusters
# cl <- cluster_dimer(dihedral = pi/4,  a=35, b=20)
cl <- cluster_helix(5, R0=20, pitch=30, 
                    delta=pi/2, delta0=0, right=TRUE,
                    a=15/2, b=15/2, c=15,
                    angles="helix")

nn <- 19
Incidence <- rep(seq(0,pi,length=nn), 3)
Incidence * 180/pi
Axes <- rep(c("x", "y", "z"), each=nn)

## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
d <- spectrum_dispersion(cl, gold, medium=1.33, 
                                     polarisation = "circular",
                                     Incidence=Incidence, Axes=Axes,
                                     method = 'solve')

dr <- spectrum_oa(cl, gold, medium=1.33, quadrature = "qmc",
                             Nq = 300, method = 'solve')

dd <- subset(d, type=="dichroism" & !(variable == "absorption"))
ddr <- subset(dr, type=="dichroism" & !(variable == "absorption"))

p <- 
ggplot(dd, aes(wavelength, value)) + 
  facet_grid(variable~Axes, scales="free") +
  geom_line(aes(colour=factor(Incidence*180/pi))) +
  geom_line(data=ddr, lty=2, col="black",lwd=1) +
  guides(colour=guide_legend(ncol=2))+
  labs(y=expression(sigma[CD]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="incident angle")

p
 

