library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


# dielectric function
gold <- epsAu(seq(400, 700))

# define a cluster of particles

dimer_model <- function(d, ...){
  cl <- cluster_dimer(d = d, a=30, b=30, c=30, dihedral = 0)
  spectrum_dispersion(cl, Incidence = pi/2, Axes="x", ...)
}

## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

single <- cluster_single(30,30,30)
ref <- spectrum_dispersion(single, Incidence = pi/2, Axes="y",  material = gold)
linear <- mdply(data.frame(d=c(seq(50, 100,by=10))), dimer_model, material = gold)
d1 <- subset(linear, type=="cross-section" & variable =="extinction")
p <- ggplot(d1, aes(wavelength, value)) +
  facet_grid(polarisation~., scales="free") + 
  geom_line(aes(colour=factor(d)))+ 
  geom_line(data=subset(ref, type=="cross-section" & variable =="extinction"), lty=2)+ 
  labs(x = "wavelength /nm", y = expression(sigma[ext]/nm^2), 
       colour = "d /nm") +
  # theme(legend.position=c(0.8,0.2)) +
  theme()


p

