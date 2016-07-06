library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


# dielectric function
gold <- epsAu(seq(500, 800))

# define a cluster of particles

dimer_model <- function(d, orientation = c("head-to-tail", "side-by-side"), ...){
  if(orientation == "head-to-tail") {
    cl <- cluster_dimer(d = d, a=50, b=20, c=20, dihedral = 0)
    res <- spectrum_dispersion(cl, Incidence = pi/2, Axes="x", ...)
    return(res[res$polarisation == "p", ])
  }
  if(orientation == "side-by-side"){
    cl <- cluster_dimer(d = d, a=20, b=20, c=50, dihedral = 0)
    res <- spectrum_dispersion(cl, Incidence = pi/2, Axes="x", ...)
    return(res[res$polarisation == "s", ])
  }
}

## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

single <- cluster_single(50,20,20)
ref <- spectrum_dispersion(single, Incidence = pi/2, Axes="x",  material = gold)
ref <- subset(ref, polarisation=="p")
ref$polarisation <- NULL
linear <- mdply(expand.grid(d=c(seq(100, 500,by=50)), orientation = c("head-to-tail", "side-by-side"), stringsAsFactors = FALSE),
                dimer_model, material = gold)
d1 <- subset(linear, type=="cross-section" & variable =="extinction")

p1 <- ggplot(d1, aes(wavelength, value)) +
  facet_grid(polarisation~., scales="free") + 
  geom_line(aes(colour=factor(d), group=d))+ 
  scale_x_continuous(expand=c(0,0)) +
  scale_color_hue(h=c(0,230)+15)+
  geom_line(data=subset(ref, type=="cross-section" & variable =="extinction"), lty=2,lwd=1)+ 
  labs(x = "Wavelength /nm", y = expression(sigma[ext]/nm^2), 
       colour = "d /nm")+
  expand_limits(y=0)+
  theme_bw(12) + theme(strip.background=element_blank(), 
                       strip.text=element_blank(),
                       legend.position="right",
                       panel.margin=unit(2,"mm"), panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank())
p1