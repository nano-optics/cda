## ----start, message=FALSE, echo=FALSE, results='hide'--------------------
knitr::read_demo("dispersion", package="cda")

## ----setup, results='hide', echo=FALSE-----------------------------------
library(cda)
library(rgl)
library(dielectric)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))

## ----demo, results='hide'------------------------------------------------
# dielectric function
wvl <- seq(400, 800, length=200)
gold <- epsAu(wvl)

# cluster, and single particle for reference
cl <- cluster_dimer(dihedral = 0,  a=35, b=20)
clr <- cluster_single(a=35, b=20)

nn <- 7
Incidence <- rep(seq(0,pi/2,length=nn), 3)
Axes <- rep(c("x", "y", "z"), each=nn)

d <- spectrum_dispersion(cl, gold, medium=1.33, 
                                     polarisation = "linear",
                                     Incidence=Incidence, Axes=Axes,
                                     method = 'solve')

id <- c(1, 1+nn, 1+2*nn)
dr <- spectrum_dispersion(clr, gold, medium=1.33, 
                                    polarisation = "linear",
                                    Incidence=Incidence[id], Axes=Axes[id],
                                    method = 'solve')

## ----plot, echo=FALSE----------------------------------------------------
dd <- subset(d, type=="cross-section" & variable == "extinction")
ddr <- subset(dr, type=="cross-section" & variable == "extinction")

p <- 
ggplot(dd, aes(wavelength, value, colour=factor(Incidence*180/pi))) + 
  facet_grid(polarisation~Axes, scales="free") +
  geom_line() +
  geom_line(data=ddr, lty=2, col="black") +
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="angle")

p

