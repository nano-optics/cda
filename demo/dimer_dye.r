library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))

wavelength = seq(450,650)
medium <- 1.33

dye = alpha_bare(wavelength)

model_dimer <- function(d, orientation="side", ...){
  if(orientation=="side-by-side")
  cl =  cluster_dimer(d, 1, 0, 0, 0, 0, 0) else
    cl =  cluster_dimer(d, 0, 0, 1, 0, 0, 0)
    
  xsec = spectrum_oa(cl, dye, ..., quadrature = 'qmc', Nq = 100, method = "solve")
}

dimers <- mdply(expand.grid(d=c(seq(0.7, 2, by=0.2)), 
                            orientation=c("head-to-head", "side-by-side"), 
                            stringsAsFactors = FALSE), 
                model_dimer, medium=medium)

cl <- cluster_single(1,0,0)
ref <- spectrum_oa(cl, dye, medium=medium, quadrature = 'qmc', Nq = 100, method = "solve")


p <- ggplot(subset(dimers, type=="cross-section" & variable =="extinction"), 
            aes(wavelength, value)) +
  geom_line(aes(colour=factor(d))) +
  facet_grid(.~orientation)+
  # geom_line(data=manual) +
  geom_line(data=subset(ref, type=="cross-section" & variable =="extinction"),
            map=aes(wavelength, value),lty=2) +
  labs(x = "wavelength /nm", y = expression(sigma[abs]/nm^2), 
       colour = "d /nm") 

p

