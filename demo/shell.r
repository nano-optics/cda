## ---- setup,echo=FALSE ---------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


## ---- demo ---------------
wavelength = seq(450,600, length=50)
medium <- 1.33

dye <- alpha_bare(wavelength)

model_shell <- function(rho=1, R0=3, d=0.5, ...){
  
  N <- dye_coverage(rho, R0+d)
  message(N)
  cl =  cluster_shell(N, R0=R0, d=d, 1, 1, 1, orientation = "radial", position = 'fibonacci');
  xsec = spectrum_oa(cl, dye, ..., quadrature = 'cheap',  method = "solve")
}

shells <- mdply(data.frame(rho=c( 0.1, 0.5, 1)), 
                model_shell, medium=medium, .progress='text')

ref <- spectrum_oa(cluster_single(1, 1, 1), dye, medium=medium, 
                   quadrature = 'qmc', Nq = 100, method = "solve")


## ---- plot,echo=TRUE,fig.width=8 ---------------

p <- ggplot(subset(shells, type=="cross-section" & variable =="extinction"), 
            aes(wavelength, value)) +
  geom_line(aes(colour=factor(rho))) +
  geom_line(data=subset(ref, type=="cross-section" & variable =="extinction"),
            map=aes(wavelength, value),lty=2) +
  labs(x = "wavelength /nm", y = expression(sigma[abs]/nm^2), 
       colour = "d /nm") 

p


