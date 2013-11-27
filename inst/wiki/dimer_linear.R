
## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)
library(knitr)


## ----setup,echo=FALSE----------------------------------------------------
knit_hooks$set(rgl = function(before, options, envir) {
  # if a device was opened before this chunk, close it
  if (before && rgl.cur() > 0) rgl.close()
  hook_rgl(before, options, envir)
})
rgl_annotate = function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
}
theme_set(theme_minimal())


## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="dimerlinear-"----

# dielectric function
wvl <- seq(500, 800)
gold <- epsAu(wvl)



## ----comparison,echo=TRUE,tidy=FALSE,fig.path="dimerlinear-",fig.width=8----
  
dimer <- function(d=100, ...){
  
  r <- cbind(c(0,0), c(-d/2, d/2), c(0, 0))
  sizes <- equal_sizes(a=50, b=20, c=20, N=2)  
  angles <- matrix(0, ncol=3, nrow=2)

  cl <- list(r=r, sizes=sizes, angles=angles)
  dispersion_spectrum(cl, material = gold)
  
}
params <- data.frame(d=seq(50, 500, by=10))
comparison <- mdply(params, dimer)

## compare with the single-particle response
single <- dispersion_spectrum(list(r=cbind(0,0,0), angles = cbind(0,0,0),
                                   sizes=cbind(a=50, b=20, c=20)), material = gold)

p <- 
  ggplot(data=subset(comparison, polarisation == "p")) + 
  facet_grid(type~., scales="free") +
  geom_line(aes(wavelength, value/2, 
                colour=d, group=d)) +
  geom_line(aes(wavelength, value), colour="red", linetype="dashed",
            data=subset(single, polarisation == "p")) +
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="distance /nm") 

p


