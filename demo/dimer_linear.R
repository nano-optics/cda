
## @knitr load
library(cda)
library(rgl)
library(ggplot2)


## @knitr setup

rgl_annotate = function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
}
theme_set(theme_minimal())


## @knitr cluster

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)



## @knitr comparison
  
dimer <- function(d=100, ...){
  
  r <- cbind(c(0,0), c(-d/2, d/2), c(0, 0))
  sizes <- equal_sizes(a=50, b=20, c=20, N=2)  
  angles <- matrix(0, ncol=3, nrow=2)

  cl <- list(r=r, sizes=sizes, angles=angles)
  linear_extinction_spectrum(cl, material = gold)
  
}
params <- data.frame(d=seq(100, 500, length=50))
comparison <- mdply(params, dimer)

p <- 
  ggplot(data=comparison) + 
  facet_grid(variable~., scales="free") +
  geom_line(aes(wavelength, value, 
                colour=d, group=d)) +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="distance /nm") 

p


