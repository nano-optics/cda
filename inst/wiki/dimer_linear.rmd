library(cda)
library(ggplot2)

wvl <- seq(400,900)
gold <- epsAu(wvl)


if(interactive() && require(rgl)){ # display RGL window
  cl <- makeRodChain(N=2, pitch=0.2)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
}

onedimer <- function(d = 200){
  clust <- makeRodChain(N=2, pitch=d*1e-3)
  m <- linear_extinction_spectrum(clust, gold, n=1.5, progress=FALSE)
  invisible(m)
}

params <- data.frame(d=seq(200, 400, by=10))
comparison <- mdply(params, onedimer, .progress="none")
## str(comparison)

p <- 
  ggplot(data=comparison)+
  labs(y=expression(sigma[ext]*" /"*mu*m^2),
       x=expression(wavelength*" /"*mu*m),
       colour = expression(d/nm), linetype=expression(polarisation))+
  geom_path(aes(wavelength, value, linetype=variable,
                colour=d,
                group=interaction(d,variable)))+ theme_bw()

p
