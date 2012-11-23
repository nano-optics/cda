library(cda)
library(ggplot2)

wvl <- seq(400,900)
gold <- epsAu(wvl)

if(interactive() && require(rgl)){ # display RGL window
  cl <- makeRodChain(N=10, pitch=0.5)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
}

onechain <- function(N, n=1.5){
  clust <- makeRodChain(N=N, pitch=0.5)
  m <- linear_extinction_spectrum(clust, gold, n=n, progress=FALSE)
  invisible(m)
}

params <- data.frame(N=c(1, 10, 50))
comparison <- mdply(params, onechain, .progress="none")
## str(comparison)

p <- 
  ggplot(data=comparison)+
labs(y=expression(sigma[ext]*" /"*mu*m^2),
       x=expression(wavelength*" /"*mu*m),
       colour = expression(N), linetype=expression(polarisation))+
  geom_path(aes(wavelength, value, linetype=variable,
                colour=factor(N),
                group=interaction(N,variable)))+ theme_bw()

p
