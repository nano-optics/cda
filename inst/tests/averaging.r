
library(cda)

gold <- cda::epsAu(seq(400, 800))

ellipsoid <- function(a0=0.05, ar=2.7){

  V <- 4/3*pi * a0^3
  b <- c <- a0/(ar^(1/3))
  a <- ar*b
  data.frame(a=a, b=b, c=c)
}


single_spectrum <- function( a0=0.01, ar=1.9, n=1.5){
    
  ell <- ellipsoid(a0, ar)
  kn <- 2*pi/gold$wavelength * n
  alpha <- polarizability.ellipsoid(gold$wavelength, gold$epsilon, ell$a, ell$b, ell$c, n)
  rowMeans( 4*pi*kn*Im(alpha))

}


onecluster <- function(N=50, averaging="grid", a0=0.01, ar=1.9, n=1.5){
  r <- cbind(0, 0, 0)
  sizes <- ellipsoid(a0, ar)
  angles <- cbind(0, 0, 0)
  clust <- list(r = r, sizes = sizes, angles = angles)
  
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=N, averaging=averaging, progress=FALSE)
  invisible(subset(m, variable == "extinction" & type == "cross section"))
  
}

params <- expand.grid(N=seq(5, 100, by=5),
                      averaging=c("grid", "GL", "QMC"),
                      stringsAsFactors=FALSE)

comparison <- mdply(params, onecluster)

single <- data.frame(wavelength=gold$wavelength, value=single_spectrum())

p <- 
ggplot(comparison) + facet_wrap(~averaging)+
  geom_path(aes(wavelength, value, colour=N, group=N))+
  geom_path(aes(wavelength, value), data=single, linetype=2) +
  labs(y=expression(sigma[ext]/mu*m^2),
       x=expression(wavelength/mu*m), colour=expression(N))+
  theme_minimal()

p


## dimer

onecluster2 <- function(N=50, averaging="grid",
                       d=100e-3,scale=1, dihedral=45, alpha1=0, alpha2=0,
                       a0=0.01, ar=1.9, n=1.5){
  
  ell <- ellipsoid(a0, ar)
  
  clust <- makeDimerDihedral(d=scale*d,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=scale*ell$a, b=scale*ell$b)
  
  m <- circular_dichroism_spectrum(clust, gold, n=n,  N=N, averaging=averaging,
                                   progress=FALSE, result.mat=FALSE, full=TRUE)

  invisible(subset(m, type=="cross section"))
}

params2 <- expand.grid(N=seq(5, 100, by=10),
                      averaging=c("grid", "GL", "QMC"),
                      stringsAsFactors=FALSE)

comparison2 <- mdply(params2, onecluster2)
single$variable <- "extinction"

p2 <- 
ggplot(comparison2) + facet_grid(variable~averaging)+
  geom_path(aes(wavelength, value, colour=N, group=N))+
  geom_path(aes(wavelength, value*2), data=single, linetype=2) +
  labs(y=expression(sigma[ext]/mu*m^2),
       x=expression(wavelength/mu*m), colour=expression(N))+
  theme_minimal()

p2


