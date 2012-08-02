library(cda)

library(cda)
library(ggplot2)

wvl <- seq(400,900, by=1)
gold <- epsAu(wvl)


onedimer <- function(d = 200){
  clust <- makeRodChain(N=2, pitch=d*1e-3)
  m <- linear_extinction_spectrum(clust, gold, n=1.5, progress=FALSE)
  invisible(m)
}

## params <- data.frame(d=seq(200, 400, by=10))
## comparison <- mdply(params, onedimer, .progress="none")
## str(comparison)

## p <- 
##   ggplot(data=comparison)+
##   labs(y=expression(sigma[ext]*" /"*mu*m^2),
##        x=expression(wavelength*" /"*mu*m),
##        colour = expression(d/nm), linetype=expression(polarisation))+
##   geom_path(aes(wavelength, value, linetype=variable,
##                 colour=d,
##                 group=interaction(d,variable)))+ theme_bw()

## p

 ## dispersion_spectrum(const arma::colvec kn, const arma::cx_mat& Beta, const arma::mat& R, \
 ##        			      const arma::mat& Euler, const arma::mat& Angles, \
 ##        			      const int full, const int progress)



clust <- makeRodChain(N=2, pitch=200e-3)

n <- 1.5
k0 <- 2 * pi/gold$wavelength
kn <- k0 * n
invalpha <- make.invalpha(clust, gold, polarizability.fun = polarizability.ellipsoid, 
                          n = n, kuwata = TRUE)

Angles <- as.matrix(expand.grid(phi=seq(0, pi, length=100), theta=0, psi=0))

## Angles <- rbind(c(pi, 1, 1),
##                 c(0, pi, 1))

## res <- linear$linear_extinction_spectrum(kn, invalpha, 
##                                          clust$r, clust$angles, as.integer(TRUE))



dispersion_spectrum <- function (cluster,  angles, material, n = 1.5, progress = FALSE, full=TRUE) 
{
    k0 <- 2 * pi/material$wavelength
    kn <- k0 * n

    Nwavelengths <- length(k0)
    Nparticles <- nrow(cluster$r)
    Nangles <- nrow(angles)
    
    invalpha <- make.invalpha(cluster, material, polarizability.fun = polarizability.ellipsoid, 
        n = n, kuwata = TRUE)
    
    stopifnot(is.matrix(invalpha),
              is.matrix(angles),
              is.matrix(clust$r), 
              is.matrix(clust$angles))

    stopifnot(ncol(invalpha)/3 == Nparticles,
              nrow(invalpha) == Nwavelengths)

       
    res <- linear$dispersion_spectrum(kn, invalpha, 
                                      clust$r, clust$angles, angles,
                                      as.integer(full), as.integer(progress))
    rbind(
          data.frame(wavelength = rep(material$wavelength, each=Nangles),
                     angles=angles[rep(seq.int(Nangles), Nwavelengths), ],
                     value = c(res[, 1, , drop=TRUE]),
                     type="extinction", polarisation="p"),
          data.frame(wavelength = rep(material$wavelength, each=Nangles),
                     angles=angles[rep(seq.int(Nangles), Nwavelengths), ],
                     value = c(res[, 2, , drop=TRUE]),
                     type="absorption", polarisation="p"),
          data.frame(wavelength = rep(material$wavelength, each=Nangles),
                     angles=angles[rep(seq.int(Nangles), Nwavelengths), ],
                     value = c(res[, 3, , drop=TRUE]),
                     type="extinction", polarisation="s"),
          data.frame(wavelength = rep(material$wavelength, each=Nangles),
                     angles=angles[rep(seq.int(Nangles), Nwavelengths), ],
                     value = c(res[, 4, , drop=TRUE]),
                     type="absorption", polarisation="s")
          )
}


clust <- makeRodChain(N=50, pitch=500e-3)
results <- dispersion_spectrum(clust, Angles, gold, progress=TRUE)

test <- melt(results, meas="value")

ggplot(test, aes(wavelength, value, colour=angles.phi, group=angles.phi)) +
  facet_grid(type ~ polarisation, scales="free") +
  geom_path()
