library(cda)
library(ggplot2)

wvl <- seq(500,900)
gold <- epsAu(wvl)

if(interactive() && require(rgl)){ # display RGL window
  cl <- makeDimerCluster(d = 200e-3, phi=pi/4)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
}

## 
## Functions for circular dichroism (CD)
## 

##' Simulate a CD spectrum
##'
##' CD spectrum
##' @title circular_dichroism_spectrum 
##' @param cluster cluster (list)
##' @param material material
##' @param n refractive index medium
##' @param N number of integration points
##' @param averaging averaging method, using either Quasi Monte Carlo, Gauss Legendre, or regular grid
##' @param full logical use full (retarded) dipolar field
##' @param progress print progress lines
##' @param result.matrix logical return the results as a matrix
##' @export
##' @family user_level circular_dichroism
##' @author baptiste Auguie
circular_dichroism_spectrum <- function(cluster, material, n=1.33, N=100, averaging = c("QMC","GL","grid"),
                                         full=TRUE, progress=FALSE, result.matrix=FALSE){

  averaging <- match.arg(averaging)

  wavelength <- material[["wavelength"]]
  k0 <- 2*pi/wavelength
  kn <- k0*n
  
  invalpha <- make.invalpha(cluster, material, polarizability.fun=polarizability.ellipsoid, n=n, kuwata=TRUE)
  
  if(averaging == "QMC") # Quasi Monte Carlo, using Halton sequence from randtoolbox
    {
      nodes <- halton(N, dim = 2, normal=FALSE)
  
      res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, as.matrix(nodes),
                                             as.integer(full), as.integer(progress))
  
    }
  
  if(averaging == "GL") # Gauss Legendre quadrature, using nodes and weights from statmod
    {
      ## divide phi[0-360] in 2*sqrt(N/2), and psi[-90-90] in sqrt(N/2) -> N points
      rndN <- ceiling(sqrt(N/2))
      GL <- gauss.quad(2*rndN)
      GL2 <- gauss.quad(rndN)
      
      invalpha <- make.invalpha(cluster, material, polarizability.fun=polarizability.ellipsoid, n=n, kuwata=TRUE)
      res <- cd$circular_dichroism_spectrum(kn, invalpha, cluster$r, cluster$angles,
                                            as.matrix(cbind(GL$nodes, GL$weights)), as.matrix(cbind(GL2$nodes, GL2$weights)),
                                            as.integer(full), as.integer(progress))
      
    }
  if(averaging == "grid") # regular grid
    {
      a <- seq(0 + 0.01/sqrt(N),1-0.01/sqrt(N), length=round(sqrt(N))) # remove end points that cause problems
      nodes <- expand.grid(x = a, y=a)
      
      res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, as.matrix(nodes),
                                             as.integer(full), as.integer(progress))
      
    }
  
  if(result.matrix){
    ## extinction, absorption, CD ext, CD abs
    return(cbind(wavelength = wavelength, extinction=res[, 1], absorption=res[, 2], 
                 CDext=res[,3], CDabs=res[,4]))
  } else {

    d <- data.frame(wavelength = wavelength,
                    extinction=res[, 1], absorption=res[, 2], scattering=res[, 1]-res[, 2],
                    CDext=res[,3], CDabs=res[,4], CDsca=res[,3]-res[,4])

    L2eV <- 6.62606896e-34 * 299792458/1.602176487e-19
    m <- melt(transform(d, energy = L2eV / wavelength * 1e6), id=c("wavelength", "energy"))
  
    m$type <- m$variable

    ## alternative using car::recode
    ## m$type <- recode(m$variable, "c('extinction', 'absorption', 'scattering')='cross-section'; else='CD'")
    ## m$variable <- recode(m$variable, "c('extinction', 'CDext')='extinction';c('absorption', 'CDabs')='absorption';else='scattering'")
    
    levels(m$type) <- list(CD="CDext",CD="CDabs",CD="CDsca",
                           `cross section`="extinction",
                           `cross section`="absorption",
                           `cross section`="scattering")
    
    levels(m$variable) <- list(extinction="extinction",
                               absorption="absorption",
                               scattering="scattering",
                               extinction="CDext",
                               absorption="CDabs",
                               scattering="CDsca")
    
    return(m)
    
  }
  
}


onecluster <- function(handedness="right"){
  test <- handedness == "right"
  clust <- makeDimerCluster(d = 200e-3, phi=pi/4, right = test)
  m <- circular_dichroism_spectrum(clust, gold, n=1.33, N=36, progress=FALSE)
  invisible(m)
}

test <- onecluster()
params <- data.frame(handedness=c("right", "left"))
comparison <- mdply(params, onecluster, .progress="none")
## str(comparison)

p <- ggplot(comparison) + facet_grid(type~., scales="free") +
  geom_hline(yintercept=0) +
  geom_path(aes(energy, value, colour=handedness, linetype=variable,
                group=interaction(handedness, variable))) +
  labs(x="energy / eV", y=expression(sigma / mu*m^2))

p

