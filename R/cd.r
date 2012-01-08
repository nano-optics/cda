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
circular_dichroism_spectrum <- function(cluster, material, n=1.33, N=100, averaging = c("QMC","GL","grid"),
                                         full=TRUE, progress=FALSE, result.matrix=FALSE){

  averaging <- match.arg(averaging)

  k0 <- 2*pi/material$wavelength
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
    return(cbind(wavelength = material$wavelength, extinction=res[, 1], absorption=res[, 2], 
                 CDext=res[,3], CDabs=res[,4]))
  } else {

    d <- data.frame(wavelength = material$wavelength,
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
