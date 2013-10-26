## 
## Functions for circular dichroism (CD)
## 

##' Simulate a CD spectrum
##'
##' CD spectrum
##' @title circular_dichroism_spectrum 
##' @param cluster cluster (list)
##' @param material material
##' @param medium refractive index medium
##' @param N number of integration points
##' @param averaging averaging method, using either Gauss Legendre quadrature (default), Quasi Monte Carlo, regular grid, or "cheap" (3 axes)
##' @param iterative logical, increase N until convergence (QMC only)
##' @param precision relative diff between two runs (QMC only)
##' @param Nmax maximum N if convergence not attained (QMC only)
##' @param full logical use full (retarded) dipolar field
##' @param progress print progress lines
##' @param result.matrix logical return the results as a matrix
##' @export
##' @family user_level circular_dichroism
##' @author baptiste Auguie
circular_dichroism_spectrum <- function(cluster, material, medium=1.33, N=100, 
                                        averaging = c("GL","QMC","grid", "cheap"),
                                        iterative=FALSE, precision=1e-3, Nmax=1e4,
                                        full=TRUE, progress=FALSE, 
                                        result.matrix=FALSE){

  averaging <- match.arg(averaging)

  wavelength <- material[["wavelength"]]
  k0 <- 2*pi/wavelength
  kn <- k0*medium
  
  invalpha <- inverse_polarizability(cluster, material, 
                            polarizability_fun=polarizability_ellipsoid, 
                            medium=medium, kuwata=TRUE)
  
  if(averaging == "QMC") # Quasi Monte Carlo, using Halton sequence from randtoolbox
    {
      nodes <- halton(N, dim = 2, normal=FALSE, init=TRUE)
      res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, 
                                             as.matrix(nodes),
                                             as.integer(full), as.integer(progress))
      
      ## iterative improvement: add new points until convergence or Nmax reached
      if(iterative){
        converged <- FALSE
        Ntot <- N
        while(Ntot < Nmax && !converged){
          oldN <- Ntot
          old <- res[,1]
          Ntot <- Ntot + 100
          nodes <- halton(100, dim = 2, normal=FALSE, init=FALSE)
          ## xsec at new points
          newres <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, 
                                                    as.matrix(nodes),
                                                    as.integer(full), as.integer(progress))
          
          ## average of the two results
          res <- (oldN * res + 100 * newres) / (oldN + 100)
          
          test <- max(abs(old - res[,1]) / res[,1]) # max relative difference in extinction cross section
          
          message("N:", Ntot, "; relative error: " , test)
          converged <- test < precision
        }
      }
    }
  
  if(averaging == "GL") # Gauss Legendre quadrature, using nodes and weights from statmod
    {
      ## divide phi[0-360] in 2*sqrt(N/2), and psi[-90-90] in sqrt(N/2) -> N points
      rndN <- ceiling(sqrt(N/2))
      GL <- gauss.quad(2*rndN)
      GL2 <- gauss.quad(rndN)
      
      res <- cd$circular_dichroism_spectrum(kn, invalpha, cluster$r, cluster$angles,
                                            as.matrix(cbind(GL$nodes, GL$weights)), 
                                            as.matrix(cbind(GL2$nodes, GL2$weights)),
                                            as.integer(full), as.integer(progress))
      
    }
  if(averaging == "grid") # regular grid
    {
      a <- seq(0 + 0.01/sqrt(N),1-0.01/sqrt(N), length=round(sqrt(N))) # remove end points that cause problems
      nodes <- expand.grid(x = a, y=a)
      
      res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, 
                                             as.matrix(nodes),
                                             as.integer(full), as.integer(progress))
      
    }
 
  if(averaging == "cheap") # 3 directions
  {
    
#     6 directions
#     nodes <- rbind(c(1/2, 0), # +x is phi=0, psi=0
#                    c(1/2, 1/2), # -x is phi=pi, psi=0
#                    c(1/2, 1/4), # +y is phi=pi/2, psi=0
#                    c(1/2, 3/4), # -y is phi=3pi/2, psi=0
#                    c(1, 1/4), # +z is phi=pi/2, psi=pi/2
#                    c(0, 1/4)) # -z is phi=pi/2, psi=-pi/2
    
    # psi, phi
    nodes <- rbind(c(1/2, 0), # +x is phi=0, psi=0
                   c(1/2, 1/4), # +y is phi=pi/2, psi=0
                   c(1, 1/4)) # +z is phi=pi/2, psi=pi/2
                    
    res <- cd$circular_dichroism_spectrum2(kn, invalpha, cluster$r, cluster$angles, 
                                           as.matrix(nodes),
                                           as.integer(full), as.integer(progress))
    
  }
  
  if(result.matrix){
    ## extinction, absorption, CD ext, CD abs
    return(cbind(wavelength = wavelength, extinction=res[, 1], absorption=res[, 2], 
                 CDext=res[,3], CDabs=res[,4]))
  } else {

    d <- data.frame(wavelength = wavelength,
                    extinction=res[, 1], absorption=res[, 2], scattering=res[, 1]-res[, 2],
                    CDext=res[,3], CDabs=res[,4], CDsca=res[,3]-res[,4]) # L - R

    L2eV <- 6.62606896e-34 * 299792458/1.602176487e-19
    m <- melt(transform(d, energy = L2eV / wavelength * 1e9), id=c("wavelength", "energy"))
  
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

if (getRversion() >= "2.15.1")
 utils::globalVariables("wavelength")

