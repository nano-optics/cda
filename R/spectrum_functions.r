## 
## High-level functions for full spectra simulations
## 


##' Orientation-averaged spectrum
##'
##' OA spectrum
##' @title spectrum_oa
##' @param cluster cluster (list)
##' @param material material
##' @param medium refractive index medium
##' @param Nq number of integration points
##' @param quadrature quadrature method, using either Gauss Legendre quadrature (default), Quasi Monte Carlo, random, or "cheap" (3 Axes)
##' @param iterative logical, increase N until convergence (QMC only)
##' @param precision relative diff between two runs (QMC only)
##' @param Qmax maximum N if convergence not attained (QMC only)
##' @param method: linear system (ls), or order-of-scattering (oos)
##' @param dN iterative increase in N (QMC only)
##' @param full logical use full (retarded) dipolar field
##' @param Nsca quadrature points for scattering cross-section
##' @param cg logical use conjugate gradient solver
##' @param born logical first Born approx as cg guess
##' @param maxiter integer termination of conjugate gradient solver
##' @param tol double, tolerance of conjugate gradient solver
##' @param progress print progress lines
##' @param verbose display messages
##' @param result.matrix logical return the results as a matrix
##' @importFrom reshape2 melt
##' @export
##' @family user_level circular_dichroism
##' @author baptiste Auguie
##' @references
##' Y. Okada, Efficient numerical orientation quadrature of light scattering properties 
##' with a quasi-Monte-Carlo method, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 9, June 2008, Pages 1719-1742.
spectrum_oa <- function(cluster, material, medium=1.33, 
                        quadrature = c("gl","qmc","random", "cheap"), Nq=100, 
                        iterative=FALSE, precision=1e-3, Qmax=1e4, dN=Nq, 
                        method = c("ls", "oos"),
                        Nsca = 50,
                        cg = FALSE, born=FALSE, 
                        maxiter = 30, tol=1e-4,
                        full=TRUE, progress=FALSE, verbose=TRUE,
                        result.matrix=FALSE){
  
  quadrature <- match.arg(quadrature)
  method <- match.arg(method)
  
  ## check whether material parameters correspond to 
  ## epsilon (NPs) pr alpha (dyes)
  isMolecular <- "alpha" %in% names(material)
  
  if(isMolecular){
    alphabar <- material[["alpha"]]
    alphaeff = alpha_embedded(alphabar, medium)
    Alpha = alpha_rescale(alphaeff, cluster[["sizes"]])
    
  } else {
    Alpha <- alpha_ellipsoid(material[["wavelength"]], 
                             material[["epsilon"]], medium, 
                             cluster[["sizes"]])
  }
  
  # incident field
  kn <- 2*pi/material[["wavelength"]]*medium
  
  Incidence <- quadrature_sphere(Nq, quadrature)
  
  
  Scattering <- quadrature_sphere(Nquad=Nsca, "gl")
  
  
  if(method == "ls") {
    
    results <- cpp_oa_spectrum(kn, medium, cluster$positions, Alpha, 
                               cluster$angles, 
                               Incidence$nodes, 
                               Incidence$weights,
                               Scattering$nodes,
                               Scattering$weights,
                               full, cg, born, maxiter, tol, progress)
    
    
  } else if(method == "oos") {
    results <- cpp_oa_spectrum_oos(kn, medium, cluster$positions, Alpha, 
                                   cluster$angles, 
                                   Incidence$nodes, 
                                   Incidence$weights,
                                   Scattering$nodes,
                                   Scattering$weights,
                                   full, maxiter, tol, progress)
  }
  ## iterative improvement: add new points 
  ## until convergence or Qmax reached
  if(iterative && quadrature == "qmc"){
    converged <- FALSE
    Ntot <- Nq
    while(Ntot < Qmax && !converged){
      oldN <- Ntot
      old <- results[,1]
      Ntot <- Ntot + dN
      Incidence <- quadrature_sphere(dN, quadrature, FALSE)
      ## xsec at new points
      
      if(method == "ls") {
        newres <- cpp_oa_spectrum(kn, cluster$positions, Alpha, 
                                  cluster$angles, 
                                  Incidence$nodes, 
                                  Incidence$weights,
                                  Scattering$nodes,
                                  Scattering$weights,
                                  full, cg, born, maxiter, 
                                  tol, progress)
      } else if(method == "oos") {
        newres <- cpp_oa_spectrum_oos(kn, cluster$positions, Alpha, 
                                      cluster$angles, 
                                      Incidence$nodes, 
                                      Incidence$weights,
                                      Scattering$nodes,
                                      Scattering$weights,
                                      full, maxiter, tol, progress)
      }
      
      ## average of the two results
      results <- (oldN * results + dN * newres) / (oldN + dN)
      
      test <- max(abs(old - results[,1]) / results[,1]) 
      ## max relative difference in extinction cross section
      if(verbose)
        message("N:", Ntot, "; relative error: " , test)
      converged <- test < precision
    }
  }
  
  
  if(result.matrix){
    ## lambda, cext, cabs, csca, dext, dabs, dsca
    return(cbind(wavelength = material[["wavelength"]], results))
    
  } else {
    
    d <- data.frame(material[["wavelength"]], results) # L - R
    names(d) <- c("wavelength", 'cext', 'cabs', 'csca',
                  "dext", "dabs", "dsca")
    
    ## checking consistency of csca
    sca2 <- d[,"cext"] - d[,"cabs"]
    error_csca <- max(abs(sca2 - d[,"csca"])/ abs(sca2 + d[,"csca"]))
    if(error_csca > 1e-3)
      warning(sprintf('consistency check: cext - cabs differs from csca by %.1f %%, try increasing Nsca [%i]',
                      100*error_csca, Nsca))
    
    L2eV <- 6.62606896e-34 * 299792458/1.602176487e-19
    m <- melt(transform(d, energy = L2eV / wavelength * 1e9), 
              id=c("wavelength", "energy"))
    
    m$type <- m$variable
    
    levels(m$type) <- list(`cross-section`="cext",
                           `cross-section`="cabs",
                           `cross-section`="csca",
                           `dichroism`="dext",
                           `dichroism`="dabs",
                           `dichroism`="dsca")
    
    levels(m$variable) <- list(extinction="cext",
                               absorption="cabs",
                               scattering="csca",
                               extinction="dext",
                               absorption="dabs",
                               scattering="dsca")
    
    return(m)
    
  }
  
}



##' dispersion spectrum
##'
##' dispersion spectrum
##' @title spectrum_dispersion
##' @param cluster list describing a cluster
##' @param material list
##' @param medium medium refractive index
##' @param angles of incident field in radians
##' @param Axes of incident field rotation character vector from ('x', 'y', 'z')
##' @param polarisation linear or circular polarisation
##' @param method: linear system (ls), or order-of-scattering (oos)
##' @param Nscat number of quadrature points in calculation of csca
##' @param cg logical, use conjugate gradient solver
##' @param born logical, use first Born approx as cg guess
##' @param maxiter integer termination of conjugate gradient solver
##' @param tol double, tolerance of conjugate gradient solver
##' @param progress logical, display progress bar
##' @return data.frame
##' @note The incident wavevector is along the z direction.
##' @export
##' @family user_level cda
##' @author baptiste Auguie
spectrum_dispersion <- function (cluster, material, medium = 1.33,
                                 Incidence=0, Axes='z', 
                                 polarisation=c("linear", "circular"), 
                                 method = c("ls", "oos"),
                                 cg = FALSE, born = FALSE, 
                                 Nsca = 50,
                                 maxiter=30, tol=1e-4,
                                 full = TRUE,
                                 progress = FALSE) 
{
  
  ## check whether material parameters correspond to 
  ## epsilon (NPs) pr alpha (dyes)
  isMolecular <- "alpha" %in% names(material)
  
  if(isMolecular){
    message("seems we are dealing with molecular dipoles...")
    alphabar <- material[["alpha"]]
    alphaeff = alpha_embedded(alphabar, medium)
    Alpha = alpha_rescale(alphaeff, cluster[["sizes"]])
    
  } else {
    Alpha <- alpha_ellipsoid(material[["wavelength"]], 
                             material[["epsilon"]], medium, 
                             cluster[["sizes"]])
  }
  
  # incident field
  kn <- 2*pi/material[["wavelength"]]*medium
  
  polarisation <- match.arg(polarisation)  
  method <- match.arg(method)
  
  if(polarisation == "linear") polarisation <- 0L else 
    if(polarisation == "circular") polarisation <- 1L 
  
  
  Nl <- length(kn)
  Nr <- ncol(cluster$positions)
  Ni <- length(Incidence)
  
  if(length(Axes) == 1) Axes <- rep(Axes, length.out=Ni)
  Axeso <- Axes # original codes
  Axes <- as.integer(factor(Axes, levels=c('x','y','z')))-1L
  stopifnot(all(Axes %in% c(0L, 1L, 2L)), !any(is.na(Axes)))
  
  stopifnot(Ni == length(Axes))
  stopifnot(is.matrix(Alpha), is.vector(Incidence), 
            is.matrix(cluster$positions), 
            is.matrix(cluster$angles))
  
  stopifnot(nrow(Alpha)/3 == Nr, 
            ncol(Alpha) == Nl)
  
  Scattering <- quadrature_sphere(Nquad=Nsca, "gl")
  
  if(method == "ls") {
    res <- cpp_dispersion_spectrum(kn, medium, cluster$positions, Alpha, 
                                   cluster$angles, Incidence, Axes, 
                                   Scattering$nodes,
                                   Scattering$weights,
                                   polarisation, full, cg, born, 
                                   maxiter, tol, progress)
    
  } else if(method == "oos") {
    res <- cpp_dispersion_spectrum_oos(kn, medium, cluster$positions, Alpha, 
                                       cluster$angles, Incidence, Axes, 
                                       Scattering$nodes,
                                       Scattering$weights,
                                       polarisation, full, 
                                       maxiter, tol, progress)
  }
  
  
  Incidence <- Incidence[rep(seq_len(Ni), Nl)]
  Axes <- Axeso[rep(seq_len(Ni), Nl)]
  wavelength <- rep(material$wavelength, each = Ni)
  
  val <- list(cext1 = as.vector(res[, 1, , drop = TRUE]),
              cabs1 = as.vector(res[, 2, , drop = TRUE]),
              csca1 = as.vector(res[, 3, , drop = TRUE]),
              cext2 = as.vector(res[, 4, , drop = TRUE]),
              cabs2 = as.vector(res[, 5, , drop = TRUE]),
              csca2 = as.vector(res[, 6, , drop = TRUE]))
  
  
  ## checking consistency of csca
  sca1 <- c(val[["csca1"]], val[["csca2"]])
  sca2 <- c(val[["cext1"]], val[["cext2"]]) - c(val[["cabs1"]], val[["cabs2"]])
  error_csca <- max(abs(sca2 - sca1)/ abs(sca2 + sca1))
  if(error_csca > 1e-3)
    warning(sprintf('consistency check: cext - cabs differs from csca by %.1f %%, try increasing Nsca [%i]',
                    100*error_csca, Nsca))
  
  val <- transform(val, 
                   dext = cext1 - cext2,
                   dabs = cabs1 - cabs2,
                   dsca = csca1 - csca2)
  
  
  
  results <- 
    rbind(data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["cext1"]],
                     variable = "extinction", type = "cross-section",
                     polarisation = "1"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["cabs1"]],
                     variable = "absorption", type = "cross-section",
                     polarisation = "1"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["csca1"]],
                     variable = "scattering", type = "cross-section",
                     polarisation = "1"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["cext2"]],
                     variable = "extinction", type = "cross-section",
                     polarisation = "2"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["cabs2"]],
                     variable = "absorption", type = "cross-section",
                     polarisation = "2"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["csca2"]],
                     variable = "scattering", type = "cross-section",
                     polarisation = "2"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["dext"]],
                     variable = "extinction", type = "dichroism",
                     polarisation = "3"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["dabs"]],
                     variable = "absorption", type = "dichroism",
                     polarisation = "3"),
          data.frame(wavelength = wavelength, Incidence = Incidence,
                     Axes=Axes,
                     value = val[["dsca"]],
                     variable = "scattering", type = "dichroism",
                     polarisation = "3"))
  
  if(polarisation == 0L)
    results$polarisation <- factor(results$polarisation, labels = c("p", "s", "p-s"))
  
  if(polarisation == 1L)
    results$polarisation <- factor(results$polarisation, labels = c("R", "L", "R-L"))
  
  ## reorder levels
  results$type <- factor(results$type, levels=c("cross-section", "dichroism"))
  invisible(results)
}


