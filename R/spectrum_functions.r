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
##' @param dN iterative increase in N (QMC only)
##' @param Nsca quadrature points for scattering cross-section
##' @param method linear system (solve), conjugate-gradient (cg), or order-of-scattering (oos)
##' @param maxiter integer termination of iterative solver
##' @param tol double, tolerance of iterative solver
##' @param progress print progress lines
##' @param verbose display messages
##' @importFrom reshape2 melt
##' @export
##' @family user_level circular_dichroism spectrum
##' @author baptiste Auguie
##' @references
##' Y. Okada, Efficient numerical orientation quadrature of light scattering properties
##' with a quasi-Monte-Carlo method, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 9, June 2008, Pages 1719-1742.
spectrum_oa <- function(cluster, material, medium=1.33,
                        quadrature = c("gl","qmc","random", "cheap"), Nq=100,
                        iterative=FALSE, precision=1e-3, Qmax=1e4, dN=Nq,
                        method = c("solve", "cg", "oos"),
                        Nsca = 50,
                        maxiter = 30, tol=1e-4,
                        progress=FALSE, verbose=TRUE){
  
  quadrature <- match.arg(quadrature)
  method <- match.arg(method)
  inversion <- switch(method, "solve" = 0L, "cg" = 1L, "oos" = 2L)
  
  ## check whether material parameters correspond to
  ## epsilon (NPs) pr alpha (dyes)
  isMolecular <- "alpha" %in% names(material)
  
  if(isMolecular){
    alphabar <- material[["alpha"]]
    alphaeff = alpha_embedded(alphabar, medium)
    Alpha = alpha_rescale(alphaeff, cluster[["sizes"]])
    
  } else {  ## other alternative must be dielectric function
    stopifnot("epsilon" %in% names(material))
    Alpha <- alpha_ellipsoid(cluster[["sizes"]], material, medium)
  }
  
  # incident field
  kn <- 2*pi/material[["wavelength"]]*medium
  
  Incidence <- quadrature_sphere(Nq=Nq, quadrature)
  
  Scattering <- quadrature_sphere(Nq=Nsca, "gl")
  
  
  results <- cpp_oa_spectrum(kn, medium, cluster$positions, Alpha,
                             cluster$angles,
                             Incidence$nodes,
                             Incidence$weights,
                             Scattering$nodes,
                             Scattering$weights,
                             inversion, maxiter, tol, progress)
  
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
      
      newres <- cpp_oa_spectrum(kn, cluster$positions, Alpha,
                                cluster$angles,
                                Incidence$nodes,
                                Incidence$weights,
                                Scattering$nodes,
                                Scattering$weights,
                                inversion, maxiter,
                                tol, progress)
      
      ## average of the two results
      results <- (oldN * results + dN * newres) / (oldN + dN)
      
      test <- max(abs(old - results[,1]) / results[,1])
      ## max relative difference in extinction cross section
      if(verbose)
        message("N:", Ntot, "; relative error: " , test)
      converged <- test < precision
    }
  }
  
  
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
  d[["energy"]] <- L2eV / d[["wavelength"]] * 1e9
  m <- melt(d, id=c("wavelength", "energy"))
  
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



##' dispersion spectrum
##'
##' dispersion spectrum
##' @title spectrum_dispersion
##' @param cluster list describing a cluster
##' @param material list
##' @param medium medium refractive index
##' @param Incidence angular directions of incident field
##' @param Axes incident field rotation axis
##' @param polarisation linear or circular polarisation
##' @param Nsca number of quadrature points in calculation of csca
##' @param method linear system (solve), conjugate-gradient (cg), or order-of-scattering (oos)
##' @param maxiter integer termination of iterative solver
##' @param tol double, tolerance of iterative solver
##' @param progress logical, display progress bar
##' @return data.frame
##' @note The incident wavevector is along the z direction.
##' @export
##' @family user_level cda spectrum
##' @author baptiste Auguie
spectrum_dispersion <- function (cluster, material, medium = 1.33,
                                 Incidence=0, Axes='z',
                                 polarisation=c("linear", "circular"),
                                 method = c("solve", "cg", "oos"),
                                 Nsca = 50,
                                 maxiter=30, tol=1e-4,
                                 progress = FALSE)
{
  
  method <- match.arg(method)
  inversion <- switch(method, "solve" = 0L, "cg" = 1L, "oos" = 2L)
  
  ## check whether material parameters correspond to
  ## epsilon (NPs) pr alpha (dyes)
  isMolecular <- "alpha" %in% names(material)
  
  if(isMolecular){
    message("seems we are dealing with molecular dipoles...")
    alphabar <- material[["alpha"]]
    alphaeff = alpha_embedded(alphabar, medium)
    Alpha = alpha_rescale(alphaeff, cluster[["sizes"]])
    
  } else {
    Alpha <- alpha_ellipsoid(cluster[["sizes"]], material, medium)
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
  
  Scattering <- quadrature_sphere(Nq=Nsca, "gl")
  
  res <- cpp_dispersion_spectrum(kn, medium, cluster$positions, Alpha,
                                 cluster$angles, Incidence, Axes,
                                 Scattering$nodes,
                                 Scattering$weights,
                                 polarisation, inversion,
                                 maxiter, tol, progress)
  
  
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
  
  ## dichroism
  val[["dext"]] <- val[["cext1"]] - val[["cext2"]]
  val[["dabs"]] <- val[["cabs1"]] - val[["cabs2"]]
  val[["dsca"]] <- val[["csca1"]] - val[["csca2"]]
  
  ## combine results in long format
  ## a bit verbose, but doesn't rely on extra packages
  
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
