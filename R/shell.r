
##' shell spectrum
##'
##' shell spectrum
##' @title spectrum_shell
##' @param cluster list describing a cluster
##' @param material list
##' @param medium medium refractive index
##' @param angles of incident field in radians
##' @param axes of incident field rotation character vector from ('x', 'y', 'z')
##' @param polarisation linear or circular polarisation
##' @param method: linear system (ls), or order-of-scattering (oos)
##' @param cg logical, use conjugate gradient solver
##' @param born logical, use first Born approx as cg guess
##' @param nmax integer termination of conjugate gradient solver
##' @param tol double, tolerance of conjugate gradient solver
##' @param progress logical, display progress bar
##' @return data.frame
##' @note The incident wavevector is along the z direction.
##' @export
##' @family user_level cda
##' @author baptiste Auguie
spectrum_shell <- function (cluster, wavelength, 
                            fun=polarizability_dye, ...,
                            medium = 1.33, core=TRUE,
                            Nquad = 100, 
                            averaging = c("cheap", "QMC", "GL"),  
                            method = c("ls", "oos"),
                            iterative = FALSE, precision = 0.001, 
                            Qmax = 10000, dN = Nquad, 
                            cg = TRUE, born = TRUE, nmax = 10, 
                            tol = 1e-02, full = TRUE, 
                            progress = FALSE, verbose = TRUE, 
                            result.matrix = FALSE) 
{
  averaging <- match.arg(averaging)
  method <- match.arg(method)
  
  silver <- epsAg(wavelength)
  epsilon <- medium^2
  prefact <- ((epsilon + 2 /3)^2) / sqrt(epsilon)
  
  k0 <- 2 * pi/wavelength
  kn <- k0 * medium
  polar <- apply(cluster[["sizes"]], MARGIN = 1, FUN = fun, 
                 wavelength = wavelength, medium=medium, ...)
  if(core){
    core <- polarizability_ellipsoid(wavelength, 
                                     epsilon = silver$epsilon,
                                     a=cluster$R0,b=cluster$R0, 
                                     medium=medium)
    Alpha <- rbind(t(do.call(cbind, polar)) * prefact, t(core))
    cluster$r <- rbind(cluster$r, c(0,0,0))
    cluster$angles <- rbind(cluster$angles, c(0,0,0))
  } else {
    Alpha <- t(do.call(cbind, polar)) * prefact
  }
  # browser()
  quadrature <- integration_points(averaging, Nquad)
  
  if(method == "ls") {
    results <- cpp_spectrum_averaging(kn, cluster$r, Alpha, 
                                      cluster$angles, 
                                      as.matrix(quadrature$angles), 
                                      quadrature$weights,
                                      full, cg, born, nmax, tol, progress)
  } else if(method == "oos") {
    results <- cpp_spectrum_averaging_os(kn, cluster$r, Alpha, 
                                         cluster$angles, 
                                         as.matrix(quadrature$angles), 
                                         quadrature$weights,
                                         full, nmax, tol, progress)
  }
  
  if (iterative && averaging == "QMC") {
    converged <- FALSE
    Ntot <- Nquad
    while (Ntot < Qmax && !converged) {
      oldN <- Ntot
      old <- results[, 1]
      Ntot <- Ntot + dN
      quadrature <- integration_points(averaging, dN, FALSE)
      
      if(method == "ls") {
        newres <- cpp_spectrum_averaging(kn, cluster$r, Alpha, 
                                         cluster$angles, 
                                         as.matrix(quadrature$angles), 
                                         quadrature$weights,
                                         full, cg, born, nmax, 
                                         tol, progress)
      } else if(method == "oos") {
        newres <- cpp_spectrum_averaging_os(kn, cluster$r, Alpha, 
                                            cluster$angles, 
                                            as.matrix(quadrature$angles), 
                                            quadrature$weights,
                                            full, nmax, tol, progress)
      }
      
      results <- (oldN * results + dN * newres)/(oldN +  dN)
      test <- max(abs(old - results[, 1])/results[, 1])
      if (verbose) 
        message("N:", Ntot, "; relative error: ", test)
      converged <- test < precision
    }
  }
  
  if(core){
  cl <- list(r=matrix(c(0,0,0), nrow=1),
             sizes = cluster$R0*equal_sizes(1,1,1,1),
             angles=equal_angles(0,0,0,1))
  
  ref <- spectrum_oa(cluster = cl,material = silver,   
                      method = "ls", 
                      averaging = "QMC", Nquad = 10,
                      medium=medium, progress = FALSE, 
                      result.matrix = TRUE)
  results <- results - ref[,-1]
} 
  
  if(result.matrix){
    ## extinction, absorption, scattering, CD ext, CD abs, CD sca
    return(cbind(wavelength = wavelength, results))
    
  } else {
    d <- data.frame(wavelength, results[,1:3])
    names(d) <- c("wavelength", "extinction", "absorption", "scattering")
    m <- reshape2::melt(d, id = c("wavelength"))
  }
  return(m)
  
}
