##' extinction spectrum
##'
##' extinction spectrum
##' @title extinction_spectrum
##' @param cluster list describing a cluster
##' @param material list
##' @param medium medium refractive index
##' @param progress logical, display progress bar
##' @return data.frame
##' @export
##' @family user_level linear
##' @author baptiste Auguie
linear_extinction_spectrum <- function(cluster, material, medium=1.33, progress=FALSE){
  
  k0 <- 2*pi/material$wavelength
  kn <- k0*medium
  invalpha <- inverse_polarizability(cluster, material, 
                                     polarizability_fun=polarizability_ellipsoid, 
                                     medium=medium, kuwata=TRUE)
  
    Nwavelengths <- length(k0)
    Nparticles <- nrow(cluster$r)
    
    stopifnot(is.matrix(invalpha),
              is.matrix(cluster$r), 
              is.matrix(cluster$angles))

    stopifnot(ncol(invalpha)/3 == Nparticles,
              nrow(invalpha) == Nwavelengths)
  
  res <- data.frame(linear$linear_extinction_spectrum(kn, invalpha, cluster$r,
                          cluster$angles, as.integer(progress)))
  names(res) <- c('s', 'p')
  res$wavelength <- material$wavelength
  melt(res, meas=c("s", "p"))

}
