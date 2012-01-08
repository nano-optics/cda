## 
## Functions for calculations with fixed incidence and linear polarisation
## 

##' extinction spectrum
##'
##' extinction spectrum
##' @title extinction_spectrum
##' @param cluster list describing a cluster
##' @param material list
##' @param n medium refractive index
##' @param progress logical, display progress bar
##' @return data.frame
##' @author baptiste Auguié
linear_extinction_spectrum <- function(cluster, material, n=1.33, progress=FALSE){
  
  k0 <- 2*pi/material$wavelength
  kn <- k0*n
  invalpha <- make.invalpha(cluster, material, polarizability.fun=polarizability.ellipsoid, n=n, kuwata=TRUE)
  
  res <- data.frame(linear$linear_extinction_spectrum(kn, invalpha, cluster$r,
                          cluster$angles, as.integer(progress)))
  names(res) <- c('s', 'p')
  res$wavelength <- material$wavelength
  melt(res, meas=c("s", "p"))

}


