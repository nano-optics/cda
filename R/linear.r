##' extinction spectrum
##'
##' extinction spectrum
##' @title extinction_spectrum
##' @param cluster list describing a cluster
##' @param material list
##' @param n medium refractive index
##' @param progress logical, display progress bar
##' @return data.frame
##' @export
##' @family user_level linear
##' @author baptiste Auguie
linear_extinction_spectrum <- function(cluster, material, n=1.33, progress=FALSE){
  
  k0 <- 2*pi/material$wavelength
  kn <- k0*n
  invalpha <- make.invalpha(cluster, material, polarizability.fun=polarizability.ellipsoid, n=n, kuwata=TRUE)
  
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

##' dispersion spectrum
##'
##' dispersion spectrum
##' @title dispersion_spectrum
##' @param cluster list describing a cluster
##' @param angles matrix of incident Euler angles
##' @param material list
##' @param n medium refractive index
##' @param invert logical
##' @param progress logical, display progress bar
##' @return data.frame
##' @export
##' @family user_level cda
##' @author baptiste Auguie
dispersion_spectrum <- function (cluster, angles, material, n = 1.33,
                                 polarisation=c("linear", "circular"), 
                                 invert = FALSE, progress = FALSE) 
{
    k0 <- 2 * pi/material$wavelength
    kn <- k0 * n
    polarisation <- match.arg(polarisation)
    polarisation <- if(polarisation == "linear") 0L else if(polarisation == "circular") 1L 
    
    invalpha <- make.invalpha(cluster, material, polarizability.fun = polarizability.ellipsoid, 
        n = n, kuwata = TRUE)
    Nwavelengths <- length(k0)
    Nparticles <- nrow(cluster$r)
    Nangles <- nrow(angles)
    stopifnot(is.matrix(invalpha), is.matrix(angles), is.matrix(cluster$r), 
        is.matrix(cluster$angles))
    stopifnot(ncol(invalpha)/3 == Nparticles, nrow(invalpha) == 
        Nwavelengths)
    res <- cda$dispersion_spectrum(kn, invalpha, cluster$r, 
        cluster$angles, angles, as.integer(polarisation), as.integer(invert), as.integer(progress))

    angles <- angles[rep(seq.int(Nangles), Nwavelengths), ]
    wavelength <- rep(material$wavelength, each = Nangles)
    
    ## print(str(res))
    results <- 
    rbind(data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 1, , drop = TRUE]),
                     type = "extinction", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 2, , drop = TRUE]),
                     type = "absorption", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 3, , drop = TRUE]),
                     type = "scattering", polarisation = "1"),
          data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 4, , drop = TRUE]),
                     type = "extinction", polarisation = "2"),
          data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 5, , drop = TRUE]),
                     type = "absorption", polarisation = "2"),
          data.frame(wavelength = wavelength, angles = angles,
                     value = c(res[, 6, , drop = TRUE]),
                     type = "scattering", polarisation = "2"))
    if(polarisation == 0L)
      results$polarisation <- factor(results$polarisation, labels = c("p", "s"))
    
    if(polarisation == 1L){
    
      results <- rbind(results,  data.frame(wavelength = wavelength, angles = angles,
                                            value = c(res[, 1, , drop = TRUE]  - res[, 4, , drop = TRUE]),
                                            type = "extinction", polarisation = "CD"))
      ## results$polarisation <- factor(results$polarisation, labels = c("R", "L"))
    }
    
    invisible(results)
}
