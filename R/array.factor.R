

##' C++ calculation of the array factor
##'
##' C++ calculation of the array factor
##' @title array factor
##' @param wavelength wavelength in microns
##' @param N half number of dipoles along one side
##' @param pitch pitch in microns
##' @return S
##' @export
##' @family user_level array
##' @author baptiste Auguie
array_factor <- function(wavelength, N, pitch){

  k <- 2*pi/wavelength
  rj <- expand.grid(-N:N, -N:N) * pitch
  rj <- rj[-((dim(rj)[1]-1)/2 + 1),1:2] # remove rj=(0,0)
  S <- array$array_factor(k, as.matrix(rj))

  invisible(data.frame(wavelength=wavelength, S=S))
  
}
