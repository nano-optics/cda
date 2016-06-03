##' C++ calculation of the array factor
##'
##' Brute-force numerical evaluation of the truncated 2D sum of dipole fields in a finite square array
##' @title array factor
##' @param wavelength wavelength in nm
##' @param N half the number of dipoles along one side
##' @param pitch pitch in nm
##' @return complex array factor
##' @export
##' @family user_level array
##' @author baptiste Auguie
##' @examples
##' S <- array_factor(seq(400, 600),  10,  500)
##' str(S)
array_factor <- function(wavelength, N, pitch){

  k <- 2*pi/wavelength
  rj <- expand.grid(-N:N, -N:N) * pitch
  rj <- rj[-((dim(rj)[1]-1)/2 + 1),1:2] # remove rj=(0,0)
  S <- cpp_array_factor(k, as.matrix(rj))

  data.frame(wavelength=wavelength, S=S)
  
}

#' Precomputed array factor for a square lattice at normal incidence
#'
#' Exact calculation of the array factor using code from Javier Garcia de Abajo 
#' (part of the pxtal program for multiple scattering calculations 
#' in infinite layered 2D arrays)
#'
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'   \item{wavelength}{normalised wavelength lambda/pitch}
#'   \item{Qx}{in-plane component of the wavevector (0, since normal incidence)}
#'   \item{Gxx}{complex value of the array factor}
#' }
#' @source Javier Garcia de Abajo
"G0"

#' Precomputed array factor for a square lattice at normal incidence
#'
#' Exact calculation of the array factor using code from Javier Garcia de Abajo 
#' (part of the pxtal program for multiple scattering calculations 
#' in infinite layered 2D arrays)
#'
#' @format A list of two interpolation functions:
#' \describe{
#'   \item{re}{real part of G0}
#'   \item{im}{imaginary part of G0}
#' }
#' @source Javier Garcia de Abajo
"gfun"
