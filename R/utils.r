
##' dye_coverage
##'
##' @title dye_coverage
##' @export
##' @family user_level cda utility
##' @author baptiste Auguie
dye_coverage <- function (rho, R) 
{
  
  area = 4*pi*R^2
  N = ceiling(area * rho)
  
}

##' Spheroid described by effective radius and aspect ratio
##'
##' Describe a spheroid by the aspect ratio and effective radius of an equi-volume sphere
##' V = 4/3 pi rv^3 = 4/3 pi a^2 c
##' c = h * a
##' @title spheroid_ar
##' @export
##' @family user_level cda utility
##' @author baptiste Auguie
spheroid_ar <- function(rv, h, type=c("prolate","oblate")){
  type <- match.arg(type)
  if(type == "prolate"){
  a = rv * h^(-1/3)
  c = h * a
  } else if(type == "oblate"){
    a = rv * h^(-1/3)
    c = h * a
  } 
  list(a=a, b=a, c=c)
}


## small undocumented utility functions

##' @noRd
##' @family user_level cda utility
##' @export
equal_sizes <- function(a, b, c, N){
  
  rbind(a = rep(a,N),
        b = rep(b,N),
        c = rep(c,N))
}

##' @noRd
##' @family user_level cda utility
##' @export
equal_angles <- function(phi, theta, gamma, N){
  
  rbind(phi = rep(phi,N),
        theta = rep(theta,N),
        gamma = rep(gamma,N))
}
