
##' dye_coverage
##'
##' @title dye_coverage
##' @export
##' @param rho surface demsity
##' @param R radius
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
##' @param rv equivolume sphere radius
##' @param h aspect ratio
##' @param type class of spheroid
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


## small utility functions


##' Utility function to create clusters 
##'
##' Identical particles
##' @title equal_sizes
##' @param a semi-axis along x
##' @param b semi-axis along y
##' @param c semi-axis along z
##' @param N number of particles
##' @return 3xN matrix
##' @author baptiste Auguie
##' @export
##' @family user_level cda utility
equal_sizes <- function(a, b, c, N){
  
  rbind(a = rep(a,N),
        b = rep(b,N),
        c = rep(c,N))
}


##' Utility function to create clusters 
##'
##' Identical particles
##' @title equal_angles
##' @param phi Euler angle
##' @param theta Euler angle
##' @param gamma Euler angle
##' @param N number of particles
##' @return 3xN matrix
##' @author baptiste Auguie
##' @export
##' @family user_level cda utility
equal_angles <- function(phi, theta, gamma, N){
  
  rbind(phi = rep(phi,N),
        theta = rep(theta,N),
        gamma = rep(gamma,N))
}

