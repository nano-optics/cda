
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


compose_euler <- function( phi1,theta1,psi1, phi2,theta2,psi2 ){
  
  cosphi = cos(phi1) ; cospsi = cos(psi1) ; costheta = cos(theta1) 
  sinphi = sin(phi1) ; sinpsi = sin(psi1) ; sintheta = sin(theta1) 
  
  R11 = cosphi*costheta*cospsi - sinphi*sinpsi 
  R21 = sinphi*costheta*cospsi + cosphi*sinpsi 
  R31 = -sintheta*cospsi 
  
  R12 = -cosphi*costheta*sinpsi - sinphi*cospsi 
  R22 = -sinphi*costheta*sinpsi + cosphi*cospsi 
  R32 = sintheta*sinpsi 
  
  R13 = cosphi*sintheta 
  R23 = sinphi*sintheta 
  R33 = costheta 
  
  
  cosphi = cos(phi2) ; cospsi = cos(psi2) ; costheta = cos(theta2) 
  sinphi = sin(phi2) ; sinpsi = sin(psi2) ; sintheta = sin(theta2) 
  
  S11 = cosphi*costheta*cospsi - sinphi*sinpsi 
  S21 = sinphi*costheta*cospsi + cosphi*sinpsi 
  S31 = -sintheta*cospsi 
  
  S12 = -cosphi*costheta*sinpsi - sinphi*cospsi 
  S22 = -sinphi*costheta*sinpsi + cosphi*cospsi 
  S32 = sintheta*sinpsi 
  
  S13 = cosphi*sintheta 
  S23 = sinphi*sintheta 
  S33 = costheta 
  
  # combined rotation
  T11 = R11 * S11 + R12 * S21 + R13 * S31 
  T12 = R11 * S12 + R12 * S22 + R13 * S32 
  T13 = R11 * S13 + R12 * S23 + R13 * S33 
  
  T21 = R21 * S11 + R22 * S21 + R23 * S31 
  T22 = R21 * S12 + R22 * S22 + R23 * S32 
  T23 = R21 * S13 + R22 * S23 + R23 * S33 
  
  T31 = R31 * S11 + R32 * S21 + R33 * S31 
  T32 = R31 * S12 + R32 * S22 + R33 * S32 
  T33 = R31 * S13 + R32 * S23 + R33 * S33 
  
  list(theta = acos(T33) ,
  phi = atan2(T23, T13) ,
  psi = atan2(T32, -T31) )
  
}

