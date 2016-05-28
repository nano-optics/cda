
##' dye_coverage
##'
##' @title dye_coverage
##' @export
##' @family user_level cda util
##' @author baptiste Auguie
dye_coverage <- function (rho, R) 
{
  
  area = 4*pi*R^2
  N = ceiling(area * rho)
  
}

##' spheroid_ar
##'
##' Describe a spheroid by the aspect ratio and effective radius of an equi-volume sphere
##' V = 4/3 pi rv^3 = 4/3 pi a^2 c
##' c = h * a
##' @title spheroid_ar
##' @export
##' @family user_level cda util
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

