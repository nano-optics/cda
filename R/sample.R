##
## Functions for creating specific geometries
##

##' Fibonacci coverage of a sphere
##'
##' Produces a set of points that covers rather uniformly the unit sphere with N points
##' with a spiral-like pattern based on a Fibonacci sequence
##' @title sample_fibonacci
##' @param N number of points
##' @describeIn sample_random 
##' @export
##' @family low_level sample fibonacci sampling of a sphere
##' @author baptiste Auguie
##' @export
sample_fibonacci <- function(N=301){
  N0 <- N
  if(N%%2 == 1) N0 <- N+1
  P <- (N0-1)/2
  ii <- seq(-P,P,by=1)
  # note: uses latitude (internally), not colatitude
  # but we don't use angles, just xyz in the end
  lat <- asin(2*ii/N0)
  Phi <- (1+sqrt(5))/2
  long <- 2*pi*ii/Phi
  
  # return N xyz positions
  rbind(x = cos(long)*cos(lat), 
        y = sin(long)*cos(lat), 
        z = sin(lat))[,1:N]
}

