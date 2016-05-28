

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

