## 
## Functions to calculate the dipole polarizability
## 

##' principal polarizability components for an ellipsoidal particle
##'
##' uses the Kuwata prescription
##' @title polarizability.ellipsoid
##' @param lambda wavelength in microns
##' @param epsilon complex permittivity
##' @param a semi-axis in um
##' @param b semi-axis in um
##' @param c semi-axis in um
##' @param n surrounding medium
##' @param kuwata logical, use Kuwata or Clausius Mossotti
##' @return matrix of polarizability
##' @export
##' @family user_level polarizability
##' @author baptiste Auguie
polarizability.ellipsoid <- function(lambda, epsilon, a=0.05, b=0.03 , c=b, n = 1.33, kuwata= TRUE) 
{
  if(kuwata){
    V <- 4 * pi/3 * a * b * c
    chi.a <- La(a, b, c)
    chi.b <- La(b, a, c)
    chi.c <- La(c, a, b)

    alpha.kuwata.a <- alpha_kuwata(lambda, epsilon, V, a, chi.a, 
        n^2)
    alpha.kuwata.b <- alpha_kuwata(lambda, epsilon, V, b, chi.b, 
        n^2)
    alpha.kuwata.c <- alpha_kuwata(lambda, epsilon, V, c, chi.c, 
        n^2)
    return(cbind(alpha.kuwata.a, alpha.kuwata.b, alpha.kuwata.c))
  } else {
    cm <- a^3*(epsilon - n^2) / (epsilon +2*n^2)
    return(cbind(cm,cm,cm))
  }
}

##' Shape factor for an ellipsoid
##'
##' calculates the shape factor for a general ellipsoid
##' @title La
##' @param a semi-axis in um
##' @param b semi-axis in um
##' @param c semi-axis in um
##' @return shape factor along a
##' @author baptiste Auguie
##' @export
##' @family user_level polarizability
La <- function (a = 0.05, b = a, c = a) 
{
  ## scaled version to help convergence
  V <- a*b*c
  b <- b/a
  c <- c/a
    integrand <- function(q) {
      fq <- (q + 1) * (q + b^2) * (q + c^2)
      1/((1 + q) * sqrt(fq))
    }
  I1 <- integrate(integrand, lower = 0, upper = Inf)$value
  V/2 * I1 / a^3
  
}
##' polarizability
##'
##' prescription from Kuwata
##' @title alpha_kuwata
##' @aliases alpha_kuwata Kuwata.A Kuwata.B
##' @param lambda wavelength
##' @param epsilon permittivity
##' @param V volume
##' @param axis semi-axis along incident field
##' @param L shape factor
##' @param epsilon.medium medium permittivity
##' @return polarizability
##' @export
##' @family user_level polarizability
##' @author baptiste Auguie
alpha_kuwata <-
function (lambda, epsilon, V, axis, L, epsilon.medium = 1.33^2) 
{
    A <- Kuwata.A(L)
    B <- Kuwata.B(L)
    x <- 2 * pi * axis/lambda
    denom <- (L + epsilon.medium/(epsilon - epsilon.medium)) + A * epsilon.medium * 
        x^2 + B * epsilon.medium^2 * x^4 - (0+1i)/3 * 4 * pi^2 * epsilon.medium^(3/2) * 
        V/lambda^3
    V/denom/(4 * pi)
}

Kuwata.A <- function(L){
  -0.4865*L - 1.046*L^2 + 0.8481*L^3
}

Kuwata.B <- function(L){
  0.01909*L + 0.19999 * L^2 + 0.6077 * L^3
}	
	

##' inverse polarizability tensors
##'
##' calculates and formats the principal polarizability of several particles
##' @title make.invalpha
##' @param cluster cluster
##' @param material material
##' @param polarizability.fun polarizability function
##' @param ... additional arguments passed to polarizability.fun
##' @return  matrix with each row being the 3 principal values of each polarizability tensor
##' @export
##' @family user_level polarizability
##' @author Baptiste Auguie
make.invalpha <- function(cluster, material, polarizability.fun=polarizability.ellipsoid, ...){
  polar <- mlply(cluster$sizes, polarizability.fun,
                 lambda=material$wavelength,
                 epsilon=material$epsilon, ...)
  invalpha <- as.matrix(1/do.call(cbind, polar))
  
}
