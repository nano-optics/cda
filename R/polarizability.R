##
## Dipole polarizability
##

##' Principal polarisability components for a dye molecule
##'
##' The dye is modelled as a sum of Lorentz oscillators
##' @title alpha_dye
##' @param sizes matrix of particle sizes (scaling factors for polarisability tensor)
##' @param wavelength wavelength in nm
##' @param medium refractive index of incident medium
##' @param ... further parameters passed to the Lorentzian function
##' @return matrix of polarisability
##' @export
##' @family user_level polarisability
##' @author baptiste Auguie
alpha_dye <- function(sizes, wavelength, medium, ...){
  
  # scalar polar, wavelength dependent
  alphabar = alpha_bare(wavelength, ...)
  
  # effective macroscopic polarisability
  alphaeff = alpha_embedded(alphabar[["alpha"]], medium)
  
  # rescaling with sizes
  Alpha = alpha_rescale(alphaeff, sizes)
  
}

##' Principal polarisability components for an ellipsoidal particle
##'
##' This long-wavelength polarisability approximation uses the Kuwata prescription
##' @title alpha_ellipsoid
##' @param sizes matrix of cluster sizes in nm
##' @param material data.frame with wavelength and epsilon
##' @param medium refractive index of surrounding medium
##' @return matrix of polarisability
##' @export
##' @family user_level polarisability
##' @author baptiste Auguie
##' @references
##' Kuwata et al. Resonant light scattering from metal nanoparticles: Practical analysis beyond Rayleigh approximation Appl. Phys. Lett. 83, 22 (2003)
##' @details
##' The Kuwata prescription includes semi-empirical terms of radiative correction and dynamic depolarisation to better match the fully retarded dipolar response in a reasonable range of (subwavelength) sizes and aspect ratios.
## NOTE: implementation is neither clear nor efficient, should probably do it in c++, or vectorise everything
alpha_ellipsoid <- function(sizes, material, medium){
  
  Nr <- ncol(sizes)
  Nl = length(material$wavelength)
  Alpha = matrix(NA, 3*Nr, Nl)
  
  for(jj in seq_len(Nr)){
    a = sizes[1,jj]
    b = sizes[2,jj]
    c = sizes[3,jj]
    V <- 4 * pi/3 * a * b * c
    chi <- depolarisation(a, b, c)
    ind = (jj-1)*3 # 0-index corresponding to particle jj
    Alpha[ind+1, ] = alpha_kuwata(material$wavelength, material$epsilon, V, a, chi[1], medium)
    Alpha[ind+2, ] = alpha_kuwata(material$wavelength, material$epsilon, V, b, chi[2], medium)
    Alpha[ind+3, ] = alpha_kuwata(material$wavelength, material$epsilon, V, c, chi[3], medium)
  }
  
  Alpha
}

## ------------- mid-level ------------- 

##' Bare (intrinsic) polarizability of a dye in vacuum
##'
##' Sum of lorentz oscillators
##' @title alpha_bare
##' @param wavelength wavelength in nm
##' @param alpha_inf scalar real offset
##' @param alpha_k vector of oscillator strengths
##' @param lambda_k vector of oscillator wavelengths 
##' @param mu_k vector of oscillator damping terms
##' @return data.frame
##' @export
##' @family user_level polarisability
##' @author baptiste Auguie
alpha_bare <- function(wavelength=seq(300,800),
                       alpha_inf=9.6e-39,
                       alpha_k=c(5.76e-38),
                       lambda_k=c(526),
                       mu_k=c(10000)){
  
  eps0 <- 8.854e-12
  nm3 <- 1e27
  prefact <- nm3/(4*pi*eps0)
  
  lorentz <- function(alpha, lambda, mu) {
    alpha * wavelength/mu * (1 - 1 / (1 - lambda^2/wavelength^2 - 1i*lambda^2/(wavelength*mu)))
  }
  
  sums <- mapply(lorentz, alpha = alpha_k, lambda = lambda_k,
                 mu = mu_k, SIMPLIFY=TRUE)
  alpha <- alpha_inf - rowSums(sums)
  
 data.frame(wavelength=wavelength, alpha = prefact * alpha)
}


alpha_rescale <- function(alpha, sizes){
  
  # a, b, c are normalised and stacked in a single row
  if(ncol(sizes)==1) scaling <- matrix(sizes/sum(sizes), nrow=3, ncol=1) else
    scaling = sweep(sizes, 2, colSums(sizes), `/`)
  # alternative
  # abc = colSums(sizes);
  # scaling = sizes %*% diag(1/abc)
  
  # assuming all particles share same polarizability function
  # and multiplied by a given value of alpha
  tcrossprod(as.vector(scaling), alpha)
}


alpha_embedded <- function(alphabar, medium){
  
  eps_m = medium^2
  Lfact =  (eps_m + 2) /3
  1/eps_m * Lfact^2 * alphabar
}


## ------------- low-level ------------- 

##' Depolarisation factor for an ellipsoid
##'
##' calculates the 3 depolarisation factors for a general ellipsoid
##' @title depolarisation
##' @param x1 semi-axis in nm
##' @param x2 semi-axis in nm
##' @param x3 semi-axis in nm
##' @importFrom stats integrate
##' @return shape factor along x1
##' @author baptiste Auguie
##' @family user_level polarizability
depolarisation <- function (x1, x2 = x1, x3 = x2)
{
  ## scaled version to help integration
  V = x1*x2*x3
  
  integrand <- function(q, r, s) {
    1/((1 + q) * sqrt((q + 1) * (q + r^2) * (q + s^2)))
  }
  
  I1 <- integrate(integrand, r=x2/x1, s=x3/x1, lower = 0, upper = Inf)
  I2 <- integrate(integrand, r=x1/x2, s=x3/x2, lower = 0, upper = Inf)
  I3 <- integrate(integrand, r=x1/x3, s=x2/x3, lower = 0, upper = Inf)
  
  V/2 * c(I1$value / x1^3,
          I2$value / x2^3,
          I3$value / x3^3)
}

##' polarizability
##'
##' prescription from Kuwata
##' @title alpha_kuwata
##' @aliases alpha_kuwata Kuwata.A Kuwata.B
##' @param wavelength wavelength
##' @param epsilon permittivity
##' @param V volume
##' @param axis semi-axis along incident field
##' @param L shape factor
##' @param medium refractive index
##' @return polarizability
##' @family user_level polarizability
##' @author baptiste Auguie
##' @references
##' Kuwata et al. Resonant light scattering from metal nanoparticles: Practical analysis beyond Rayleigh approximation Appl. Phys. Lett. 83, 22 (2003)
alpha_kuwata <- function (wavelength, epsilon, V, axis, L, medium = 1.33)
{
  A <- Kuwata_A(L)
  B <- Kuwata_B(L)
  x0 <- 2 * pi * axis/wavelength
  eps_m <- medium^2
  denom <- (L + eps_m/(epsilon - eps_m)) +
    A * eps_m * x0^2 +
    B * eps_m^2 * x0^4 -
    (0+1i)/3 * 4 * pi^2 * eps_m^(3/2) * V/wavelength^3
  
  V/denom/(4 * pi)
}

Kuwata_A <- function(L){
  -0.4865*L - 1.046*L^2 + 0.8481*L^3
}

Kuwata_B <- function(L){
  0.01909*L + 0.19999 * L^2 + 0.6077 * L^3
}
