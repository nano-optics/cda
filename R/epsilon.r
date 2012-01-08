##' permittivity silver
##'
##' analytical dielectric function of Silver (Drude model)
##' @title epsAg
##' @param lambda wavelength in nm
##' @param epsilon.infty background dielectric constant
##' @param lambda.p plasma wavelength
##' @param mu.p damping constant
##' @return data.frame
##' @author baptiste Auguié
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
##' @examples
##' require(dielectric) ; data(AgPalik)
##' wvl <- seq(300, 900)
##' silver <- epsAg(wvl)
##' 
##' matplot(silver$wavelength, cbind(Re(silver$epsilon), Im(silver$epsilon)), t="l", lty=1)
##' matpoints(AgPalik$wavelength, cbind(Re(AgPalik$epsilon), Im(AgPalik$epsilon)))

epsAg <- function(lambda, epsilon.infty = 4,
                  lambda.p = 282, mu.p = 17000){
  
  data.frame(wavelength=lambda*1e-3, epsilon=
             epsilon.infty*(1 - 1 / (lambda.p^2*(1/lambda^2 + 1i / (mu.p*lambda)))))
}

##' permittivity gold
##'
##' analytical dielectric function of Au (Drude model + interband transitions)
##' @title epsAu
##' @param lambda wavelength in nm
##' @param epsilon.infty background dielectric constant
##' @param lambda.p plasma wavelength
##' @param mu.p damping constant
##' @param A1 A1
##' @param phi1 phi1
##' @param lambda1 lambda1
##' @param mu1 mu1
##' @param A2 A2
##' @param phi2 phi2
##' @param lambda2 lambda2
##' @param mu2 mu2
##' @return data.frame
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
##' @author baptiste Auguié
##' require(dielectric) ; data(AuJC)
##' wvl <- seq(300, 900)
##' gold <- epsAu(wvl)
##' 
##' matplot(gold$wavelength, cbind(Re(gold$epsilon), Im(gold$epsilon)), t="l", lty=1)
##' matpoints(AuJC$wavelength, cbind(Re(AuJC$epsilon), Im(AuJC$epsilon)))

epsAu <- function(lambda, epsilon.infty = 1.54,
                  lambda.p = 177.5, mu.p = 14500,
                  A1 = 1.27, phi1 = -pi/4, lambda1 = 470, mu1 = 1900,
                  A2 = 1.1, phi2 = -pi/4, lambda2 = 325, mu2 = 1060){
  eps.drude <- 
    epsilon.infty*(1 - 1 / (lambda.p^2*(1/lambda^2 + 1i / (mu.p*lambda))))

  data.frame(wavelength=lambda*1e-3, epsilon=
  eps.drude + A1 / lambda1 * (exp(1i*phi1)  / (1/lambda1 - 1/lambda - 1i/mu1) +
                              exp(-1i*phi1) / (1/lambda1 + 1/lambda + 1i/mu1)) +
              A2 / lambda2 * (exp(1i*phi2)  / (1/lambda2 - 1/lambda - 1i/mu2) +
                              exp(-1i*phi2) / (1/lambda2 + 1/lambda + 1i/mu2))
             )
}

