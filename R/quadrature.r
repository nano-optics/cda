##
## Functions for orientation averaging
##


##' Quadrature points on a sphere
##'
##' Numerical integration points for angular averaging
##' @title quadrature_sphere
##' @param Nq number of integration points
##' @param quadrature quadrature method, using either Gauss Legendre quadrature (default), Quasi Monte Carlo, regular grid, or "cheap" (3 axes)
##' @param init (qmc method only) logical: restart, or continue from previous call
##' @importFrom randtoolbox halton
##' @importFrom statmod gauss.quad
##' @export
##' @family low_level quadrature
##' @author baptiste Auguie
quadrature_sphere <- function(Nq = 30,
    quadrature = c("qmc", "gl", "cheap", "random", 'grid','grid2', 'fibonacci', 'fibonacci2'),
    init = TRUE){

  quadrature <- match.arg(quadrature)

  if(quadrature == "cheap"){

    # nodes <- cbind(c(0, pi/2, 0), # along z
    #                 c(pi/2, 0, 0), # along x
    #                 c(pi/2, pi/2, 0)) # along y

    nodes <- rbind(alpha = c(0, pi/2, pi/2),
                        beta  = c(pi/2, 0, pi/2),
                        gamma = c(0, 0, 0))

    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }

  if(quadrature == "qmc"){ # quasi monte-carlo

    p <- randtoolbox::halton(Nq, dim = 2, normal=FALSE, init=init)

    alpha <- p[,1]*2*pi
    beta <- acos(2*p[,2] - 1) # cos(beta) in [-1,1]
    nodes <- rbind(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }

  if(quadrature == "random"){ # monte-carlo with random points

    alpha <- runif(Nq, 0, 2*pi) # uniform [-pi,pi]
    beta <- acos(runif(Nq, -1, 1)) # cos-uniform [-1,1]
    nodes <- rbind(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }
  
  if(quadrature == "grid"){ # grid in acos beta and alpha
    # might have slightly more than N total points
    # N = 2P^2 -> P = [sqrt(N/2)]
    # now check +/-1 on each one and pick best combo
    P <- sqrt(Nq/2)
    opts <- expand.grid(N1 = round(P) + seq(-1,1), N2=2*round(P) + seq(-1,1))
    prod <- opts$N1 * opts$N2
    best <- which.min(abs(Nq - (prod+2))) # add two poles
    Nbeta <- opts$N1[best]
    Nalpha <- opts$N2[best]
    
    alpha <- seq(0, 2*pi *(1 - 1/Nalpha), by = 2*pi/Nalpha)
    beta <- acos(seq(-1+1/(2*Nbeta+1),1-1/(2*Nbeta+1), length.out = Nbeta)) 
    # cos-uniform ]-1,1[ exclude poles as otherwise many points there
    # nodes <- rbind(expand.grid(alpha=alpha, beta=beta, gamma=0),
    #                c(0,0,0), c(0, pi, 0)) # add two poles
    nodes <- expand.grid(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/nrow(nodes), nrow(nodes))
    
    return(list(nodes=t(nodes), weights=weights))
  }
  if(quadrature == "grid2"){ # grid in acos beta and alpha
    # might have slightly more than N total points
    # N = 2P^2 -> P = [sqrt(N/2)]
    # now check +/-1 on each one and pick best combo
    P <- sqrt(Nq/2)
    opts <- expand.grid(N1 = round(P) + seq(-1,1), N2=2*round(P) + seq(-1,1))
    prod <- opts$N1 * opts$N2
    best <- which.min(abs(Nq - (prod+2))) # add two poles
    Nbeta <- opts$N1[best]
    Nalpha <- opts$N2[best]
    
    alpha <- seq(0, 2*pi *(1 - 1/Nalpha), by = 2*pi/Nalpha)
    # beta <- acos(seq(-1+1/(2*Nbeta+1),1-1/(2*Nbeta+1), length.out = Nbeta)) 
    step <- 2/Nbeta
    beta <- acos(seq(-1+step,1-step, length.out = Nbeta)) 
    # cos-uniform ]-1,1[ exclude poles as otherwise many points there
    # nodes <- rbind(expand.grid(alpha=alpha, beta=beta, gamma=0),
    #                c(0,0,0), c(0, pi, 0)) # add two poles
    nodes <- expand.grid(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/nrow(nodes), nrow(nodes))
    
    return(list(nodes=t(nodes), weights=weights))
  }
  
  
  if(quadrature == "fibonacci"){ # spiral
    
    N0 <- Nq
    if(Nq%%2 == 1) N0 <- Nq+1
    P <- (N0-1)/2
    ii <- seq(-P,P,by=1)
    # asin(2*ii/N0) latitude = asin(x) -> colatitude = pi/2 - asin(x)
    beta <- pi/2 - asin(2*ii/N0)
    golden <- (1+sqrt(5))/2
    alpha <- (2*pi*ii/golden) %% (2*pi)
    
    nodes <- rbind(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }
  
  
  if(quadrature == "fibonacci2"){ # spiral
    
    N0 <- Nq
    if(Nq%%2 == 1) N0 <- Nq+1
    
    ii <- seq(0,N0,by=1)
    
    beta <- acos(1 - (2*ii+1)/N0)
    golden <- (1+sqrt(5))/2
    alpha <- (2*pi*ii/golden) %% (2*pi)
    
    nodes <- rbind(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }
  

  if(quadrature == "gl"){ #  gauss legendre along cos beta, grid alpha
    # might have slightly more than N total points
    # N = 2P^2 -> P = [sqrt(N/2)]
    # now check +/-1 on each one and pick best combo
    P <- sqrt(Nq/2)
    opts <- expand.grid(N1 = round(P) + seq(-1,1), N2=2*round(P) + seq(-1,1))
    prod <- opts$N1 * opts$N2
    best <- which.min(abs(Nq - prod))
    Nbeta <- opts$N1[best]
    Nalpha <- opts$N2[best]

    # scale the coordinates from (-1, 1) to (0, 2pi)
    # and cos beta in (-1, 1) resp. (that one unnecessary btw)
    psi1=-1; psi2=1;
    C2 = (psi2 - psi1) / 2; D2 = (psi2 + psi1) / 2;
    
    GL_alpha <- seq(0, 2*pi *(1 - 1/Nalpha), by = 2*pi/Nalpha)
    GL_beta <- statmod::gauss.quad(Nbeta)

    alpha = GL_alpha
    beta  = acos(GL_beta$nodes*C2 + D2)

    # grid of angles, gamma is constant here
    nodes <- expand.grid(alpha=alpha, beta=beta, gamma=0)
    # corresponding weights for 1D quadrature
    weights <- expand.grid(alpha=rep(1/Nalpha, Nalpha),
                           beta=GL_beta$weights)
    # combine the weigths
    weights <- C2 / 2 * weights$alpha * weights$beta

    return(list(nodes=t(as.matrix(nodes)), weights=weights))
  }


}
