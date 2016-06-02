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
    quadrature = c("qmc", "gl", "cheap", "random"),
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

    alpha <- runif(Nq, -pi, pi) # uniform [-pi,pi]
    beta <- acos(runif(Nq, -1, 1)) # cos-uniform [-1,1]
    nodes <- rbind(alpha=alpha, beta=beta, gamma=0)
    weights <- rep(1/ncol(nodes), ncol(nodes))
    return(list(nodes=nodes, weights=weights))
  }

  if(quadrature == "gl"){ #  gauss legendre
    # might have slightly more than N total points
    rndN <- ceiling(sqrt(Nq/2))

    # scale the coordinates from (-1, 1) to (0, 2pi)
    # and cos beta in (-1, 1) resp. (that one unnecessary btw)
    phi1=0; phi2=2*pi; psi1=-1; psi2=1;
    C1 = (phi2 - phi1) / 2;  D1 = (phi2 + phi1) / 2;
    C2 = (psi2 - psi1) / 2; D2 = (psi2 + psi1) / 2;

    GL_alpha <- statmod::gauss.quad(2*rndN)
    GL_beta <- statmod::gauss.quad(rndN)

    alpha = GL_alpha$nodes*C1 + D1
    beta  = acos(GL_beta$nodes*C2 + D2)

    # grid of angles, gamma is constant here
    nodes <- expand.grid(alpha=alpha, beta=beta, gamma=0)
    # corresponding weights for 1D quadrature
    weights <- expand.grid(alpha=GL_alpha$weights,
                           beta=GL_beta$weights)
    # combine the weigths
    weights <- C1 * C2 / (4*pi) * weights$alpha * weights$beta

    return(list(nodes=t(as.matrix(nodes)), weights=weights))
  }


}
