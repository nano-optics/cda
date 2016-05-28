##
## Clusters are lists with fields positions, angles and sizes
## which encode the geometric information about a cluster of particles.
## Basic S3 methods are provided to get an overview of a cluster's definition
##

##' @noRd
##' @importFrom utils str
##' @export
print.cluster <- function(x, ...)  str(unclass(x))

##' @noRd
##' @export
length.cluster <- function(x) nrow(x[["sizes"]])

##' @noRd
##' @export
plot.cluster <- function(x, ...){
  visualise_rgl(x, ...)
}

##' @noRd
##' @export
visualise <- function (x, type, outfile=NULL, ...)
  UseMethod("visualise")

##' @noRd
##' @export
visualise.cluster <- function(x, type=c("rgl", "povray"), outfile=NULL, ...){
  type <- match.arg(type)
  if(type == "rgl")
    visualise_rgl(x, outfile=outfile, ...) else if(type == "povray") 
      visualise_povray(x, outfile=outfile, ...)
}

##' @noRd
##' @export
equal_sizes <- function(a, b, c, N){
  
  rbind(a = rep(a,N),
        b = rep(b,N),
        c = rep(c,N))
}

##' @noRd
##' @export
equal_angles <- function(phi, theta, gamma, N){
  
  rbind(phi = rep(phi,N),
        theta = rep(theta,N),
        gamma = rep(gamma,N))
}


##' cluster_single
##'
##' a single particle
##' @title cluster_single
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @param phi angle
##' @param theta angle
##' @param psi angle
##' @return list of class cluster with fields: positions, sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
##' @examples 
##' b = cluster_ball(1000); 
##' library(rgl)
##' open3d()
##' rgl.spheres(t(b$positions))
##' rgl.spheres(0,0, 0, radius = 15, alpha=0.2, col="red")
cluster_single <- function(a, b=a, c=b, phi=0, theta=0, psi=0)
  structure(list(positions=matrix(c(0, 0, 0),3,1), 
                 angles=matrix(c(phi, theta, psi),3,1), 
                 sizes=matrix(c(a,b,c),3,1)),
            class="cluster")



##' cluster_ball
##'
##' a ball of dipoles on a cubic lattice
##' @title cluster_ball
##' @param N number of particles
##' @param R0 ball radius
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @return list of class cluster with fields: positions, sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
##' @examples 
##' b = cluster_ball(1000); 
##' library(rgl)
##' open3d()
##' rgl.spheres(t(b$positions))
##' rgl.spheres(0,0, 0, radius = 15, alpha=0.2, col="red")
cluster_ball <- function(N, R0=15, a=1, b=1, c=b){
  
  Nc <- 2*N # approx 1/2 dipoles are unoccupied in the cube
  Nl <- (ceiling(Nc^(1/3)) - 1) # linear dimension
  rr <- seq(-1, 1, length=Nl)
  m <- as.matrix(expand.grid(x=rr, y=rr, z=rr))
  distances <- m[,1]^2 + m[,2]^2 + m[,3]^2
  positions <- as.matrix(t(m[distances <= 1,]))
  N <- ncol(positions) # true N
  message(N)
  sizes <- equal_sizes(a=a, b=b, c=c, N=N)
  angles <- equal_angles(0,0,0, N=N)
  
  structure(list(positions = R0*positions,
                 sizes = sizes,
                 angles = angles,
                 N = N, R0=R0),
            class="cluster")
  
}


##' cluster_array
##'
##' cluster with a square of nanorods
##' @title cluster_array
##' @param N number of particles
##' @param pitch center-to-center distance
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @return list of class cluster with fields: positions, sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_array <- function(N, pitch=500, a=50, b=30, c=b){
  
  sqN <- round(sqrt(N))
  xyz <- expand.grid(x=seq_len(sqN)* pitch, y=seq_len(sqN)* pitch, z=0)
  positions <- t(as.matrix(xyz))
  sizes <- equal_sizes(a=a, b=b, c=c, N=sqN^2)
  angles <- equal_angles(0,0,0, N=sqN^2)
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = angles),
            class="cluster")
  
}

##' cluster_chain
##'
##' cluster with a chain of nanorods, chain axis along x
##' @title cluster_chain
##' @param N number of particles
##' @param pitch center-to-center distance
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @return list of class cluster with fields: positions, sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_chain <- function(N, pitch=500, a=50, b=30, c=b){
  
  positions <- t(as.matrix(expand.grid(x=seq_len(N) * pitch, y=0, z=0)))
  sizes <- equal_sizes(a=a, b=b, c=c, N=N)
  angles <- equal_angles(0,0,0, N=N)
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = angles),
            class="cluster")
  
}


##' cluster_dimer
##'
##' cluster with two nanorods, dimer axis along z
##' first rod along x at (0, 0, -d/2)
##' second rod at (0, 0, d/2)
##' @title cluster_dimer
##' @param d center-to-center distance
##' @param dihedral dihedral angle
##' @param alpha1 angle first rod
##' @param alpha2 angle second rod
##' @param a semi axis
##' @param b semi axis
##' @return list of class cluster with fields: positions, sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_dimer <- function(d=100,
                          a=35, b=12, c=b,
                          dihedral=pi/4, alpha1=0, alpha2=0){
  
  positions <- rbind(c(0, 0),
                     c(0, 0),
                     c(-d/2, d/2))
  angles <- rbind(c(dihedral, 0),
                  c(alpha1, alpha2),
                  c(0, 0))
  sizes <- equal_sizes(a, b, c, 2)
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = angles),
            class="cluster")
  
}


##' cluster_shell
##'
##' shell cluster
##' @title cluster_shell
##' @param N number of particles
##' @param R0 radius of core
##' @param d distance
##' @param a ellipsoid semi-axis
##' @param b ellipsoid semi-axis
##' @param c ellipsoid semi-axis
##' @param position type of random coverage
##' @param position type of angular orientation
##' @param seed random seed for reproducibility
##' @param ... extra arguments (ignored)
##' @importFrom stats runif
##' @return list
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_shell <- function(N=50, R0=30, d=1,
                          a=1, b=1, c=1, # a way to select dipole orientation
                          orientation=c("radial", "flat", "random"),
                          position = c("fibonacci", "hc", "random"),
                          exclusion = 5*N^(-1/4), seed=123, ...){
  
  ## argument check
  position <- match.arg(position)
  orientation <- match.arg(orientation)
  
  set.seed(seed) # reproducible randomness
  
  ## the shell is at a distance d from R0
  R <- R0 + d
  
  ## point picking
  if(position == "random"){
    positions <- R * sample_random(N)
  } else if(position == "hc"){
    positions <- R * sample_hc(N, exclusion/R, ...)
  } else if(position == "fibonacci"){
    positions <- R * sample_fibonacci(N)
  }
  
  sizes <- equal_sizes(a, b, c, N)
  
  ## dipole orientation
  if(orientation == "random"){
    
    angles <- rbind(phi = runif(N,0,2*pi),
                    theta = runif(N,0,2*pi),
                    psi = rep(0, N))
    
  } else if(orientation == "flat"){
    
    
    phi <- atan2(positions[2,], positions[1,])
    theta <- acos(positions[3,]/R)
    tangent1 <- rbind(-sin(phi), cos(phi), 0)
    tangent2 <- rbind(cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta))
    
    tangent <- matrix(runif(N, -1, 1), nrow=3, ncol=N, byrow=TRUE) * tangent1 +
      matrix(runif(N, -1, 1), nrow=3, ncol=N, byrow=TRUE) * tangent2
    
    phi <- atan2(tangent[2,], tangent[1,])
    theta <- acos(tangent[3,]/sqrt(colSums(tangent*tangent)))
    angles <- rbind(phi=phi, theta=theta, psi=0)
    
    
  } else if(orientation == "radial"){
    
    
    phi <- atan2(positions[2,], positions[1,])
    theta <- acos(positions[3,]/R)
    angles <- rbind(phi=phi, theta=theta, psi=0)
    
  }
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = angles,
                 R0 = R0, d = d),
            class="cluster")
  
}

##' @noRd
##' @export
helix <- function(R0=500, pitch=600, N=5,
                  delta=pi/8, delta0=pi/2, n.smooth=100*N, right=TRUE){
  
  handedness <- (-1)^(!right)
  phase <- seq(from=delta0, by=delta, length=N)
  phase2 <- seq(from=delta0, max(phase), length=n.smooth)
  
  xyz <- function(ph){
    x = R0*cos(ph)
    y = R0*sin(ph)
    z <- handedness * ph * pitch / (2*pi)
    
    rbind(x, y, z)
  }
  positions <- xyz(phase)
  centering <- mean(positions[3,])
  positions[3,] <- positions[3,] - centering
  positions2 <- xyz(phase2)
  positions2[3,] <- positions2[3,] - centering
  xp <-  - positions[2,]
  yp <-   positions[1,]
  zp <-  handedness * pitch / (2*pi)
  n <- sqrt(xp^2+yp^2+zp^2)
  
  phi <-  atan2(yp, xp)
  theta <-  acos(zp/n)
  
  list(positions=as.matrix(positions),
       angles = rbind(phi, theta, 0),
       R0=R0, smooth=as.matrix(positions2))
}


##' cluster_helix
##'
##' helical cluster of ellipsoids
##' @title cluster_helix
##' @param N number of particles
##' @param R0 radius of helix
##' @param pitch pitch of helix
##' @param delta angle between particles
##' @param delta0 initial angle
##' @param right logical, helicity
##' @param a ellipsoid semi-axis
##' @param b ellipsoid semi-axis
##' @param c ellipsoid semi-axis
##' @param angles type of angular orientation
##' @param seed random seed for reproducibility
##' @param ... extra arguments (ignored)
##' @return list
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_helix <- function(N=5,
                          a=10, b=10, c=20,
                          R0=100, pitch=200,
                          delta=pi/5, delta0=0, right=TRUE,
                          angles=c("helix", "random", "fixed"),
                          seed=123, ...){
  
  hel <- helix(R0=R0, pitch=pitch, N=N, delta=delta, delta0=delta0, right=right)
  nr <- NROW(hel$angles)
  positions <- as.matrix(hel$positions)
  set.seed(seed) # reproducible
  angles <- switch(match.arg(angles),
                   "helix" = hel$angles,
                   "fixed" = cbind(rep(0, nr),
                                   rep(0, nr),
                                   rep(0, nr)),
                   "random" = matrix(cbind(runif(nr, 0, 2*pi),
                                           runif(nr, 0, 2*pi),
                                           runif(nr, 0, 2*pi)), ncol=3, byrow=T))
  
  sizes <- equal_sizes(a, b, c, N)
  
  structure(list(positions = positions,
                 sizes = sizes,
                 angles = as.matrix(angles),
                 R0 = R0),
            class="cluster")
}
