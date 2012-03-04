## 
## Functions for the creation of special clusters
## 


##' clust.equalsizes
##'
##' generate a matrix of equal particle sizes
##' @title clust.equalsizes
##' @param a a
##' @param b b
##' @param c c
##' @param N N
##' @return matrix Nx3
##' @author baptiste Auguie
clust.equalsizes <- function(a, b, c, N)
ldply(1:N, function(.) data.frame(a=a, b=b, c=c))


##' clust.equalangles
##'
##' generate a matrix of equal angles
##' @title clust.equalangles
##' @param phi phi
##' @param theta theta
##' @param psi psi
##' @param N N
##' @return matrix Nx3
##' @author baptiste Auguie
clust.equalangles <- function(phi=0, theta=0, psi=0, N)
cbind(phi=rep(phi, N), theta=rep(theta, N), psi=rep(psi, N))

##' helix curve
##'
##' add particles on an helix
##' @title helix
##' @param R0 radius in microns
##' @param pitch pitch in microns
##' @param N number of particles
##' @param delta twist angle between two particles
##' @param delta0 angle shift at origin
##' @param n.smooth  number of interpolation points (for plotting purposes)
##' @param right logical, handedness
##' @return list with r,  sizes,  invalpha,  angles, R0 and smooth interp. points
##' @author baptiste Auguie
##' @examples 
##' cl <- helix(0.5, 1, 36, delta=pi/6, n.smooth=1e3) ; str(cl)
##' \dontrun{require(rgl)
##' open3d()
##' spheres3d(cl$smooth, radius=0.01,col=2)
##' spheres3d(cl$positions, radius=0.1, col="gold")
##' ## ellipsoids are oriented following the helix
##' sizes <- clust.equalsizes(0.04,0.02,0.02,NROW(cl$positions))
##' rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold") }

helix <- function(R0=0.5, pitch=0.6, N=5, 
                  delta=pi/8, delta0=pi/2, n.smooth=100*N, right=TRUE){

  handedness <- (-1)^(!right)
  phase <- seq(from=delta0, by=delta, length=N)
  phase2 <- seq(from=delta0, max(phase), length=n.smooth)
  
  xyz <- function(ph){
    x = R0*cos(ph)
    y = R0*sin(ph)
    z <- handedness * ph * pitch / (2*pi)
    
    data.frame(x=x, y=y, z=z)
}
  positions <- xyz(phase)
  centering <- mean(positions$z)
  positions$z <- positions$z - centering
  positions2 <- xyz(phase2)
  positions2$z <- positions2$z - centering
  xp <-  - positions$y 
  yp <-   positions$x 
  zp <-  handedness * pitch / (2*pi)
  n <- sqrt(xp^2+yp^2+zp^2) 
  
  phi <-  atan2(yp, xp)
  psi <-  asin(zp/n)

  list(positions=as.matrix(positions),
       angles = cbind(phi, pi/2, psi),
       radius=R0, smooth=as.matrix(positions2))
}

##' helix.zt
##'
##' create a matrix of xyz positions for ellipsoids arranged on a cylinder with an helical twist
##' @title helix2
##' @param R0 radius in um
##' @param pitch pitch in um
##' @param z positions of particles along cylinder axis
##' @param theta angular position along cylinder
##' @param right logical, handedness
##' @return matrix N x 3
##' @author baptiste Auguie
helix.zt <- function (R0 = 0.5, pitch = 0.6, z=runif(100,-0.5,0.5), theta=runif(100,0,2*pi), right = TRUE) {
  x = R0 * cos(theta)
  y = R0 * sin(theta)
  handedness <- (-1)^(!right)
  positions <- data.frame(x = x, y = y, z = z)
    
  xp <- -positions$y
  yp <- positions$x
  zp <- handedness * pitch/(2 * pi)
  n <- sqrt(xp^2 + yp^2 + zp^2)
  phi <- atan2(yp, xp)
  psi <- asin(zp/n)
  
  list(positions = as.matrix(positions), angles = cbind(phi, pi/2, psi), radius = R0)
}




##' make a cluster of spheres
##'
##' make a cluster of spheres
##' @title makeSpheresCluster
##' @param radius radius
##' @param N number of spheres
##' @param R0 radius helix
##' @param pitch pitch helix
##' @param delta phase helix
##' @param delta0 initial phase helix
##' @param right logical, handedness
##' @return list of positions, sizes and angles
##' @author baptiste Auguie
makeSpheresCluster <- function(radius = 0.005, N, R0=0.5, pitch=0.6, 
                               delta=pi/8, delta0=0, right=TRUE){

  hel <- helix(N=N, R0=R0, pitch=pitch, delta=delta, delta0=delta0, right=right)
  sizes <- clust.equalsizes(a=radius, b=radius, c=radius, N=N)
  angles <- clust.equalangles(N=N)
  list(r=hel$positions, sizes=sizes, angles=angles)
  
}
##' makeDimerCluster
##'
##' makeDimerCluster
##' first rod along x at (0, 0, -d/2)
##' second rod at (0, 0, d/2),  rotated by theta (z) in [0, pi], psi (z') in [-pi/2, pi/2]
##' @title makeDimerCluster
##' @param d center-to-center distance
##' @param phi longitude
##' @param psi latitude
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @param right logical,  handedness
##' @return list with r,  sizes,  angles
##' @author baptiste Auguie
makeDimerCluster <- function(d=a, 
                             phi=pi/4, psi=0,
                             a=35e-3, b=12e-3, c=b,
                             right=TRUE){

  r <- cbind(c(0,0), c(0, 0), c(-d/2, d/2))
  sizes <- clust.equalsizes(a=a, b=b, c=c, N=2)
  phi <- if(right) phi else -phi
  angles <- cbind(c(0, phi), c(0, pi/2), c(0, psi))
  list(r=r, sizes=sizes, angles=angles)
  
}

##' makeDimerDihedral
##'
##' makeDimerDihedral
##' first rod along x at (0, 0, -d/2)
##' second rod at (0, 0, d/2)
##' @title makeDimerDihedral
##' @param d center-to-center distance
##' @param dihedral dihedral angle
##' @param alpha1 angle first rod
##' @param alpha2 angle second rod
##' @param a semi axis
##' @param b semi axis
##' @param right handedness
##' @return list with r,  sizes,  angles
##' @author baptiste Auguie
makeDimerDihedral <- function(d=a, 
                             dihedral=0, alpha1=0, alpha2=0,
                             a=35e-3, b=12e-3, 
                              right=TRUE){

  r <- cbind(c(0,0), c(0, 0), c(-d/2, d/2))
  sizes <- clust.equalsizes(a=a, b=b, c=b, N=2)  
  angles <- cbind(c(dihedral, 0), c(pi/2, pi/2), c(alpha1, alpha2))
  list(r=r, sizes=sizes, angles=angles)
  
}


##' makeRodChain
##'
##' makeRodChain
##' @title makeRodChain
##' @param N number of rods
##' @param pitch pitch
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @return list with r,  sizes, angles
##' @author baptiste Auguie
makeRodChain <- function(N, pitch=0.5, a=50e-3, b=30e-3,c=b){

  r <- as.matrix(expand.grid(x=0, y=seq(1,N) * pitch, z=0))
  sizes <- clust.equalsizes(a=a, b=b, c=c, N=N)
  angles <- clust.equalangles(N=N)
  list(r=r, sizes=sizes, angles=angles)
  
}

##' makeHelixCluster
##'
##' makeHelixCluster
##' @title makeHelixCluster
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
##' @param ... ignored arguments
##' @return list 
##' @author baptiste Auguie
makeHelixCluster <- function(N=5, R0=12e-3, pitch=15e-3, 
                             delta=pi/2, delta0=0, right=TRUE,
                             a=0.005, b=a, c=a,
                             angles=c("helix", "random", "fixed"),
                             seed=123, ...){

 hel <- helix(R0=R0, pitch=pitch, N=N, delta=delta, delta0=delta0, right=right)
 nr <- NROW(hel$angles)
 r <- hel$positions
 set.seed(seed) # always have same first particles
 angles <- switch(match.arg(angles),
                  "helix" = hel$angles,
                  "fixed" = cbind(rep(0, nr),
                    rep(0, nr),
                    rep(0, nr)),
                  "random" = matrix(cbind(runif(nr, 0, 2*pi),
                    runif(nr, 0, 2*pi),
                    runif(nr, 0, 2*pi)), ncol=3, byrow=T))
 
 sizes <- clust.equalsizes(a, b, c, N)
 list(r=r, sizes=sizes, angles=angles)
}


##' helical cluster z theta
##'
##' helical cluster z theta
##' @title makeHelixCluster.zt
##' @param z z
##' @param theta theta
##' @param R0 R0
##' @param pitch pitch 
##' @param delta0 delta0
##' @param right right
##' @param a a
##' @param b b
##' @param c c
##' @return cluster
##' @author baptiste Auguie
makeHelixCluster.zt <- function(z = runif(10, -0.5, 0.5), theta = runif(10, 0, 2 * pi),
                             R0 = 0.5, pitch = 6, delta0=0, 
                             right = TRUE, 
                             a=0.035, b=0.012, c=0.012){
  hel <- helix.zt(R0=R0, pitch=pitch, z=z, theta=theta, right=right)
 
  r <- hel$positions
  angles <- hel$angles
  N <- nrow(r)
  
  sizes <- clust.equalsizes(a, b, c, N)
   list(r=r, sizes=sizes, angles=angles)
  
}
