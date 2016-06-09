euler_xyz <- function(phi, theta, psi){
  
  Rot = matrix(NA, 3,3)
  cosphi = cos(phi); cospsi = cos(psi); costheta = cos(theta);
  sinphi = sin(phi); sinpsi = sin(psi); sintheta = sin(theta);
  
  Rot[1,1] = cosphi*costheta*cospsi - sinphi*sinpsi
  Rot[1,2] = sinphi*costheta*cospsi + cosphi*sinpsi
  Rot[1,3] = -sintheta*cospsi
  
  Rot[2,1] = -cosphi*costheta*sinpsi - sinphi*cospsi
  Rot[2,2] = -sinphi*costheta*sinpsi + cosphi*cospsi
  Rot[2,3] = sintheta*sinpsi
  
  Rot[3,1] = cosphi*sintheta
  Rot[3,2] = sinphi*sintheta
  Rot[3,3] = costheta
  
  Rot
  
}

library(rgl)
library(cd)


rgl.ellipsoid <- function (x=0,y=0,z=0, a = 1,b=1,c=1, phi=0,theta=0,psi=0,
                           subdivide = 3, smooth = TRUE, ...) 
{
  
  sphere <- rgl::subdivision3d(cube3d(...), subdivide)
  class(sphere) <- c("mesh3d","shape3d")
  
  norm <- sqrt(sphere$vb[1, ]^2 + sphere$vb[2, ]^2 + sphere$vb[3, 
                                                               ]^2)
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
  sphere$vb[4, ] <- 1
  sphere$normals <- sphere$vb
  result <- rgl::scale3d(sphere, a,b,c)
  rotM <- euler_xyz(psi, phi, 0*theta)
  result <- rgl::rotate3d(result,matrix=rotM)
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}

rgl.ellipsoids <- function(positions, sizes, angles,...){
  
  N <- NROW(positions)
  ll <- lapply(seq(1,N), function(ii)
    rgl.ellipsoid(positions[ii,1],positions[ii,2],positions[ii,3],
                  sizes[ii,1],sizes[ii,2],sizes[ii,3],
                  angles[ii,1],angles[ii,2],angles[ii,3], ...))
  
  rgl::shapelist3d(ll,...)
  
}

##' Add axes to a rgl scene
##'
##' x, y, z axes
##' @title rgl_annotate
##' @return draw axes
##' @author baptiste Auguie
##' @export
##' @family user_level rgl
rgl_annotate <- function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
  axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
  axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
  axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
  title3d('','','x axis','y axis','z axis')
}







cl1 <- cluster_shell(N=31, R0=10, a=1,b=0.2,c=0.2, 
                     orientation = "flat", position="fibonacci")
# open3d(FOV=2, windowRect=c(0,0,900,900), zoom=0.9)
rgl.spheres(0,0,0, radius=cl1$R0)
rgl.ellipsoids(cl1$r, cl1$sizes, cl1$angles, col="red")

