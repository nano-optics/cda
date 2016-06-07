visualise_rgl <- function(cl, outfile=NULL, show_core=TRUE, ...){
  rgl.ellipsoids(cl$positions, cl$sizes, cl$angles, ...)
  
  if("R0" %in% names(cl) && show_core) 
    rgl::rgl.spheres(0,0,0, radius=cl$R0, col="grey", alpha=0.9)
  
  if(!is.null(outfile))
    rgl::rgl.snapshot( outfile, fmt = "png", top = TRUE )
}


##' creates an rgl ellipsoid
##'
##' deforms, rotate, and translate a sphere
##' @title rgl.ellipsoid
##' @param x x
##' @param y y
##' @param z z
##' @param a axis
##' @param b axis
##' @param c axis
##' @param phi phi
##' @param theta theta
##' @param psi psi
##' @param subdivide subdivision
##' @param smooth smoothing
##' @param ... additional params
##' @export
##' @return an rgl mesh
##' @author baptiste Auguie
##' @family user_level rgl
##' @examples
##' \dontrun{ require(rgl) ;  ee <- rgl.ellipsoid()
##' shapelist3d(ee) }
rgl.ellipsoid <- function (x=0, y=0, z=0, a = 1, b=1, c=1, phi=0, theta=0, psi=0,
                           subdivide = 3, smooth = TRUE, ...)
{
  
  sphere <- rgl::subdivision3d(rgl::cube3d(...), subdivide)
  class(sphere) <- c("mesh3d","shape3d")
  
  norm <- sqrt(sphere$vb[1, ]^2 + 
                 sphere$vb[2, ]^2 + 
                 sphere$vb[3, ]^2 )
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
  sphere$vb[4, ] <- 1
  sphere$normals <- sphere$vb
  result <- rgl::scale3d(sphere, a,b,c)
  rotM <- cpp_euler_passive(phi,theta,psi)
  result <- rgl::rotate3d(result,matrix=rotM)
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}

##' Create a list of rgl ellipsoids oriented in space
##'
##' each ellipsoid is specified by its position, dimensions, and Euler angles
##' @title rgl.ellipsoids
##' @param positions matrix of positions
##' @param sizes matrix of axis lengths
##' @param angles matrix of Euler angles
##' @param colour colour of each ellipsoid
##' @param ... additional params
##' @return rgl mesh
##' @author baptiste Auguie
##' @export
##' @family user_level rgl
##' @examples
##' cl <- helix(0.5, 1, 36, delta=pi/6, n.smooth=1e3)
##' sizes <- equal_sizes(0.04,0.02,0.02,NROW(cl$positions))
##' \dontrun{ require(rgl) ; rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold") }
rgl.ellipsoids <- function(positions, sizes, angles, colour = "red", ...){
  
  N <- NCOL(positions)
  colours <- rep(colour, length.out=N)
  ll <- lapply(seq(1,N), function(ii)
    rgl.ellipsoid(positions[1,ii],positions[2,ii],positions[3,ii],
                  sizes[1,ii],sizes[2,ii],sizes[3,ii],
                  angles[1,ii],angles[2,ii],angles[3,ii], col = colours[ii], ...))
  
  rgl::shapelist3d(ll)
  
}

##' Add axes to a rgl scene
##'
##' x, y, z axes
##' @title rgl_annotate
##' @return draw axes
##' @author baptiste Auguie
##' @noRd
##' @family user_level rgl
rgl_annotate <- function(){
  rgl::axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
  rgl::axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
  rgl::axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
  rgl::axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
  rgl::title3d('','','x axis','y axis','z axis')
}
