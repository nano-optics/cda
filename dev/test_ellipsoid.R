
.euler <- function(phi, theta, psi){
  alpha <- phi
  beta <- theta
  gamma <- psi
  Ra <- matrix(c(cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, 0, 0, 1), ncol=3, byrow=T)
  Rb <- matrix(c(1, 0, 0, 0, cos(beta), sin(beta), 0, -sin(beta), cos(beta)), ncol=3, byrow=T)
  Rc <- matrix(c(cos(gamma), sin(gamma), 0, -sin(gamma), cos(gamma), 0, 0, 0, 1), ncol=3, byrow=T)

  return(Rc%*%Rb%*%Ra)

}


test_ellipsoid <- function (x=0,y=0,z=0, a = 3,b=1,c=1, phi=0,theta=0,psi=0,
                           subdivide = 3, smooth = TRUE, ..., manual=TRUE)
{
  
  rotM <- .euler(phi,theta,psi)
  print(rotM)

  sphere <- rgl::subdivision3d(cube3d(...), subdivide)
  class(sphere) <- c("mesh3d","shape3d")

  norm <- sqrt(sphere$vb[1, ]^2 + sphere$vb[2, ]^2 + sphere$vb[3,
                                                               ]^2)
  sphere$vb[1, ] <- sphere$vb[1, ]/norm
  sphere$vb[2, ] <- sphere$vb[2, ]/norm
  sphere$vb[3, ] <- sphere$vb[3, ]/norm
  sphere$vb[4, ] <- 1

  if(manual){
  sphere$vb[1, ] <- a*sphere$vb[1, ]
  sphere$vb[2, ] <- b*sphere$vb[2, ]
  sphere$vb[3, ] <- c*sphere$vb[3, ]
  sphere$vb[4, ] <- 1
#  sphere$vb[1:3,] <- t(t(sphere$vb[1:3,]) %*% rotM)
#    sphere$vb[1:3,] <- t(rotM) %*%  sphere$vb[1:3,]
  sphere$vb[1:3,] <- rotM %*%  sphere$vb[1:3,]
  sphere$normals <- sphere$vb
  result <- sphere
  } else {
    result <- rgl::scale3d(sphere, a,b,c)
    result <- rgl::rotate3d(result,matrix=rotM)
  }

print(str(result))
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}

rgl.close()
rgl.open()
e1 <- test_ellipsoid(phi=0,theta=0,psi=0, color='red')
e2 <- test_ellipsoid(phi=pi/4,theta=0,psi=0, color='green')
e3 <- test_ellipsoid(phi=0,theta=pi/2,psi=pi/4, color='blue')
view3d( theta = 0, phi = 90)
rgl::shapelist3d(list(e1,e2,e3))
