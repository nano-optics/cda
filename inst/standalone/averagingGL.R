

## norm1 <- function(x) Matrix::norm(matrix(x),'f')
require(Matrix)
require(cda)
require(plyr)
norm <- function(x) as.numeric(sqrt(t.default(x) %*% x))
# norm(c(1,1,0))
# 
# E0 <- c(1, 0, 0)
# k <- 6/600 *c(0,0,1)


cext <- function(k, E0, Ei, P){
  Re(4*pi* norm(k) /(E0%*%Conj(E0)) * Im(P%*%Conj(Ei)))
}

cext_avg <- function(k, E0, Ei, P){
  
  res <- rep(0, ncol(P))
  for (ii in seq.int(ncol(P)))
    res[ii] <- cext(k, E0, Ei[,ii], P[,ii])
  res
}


incident_field <- function(E0, k, r){
  kr <- crossprod(k, r)
  
  expikr <- exp(1i*kr)
  as.vector(matrix(c(E0[1]*expikr, E0[2]*expikr, 
                     E0[3]*expikr), nrow=3, byrow=T))  
}

incident_fields <- function(E0, k, r, angles){
  
  Ei <- matrix(0, ncol=ncol(angles), nrow=3*ncol(r))
  for (ii in seq(1, ncol(angles))){
    
    rot <- rotation(angles[,ii])
    k_r  <- t(rot) %*% k
    E0_r  <- t(rot) %*% E0
    Ei[,ii] <- incident_field(E0_r, k_r, r)
    
  }
  
  Ei
}

rotation <- function(x){
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  Ra <- matrix(c(cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, 0, 0, 1), ncol=3, byrow=T)
  Rb <- matrix(c(1, 0, 0, 0, cos(beta), sin(beta), 0, -sin(beta), cos(beta)), ncol=3, byrow=T)
  Rc <- matrix(c(cos(gamma), sin(gamma), 0, -sin(gamma), cos(gamma), 0, 0, 0, 1), ncol=3, byrow=T)
  
  return(Rc%*%Rb%*%Ra)
  
}

interaction_matrix <-  function(r, kn, beta, euler){
  
  N <- ncol(r)
  A <- matrix(0, 3*N, 3*N)
  
  I <- diag(3)
  
  for (jj in 1:N){
    for (kk in 1:N){
      
      if (jj != kk){          
        
        ind1 <- ((jj-1)*3+1):(jj*3)
        ind2 <- ((kk-1)*3+1):(kk*3)
        
        rk_to_rj <- r[,jj] - r[,kk]
        rjk <- norm(rk_to_rj)
        rjk_hat <- rk_to_rj / rjk
        rjkrjk <- outer(rjk_hat, rjk_hat)
        
        Ajk <- exp(1i*kn*rjk) / rjk * (kn^2*(rjkrjk - I) + (1i*kn*rjk-1) / rjk^2 * (3*rjkrjk - I))           
        A[ind1, ind2] <- Ajk
        
      } else {
        
        ind1 <- ((jj-1)*3+1):(jj*3)
        rot <- rotation(euler[,jj])
        A[ind1, ind1]  <- t(rot) %*% diag(beta[ind1]) %*% rot
        
      }
    }
  }
  invisible(A)
}



# 
# r <- t(expand.grid(x=seq(-2, 2),
#                    y=seq(0, 1),
#                    z=seq(0, 1)))
# N <- ncol(r)
# kn <- 6/600
# beta <- rep(1:3, N)
# euler <- matrix(c(theta=0, phi=0, psi=0), ncol=N, nrow=3, byrow=FALSE)

# A <- interaction_matrix(r, kn, beta, euler)

# Ei <- incident_field(E0, k, r)
# P <- jitter(Re(Ei)) + Ei
# cext(k, E0, Ei, P)



simulation <- function(wavelength, cluster, N=10){
  kn <- 2*pi/wavelength
  N <- ncol(cluster$r)
  alpha <- do.call(cbind, mapply(cda::polarizability_ellipsoid, 
                                 a=cluster$sizes[1,], b=cluster$sizes[2,], c=cluster$sizes[3,],
                                 MoreArgs=list(wavelength=wavelength, 
                                               epsilon=dielectric::epsAg(wavelength)$epsilon, 
                                               medium=1.0), SIMPLIFY=FALSE))
  beta <- 1/alpha
  
  angles <- cbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                 c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                 c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2
  
  rndN <- ceiling(sqrt(N/2))
  GL_phi <- statmod::gauss.quad(2*rndN)
  GL_psi <- statmod::gauss.quad(rndN)
#   
# scale the coordinates from (-1, 1) to (0, 2pi) and (-pi/2, pi/2) resp.
  phi1=0; phi2=2*pi; psi1=-pi/2; psi2=pi/2; 
  C1 = (phi2 - phi1) / 2;  D1 = (phi2 + phi1) / 2;
  C2 = (psi2 - psi1) / 2; D2 = (psi2 + psi1) / 2;
  
  phi = GL_phi$nodes*C1 + D1  
  psi = GL_psi$nodes*C2 + D2
  
  # grid of angles, theta is constant here
  grid <- expand.grid(phi=phi, theta=pi/2, psi=psi)
  # corresponding weights for 1D quadrature
  weights <- expand.grid(phi=GL_phi$weights, psi=GL_psi$weights)
  # combine the weigths for each point; cos(psi) comes from the Jacobian in the integral
  grid$weights <- C1 * C2 /(4*pi) * cos(grid$psi) * weights$phi * weights$psi
    
#   phi = QMC(ll, 1)*2*pi
#   psi = asin(2*QMC(ll, 0) - 1);
  
  res <- cd$circular_dichroism_spectrum(kn, invalpha, cluster$r, cluster$angles,
                                        as.matrix(cbind(GL$nodes, GL$weights)), 
                                        as.matrix(cbind(GL2$nodes, GL2$weights)),
                                        as.integer(full), as.integer(progress))
  
  
  A <- interaction_matrix(cluster$r, kn, beta, cluster$angles)
  E0L=c(0,1,1i)
  E0R=c(0,1i,1)
  k0=c(1,0,0)
  
  EiL <- incident_fields(E0L, k=kn*k0, r=cluster$r, angles)
  EiR <- incident_fields(E0R, k=kn*k0, r=cluster$r, angles)
 
  PL <- solve(A, EiL)
  PR <- solve(A, EiR)
#   browser()
  cext_avg(kn*k0, E0L, EiL, PL) - cext_avg(kn*k0, E0R, EiR, PR)
}

cluster <- list(r = cbind(c(0, 0, 0),
                          c(0, 0, 400)),
                angles = cbind(c(0, 0, 0),
                               c(pi/4, pi/2, 0)),
                sizes = cbind(c(40, 20, 20),
                              c(40, 20, 20)))
# simulation(500, cluster)
wavelength <- seq(200, 900)
# material <- epsAu(wavelength)

test <- sapply(wavelength, simulation, cluster=cluster)

plot(wavelength,test, t="l")

cluster2 <- lapply(cluster, t)
require(cda)
# test2 <- dispersion_spectrum(cluster2, c(0), "x", medium=1, 
#                              material=epsAu(wavelength))

test2 <- circular_dichroism_spectrum(cluster2, averaging="cheap", medium=1, 
                             material=epsAg(wavelength))

with(subset(test2, type == "CD" & variable == "extinction"), 
     lines(wavelength, value, col="red"))



alpha <- cda::polarizability_ellipsoid(wavelength, a=40, b=20, c=20,
                                       medium=1.0,
                                       epsilon = dielectric::epsAu(wavelength)$epsilon)


