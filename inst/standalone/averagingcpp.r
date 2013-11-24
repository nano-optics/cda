

require(cda)
## norm1 <- function(x) Matrix::norm(matrix(x),'f')
require(Matrix)
require(plyr)
norm <- function(x) as.numeric(sqrt(t.default(x) %*% x))
# norm(c(1,1,0))
# 
# E0 <- c(1, 0, 0)
# k <- 6/600 *c(0,0,1)

incident_field <- function(E0, k, r){
  kr <- crossprod(k, r)
  
  expikr <- exp(1i*kr)
  as.vector(matrix(c(E0[1]*expikr, E0[2]*expikr, 
                     E0[3]*expikr), nrow=3, byrow=T))  
}

cext <- function(k, E0, Ei, P){
  Re(4*pi* norm(k) /(E0%*%Conj(E0)) * Im(P%*%Conj(Ei)))
}

# cext_avg <- function(k, E0, Ei, P){
#   
#   res <- 0
#   for (ii in seq.int(ncol(P)))
#     res <- res + cext(k, E0, Ei[,ii], P[,ii])
#   res
# }

# double absorption(const double kn, const arma::cx_colvec& P, const arma::cx_mat& diagBeta)
# {
#   const double c = 4*arma::math::pi()*kn*(imag(cdot(diagBeta * P, P)) -  \
#                                           kn*kn*kn* 2/3 * real(cdot(P, P))); 
#   return c;
# }

cabs_avg <- function(kn, Alpha, P){
  Eexc <- c(Alpha %*% P)
  
  4*pi* kn * Re(Im(Conj(Eexc) %*% c(P)) - kn^3*2/3* Conj(c(P))%*%c(P))/ ncol(P)
}

cext_avg <- function(k, E0, Ei, P){
  
  Re(4*pi* norm(k) /(E0%*%Conj(E0)) * Im(c(P)%*%Conj(c(Ei)))) / ncol(P)
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
cda$euler(0,pi/2,pi/3)
rotation(c(0,pi/2,pi/3))

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

diagonal_blocks <-  function(r, kn, beta, euler){
  
  N <- ncol(r)
  A <- matrix(0, 3*N, 3*N)
  
  
  for (jj in 1:N){
    for (kk in 1:N){
      
      if (jj == kk){          
        
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



simulation <- function(wavelength, cluster){
  kn <- 2*pi/wavelength
  N <- ncol(cluster$r)
  alpha <- do.call(cbind, mapply(cda::polarizability_ellipsoid, a=cluster$sizes[1,], b=cluster$sizes[2,], c=cluster$sizes[3,], MoreArgs=list(wavelength=wavelength, epsilon=dielectric::epsAg(wavelength)$epsilon, medium=1.0), SIMPLIFY=FALSE))
  beta <- 1/alpha
  
  angles <- cbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                 c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                 c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2
#   angles <- angles*0
  A <- interaction_matrix(cluster$r, kn, beta, cluster$angles)
  Adiag <- diagonal_blocks(cluster$r, kn, beta, cluster$angles)
  E0L=1/sqrt(2)*c(0,1,1i)
  E0R=1/sqrt(2)*c(0,1i,1)
  k0=c(1,0,0)
  
  EiL <- incident_fields(E0L, k=kn*k0, r=cluster$r, angles)
  EiR <- incident_fields(E0R, k=kn*k0, r=cluster$r, angles)
 
  PL <- solve(A, EiL)
  PR <- solve(A, EiR)
  
  eL <- cext_avg(kn*k0, E0L, EiL, PL)
  eR <- cext_avg(kn*k0, E0R, EiR, PR)
  aL <- cabs_avg(kn, Adiag, PL)
  aR <- cabs_avg(kn, Adiag, PR)  
  
  c(extinction = 0.5*(eL+eR),
    absorption = 0.5*(aL+aR),
    CDext=eL-eR,
    CDa=aL-aR)
}

simulationcpp <- function(wavelength, cluster){
  kn <- 2*pi/wavelength
  N <- ncol(cluster$r)
  alpha <- do.call(cbind, mapply(cda::polarizability_ellipsoid, a=cluster$sizes[1,], b=cluster$sizes[2,], c=cluster$sizes[3,], MoreArgs=list(wavelength=wavelength, epsilon=dielectric::epsAg(wavelength)$epsilon, medium=1.0), SIMPLIFY=FALSE))
  beta <- 1/alpha
  
  angles <- cbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                  c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                  c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2
# angles <- angles*0
  A <- cda$interaction_matrix(cluster$r, kn, beta, cluster$angles, 1L)
  Adiag <- diagonal_blocks(cluster$r, kn, beta, cluster$angles)
  E0L=1/sqrt(2)*c(0,1,1i)
  E0R=1/sqrt(2)*c(0,1i,1)
  k0=c(1,0,0)
#   
#   EiL <- incident_fields(E0L, k=kn*k0, r=cluster$r, angles)
#   EiR <- incident_fields(E0R, k=kn*k0, r=cluster$r, angles)
#   
#   PL <- solve(A, EiL)
#   PR <- solve(A, EiR)
#   
#   eL <- cda$extinction(kn, PL, EiL)
#   eR <- cda$extinction(kn, PR, EiR)
#   aL <- cda$absorption(kn, PL, Adiag)
#   aR <- cda$absorption(kn, PR, Adiag)
  
  axis_angles <- c(0, pi/2, pi/2)
  axis_axes <- c(0L, 1L, 2L)
  
  test <- dispersion$dispersion(cluster$r, A, beta, kn, 
                                axis_angles, axis_axes, 
                                cluster$angles, 1L)
  aa <- colMeans(test)
  weights <- rep(1, length=ncol(angles))
  tmp <- cd$averaging2(cluster$r, A, Adiag, kn, angles, weights)
#   eL <- tmp[1,1]
  eL <- aa[1]
  eR <- tmp[3,1]
  aL <- tmp[2,1]
  aR <- tmp[4,1]
  
  c(extinction = 0.5*(eL+eR),
    absorption = 0.5*(aL+aR),
    CDext=eL-eR,
    CDa=aL-aR)
  
}

cluster <- list(r = cbind(c(0, 0, 0),
                          c(0, 0, 400)),
                angles = cbind(c(0, 0, 0),
                               c(pi/6, pi/2, 0)),
                sizes = cbind(c(40, 20, 20),
                              c(40, 20, 20)))

cluster$sizes <- cluster$sizes*1
# simulation(500, cluster)
wavelength <- data.frame(wavelength=seq(300, 500))
# material <- epsAu(wavelength)

test <- mdply(wavelength, simulation, cluster=cluster)
test2 <- mdply(wavelength, simulationcpp, cluster=cluster)
library(ggplot2)
library(reshape2)
m <- melt(test, id="wavelength")
m2 <- melt(test2, id="wavelength")

p <- ggplot(m, aes(wavelength, value, colour=variable))+
  facet_grid(variable~., scales="free")+
  geom_line(size=2, alpha=0.5) +
  geom_line(data=m2, linetype=2, size=2)

print(p)
# 
# cluster2 <- lapply(cluster, t)
# # test2 <- dispersion_spectrum(cluster2, c(0), "x", medium=1, 
# #                              material=epsAu(wavelength))
# 
# test2 <- circular_dichroism_spectrum(cluster2, averaging="cheap", medium=1, 
#                              material=epsAg(wavelength))
# 
# with(subset(test2, type == "CD" & variable == "extinction"), 
#      lines(wavelength, value, col="red"))
# 
# 
# 
# alpha <- cda::polarizability_ellipsoid(wavelength, a=40, b=20, c=20,
#                                        medium=1.0,
#                                        epsilon = dielectric::epsAu(wavelength)$epsilon)


