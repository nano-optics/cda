

## norm1 <- function(x) Matrix::norm(matrix(x),'f')
require(Matrix)
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
  Re(4*pi* norm(k) /(E0%*%Conj(E0)) * Im(P%*%(Conj(Ei))))
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

cluster <- list(r = cbind(c(0, 0, 0),
                          c(0, 0, 200)),
                angles = cbind(c(0, 0, 0),
                               c(pi/5, pi/3, 0)),
                sizes = cbind(c(40, 30, 30),
                              c(50, 20, 20)))


simulation <- function(wavelength, cluster, E0=c(1,0,0), k0=c(0,0,1)){
  kn <- 2*pi/wavelength
  N <- ncol(cluster$r)
  alpha <- do.call(cbind, mapply(cda::polarizability_ellipsoid, a=cluster$sizes[1,], b=cluster$sizes[2,], c=cluster$sizes[3,], MoreArgs=list(wavelength=wavelength, epsilon=dielectric::epsAu(wavelength)$epsilon, medium=1.0), SIMPLIFY=FALSE))
  beta <- 1/alpha
  A <- interaction_matrix(cluster$r, kn, beta, cluster$angles)
  Ei <- incident_field(E0, k=kn*k0, r=cluster$r)
  P <- solve(A, Ei)
  cext(kn*k0, E0, Ei, P)
}

# simulation(500, cluster)
wavelength <- seq(300, 900)
# material <- epsAu(wavelength)

test <- sapply(wavelength, simulation, cluster=cluster, E0=c(0,1,0))

plot(wavelength,test, t="l")

cluster2 <- lapply(cluster, t)
require(cda)
test2 <- dispersion_spectrum(cluster2, c(0), "x", medium=1, 
                             material=epsAu(wavelength))

with(subset(test2, type == "extinction" & polarisation == "s"), 
     lines(wavelength, value, col="red"))



alpha <- cda::polarizability_ellipsoid(wavelength, a=40, b=20, c=20,
                                       medium=1.0,
                                       epsilon = dielectric::epsAu(wavelength)$epsilon)


