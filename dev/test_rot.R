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


euler_active <- function(phi, theta, psi){
  euler_xyz(-phi,-theta,-psi)
}

vi <- c(1,0,0)
vj <- c(0,1,0)
vk <- c(0,0,1)

euler_active(0,0,0)

## one axis
R1 <- euler_active(pi/3,0,0) # rotate along z
R2 <- euler_active(0, pi/3, 0)# rotate along x
R3 <- euler_active(0,0,pi/3)# rotate along z (same as R1)

R1 %*% vi ; matrix(c(1/2, sqrt(3)/2, 0))
R1 %*% vj ; matrix(c(-sqrt(3)/2, 1/2, 0))
R1 %*% vk ; matrix(vk)

R2 %*% vi ; matrix(vi)
R2 %*% vj ; matrix(c(0, 1/2, sqrt(3)/2))
R2 %*% vk ; matrix(c(0,0,1))

R3 %*% vi ; matrix(c(1/2, sqrt(3)/2, 0))
R3 %*% vj ; matrix(c(-sqrt(3)/2, 1/2, 0))
R3 %*% vk ; matrix(c(0,0,1))


## three axes
R123 <- euler_active(pi/3,pi/2,pi/3) # rotate along z then along x' then along z''

R123 %*% vi
R123 %*% vj 
R123 %*% vk 







