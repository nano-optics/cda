library(cda)


r <- t(expand.grid(x=seq(-2, 2),
                   y=seq(0, 1),
                   z=seq(0, 1)))
N <- ncol(r)
kn <- 6/600
beta <- rep(1:3, N)
euler <- matrix(c(theta=0, phi=0, psi=0), ncol=N, nrow=3, byrow=FALSE)

A <- cda$interaction_matrix(r, kn, beta, euler, 0L)
D <- cda$block_diagonal(beta, euler)

source("~/Documents/github/cda/inst/standalone/functions.r")

A2 <- interaction_matrix(r, kn, beta, euler)
all.equal(A, A2)

all.equal(A[1:3,1:3], D[1:3,1:3])
all.equal(A2[4:6,4:6], D[4:6,4:6])

angles <- matrix(c(theta=0, phi=0, psi=0), ncol=4, nrow=3, byrow=FALSE)
Ei <- incident_field(c(1,0,0),k=c(0,0,1), r)
Ei2 <- cda$incident_field(c(1,0,0),k=c(0,0,1), r, angles)



