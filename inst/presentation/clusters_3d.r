## generate 3D views of some predefined cluster shapes using RGL rendering or povray
library(cda)

radius <- 0.3
pitch <- 0.4
hel <- helix(N=35, R0=radius, pitch=pitch, delta=pi/9,delta0=2.5*pi/6,right=FALSE)

positions <- hel$positions
angles <- hel$angles 
sizes <- matrix(c(0.035, 0.012, 0.012), ncol=3,
                        nrow=nrow(positions), byrow=T)

particles.povray(positions,
                 angles, sizes, out="pos-helix-ell.pov")
curve.povray(hel$smooth, 0.5, radius=radius,out="pos-helix-ell.pov")


