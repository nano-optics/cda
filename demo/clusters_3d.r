## generate 3D views of some predefined cluster shapes using RGL rendering or povray
library(cda)

radius <- 0.35
pitch <- 1.0
hel <- helix(N=20, R0=radius, pitch=pitch, delta=pi/9,delta0=2.5*pi/6,right=FALSE)

positions <- hel$positions
angles <- hel$angles 
sizes <- matrix(c(0.035, 0.012, 0.012), ncol=3,
                        nrow=nrow(positions), byrow=T)

particles.povray(positions,
                 angles, sizes, out="pos-helix-ell.pov")
curve.povray(hel$smooth, 0.5, radius=radius,out="pos-helix-ell.pov")

  euler <- system.file("povray", "euler.pov", package = "cda")
  helix.template.ell <- system.file("povray", "template-helix-ell.pov", package = "cda")
  file.copy(helix.template.ell, ".")
  file.copy(euler, ".")
  

if(interactive() && require(rgl)){ # display RGL window
  
  open3d()
  cl <- makeSpheresCluster(N=7, radius=5e-3, R0=12e-3, pitch=15e-3,
                           delta=pi/2, right=TRUE)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
  hel <- helix(N = 7, R0 = 12e-3, pitch = 15e-3, delta = pi/2, 
                     delta0 = 0, right = TRUE)
  lines3d(hel$smooth, lwd=3, col="blue")
  rgl.snapshot( "helix-rgl.png", fmt="png", top=TRUE )
  particles.povray(cl$r,
                   cl$angles,
                   cl$sizes, out="pos-helix.pov")
  curve.povray(hel$smooth, 0.5, radius=12e-3,out="pos-helix.pov")
  
  euler <- system.file("povray", "euler.pov", package = "cda")
  helix.template <- system.file("povray", "template-helix.pov", package = "cda")
  file.copy(helix.template, ".")
  file.copy(euler, ".")
  
  open3d()
  cl <- makeRodChain(N=10, pitch=0.15)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
  rgl.snapshot( "chain-rgl.png", fmt="png", top=TRUE )

  open3d()
  cl <- makeDimerCluster(d = 200e-3, phi=pi/4)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
  rgl.snapshot( "dimer-rgl.png", fmt="png", top=TRUE )
  particles.povray(cl$r,
                   cl$angles,
                   cl$sizes, out="pos-dimer.pov")
  
  axes <- system.file("povray", "axes.pov", package = "cda")
  euler <- system.file("povray", "euler.pov", package = "cda")
  dimer.template <- system.file("povray", "template-dimer.pov", package = "cda")
  file.copy(dimer.template, ".")
  file.copy(axes, ".")
  file.copy(euler, ".")
  
  print("now run povray template-helix.pov")
  print("now run povray template-helix-ell.pov")
  print("now run povray template-dimer.pov")

  }
