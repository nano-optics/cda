## generate 3D views of some predefined cluster shapes using RGL rendering or povray
library(cda)  
require(rgl)

cl1 <- cluster_dimer(d=100, 
                          dihedral=45*pi/180, alpha1=10*pi/180, alpha2=0,
                          a=35, b=12, 
                          right=TRUE)

open3d()
rgl.ellipsoids(cl1$r, cl1$sizes, cl1$angles, col="gold")
axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
  
cl2 <- cluster_chain(10, pitch=100, a=50, b=30)
open3d()
rgl.ellipsoids(cl2$r, cl2$sizes, cl2$angles, col="gold")
axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x','y','z')
  
  
cl3 <- cluster_helix(N=20, R0=500, pitch=1000, 
                          delta=pi/7, delta0=0, right=TRUE,
                          a=50, b=20, c=20,
                          angles=c("helix", "random", "fixed"),
                          seed=123)
hel3 <- helix(N = 20, R0 = 500, pitch = 1000, delta = pi/7, 
             delta0 = 0, right = TRUE)
  
open3d()
rgl.ellipsoids(cl3$r, cl3$sizes, cl3$angles, col="gold")
lines3d(hel3$smooth, lwd=1, col="red")
axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x','y','z')



  