
## @knitr load
library(cda)
library(rgl)


## @knitr setup

rgl_annotate = function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
}


## @knitr manual

cl1 <- list(r = rbind(c(0, 0, 0),
                      c(0, 500, 0)),
            angles = rbind(c(0, 0, 0),
                           c(0, pi/2, pi/3)),
            sizes = rbind(c(50, 20, 20),
                          c(40, 30, 30)))
rgl.ellipsoids(cl1$r, cl1$sizes, cl1$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()
                      


## @knitr dimer
cl2 <- cluster_dimer(d=100, 
                          dihedral=45*pi/180, alpha1=10*pi/180, alpha2=0,
                          a=35, b=12)

rgl.ellipsoids(cl2$r, cl2$sizes, cl2$angles, col="gold")
rgl_annotate()


## @knitr chain
cl2 <- cluster_chain(10, pitch=100, a=50, b=30)
rgl.ellipsoids(cl2$r, cl2$sizes, cl2$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 30, fov = 70, zoom = 1)
rgl_annotate()


## @knitr helix
cl3 <- cluster_helix(N=20, R0=500, pitch=1000, 
                          delta=pi/7, delta0=0, right=TRUE,
                          a=100, b=50, c=20,
                          angles="helix")
hel3 <- helix(N = 20, R0 = 500, pitch = 1000, delta = pi/7, 
             delta0 = 0, right = TRUE)
  
rgl.ellipsoids(cl3$r, cl3$sizes, cl3$angles, col="gold")
lines3d(hel3$smooth, lwd=1, col="red")
rgl_annotate()


