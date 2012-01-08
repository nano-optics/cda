library(cda)
library(ggplot2)

wvl <- seq(400,700)
gold <- epsAu(wvl)

if(interactive() && require(rgl)){ # display RGL window
  open3d()
  cl <- makeSpheresCluster(N=7, radius=5e-3, R0=12e-3, pitch=15e-3,
                           delta=pi/2, right=TRUE)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
  hel <- helix(N = 7, R0 = 12e-3, pitch = 15e-3, delta = pi/2, 
                     delta0 = 0, right = TRUE)
  lines3d(hel$smooth, lwd=3, col="blue")
}


onecluster <- function(N,Nq=30){
  clust <- makeSpheresCluster(N=N, radius=5e-3, R0=12e-3, pitch=15e-3,
                              delta=pi/2, right=TRUE)
  m <- circular_dichroism_spectrum(clust, gold, n=1.33, N=Nq, progress=FALSE)
  m$N <- N
  invisible(m)
}

params <- data.frame(N=seq(2,7))
comparison <- mdply(params, onecluster, .progress="none")

p <- 
ggplot(comparison) + facet_grid(type~variable, scales="free") +
  geom_hline(yintercept=0) +
  geom_path(aes(energy, value*1e6, colour=factor(N))) + 
  labs(x="energy / eV", y=expression(sigma / mu*m^2), colour="N")

p


