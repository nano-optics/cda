library(cda)
library(ggplot2)

wvl <- seq(500,900)
gold <- epsAu(wvl)

if(interactive() && require(rgl)){ # display RGL window
  cl <- makeDimerCluster(d = 200e-3, phi=pi/4)
  rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
}


onecluster <- function(handedness="right"){
  test <- handedness == "right"
  clust <- makeDimerCluster(d = 200e-3, phi=pi/4, right = test)
  m <- circular_dichroism_spectrum(clust, gold, n=1.33, N=36, progress=FALSE)
  invisible(m)
}

test <- onecluster()
params <- data.frame(handedness=c("right", "left"))
comparison <- mdply(params, onecluster, .progress="none")
## str(comparison)

p <- ggplot(comparison) + facet_grid(type~., scales="free") +
  geom_hline(yintercept=0) +
  geom_path(aes(energy, value, colour=handedness, linetype=variable,
                group=interaction(handedness, variable))) +
  labs(x="energy / eV", y=expression(sigma / mu*m^2))

p

