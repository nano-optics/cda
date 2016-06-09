library(cd)
library(rgl)

wavelength <- seq(300, 800)
test <- polarizability_dye(wavelength, a=1, b=0.5, c=0.2)

test2 <- polarizability_ellipsoid(wavelength, epsAg(wavelength)$epsilon, 
                                  a = 0.8, b=0.8, c=0.8)

matplot(wavelength, Im(test2), t="l", col=1)
matlines(wavelength, Im(test), t="l", col=2)


wavelength <- seq(300, 600, length=50)
N <- 50
cl <- cluster_shell(N = N, d=1, R0 = 15,  a=1, b=1, c=1, 
                    position = "quasi-random", orientation = "radial")

cl1 <- cluster_shell(N = 2, d=1, R0 = 50, a=1, b=0, c=0, 
                    position = "quasi-random", orientation = "radial")
# cl1 <- cl
rgl.spheres(0,0,0,15)
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="red")

test1 <- spectrum_shell(cluster = cl1, core=FALSE, 
                        wavelength = wavelength, 
                        method = "ls",
                       averaging = "QMC", Nquad = 50,
                       medium=1, progress = TRUE)

test <- spectrum_shell(cluster = cl, core=TRUE, wavelength = wavelength, 
                       method = "ls", cg=FALSE,
                       averaging = "QMC", Nquad = 50,tol=1e-4, nmax=20,
                       medium=1.33, progress = TRUE)

library(ggplot2)
ggplot(test, aes(wavelength, value/N, colour=variable,group=variable)) + 
  geom_line() + 
  # geom_line(aes(wavelength, value/N), 
  # data=subset(test2, type == "cross section"), lty=3) +
  theme()
