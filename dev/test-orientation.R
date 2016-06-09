library(cd)
library(rgl)

wavelength <- seq(300, 800)
test <- polarizability_dye(wavelength, a=1, b=0.5, c=0.2)

test2 <- polarizability_ellipsoid(wavelength, epsAg(wavelength)$epsilon, 
                                  a = 0.8, b=0.8, c=0.8)

matplot(wavelength, Im(test2), t="l", col=1)
matlines(wavelength, Im(test), t="l", col=2)



wavelength <- seq(450, 600, length=50)
N <- 50
cl <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=1, c=1, 
                    position = "quasi-random", orientation = "radial")

cl1 <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=0, c=0, 
                    position = "quasi-random", orientation = "radial")

cl2 <- cluster_shell(N = N, d=1, R0 = 10,  a=1, b=0, c=0, 
                     position = "quasi-random", orientation = "flat")


# # cl1 <- cl
# rgl.spheres(0,0,0,15)
# rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="red")


test <- spectrum_shell(cluster = cl, core=FALSE, wavelength = wavelength, 
                       method = "ls", cg=FALSE,
                       averaging = "QMC", Nquad = 50,tol=1e-4, nmax=20,
                       medium=1.0, progress = TRUE)


test1 <- spectrum_shell(cluster = cl1, core=FALSE, wavelength = wavelength, 
                       method = "ls", cg=FALSE,
                       averaging = "QMC", Nquad = 50,tol=1e-4, nmax=20,
                       medium=1.0, progress = TRUE)

test2 <- spectrum_shell(cluster = cl2, core=T, wavelength = wavelength, 
                        method = "ls", cg=FALSE,
                        averaging = "QMC", Nquad = 50,tol=1e-4, nmax=20,
                        medium=1.0, progress = TRUE)


library(ggplot2)
ggplot(test2, aes(wavelength, value/N, colour=variable,group=variable)) + 
  geom_line() + geom_line(data=test2, lty=2) + 
  geom_line(data=test2, lty=3) + 
  # geom_line(aes(wavelength, value/N), 
  # data=subset(test2, type == "cross section"), lty=3) +
  theme()
