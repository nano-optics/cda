
## @knitr load
library(cda)
library(rgl)
library(ggplot2)



theme_set(theme_minimal())


## @knitr cluster
# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

cluster_array <- function(N, pitch = 500, a = 50, b = 30, c = b, ...) 
{
  r <- as.matrix(expand.grid(x = seq(1, N) * pitch, y = seq(1, N) * pitch, z = 0))
  N2 <- NROW(r)
  sizes <- equal_sizes(a = a, b = b, c = c, N = N2)
  angles <- equal_angles(N = N2)
  list(r = r, sizes = sizes, angles = angles)
}

array <- function(N, pitch = 500, ...){
  cl <- cluster_array(N, pitch, ...)
  linear_extinction_spectrum(cluster = cl, material = gold, ...)
}
  
params <- data.frame(N=c(1, 10, 20))
#comparison <- mdply(params, array, progress=TRUE)

p <- 
  ggplot(data=comparison)+
labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm),
       colour = expression(N), linetype=expression(polarisation))+
  geom_line(aes(wavelength, value, linetype=variable,
                colour=factor(N),
                group=interaction(N,variable)))

#p

#save(comparison, params, file="cd.rda")


## @knitr comparison
data(G0)
pitch <- 500
S <- array_factor(wavelength=gold$wavelength / 1.33,
                   N=200, pitch=pitch)

S$smooth <- smooth.spline(S$wavelength, Re(S$S), df=50)$y + 
  1i * smooth.spline(S$wavelength, Im(S$S), df=50)$y

par(mfrow=c(2,1),mar=c(0,0,0,0))
plot(gold$wavelength, Re(S$S),t="l")
lines(gold$wavelength, Re(S$smooth),lty=2)
plot(gold$wavelength, Im(S$S),t="l")
lines(gold$wavelength, Im(S$smooth),lty=2)

interpolate.fun <- function(x, y){
  list(re=approxfun(x, Re(y)),
       im=approxfun(x, Im(y)))
}

gfun <- interpolate.fun(G0$wavelength, G0$Gxx)
G <- gfun$re(gold$wavelength/pitch/1.33) + 1i*gfun$im(gold$wavelength/pitch/1.33)
alpha_0 <- polarizability_ellipsoid(gold$wavelength, gold$epsilon, medium=1.33)[,1]
alpha_1 <- 1 / (1 / alpha_0 - S$S)
alpha_2 <- 1 / (1 / alpha_0 - S$smooth)
alpha_3 <- 1 / (1 / alpha_0 - G/pitch^3)
str(alpha)
k <- 2*pi * 1.33 / gold$wavelength
palette(RColorBrewer::brewer.pal(5,"Set1"))
par(mfrow=c(1,1),mar=c(2,2,1,1),lwd=2)
plot(gold$wavelength, k*Im(alpha_0),t="l", xlim=c(400,800),col="black")
#lines(gold$wavelength, k*Im(alpha_1),col=5)
lines(gold$wavelength, k*Im(alpha_2),col=1)
lines(gold$wavelength, k*Im(alpha_3),col=2)
with(subset(comparison, variable == "s" & N == 20), 
     lines(gold$wavelength, value/4/pi, lty=1,col=3))
legend("topleft", lty=c(1,1,1,1),col=c("black", 1,2,3), c("isolated", "truncated", "converged", "finite" ))
