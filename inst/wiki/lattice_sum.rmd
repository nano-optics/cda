## comparison of numerical approximation to the lattice sum for a 2D array of dipoles and its converged solution

library(cda)
library(ggplot2)

data(G0)

## number of dipoles for the numerical evaluation
N <- 100
pitch <- 0.6
lambda0 <- seq(0.45,1.0,length=300)
n <- 1.33
lambda <- lambda0 / n
lambdap <- lambda / pitch

S1 <- array_factor(wavelength=lambda,
                   N=N, pitch=pitch)


interpolate.fun <- function(x, y){
  list(re=approxfun(x, Re(y)),
       im=approxfun(x, Im(y)))
}

gfun <- interpolate.fun(G0$wavelength, G0$Gxx)


results <- data.frame(lambdap = S1$wavelength / pitch,
                      real = Re(S1$S)*pitch^3,
                      imag = Im(S1$S)*pitch^3,
                      method = "truncated")

javier <- data.frame(lambdap = lambdap,
                     real = gfun$re(lambdap),
                     imag = gfun$im(lambdap),
                     method = "converged")

p <- ggplot(melt(rbind(results, javier), id = c("lambdap", "method")),
            aes(lambdap, value, colour=variable, linetype=method)) +
  geom_vline(xintercept = c(1, sqrt(2)/2), linetype=2, colour="grey")+
  geom_path() + theme_bw() + labs(x=expression("reduced wavelength "*lambda/nh), y=expression(S/h^3))

p


