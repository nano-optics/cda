library(cda)



wavelength  =500;
epsilon=epsAg(wavelength)$epsilon;
medium=1.33;
a=80; b=50; c=30;
k0 = 2*pi/wavelength
kn = medium * k0;


V = 4 * pi/3 * a * b * c;
La = La(a, b, c);

axx = alpha_kuwata(wavelength, epsilon, V, a, La, medium);

d <- data.frame(wavelength=wavelength, re = Re(axx), im = Im(axx))
library(ggplot2)
library(reshape2)

cl <- list(r=rbind(c(0,0,0), c(150,0,0)),
           angles=rbind(c(0,0,0), c(0,0,0)),
           sizes=rbind(c(50,30,30), c(50,30,30)))

Beta <- inverse_polarizability(cl, material=epsAg(wavelength), medium=medium)
A <- cda$interaction_matrix(cl$r, kn, c(Beta), cl$angles, TRUE)
Adiag <- cda$block_diagonal(c(Beta), cl$angles)

Angles = rbind(c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0))

k0 = c(0,0,1)
E0 = c(0,1,0)
Ei <- cda$incident_field(E0, k=kn*k0, r=cl$r, Angles)
# Ei <- cda:::.incident_field(E0, k=kn*k0, r=cl$r, Angles)

P <- solve(A, Ei)
as.vector(cda$extinction(kn, P, Ei))
as.vector(cda$absorption(kn, P, Adiag))
