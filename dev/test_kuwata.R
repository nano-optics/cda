library(cda)



wavelength  =300:800;
epsilon=epsAg(wavelength)$epsilon;
medium=1.33;
a=80; b=50; c=30;


V = 4 * pi/3 * a * b * c;
La = La(a, b, c);

axx = alpha_kuwata(wavelength, epsilon, V, a, La, medium);

d <- data.frame(wavelength=wavelength, re = Re(axx), im = Im(axx))
library(ggplot2)
library(reshape2)

ggplot(melt(d, id="wavelength"), aes(wavelength, value, colour=variable)) + geom_line()


