library(cda)
library(ggplot2)

library(gridExtra)

gold <- epsAu(seq(400, 900))

## uncomment to run the code (takes a couple of minutes)
run <- FALSE

ellipsoid <- function(a0=0.05, ar=1){

  V <- 4/3*pi * a0^3
  b <- c <- a0/(ar^(1/3))
  a <- ar*b
  data.frame(a=a, b=b, c=c, V=V)
}


onecluster <- function(a0=0.04, ar=0.015, dihedral=45, alpha1=0, alpha2=0, n=1.5){
  ell <- ellipsoid(a0, ar)
  
  clust <- makeDimerDihedral(d=100e-3,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=ell$a, b=ell$b)
  
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=20,
                                   progress=FALSE, result.mat=TRUE)
  ## wavelength, extinction, absorption, CD ext, CD abs
  data.frame(wavelength=m[,1], g = m[,4] / m[,2], extinction=m[,2], CD=m[,4])
}

p <- expand.grid(ar=seq(1.0, 2.5, length=100), a0=c(0.01, 0.02, 0.03, 0.04))

refactorize <- function(d, columns=names(d)) {
  d[columns] = lapply(d[columns],
     function(.f) as.factor(.f)[, drop=TRUE])
  d
}

if(run){
all <- mdply(p, onecluster, .progress="text")

all <- refactorize(all, c("scale"))

p0 <- 
ggplot(all)+ facet_grid(a0~., scales="free")+
  geom_path(aes(wavelength, CD, colour=ar)) 

gmax <- ddply(all, .(ar, a0), summarise,
              maxi=max(g), mini=min(g), optim=max(abs(g)))
gmax <- refactorize(gmax, c("a0"))

save(all, file="all-spectra.rda")
save(gmax, file="gmax.rda")
} else load("gmax.rda")

p <- 
ggplot(gmax) + 
  geom_path(aes(ar, optim, colour=a0)) +
  theme_bw()+
  labs(x="Aspect ratio", y="Anisotropy factor")+
  opts(legend.position="none")

p

ggsave(p, file="fig4.pdf")
