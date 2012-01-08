library(cda)
library(ggplot2)

gold <- epsAu(seq(400, 900))

## model the dimer using cda
onecluster <- function(d=0.1, dihedral=45, alpha1=0, alpha2=0,
                       a=0.04, b=0.015, n=1.5){

  scale <- 1
  clust <- makeDimerDihedral(d=d,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=a, b=b)
  
  xgeom <- 2*pi*( a*scale * b*scale* b*scale)^(2/3)
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=20,
                                   progress=FALSE, result.mat=FALSE)
 
  invisible(m)
}

cd <- onecluster()

cd2plot <- with(cd, data.frame(wavelength=wavelength*1e3, value=value))
cd2plot$variable <- cd$type
levels(cd2plot$variable) <- c("cd","ext")
cd2plot$type <- cd$variable

cd2plot <- subset(cd2plot, type == "extinction")

## load T-matrix results
load("lmax6-mmax3.rda")
## conversion factor from atomic units to metric
AUM <- 0.529177208e-4

p <- 
ggplot(l6m3)+facet_grid(variable~.,scale="free")+
  geom_path(aes(wavelength, value*AUM^2))  +
  geom_path(aes(x=wavelength,y=value),data=cd2plot,linetype=2) +
  geom_hline(aes(yintercept=y),data=data.frame(y=0))+
  geom_vline(aes(xintercept=x),data=data.frame(x=-Inf)) +
  scale_x_continuous(expression(wavelength/nm), expand=c(0,0))+ 
  scale_y_continuous(expression(sigma/mu*m^2), expand=c(0,0))+ 
  theme_minimal() + opts(strip.background=theme_blank(),
                         panel.margin = unit(1, "lines"))

p

ggsave("comparison.pdf",p)
