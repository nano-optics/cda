library(cda)
library(ggplot2)
library(gridExtra)

gold <- epsAu(seq(400, 900))

onecluster <- function(full=TRUE, d=100e-3,scale=0.5, dihedral=45, alpha1=0, alpha2=0,
                       a=0.04, b=0.015, n=1.5){
  
  clust <- makeDimerDihedral(d=scale*d,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=scale*a, b=scale*b)
  
  xgeom <- 2*pi*( a*scale * b*scale* b*scale)^(2/3)
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=20,
                                   progress=FALSE, result.mat=FALSE, full=full)
  m$value <- m$value / xgeom
  invisible(m)
}

p <- expand.grid(full=c(TRUE, FALSE))
all <- mdply(p, onecluster, dihedral=10, scale=1, d=100e-3)

all1 <- subset(all, variable == "extinction" & type == "CD")

p1 <- 
  ggplot(data=all1)+
  labs(y=expression(symbol("\341")~sigma[cd]~symbol("\361") / (2*pi*a[0]^2)),
           x=expression("Energy /eV"), linetype="") + 
  geom_line(aes(energy, value, linetype=full)) +
  geom_path(aes(x,y),data=data.frame(x=c(-Inf,Inf),y=c(0,0))) +
  theme_minimal(10) + opts(legend.position="none",
                         strip.text.y=theme_blank(),
                         strip.background=theme_blank()) +
  geom_path(aes(x,y),data=data.frame(x=c(-Inf,-Inf),y=c(-Inf,Inf)))  +
  scale_linetype_manual(values=c("dashed","solid"))

p1

ggsave(p1, file="fig2.pdf")
