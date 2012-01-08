library(cda)
library(ggplot2)

library(gridExtra)
gold <- epsAu(seq(400, 1200))

onecluster <- function(scale=1, dihedral=45, alpha1=0, alpha2=0, a=0.04, b=0.015, n=1.5){
  
  clust <- makeDimerDihedral(d=100e-3,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=scale*a, b=scale*b)
  
  xgeom <- 2*pi*( a*scale * b*scale* b*scale)^(2/3)
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=20,
                                   progress=FALSE, result.mat=FALSE)
  m$value <- m$value / xgeom
  invisible(m)
}

p <- data.frame(scale=c(0.1, 1, 1.5))
all <- mdply(p, onecluster, dihedral=45)

refactorize <- function(d, columns=names(d)) {
  d[columns] = lapply(d[columns],
     function(.f) as.factor(.f)[, drop=TRUE])
  d
}

all <- refactorize(all, c("scale","variable"))

all1 <- subset(all, type == "CD")
all2 <- subset(all, type == "cross section")

p1 <- 
  ggplot(data=all1)+
  facet_grid(scale~., scale="free") +
  geom_hline(data=data.frame(yintercept=0), aes(yintercept=yintercept)) +
  labs(y=expression(symbol("\341")~sigma[cd]~symbol("\361") / (2*pi*a[0]^2)),
           x=expression("Energy /eV"), colour="") + 
  geom_line(aes(energy, value, colour=variable))+
  theme_minimal(10) + opts(legend.position="none",
                         strip.text.y=theme_blank(),
                         strip.background=theme_blank()) +
  geom_path(aes(x,y),data=data.frame(x=c(-Inf,-Inf),y=c(-Inf,Inf)))+
  xlim(1,2.6)+
  scale_x_continuous(breaks=seq(1,2.6,by=0.2),expand=rep(0,2)) +
  opts(plot.margin=unit(c(0.5,0,0,0),"lines"))

p2 <- p1 %+% all2 +
  ylab(expression(symbol("\341")~sigma~symbol("\361") / (2*pi*a[0]^2)))


g <- arrangeGrob(p1, p2, ncol=2)
if(interactive())
  grid.draw(g)

ggsave(g, file="fig3.pdf")

