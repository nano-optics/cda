library(cda)
library(ggplot2)

library(gridExtra)

gold <- epsAu(seq(400, 900))

onecluster <- function(scale=1, dihedral=10, alpha1=10, alpha2=20, a=0.04, b=0.015, n=1.5){
  
  clust <- makeDimerDihedral(d=scale*100e-3,
                             dihedral*pi/180,alpha1*pi/180,alpha2*pi/180,
                             a=scale*a, b=scale*b)
  
  xgeom <- 2*pi*( a*scale * b*scale* b*scale)^(2/3)
  m <- circular_dichroism_spectrum(clust, gold, n=n, N=20,
                                   progress=FALSE, result.mat=FALSE)
  m$value <- m$value / xgeom
  invisible(m)
}

p <- data.frame(scale=c(0.1, 1))
all <- mdply(p, onecluster, dihedral=10)

all <- subset(all, energy <= 2.6 & energy >= 1.4)
all1 <- subset(all, scale == 0.1 & !(type == "CD" & variable != "extinction"))
all2 <- subset(all, scale == 1 & !(type == "CD" & variable != "extinction"))

symmetrise <-  function(.d) {
  .dd <- subset(.d, type=="CD")
  lim <- max(abs(.dd$value))
  data.frame(x=2, y=lim*c(-1, 1),
             type=factor("CD",
               levels=c("CD", "cross section")))
}


grid.pretty <- function(range)
  pretty(range, n=5, min.n =5)

p1 <- 
  ggplot(data=all1)+
  facet_grid(type~., scale="free") +
  geom_hline(data=data.frame(yintercept=0), aes(yintercept=yintercept)) +
  geom_blank(data=symmetrise(all1), aes(x=x, y=y))+
  labs(y=expression(symbol("\341")~sigma~symbol("\361") / (2*pi*a[0]^2)),
           x=expression("Energy /eV"), colour="") + 
  geom_line(aes(energy, value, colour=variable))+
  theme_minimal(10) + opts(legend.position="none",
                         strip.text.y=theme_blank(),
                         strip.background=theme_blank()) +
  geom_path(aes(x,y),data=data.frame(x=c(-Inf,-Inf),y=c(-Inf,Inf))) +
  scale_x_continuous(expand=c(0,0),lim=c(1.4,2.6)) +
  opts(plot.margin=unit(c(0,0,0,0),"lines"))

p2 <- p1 %+% all2 + opts(legend.position=c(0.7, 0.5)) +
  geom_blank(data=symmetrise(all2), aes(x=x, y=y))


g <- arrangeGrob(p1, p2, ncol=2)
if(interactive())
  grid.draw(g)

ggsave(g, file="fig1.pdf")
