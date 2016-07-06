## ---- setup,echo=FALSE ---------------
library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
theme_set(theme_bw() + theme(strip.background=element_blank()))

## ---- demo ---------------

# dielectric function
gold <- epsAu(seq(400, 800))

# define a cluster of particles

model <- function(f=1, ...){
  
  cl <- cluster_dimer(d = 100*f, a=50*f, b=20*f, c=20*f, dihedral = pi/4)
  spectrum_oa(cl, gold, Nsca=50)
  
}

d1 <- model(0.1)
d2 <- model(1)


symmetrise_scale <- function(p, axis = c("y", "x"), panel=seq_along(facets)){
  axis <- match.arg(axis)
  gb <- ggplot_build(p)
  type <- switch(axis, "x" = "x.range", "y" = "y.range")
  fname <- as.character(p$facet[[1]])
  facets <- gb$panel$layout[[fname]]
  lims <- sapply(gb$panel$ranges, "[[", type)[,panel,drop=FALSE]
  lims2 <- as.vector(t(tcrossprod(apply(abs(lims), 2, max), c(-1,1))))
  dummy <- setNames(data.frame(rep(facets[panel], each=2), lims2), c(fname, axis))
  switch(axis, 
         "x" = p + geom_blank(data=dummy, aes(x=x, y=Inf), inherit.aes = FALSE), 
         "y" = p + geom_blank(data=dummy, aes(x=Inf, y=y), inherit.aes = FALSE))
}


p1 <- ggplot(d1, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line() + 
  labs(x = "Wavelength /nm", y = expression(sigma/nm^2), 
       colour = "") +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_brewer(palette="Set1", labels=parse(text=c('sigma[ext]','sigma[abs]','sigma[sca]'))) 

p2 <- p1 %+% d2 
pp1 <- symmetrise_scale(p1, "y", panel=2) + ggtitle("f=0.1")
pp2 <- symmetrise_scale(p2, "y", panel=2) + ggtitle("f=1")

gg <- cbind(ggplotGrob(pp1)[,-c(5,6)], ggplotGrob(pp2)[,-c(1,2)])
library(grid)
grid.newpage()
grid.draw(gg)
