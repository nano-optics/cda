## ----start, message=FALSE, echo=FALSE, results='hide'--------------------
library(rgl)
library(knitr)
knitr::knit_hooks$set(rgl = hook_rgl)
knitr::read_demo("clusters", package="cda")
knitr::read_demo("dimer_cd", package="cda")

## ----visuals, rgl=TRUE, message=FALSE, echo=FALSE, results='hide'--------
library(cda)
library(rgl)

## random cluster with colours
cl1 <- function(N=50){
  R = 50;
  s = 4;
  cl = cluster_shell(N, R, s, 1, 1, 1,'random', 'hc',
                     exclusion = 10);
  rr <- runif(N, 0.6, 1.0)
  cl$sizes = s*rbind(rr, rr, runif(N, 1.2, 2.2));
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}

cl2 <- function(N=50){
  R = 15;
  s = 4;
  cl = cluster_shell(N, R, s, 1, 1, 3,  'radial', 'fibonacci');
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}

cl3 <- function(N=50){
  R = 15;
  s = 4;
  cl = cluster_shell(N, R, s, 1, 1, 3,  'flat', 'fibonacci');
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}

cl4 <- function(N=50){
  R = 15;
  s = 4;
  cl = cluster_shell(N, R, s, 1, 1, 3,  'random', 'fibonacci');
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}

cl5 <- function(){
  
  d = 50;
  a = 50; b=20; c=20;
  dihedral = pi/4;
  cl = cluster_dimer(d, a,b,c, dihedral,0,0);
  # open3d()
  visualise(cl, col="gold")
  # rgl.close()
}

cl6 <- function(N=20){
  cl = cluster_helix(N,
                     a=10, b=10, c=20,
                     R0=100, pitch=200,
                     delta=pi/5)
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours, show_core=FALSE)
  # rgl.close()
}

cl7 <- function(N=10){
  cl = cluster_chain(N, pitch=100, a=20, b=20, c=50)
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}

cl8 <- function(N=100){
  cl = cluster_array(N, pitch=200, a=50, b=20, c=20)
  colours = scales::hue_pal()(N)
  # open3d()
  visualise(cl, col=colours)
  # rgl.close()
}


open3d()
layout3d(matrix(1:8, 2, 4, byrow=TRUE))
cl1(); next3d()
cl2(); next3d()
cl3(); next3d()
cl4(); next3d()
cl5(); next3d()
cl6(); next3d()
cl7(); next3d()
cl8()
par3d(windowRect=c(0, 100, 800, 500))

## ----setup, results='hide', echo=FALSE-----------------------------------
library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
theme_set(theme_bw() + theme(strip.background=element_blank()))

## ----demo, results='hide'------------------------------------------------
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

## ----plot, echo=FALSE----------------------------------------------------

