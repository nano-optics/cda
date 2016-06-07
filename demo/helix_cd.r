library(cda)
library(rgl)
library(dielectric)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


# dielectric function
wvl <- seq(400, 700)
gold <- epsAu(wvl)

ar <- 2
N <- 5
cl <- cluster_helix(N, R0=12, pitch=15, 
                    delta=pi/2, delta0=0, right=TRUE,
                    a=5/ar, b=5/ar, c=5,
                    angles="helix")
# visualise(cl, show_core=FALSE)

## ----comparison,echo=TRUE,tidy=FALSE,fig.path="helix-",fig.width=8-------

simulation <- function(N=3, scale=1, ar=1, ...){
  cl <- cluster_helix(N, R0=12*scale, pitch=15*scale, 
                      delta=pi/2, delta0=0, right=TRUE,
                      a=5/ar*scale, b=5/ar*scale, c=5*scale,
                      angles="helix")
  spectrum_oa(cl, material = gold, medium=1.33)
  
}


params <- expand.grid(N=seq(3, 7), ar= c(1, 1.1))
comparison <- mdply(params, simulation, .progress = "text")

symmetrise_scale <- function(p, axis = c("y", "x")){
  axis <- match.arg(axis)
  gb <- ggplot_build(p)
  type <- switch(axis, "x" = "x.range", "y" = "y.range")
  lims <- sapply(gb$panel$ranges, "[[", type)
  fname <- as.character(p$facet$facets)
  facets <- gb$panel$layout[[fname]]
  lims2 <- as.vector(t(tcrossprod(apply(abs(lims), 2, max), c(-1,1))))
  dummy <- setNames(data.frame(rep(facets, each=2), lims2), c(fname, axis))
  switch(axis, 
         "x" = p + geom_blank(data=dummy, aes(x=x, y=Inf)), 
         "y" = p + geom_blank(data=dummy, aes(x=Inf, y=y)))
}

p <- 
  ggplot(data=subset(comparison, type == "dichroism" & variable == "extinction")) + 
  facet_wrap(~ ar, scales="free" ,labeller=labeller(.cols=label_both)) +
  geom_line(aes(wavelength, value/N, 
                colour=factor(N))) +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="# particles") 

symmetrise_scale(p)


