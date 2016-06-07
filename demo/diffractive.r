library(cda)
library(rgl)
library(dielectric)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))

# dielectric function
wvl <- seq(500, 750, length=200)
gold <- epsAu(wvl)

model <- function(N, cluster, pitch = 500, ...){
  if(cluster=="chain")
  cl <- cluster_chain(N, pitch, a=30, b=50, c=30,...) else
    cl <- cluster_array(N, pitch, a=30, b=50, c=30,...) 
    
  spectrum_dispersion(cluster = cl, material = gold)
}

params <- expand.grid(N=c(1, 10, 50), 
                      cluster=c("chain","array"), stringsAsFactors = FALSE)
comparison <- mdply(params, model,.progress="text")

p <- 
  ggplot(data=subset(comparison, 
                     type == "cross-section" & variable == "extinction" & polarisation=="s"))+ 
  facet_grid(~cluster,scales="free")+
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm),
       colour = expression(N), linetype=expression(polarisation))+
  geom_line(aes(wavelength, value, 
                colour=factor(N),
                group=interaction(N,polarisation)))

p
