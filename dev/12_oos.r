setwd("~/Documents/plasmonics/cd/vignettes/userguide")

library(cd)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))

# dielectric function
wvl <- seq(500, 750, length=200)
gold <- epsAu(wvl)

model <- function(N, method="ls", maxiter=1, pitch = 500, material=gold, ...){
  
  cl <- cluster_chain(N, pitch, a=30, b=50, c=30) 
    
  tt <- system.time({
  d <- spectrum_dispersion(cluster = cl, 
                                      maxiter=maxiter, method = method, material = material, ...)
  })
  d$time <- tt[3]
  d
}

params <- rbind(expand.grid(N=c(1, 50, 100, 150, 200), 
                      method=c("ls"), 
                      maxiter = c(1),
                      stringsAsFactors = FALSE),
                expand.grid(N=c(1, 50, 100, 150, 200), 
                            method=c("oos"), 
                            maxiter = c(1:10),
                            stringsAsFactors = FALSE))


comparison <- mdply(params, model, material=gold, .progress="text",  tol = 1e-04)

p <- 
  ggplot(data=subset(comparison, method == "oos" &
                     type == "cross-section" & variable == "extinction" & polarisation=="s"))+ 
  facet_grid(N~.,scales="free")+
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm),
       colour = expression(N), linetype=expression(polarisation))+
  geom_line(aes(wavelength, value, 
                colour=factor(maxiter)))+
  geom_line(aes(wavelength, value), lty=2,
            data=subset(comparison, method == "ls" &
                        type == "cross-section" & variable == "extinction" & polarisation=="s"))

p

timings <- mdply(params, model, material = epsAu(650), .progress="text",  tol = 1e-02)

times <- ddply(timings, c("N","method", "maxiter"), summarise, t = unique(time))
p2 <- ggplot(subset(times, method=="oos"), aes(N, t)) + 
  geom_line(aes( colour=factor(maxiter))) +
  geom_line(data = subset(times, method=="ls")) 

library(gridExtra)
ggsave("12_oos.pdf", arrangeGrob(p, p2), width=8, height=4)

