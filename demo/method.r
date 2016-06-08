library(cda)
library(dielectric)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

theme_set(theme_bw() + theme(strip.background=element_blank()))


# dielectric function
gold <- epsAu(seq(400, 700))


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
cl <- cluster_dimer(d = 50, a=30, b=10, c=10, dihedral = pi/4)
gold <- epsAu(seq(500, 800))
test_fun <- function(method = "solve", ...){
  spectrum_oa(cl, gold, method=method, ...)
}

single <- spectrum_oa(cluster_single(a=30, b=10, c=10), gold, method="solve")
par <- expand.grid(method=c("cg","oos"), maxiter = seq(1,10, by=2), stringsAsFactors = FALSE)
exact <-  spectrum_oa(cl, gold, method="solve")
all <- mdply(par, test_fun, .progress="text")

p <- ggplot(subset(all, type == "cross-section" & variable == "extinction"),
            aes(wavelength, value)) + 
  facet_grid(variable~method, scales="free") +
  geom_line(aes(color=factor(maxiter))) +
  geom_line(data=subset(exact, type == "cross-section" & variable == "extinction"), lty=3) + 
  geom_line(data=subset(single, type == "cross-section" & variable == "extinction"), lty=2) + 
  labs(x = "wavelength /nm", y = expression(sigma/nm^2), 
       colour = "# iter")

p
